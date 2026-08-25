use std::error::Error;
use std::path::Path;

use clap::{Parser, Subcommand};

use hmnrandomread::{
    fusion, insert_size, BuiltProfile, Config, FastqReader, FastqRecord, FastqWriter,
    FusionConfig, FusionGenerator, Generator, InsertSizeStats, ProfileDiversity, ProfileSequencer,
    RandomGenerator, Reference,
};

const APP_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser)]
#[command(
    author = APP_NAME,
    version = VERSION,
    about = format!("{APP_NAME} CLI"),
    long_about = None,
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Simulate paired-end FASTQ reads from one or more reference genomes.
    Simulate {
        /// Reference path with an optional read count and diversity id:
        /// `path[,nb_reads[,id_diversity]]`. May be repeated.
        #[arg(long, action = clap::ArgAction::Append)]
        input_reference_fasta: Vec<String>,

        /// Diversity (SNP/indel) profile CSV.
        #[arg(long)]
        input_profile_diversity_csv: Option<String>,

        /// Sequencer error profile CSV.
        #[arg(long)]
        input_profile_sequencer_csv: Option<String>,

        /// Forward read output FASTQ (gzip-compressed).
        #[arg(long)]
        output_forward_fastq: String,

        /// Reverse read output FASTQ (gzip-compressed).
        #[arg(long)]
        output_reverse_fastq: String,

        /// Read length.
        #[arg(long, default_value_t = 150)]
        parameter_length_reads_int: usize,

        /// Mean fragment insert size.
        #[arg(long, default_value_t = 500)]
        parameter_mean_insert_int: u32,

        /// Standard deviation of the fragment insert size.
        #[arg(long, default_value_t = 50)]
        parameter_std_insert_int: u32,

        /// Identifier to select within --input-profile-sequencer-csv. Required
        /// if --input-profile-sequencer-csv is set.
        #[arg(long)]
        parameter_profile_sequencer_id_str: Option<String>,

        /// Seed for the random number generator.
        #[arg(long, default_value_t = 0)]
        parameter_seed_int: u64,
    },
    /// Build a sequencer error profile CSV (for `--input-profile-sequencer-csv`)
    /// from real paired FASTQ files or a BAM.
    BuildProfileSequencer {
        /// Identifier for the produced profile's rows (matches
        /// `--parameter-profile-sequencer-id-str` for `simulate`).
        #[arg(long)]
        parameter_id_str: String,

        /// Forward read FASTQ(s) (plain or gzip-compressed), space-separated
        /// (e.g. one per lane). Requires --input-reverse-fastq with the
        /// same count, in matching order; mutually exclusive with
        /// --input-bam.
        #[arg(long, num_args = 1..)]
        input_forward_fastq: Vec<String>,

        /// Reverse read FASTQ(s), space-separated, matching
        /// --input-forward-fastq's count and order; mutually exclusive with
        /// --input-bam.
        #[arg(long, num_args = 1..)]
        input_reverse_fastq: Vec<String>,

        /// BAM(s), space-separated; mutually exclusive with the FASTQ
        /// inputs.
        #[arg(long, num_args = 1..)]
        input_bam: Vec<String>,

        /// Output profile sequencer CSV.
        #[arg(long)]
        output_profile_sequencer_csv: String,
    },
    /// Report the mean and standard deviation of the fragment insert size
    /// (for `--parameter-mean-insert-int`/`--parameter-std-insert-int`) from
    /// one or more real BAMs.
    StatisticsInsertSize {
        /// BAM(s), space-separated; each is reported as its own row.
        #[arg(long, num_args = 1.., required = true)]
        input_bam: Vec<String>,

        /// Output CSV: `file`, `mean_insert_size`, `std_insert_size`, one row
        /// per --input-bam.
        #[arg(long)]
        output_statistics_csv: Option<String>,
    },
    /// Spike a real sample's paired FASTQ with synthetic reads supporting a
    /// gene fusion / structural-variant breakpoint, at a depth proportional
    /// to the real coverage already present at the primary breakpoint.
    FusionInSample {
        /// Reference path with an optional read count and diversity id:
        /// `path[,nb_reads[,id_diversity]]`. May be repeated. `nb_reads` is
        /// ignored — the number of fusion reads produced is derived from
        /// the real depth at --parameter-breakpoint-primary-roi times
        /// --parameter-rate-float instead.
        #[arg(long, action = clap::ArgAction::Append)]
        input_reference_fasta: Vec<String>,

        /// Diversity (SNP/indel) profile CSV, applied to the produced
        /// fusion reads.
        #[arg(long)]
        input_profile_diversity_csv: Option<String>,

        /// Sequencer error profile CSV, applied to the produced fusion reads.
        #[arg(long)]
        input_profile_sequencer_csv: Option<String>,

        /// The real sample's forward FASTQ (plain or gzip-compressed) to
        /// spike with fusion reads.
        #[arg(long)]
        input_forward_fastq: String,

        /// The real sample's reverse FASTQ (plain or gzip-compressed) to
        /// spike with fusion reads.
        #[arg(long)]
        input_reverse_fastq: String,

        /// Indexed (.bai/.csi sidecar) BAM of the same sample, used to read
        /// the pileup depth at --parameter-breakpoint-primary-roi.
        #[arg(long)]
        input_bam: String,

        /// Primary breakpoint, as `chrom:pos` (1-based). `pos` is the last
        /// reference base kept on the 5' side of the junction; the base
        /// right after it starts the 3' side.
        #[arg(long)]
        parameter_breakpoint_primary_roi: String,

        /// Secondary breakpoint, as `chrom:pos` (1-based), same convention
        /// as --parameter-breakpoint-primary-roi.
        #[arg(long)]
        parameter_breakpoint_secondary_roi: String,

        /// Fraction of the real depth at the primary breakpoint to produce
        /// as fusion-supporting read pairs. Range: 0.0-0.5.
        #[arg(long)]
        parameter_rate_float: f64,

        /// Fraction of the produced fusion read pairs assigned to the
        /// reciprocal junction orientation (secondary breakpoint's upstream
        /// sequence on the left, primary's downstream sequence on the
        /// right), rather than the primary->secondary orientation. Range:
        /// 0.0-1.0. 0.5 (default) splits evenly, as for a balanced
        /// reciprocal translocation; 0.0 produces only the
        /// primary->secondary junction (e.g. for an unbalanced fusion where
        /// only one derivative is relevant); 1.0 produces only the
        /// reciprocal junction.
        #[arg(long, default_value_t = 0.5)]
        parameter_reciprocal_rate_float: f64,

        /// Forward read output FASTQ (gzip-compressed): the input forward
        /// FASTQ plus the produced fusion reads.
        #[arg(long)]
        output_forward_fastq: String,

        /// Reverse read output FASTQ (gzip-compressed): the input reverse
        /// FASTQ plus the produced fusion reads.
        #[arg(long)]
        output_reverse_fastq: String,

        /// Read length.
        #[arg(long, default_value_t = 150)]
        parameter_length_reads_int: usize,

        /// Mean fragment insert size.
        #[arg(long, default_value_t = 500)]
        parameter_mean_insert_int: u32,

        /// Standard deviation of the fragment insert size.
        #[arg(long, default_value_t = 50)]
        parameter_std_insert_int: u32,

        /// Identifier to select within --input-profile-sequencer-csv. Required
        /// if --input-profile-sequencer-csv is set.
        #[arg(long)]
        parameter_profile_sequencer_id_str: Option<String>,

        /// Seed for the random number generator.
        #[arg(long, default_value_t = 0)]
        parameter_seed_int: u64,
    },
    /// Display the application version.
    Version,
}

fn main() {
    env_logger::init();
    let cli = Cli::parse();
    match &cli.command {
        Commands::Simulate { .. } => std::process::exit(cmd_simulate(&cli.command)),
        Commands::BuildProfileSequencer { .. } => {
            std::process::exit(cmd_build_profile_sequencer(&cli.command))
        }
        Commands::StatisticsInsertSize { .. } => {
            std::process::exit(cmd_statistics_insert_size(&cli.command))
        }
        Commands::FusionInSample { .. } => std::process::exit(cmd_fusion_in_sample(&cli.command)),
        Commands::Version => println!("{VERSION}"),
    }
}

/// Open every `--input-reference-fasta` spec into a [`Reference`].
fn load_references(specs: &[String], min_scaffold_len: u64) -> Result<Vec<Reference>, String> {
    let mut references = Vec::with_capacity(specs.len());
    for spec in specs {
        let (path, nb_reads, id_diversity) = Reference::parse_spec(spec)
            .map_err(|e| format!("invalid --input-reference-fasta '{spec}': {e}"))?;
        if !Path::new(&path).is_file() {
            return Err(format!("reference file not found: {path}"));
        }
        let reference = Reference::open(path.clone(), nb_reads, id_diversity, min_scaffold_len)
            .map_err(|e| format!("failed to open reference '{path}': {e}"))?;
        references.push(reference);
    }
    Ok(references)
}

/// Load the optional `--input-profile-diversity-csv`/`--input-profile-sequencer-csv`.
fn load_profiles(
    input_profile_diversity_csv: &Option<String>,
    input_profile_sequencer_csv: &Option<String>,
    parameter_profile_sequencer_id_str: &Option<String>,
) -> Result<(Option<ProfileDiversity>, Option<ProfileSequencer>), String> {
    let profile_diversity = match input_profile_diversity_csv {
        None => None,
        Some(path) => {
            if !Path::new(path).is_file() {
                return Err(format!("profile diversity file not found: {path}"));
            }
            let profile = ProfileDiversity::parse_csv(path)
                .map_err(|e| format!("failed to parse profile diversity '{path}': {e}"))?;
            Some(profile)
        }
    };

    let profile_sequencer = match input_profile_sequencer_csv {
        None => None,
        Some(path) => {
            let Some(id) = parameter_profile_sequencer_id_str else {
                return Err(
                    "--parameter-profile-sequencer-id-str is required when \
                     --input-profile-sequencer-csv is set"
                        .to_string(),
                );
            };
            if !Path::new(path).is_file() {
                return Err(format!("profile sequencer file not found: {path}"));
            }
            let profile = ProfileSequencer::parse_csv(path, id, true)
                .map_err(|e| format!("failed to parse profile sequencer '{path}': {e}"))?;
            Some(profile)
        }
    };

    Ok((profile_diversity, profile_sequencer))
}

fn cmd_simulate(command: &Commands) -> i32 {
    let Commands::Simulate {
        input_reference_fasta,
        input_profile_diversity_csv,
        input_profile_sequencer_csv,
        output_forward_fastq,
        output_reverse_fastq,
        parameter_length_reads_int,
        parameter_mean_insert_int,
        parameter_std_insert_int,
        parameter_profile_sequencer_id_str,
        parameter_seed_int,
    } = command
    else {
        unreachable!("cmd_simulate is only called for Commands::Simulate");
    };

    println!(
        "simulate: starting ({} reference(s))",
        input_reference_fasta.len()
    );

    if input_reference_fasta.is_empty() {
        log::error!("at least one --input-reference-fasta must be provided");
        return 1;
    }

    let min_scaffold_len = ((*parameter_length_reads_int as u64 * 2) / 3).max(1);
    let references = match load_references(input_reference_fasta, min_scaffold_len) {
        Ok(references) => references,
        Err(e) => {
            log::error!("{e}");
            return 1;
        }
    };
    let total_reads: u64 = references.iter().map(|r| r.nb_reads).sum();
    println!(
        "simulate: loaded {} reference(s), {total_reads} read pair(s) requested",
        references.len()
    );

    let (profile_diversity, profile_sequencer) = match load_profiles(
        input_profile_diversity_csv,
        input_profile_sequencer_csv,
        parameter_profile_sequencer_id_str,
    ) {
        Ok(profiles) => profiles,
        Err(e) => {
            log::error!("{e}");
            return 1;
        }
    };
    if let Some(path) = input_profile_diversity_csv {
        println!("simulate: loaded diversity profile from '{path}'");
    }
    if let Some(path) = input_profile_sequencer_csv {
        println!(
            "simulate: loaded sequencer profile from '{path}' (id={})",
            parameter_profile_sequencer_id_str.as_deref().unwrap_or("")
        );
    }

    let config = Config {
        length_reads: *parameter_length_reads_int,
        mean_insert_size: *parameter_mean_insert_int as f64,
        std_insert_size: *parameter_std_insert_int as f64,
        seed: *parameter_seed_int,
        profile_diversity,
        profile_sequencer,
    };

    let generator = match Generator::new(references, config) {
        Ok(generator) => generator,
        Err(e) => {
            log::error!("{e}");
            return 1;
        }
    };

    println!("simulate: generating reads...");
    match generator.run(output_forward_fastq, output_reverse_fastq) {
        Ok(()) => {
            println!(
                "simulate: done ({total_reads} read pair(s) requested, written to \
                 '{output_forward_fastq}' and '{output_reverse_fastq}')"
            );
            0
        }
        Err(e) => {
            log::error!("{e}");
            1
        }
    }
}

fn cmd_build_profile_sequencer(command: &Commands) -> i32 {
    let Commands::BuildProfileSequencer {
        parameter_id_str,
        input_forward_fastq,
        input_reverse_fastq,
        input_bam,
        output_profile_sequencer_csv,
    } = command
    else {
        unreachable!("cmd_build_profile_sequencer is only called for Commands::BuildProfileSequencer");
    };

    let has_fastq = !input_forward_fastq.is_empty() || !input_reverse_fastq.is_empty();
    let has_bam = !input_bam.is_empty();

    println!("build-profile-sequencer: starting (id={parameter_id_str})");

    let profile = match (has_fastq, has_bam) {
        (true, false) => {
            if input_forward_fastq.len() != input_reverse_fastq.len() {
                log::error!(
                    "--input-forward-fastq and --input-reverse-fastq must list the same \
                     number of files, in matching order ({} vs {})",
                    input_forward_fastq.len(),
                    input_reverse_fastq.len()
                );
                return 1;
            }
            for path in input_forward_fastq.iter().chain(input_reverse_fastq) {
                if !Path::new(path).is_file() {
                    log::error!("input file not found: {path}");
                    return 1;
                }
            }
            println!(
                "build-profile-sequencer: reading {} forward/reverse FASTQ pair(s)",
                input_forward_fastq.len()
            );
            BuiltProfile::from_fastq(parameter_id_str, input_forward_fastq, input_reverse_fastq)
        }
        (false, true) => {
            for path in input_bam {
                if !Path::new(path).is_file() {
                    log::error!("input file not found: {path}");
                    return 1;
                }
            }
            println!(
                "build-profile-sequencer: reading {} BAM file(s)",
                input_bam.len()
            );
            BuiltProfile::from_bam(parameter_id_str, input_bam)
        }
        (false, false) => {
            log::error!(
                "either --input-forward-fastq/--input-reverse-fastq or --input-bam must be \
                 provided"
            );
            return 1;
        }
        (true, true) => {
            log::error!(
                "--input-forward-fastq/--input-reverse-fastq and --input-bam are mutually \
                 exclusive"
            );
            return 1;
        }
    };

    let profile = match profile {
        Ok(profile) => profile,
        Err(e) => {
            log::error!("{e}");
            return 1;
        }
    };
    println!(
        "build-profile-sequencer: detected sequencer '{}'",
        profile.sequencer
    );

    match profile.write_csv(output_profile_sequencer_csv) {
        Ok(()) => {
            println!(
                "build-profile-sequencer: done (profile written to \
                 '{output_profile_sequencer_csv}')"
            );
            0
        }
        Err(e) => {
            log::error!("{e}");
            1
        }
    }
}

fn cmd_statistics_insert_size(command: &Commands) -> i32 {
    let Commands::StatisticsInsertSize {
        input_bam,
        output_statistics_csv,
    } = command
    else {
        unreachable!("cmd_statistics_insert_size is only called for Commands::StatisticsInsertSize");
    };

    println!(
        "statistics-insert-size: starting ({} BAM file(s))",
        input_bam.len()
    );

    for path in input_bam {
        if !Path::new(path).is_file() {
            log::error!("input file not found: {path}");
            return 1;
        }
    }

    let mut rows = Vec::with_capacity(input_bam.len());
    for path in input_bam {
        let stats = match InsertSizeStats::from_bam(&[path]) {
            Ok(stats) => stats,
            Err(e) => {
                log::error!("{path}: {e}");
                return 1;
            }
        };
        let basename = Path::new(path)
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.clone());
        println!(
            "statistics-insert-size: {basename} — {} read pair(s), mean={:.2}, std={:.2}",
            stats.n, stats.mean, stats.std
        );
        rows.push((basename, stats));
    }

    if let Some(output_statistics_csv) = output_statistics_csv {
        if let Err(e) = insert_size::write_csv(&rows, output_statistics_csv) {
            log::error!("{e}");
            return 1;
        }
        println!("statistics-insert-size: done (stats written to '{output_statistics_csv}')");
    }
    0
}

fn cmd_fusion_in_sample(command: &Commands) -> i32 {
    let Commands::FusionInSample {
        input_reference_fasta,
        input_profile_diversity_csv,
        input_profile_sequencer_csv,
        input_forward_fastq,
        input_reverse_fastq,
        input_bam,
        parameter_breakpoint_primary_roi,
        parameter_breakpoint_secondary_roi,
        parameter_rate_float,
        parameter_reciprocal_rate_float,
        output_forward_fastq,
        output_reverse_fastq,
        parameter_length_reads_int,
        parameter_mean_insert_int,
        parameter_std_insert_int,
        parameter_profile_sequencer_id_str,
        parameter_seed_int,
    } = command
    else {
        unreachable!("cmd_fusion_in_sample is only called for Commands::FusionInSample");
    };

    println!("fusion-in-sample: starting");

    if input_reference_fasta.is_empty() {
        log::error!("at least one --input-reference-fasta must be provided");
        return 1;
    }
    if !(0.0..=0.5).contains(parameter_rate_float) {
        log::error!(
            "--parameter-rate-float must be within [0.0, 0.5], got {parameter_rate_float}"
        );
        return 1;
    }
    if !(0.0..=1.0).contains(parameter_reciprocal_rate_float) {
        log::error!(
            "--parameter-reciprocal-rate-float must be within [0.0, 1.0], got \
             {parameter_reciprocal_rate_float}"
        );
        return 1;
    }
    for path in [input_forward_fastq, input_reverse_fastq, input_bam] {
        if !Path::new(path).is_file() {
            log::error!("input file not found: {path}");
            return 1;
        }
    }

    let primary = match fusion::Breakpoint::parse(parameter_breakpoint_primary_roi) {
        Ok(bp) => bp,
        Err(e) => {
            log::error!("invalid --parameter-breakpoint-primary-roi: {e}");
            return 1;
        }
    };
    let secondary = match fusion::Breakpoint::parse(parameter_breakpoint_secondary_roi) {
        Ok(bp) => bp,
        Err(e) => {
            log::error!("invalid --parameter-breakpoint-secondary-roi: {e}");
            return 1;
        }
    };

    let min_scaffold_len = ((*parameter_length_reads_int as u64 * 2) / 3).max(1);
    let references = match load_references(input_reference_fasta, min_scaffold_len) {
        Ok(references) => references,
        Err(e) => {
            log::error!("{e}");
            return 1;
        }
    };
    println!("fusion-in-sample: loaded {} reference(s)", references.len());

    let Some(ref_primary) = references
        .iter()
        .find(|r| r.faidx.seq_len(&primary.chrom).is_ok())
    else {
        log::error!(
            "chromosome '{}' (primary breakpoint) not found in any --input-reference-fasta",
            primary.chrom
        );
        return 1;
    };
    let Some(ref_secondary) = references
        .iter()
        .find(|r| r.faidx.seq_len(&secondary.chrom).is_ok())
    else {
        log::error!(
            "chromosome '{}' (secondary breakpoint) not found in any --input-reference-fasta",
            secondary.chrom
        );
        return 1;
    };

    let (profile_diversity, profile_sequencer) = match load_profiles(
        input_profile_diversity_csv,
        input_profile_sequencer_csv,
        parameter_profile_sequencer_id_str,
    ) {
        Ok(profiles) => profiles,
        Err(e) => {
            log::error!("{e}");
            return 1;
        }
    };
    if let Some(path) = input_profile_diversity_csv {
        println!("fusion-in-sample: loaded diversity profile from '{path}'");
    }
    if let Some(path) = input_profile_sequencer_csv {
        println!(
            "fusion-in-sample: loaded sequencer profile from '{path}' (id={})",
            parameter_profile_sequencer_id_str.as_deref().unwrap_or("")
        );
    }

    // Long enough on each side that a fragment drawn from simulate's own
    // insert-size gaussian will almost always fit entirely within one side.
    let flank_len = ((*parameter_mean_insert_int as f64 + 4.0 * *parameter_std_insert_int as f64)
        .ceil() as u64)
        .max(*parameter_length_reads_int as u64);

    let (junction_a, junction_index_a) = match fusion::build_junction(
        &ref_primary.faidx,
        &primary,
        &ref_secondary.faidx,
        &secondary,
        flank_len,
    ) {
        Ok(v) => v,
        Err(e) => {
            log::error!("failed to build the primary->secondary junction sequence: {e}");
            return 1;
        }
    };
    let (junction_b, junction_index_b) = match fusion::build_junction(
        &ref_secondary.faidx,
        &secondary,
        &ref_primary.faidx,
        &primary,
        flank_len,
    ) {
        Ok(v) => v,
        Err(e) => {
            log::error!("failed to build the secondary->primary junction sequence: {e}");
            return 1;
        }
    };

    let depth = match fusion::depth_at(input_bam, &primary.chrom, primary.pos) {
        Ok(d) => d,
        Err(e) => {
            log::error!("failed to compute depth at the primary breakpoint: {e}");
            return 1;
        }
    };
    println!(
        "fusion-in-sample: depth {depth} at primary breakpoint {}:{}",
        primary.chrom, primary.pos
    );

    let n_total = (depth as f64 * parameter_rate_float).round() as u64;
    let n_b = (n_total as f64 * parameter_reciprocal_rate_float).round() as u64;
    let n_a = n_total - n_b;
    println!(
        "fusion-in-sample: producing {n_total} fusion-supporting read pair(s) ({n_a} primary->secondary \
         + {n_b} reciprocal, --parameter-reciprocal-rate-float {parameter_reciprocal_rate_float})"
    );

    let fusion_config = FusionConfig {
        length_reads: *parameter_length_reads_int,
        mean_insert_size: *parameter_mean_insert_int as f64,
        std_insert_size: *parameter_std_insert_int as f64,
        profile_diversity,
        id_diversity: ref_primary.id_diversity.clone(),
        profile_sequencer,
    };
    let mut rng = RandomGenerator::new(*parameter_seed_int);
    let generator = FusionGenerator::new(fusion_config);

    let label_a = format!(
        "{}:{}>{}:{}",
        primary.chrom, primary.pos, secondary.chrom, secondary.pos
    );
    let label_b = format!(
        "{}:{}>{}:{}",
        secondary.chrom, secondary.pos, primary.chrom, primary.pos
    );
    let (mut forward, mut reverse) =
        generator.generate(&junction_a, junction_index_a, n_a, &label_a, &mut rng, 0);
    let (forward_b, reverse_b) =
        generator.generate(&junction_b, junction_index_b, n_b, &label_b, &mut rng, n_a);
    forward.extend(forward_b);
    reverse.extend(reverse_b);

    let produced = forward.len() as u64;
    if produced < n_total {
        log::warn!(
            "{} fusion read pair(s) out of {n_total} requested were skipped",
            n_total - produced
        );
    }

    if let Err(e) = write_spiked_fastq(input_forward_fastq, output_forward_fastq, &forward) {
        log::error!("{e}");
        return 1;
    }
    if let Err(e) = write_spiked_fastq(input_reverse_fastq, output_reverse_fastq, &reverse) {
        log::error!("{e}");
        return 1;
    }

    println!(
        "fusion-in-sample: done ({produced} fusion read pair(s) added, written to \
         '{output_forward_fastq}' and '{output_reverse_fastq}')"
    );
    0
}

/// Copy every record from `input_path` into `output_path`, then append
/// `extra` — the input FASTQ is never modified.
fn write_spiked_fastq(
    input_path: &str,
    output_path: &str,
    extra: &[FastqRecord],
) -> Result<(), Box<dyn Error>> {
    let mut reader = FastqReader::open(input_path)?;
    let mut writer = FastqWriter::create(output_path)?;
    while let Some(record) = reader.next_record()? {
        writer.write_raw(&record)?;
    }
    for record in extra {
        writer.write_record(record)?;
    }
    writer.finish()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_repeated_reference_flag() {
        let cli = Cli::parse_from([
            APP_NAME,
            "simulate",
            "--input-reference-fasta",
            "a.fa,10",
            "--input-reference-fasta",
            "b.fa,20,human",
            "--output-forward-fastq",
            "r1.fastq.gz",
            "--output-reverse-fastq",
            "r2.fastq.gz",
        ]);
        let Commands::Simulate {
            input_reference_fasta,
            parameter_length_reads_int,
            ..
        } = &cli.command
        else {
            panic!("expected Commands::Simulate");
        };
        assert_eq!(input_reference_fasta, &vec!["a.fa,10", "b.fa,20,human"]);
        assert_eq!(*parameter_length_reads_int, 150);
    }

    #[test]
    fn version_subcommand_parses() {
        let cli = Cli::parse_from([APP_NAME, "version"]);
        assert!(matches!(cli.command, Commands::Version));
    }

    #[test]
    fn build_profile_sequencer_from_fastq_parses() {
        let cli = Cli::parse_from([
            APP_NAME,
            "build-profile-sequencer",
            "--parameter-id-str",
            "n1",
            "--input-forward-fastq",
            "r1.fastq.gz",
            "--input-reverse-fastq",
            "r2.fastq.gz",
            "--output-profile-sequencer-csv",
            "profile_sequencer.csv",
        ]);
        let Commands::BuildProfileSequencer {
            parameter_id_str,
            input_forward_fastq,
            input_bam,
            ..
        } = &cli.command
        else {
            panic!("expected Commands::BuildProfileSequencer");
        };
        assert_eq!(parameter_id_str, "n1");
        assert_eq!(input_forward_fastq, &vec!["r1.fastq.gz".to_string()]);
        assert!(input_bam.is_empty());
    }

    #[test]
    fn build_profile_sequencer_accepts_space_separated_fastq_lists() {
        let cli = Cli::parse_from([
            APP_NAME,
            "build-profile-sequencer",
            "--parameter-id-str",
            "n1",
            "--input-forward-fastq",
            "lane1_r1.fastq.gz",
            "lane2_r1.fastq.gz",
            "--input-reverse-fastq",
            "lane1_r2.fastq.gz",
            "lane2_r2.fastq.gz",
            "--output-profile-sequencer-csv",
            "profile_sequencer.csv",
        ]);
        let Commands::BuildProfileSequencer {
            input_forward_fastq,
            input_reverse_fastq,
            ..
        } = &cli.command
        else {
            panic!("expected Commands::BuildProfileSequencer");
        };
        assert_eq!(
            input_forward_fastq,
            &vec!["lane1_r1.fastq.gz".to_string(), "lane2_r1.fastq.gz".to_string()]
        );
        assert_eq!(
            input_reverse_fastq,
            &vec!["lane1_r2.fastq.gz".to_string(), "lane2_r2.fastq.gz".to_string()]
        );
    }

    #[test]
    fn build_profile_sequencer_from_bam_parses() {
        let cli = Cli::parse_from([
            APP_NAME,
            "build-profile-sequencer",
            "--parameter-id-str",
            "n1",
            "--input-bam",
            "reads.bam",
            "--output-profile-sequencer-csv",
            "profile_sequencer.csv",
        ]);
        let Commands::BuildProfileSequencer { input_bam, .. } = &cli.command else {
            panic!("expected Commands::BuildProfileSequencer");
        };
        assert_eq!(input_bam, &vec!["reads.bam".to_string()]);
    }

    #[test]
    fn rejects_mixing_fastq_and_bam_inputs() {
        let cli = Cli::parse_from([
            APP_NAME,
            "build-profile-sequencer",
            "--parameter-id-str",
            "n1",
            "--input-forward-fastq",
            "r1.fastq.gz",
            "--input-reverse-fastq",
            "r2.fastq.gz",
            "--input-bam",
            "reads.bam",
            "--output-profile-sequencer-csv",
            "profile_sequencer.csv",
        ]);
        assert_eq!(cmd_build_profile_sequencer(&cli.command), 1);
    }

    #[test]
    fn rejects_no_input_provided() {
        let cli = Cli::parse_from([
            APP_NAME,
            "build-profile-sequencer",
            "--parameter-id-str",
            "n1",
            "--output-profile-sequencer-csv",
            "profile_sequencer.csv",
        ]);
        assert_eq!(cmd_build_profile_sequencer(&cli.command), 1);
    }

    #[test]
    fn statistics_insert_size_parses() {
        let cli = Cli::parse_from([
            APP_NAME,
            "statistics-insert-size",
            "--input-bam",
            "a.bam",
            "b.bam",
            "--output-statistics-csv",
            "stats.csv",
        ]);
        let Commands::StatisticsInsertSize {
            input_bam,
            output_statistics_csv,
        } = &cli.command
        else {
            panic!("expected Commands::StatisticsInsertSize");
        };
        assert_eq!(input_bam, &vec!["a.bam".to_string(), "b.bam".to_string()]);
        assert_eq!(output_statistics_csv.as_deref(), Some("stats.csv"));
    }

    #[test]
    fn statistics_insert_size_output_csv_is_optional() {
        let cli = Cli::parse_from([
            APP_NAME,
            "statistics-insert-size",
            "--input-bam",
            "a.bam",
        ]);
        let Commands::StatisticsInsertSize {
            output_statistics_csv,
            ..
        } = &cli.command
        else {
            panic!("expected Commands::StatisticsInsertSize");
        };
        assert!(output_statistics_csv.is_none());
    }

    #[test]
    fn statistics_insert_size_csv_reports_basename_not_full_path() {
        let dir = tempfile::tempdir().unwrap();
        let nested = dir.path().join("some").join("nested").join("dir");
        std::fs::create_dir_all(&nested).unwrap();
        let bam_path = nested.join("reads.sam");
        std::fs::write(
            &bam_path,
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n\
             r1\t99\tchr1\t101\t60\t50M\t=\t251\t200\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
             r1\t147\tchr1\t251\t60\t50M\t=\t101\t-200\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        )
        .unwrap();
        let csv_path = dir.path().join("stats.csv");

        let cli = Cli::parse_from([
            APP_NAME,
            "statistics-insert-size",
            "--input-bam",
            bam_path.to_str().unwrap(),
            "--output-statistics-csv",
            csv_path.to_str().unwrap(),
        ]);
        assert_eq!(cmd_statistics_insert_size(&cli.command), 0);

        let contents = std::fs::read_to_string(&csv_path).unwrap();
        let row = contents.lines().nth(1).unwrap();
        assert!(row.starts_with("reads.sam,"));
    }

    #[test]
    fn statistics_insert_size_requires_input_bam() {
        let result = Cli::try_parse_from([APP_NAME, "statistics-insert-size"]);
        assert!(result.is_err());
    }

    #[test]
    fn parses_fusion_in_sample_flags() {
        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-in-sample",
            "--input-reference-fasta",
            "genome.fa",
            "--input-forward-fastq",
            "sample_r1.fastq.gz",
            "--input-reverse-fastq",
            "sample_r2.fastq.gz",
            "--input-bam",
            "sample.bam",
            "--parameter-breakpoint-primary-roi",
            "chr9:130854064",
            "--parameter-breakpoint-secondary-roi",
            "chr22:23632600",
            "--parameter-rate-float",
            "0.1",
            "--output-forward-fastq",
            "out_r1.fastq.gz",
            "--output-reverse-fastq",
            "out_r2.fastq.gz",
        ]);
        let Commands::FusionInSample {
            input_reference_fasta,
            input_forward_fastq,
            input_bam,
            parameter_breakpoint_primary_roi,
            parameter_rate_float,
            parameter_reciprocal_rate_float,
            parameter_length_reads_int,
            ..
        } = &cli.command
        else {
            panic!("expected Commands::FusionInSample");
        };
        assert_eq!(input_reference_fasta, &vec!["genome.fa".to_string()]);
        assert_eq!(input_forward_fastq, "sample_r1.fastq.gz");
        assert_eq!(input_bam, "sample.bam");
        assert_eq!(parameter_breakpoint_primary_roi, "chr9:130854064");
        assert_eq!(*parameter_rate_float, 0.1);
        assert_eq!(*parameter_reciprocal_rate_float, 0.5);
        assert_eq!(*parameter_length_reads_int, 150);
    }

    #[test]
    fn fusion_in_sample_rejects_reciprocal_rate_out_of_range() {
        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-in-sample",
            "--input-reference-fasta",
            "genome.fa",
            "--input-forward-fastq",
            "sample_r1.fastq.gz",
            "--input-reverse-fastq",
            "sample_r2.fastq.gz",
            "--input-bam",
            "sample.bam",
            "--parameter-breakpoint-primary-roi",
            "chr9:130854064",
            "--parameter-breakpoint-secondary-roi",
            "chr22:23632600",
            "--parameter-rate-float",
            "0.1",
            "--parameter-reciprocal-rate-float",
            "1.5",
            "--output-forward-fastq",
            "out_r1.fastq.gz",
            "--output-reverse-fastq",
            "out_r2.fastq.gz",
        ]);
        assert_eq!(cmd_fusion_in_sample(&cli.command), 1);
    }

    #[test]
    fn fusion_in_sample_rejects_rate_out_of_range() {
        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-in-sample",
            "--input-reference-fasta",
            "genome.fa",
            "--input-forward-fastq",
            "sample_r1.fastq.gz",
            "--input-reverse-fastq",
            "sample_r2.fastq.gz",
            "--input-bam",
            "sample.bam",
            "--parameter-breakpoint-primary-roi",
            "chr9:130854064",
            "--parameter-breakpoint-secondary-roi",
            "chr22:23632600",
            "--parameter-rate-float",
            "0.9",
            "--output-forward-fastq",
            "out_r1.fastq.gz",
            "--output-reverse-fastq",
            "out_r2.fastq.gz",
        ]);
        assert_eq!(cmd_fusion_in_sample(&cli.command), 1);
    }

    #[test]
    fn rejects_mismatched_forward_reverse_fastq_counts() {
        let cli = Cli::parse_from([
            APP_NAME,
            "build-profile-sequencer",
            "--parameter-id-str",
            "n1",
            "--input-forward-fastq",
            "lane1_r1.fastq.gz",
            "lane2_r1.fastq.gz",
            "--input-reverse-fastq",
            "lane1_r2.fastq.gz",
            "--output-profile-sequencer-csv",
            "profile_sequencer.csv",
        ]);
        assert_eq!(cmd_build_profile_sequencer(&cli.command), 1);
    }
}
