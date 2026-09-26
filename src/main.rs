use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::path::Path;

use clap::{Parser, Subcommand};

use hmnrandomread::{
    bam_statistics, fusion, haplotype, io::FastqRawRecord, truth, FastaIndexedReader, BuiltProfile, Config, FastqReader, FastqRecord,
    FastqWriter, FusionConfig, FusionGenerator, Generator, BamStatistics, ProfileDiversity,
    ProfileSequencer, RandomGenerator, Reference,
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
        /// `path[,nb_reads[,id_diversity]]`. May be repeated. `nb_reads` may
        /// be left empty (`path,,id_diversity`) and then defaults to 0.
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
    /// Report, for each of one or more real BAMs, the mean and standard
    /// deviation of the fragment insert size (for
    /// `--parameter-mean-insert-int`/`--parameter-std-insert-int`) and the
    /// maximum depth reached, with where.
    StatisticsBam {
        /// Coordinate-sorted BAM(s), space-separated; each is reported as
        /// its own row.
        #[arg(long, num_args = 1.., required = true)]
        input_bam: Vec<String>,

        /// Output CSV: `file`, `mean_insert_size`, `std_insert_size`,
        /// `max_depth`, `max_depth_position`, one row per --input-bam.
        #[arg(long)]
        output_statistics_csv: Option<String>,
    },
    /// Build a pool of synthetic read pairs supporting a gene fusion /
    /// structural-variant breakpoint, with a truth TSV describing every pair.
    ///
    /// The pool is written on its own, not merged into a sample: dilute it
    /// into real FASTQ at a chosen allelic fraction with `fusion-spike`.
    FusionSimulate {
        /// Reference path with the number of fusion read pairs to produce and
        /// an optional diversity id: `path[,nb_reads[,id_diversity]]`. May be
        /// repeated (e.g. one file per partner chromosome); `nb_reads` and
        /// `id_diversity` are taken from whichever reference holds
        /// --parameter-breakpoint-primary-roi and must be set there.
        #[arg(long, action = clap::ArgAction::Append)]
        input_reference_fasta: Vec<String>,

        /// Diversity (SNP/indel) profile CSV, applied to the produced
        /// fusion reads.
        #[arg(long)]
        input_profile_diversity_csv: Option<String>,

        /// Sequencer error profile CSV, applied to the produced fusion reads.
        #[arg(long)]
        input_profile_sequencer_csv: Option<String>,

        /// Primary breakpoint, as `chrom:pos` (1-based): `pos` is the last
        /// reference base before the break, which falls between `pos` and
        /// `pos + 1`. The primary partner is the 5' one of the
        /// primary->secondary junction.
        #[arg(long)]
        parameter_breakpoint_primary_roi: String,

        /// Secondary breakpoint, as `chrom:pos` (1-based), same convention
        /// as --parameter-breakpoint-primary-roi. The secondary partner is
        /// the 3' one of the primary->secondary junction.
        #[arg(long)]
        parameter_breakpoint_secondary_roi: String,

        /// Probability for each partner of each read pair to be taken on
        /// the `+` strand rather than the `-` one, drawn independently for
        /// the two partners of every pair. Range: 0.0-1.0. On the `+`
        /// strand, a 5' partner contributes the bases up to the break and a
        /// 3' partner the bases after it; on the `-` strand, the reverse
        /// complement of the bases after the break (5') or up to it (3').
        /// 1.0 (default) keeps both partners on `+` (a translocation such as
        /// PML-RARA); partners on opposite strands give an inversion-type
        /// junction (e.g. EML4-ALK), so a rate of `r` mixes the four
        /// geometries at `r^2`, `r(1-r)`, `(1-r)r` and `(1-r)^2`.
        #[arg(long, default_value_t = 1.0)]
        parameter_forward_strand_rate_float: f64,

        /// Fraction of the produced fusion read pairs assigned to the
        /// reciprocal junction orientation (secondary partner's 5' part
        /// followed by the primary partner's 3' part), rather than the
        /// primary->secondary orientation. Range:
        /// 0.0-1.0. 0.5 (default) splits evenly, as for a balanced
        /// reciprocal translocation; 0.0 produces only the
        /// primary->secondary junction (e.g. for an unbalanced fusion where
        /// only one derivative is relevant, and the choice to make when
        /// quantifying a single junction, since a caller reports the two
        /// derivatives as separate events and splitting the pool halves the
        /// support behind each); 1.0 produces only the reciprocal junction.
        #[arg(long, default_value_t = 0.5)]
        parameter_reciprocal_rate_float: f64,

        /// Forward read output FASTQ (gzip-compressed) of the fusion pool.
        #[arg(long)]
        output_forward_fastq: String,

        /// Reverse read output FASTQ (gzip-compressed) of the fusion pool.
        #[arg(long)]
        output_reverse_fastq: String,

        /// Output truth TSV: one row per produced read pair (QNAME,
        /// junction, fragment bounds, junction offset, support type) under a
        /// `#key<TAB>value` metadata block. Required by `fusion-spike`.
        #[arg(long)]
        output_truth_tsv: String,

        /// Read length.
        #[arg(long, default_value_t = 150)]
        parameter_length_reads_int: usize,

        /// Mean fragment insert size.
        #[arg(long, default_value_t = 500)]
        parameter_mean_insert_int: u32,

        /// Standard deviation of the fragment insert size.
        #[arg(long, default_value_t = 50)]
        parameter_std_insert_int: u32,

        /// Minimum number of bases on each side of the junction for a read
        /// crossing it to be recorded as a `split_read` in the truth TSV;
        /// below it the pair is recorded as `short_anchor` (a real aligner
        /// typically can't call such a read as chimeric: short-read
        /// aligners need a minimum seed length — e.g. BWA-MEM's default is
        /// 19bp — and SV/fusion callers add their own minimum overhang on
        /// top). Pairs whose junction falls between the two reads are
        /// `discordant_pair`. Only classifies: every fragment crossing the
        /// junction is kept, as in a real library.
        #[arg(long, default_value_t = 20)]
        parameter_minimum_anchor_int: usize,

        /// Identifier to select within --input-profile-sequencer-csv. Required
        /// if --input-profile-sequencer-csv is set.
        #[arg(long)]
        parameter_profile_sequencer_id_str: Option<String>,

        /// Seed for the random number generator.
        #[arg(long, default_value_t = 0)]
        parameter_seed_int: u64,
    },
    /// Dilute a `fusion-simulate` pool into a real sample's paired FASTQ at a
    /// target allelic fraction, scattering the injected pairs at random
    /// positions through the output and removing the reference fragments
    /// the fusion takes the place of.
    ///
    /// `F` is the number of distinct reference fragments crossing the
    /// primary breakpoint in --input-bam — by a read or by their insert,
    /// the same class of molecule the pool holds. The fusion replaces
    /// reference molecules rather than adding to them: for every injected
    /// pair a reference fragment crossing the primary breakpoint and one
    /// crossing the secondary breakpoint are dropped from the output — for
    /// the more represented of the two derivative junctions only, since a
    /// balanced translocation's two derivatives come from the same two
    /// broken homologs. With --parameter-reciprocal-rate-float 0 (the case
    /// for quantification) that is exactly `round(f * F)` pairs injected
    /// and as many fragments removed at each breakpoint, where `f` is
    /// --parameter-allelic-fraction-float, and the molecule count at the
    /// primary breakpoint stays `F`.
    FusionSpike {
        /// The real sample's forward FASTQ (plain or gzip-compressed).
        #[arg(long)]
        input_forward_fastq: String,

        /// The real sample's reverse FASTQ (plain or gzip-compressed).
        #[arg(long)]
        input_reverse_fastq: String,

        /// Indexed (.bai/.csi sidecar), duplicate-marked BAM of the same
        /// sample, aligned from --input-forward-fastq/--input-reverse-fastq:
        /// used to count the reference molecules crossing each breakpoint
        /// (duplicates left out), to pick those to remove, by QNAME, along
        /// with every pair flagged as their duplicate, and to call the sample's own SNVs/indels
        /// around each breakpoint, which the injected reads are made to
        /// carry (phased onto one haplotype, since the fusion happened on
        /// one homolog).
        #[arg(long)]
        input_bam: String,

        /// Indexed FASTA reference(s) the pool was generated from and
        /// --input-bam aligned against, holding both breakpoints' contigs.
        /// May be repeated (e.g. one file per partner chromosome). Needed
        /// to rebuild each injected fragment with the sample's variants.
        #[arg(long, num_args = 1.., required = true, action = clap::ArgAction::Append)]
        input_reference_fasta: Vec<String>,

        /// Forward read FASTQ of the `fusion-simulate` pool.
        #[arg(long)]
        input_fusion_forward_fastq: String,

        /// Reverse read FASTQ of the `fusion-simulate` pool.
        #[arg(long)]
        input_fusion_reverse_fastq: String,

        /// Truth TSV written by `fusion-simulate` for that pool: supplies the
        /// breakpoint to measure the reference depth at, and the inventory of
        /// pairs to sample from.
        #[arg(long)]
        input_truth_tsv: String,

        /// Target allelic fraction of the fusion in the output: fusion
        /// fragments over all fragments crossing the primary breakpoint.
        /// Range: 0.0-1.0, exclusive.
        #[arg(long)]
        parameter_allelic_fraction_float: f64,

        /// Forward read output FASTQ (gzip-compressed): the sample's reads,
        /// minus the removed reference fragments, with the sampled fusion
        /// pairs scattered among them.
        #[arg(long)]
        output_forward_fastq: String,

        /// Reverse read output FASTQ (gzip-compressed): the sample's reads,
        /// minus the removed reference fragments, with the sampled fusion
        /// pairs scattered among them.
        #[arg(long)]
        output_reverse_fastq: String,

        /// Output truth TSV: the pairs actually injected, under a metadata
        /// block recording the measured depth, the fragments removed, the
        /// target fraction and the realized one.
        #[arg(long)]
        output_truth_tsv: String,

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
        Commands::StatisticsBam { .. } => {
            std::process::exit(cmd_statistics_bam(&cli.command))
        }
        Commands::FusionSimulate { .. } => std::process::exit(cmd_fusion_simulate(&cli.command)),
        Commands::FusionSpike { .. } => std::process::exit(cmd_fusion_spike(&cli.command)),
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

fn cmd_statistics_bam(command: &Commands) -> i32 {
    let Commands::StatisticsBam {
        input_bam,
        output_statistics_csv,
    } = command
    else {
        unreachable!("cmd_statistics_bam is only called for Commands::StatisticsBam");
    };

    println!(
        "statistics-bam: starting ({} BAM file(s))",
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
        let stats = match BamStatistics::from_bam(&[path]) {
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
            "statistics-bam: {basename} — {} read pair(s), mean={:.2}, std={:.2}, max depth={} \
             at {}",
            stats.n,
            stats.mean,
            stats.std,
            stats.max_depth,
            stats.max_depth_position.as_deref().unwrap_or("-")
        );
        rows.push((basename, stats));
    }

    if let Some(output_statistics_csv) = output_statistics_csv {
        if let Err(e) = bam_statistics::write_csv(&rows, output_statistics_csv) {
            log::error!("{e}");
            return 1;
        }
        println!("statistics-bam: done (stats written to '{output_statistics_csv}')");
    }
    0
}

fn cmd_fusion_simulate(command: &Commands) -> i32 {
    let Commands::FusionSimulate {
        input_reference_fasta,
        input_profile_diversity_csv,
        input_profile_sequencer_csv,
        parameter_breakpoint_primary_roi,
        parameter_breakpoint_secondary_roi,
        parameter_forward_strand_rate_float,
        parameter_reciprocal_rate_float,
        output_forward_fastq,
        output_reverse_fastq,
        output_truth_tsv,
        parameter_length_reads_int,
        parameter_mean_insert_int,
        parameter_std_insert_int,
        parameter_minimum_anchor_int,
        parameter_profile_sequencer_id_str,
        parameter_seed_int,
    } = command
    else {
        unreachable!("cmd_fusion_simulate is only called for Commands::FusionSimulate");
    };

    println!("fusion-simulate: starting");

    if input_reference_fasta.is_empty() {
        log::error!("at least one --input-reference-fasta must be provided");
        return 1;
    }
    if !(0.0..=1.0).contains(parameter_reciprocal_rate_float) {
        log::error!(
            "--parameter-reciprocal-rate-float must be within [0.0, 1.0], got \
             {parameter_reciprocal_rate_float}"
        );
        return 1;
    }
    if !(0.0..=1.0).contains(parameter_forward_strand_rate_float) {
        log::error!(
            "--parameter-forward-strand-rate-float must be within [0.0, 1.0], got \
             {parameter_forward_strand_rate_float}"
        );
        return 1;
    }
    if *parameter_length_reads_int < 2 * parameter_minimum_anchor_int {
        log::error!(
            "--parameter-length-reads-int ({parameter_length_reads_int}) must be at least twice \
             --parameter-minimum-anchor-int ({parameter_minimum_anchor_int}), otherwise no read \
             can ever fit the required anchor on both sides of the junction and no pair would \
             ever be a split_read"
        );
        return 1;
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
    println!("fusion-simulate: loaded {} reference(s)", references.len());

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

    // The pool size comes from the `nb_reads` field of the reference holding
    // the primary breakpoint, the same spec that already supplies
    // `id_diversity` — no separate count flag to keep in sync with it.
    let n_total = ref_primary.nb_reads;
    if n_total == 0 {
        log::error!(
            "no read count set for '{path}': give the number of fusion read pairs to produce in \
             its --input-reference-fasta spec, as '{path},<nb_reads>[,id_diversity]'",
            path = ref_primary.path
        );
        return 1;
    }

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
        println!("fusion-simulate: loaded diversity profile from '{path}'");
    }
    if let Some(path) = input_profile_sequencer_csv {
        println!(
            "fusion-simulate: loaded sequencer profile from '{path}' (id={})",
            parameter_profile_sequencer_id_str.as_deref().unwrap_or("")
        );
    }

    // Long enough on each side that a fragment drawn from simulate's own
    // insert-size gaussian will almost always fit entirely within one side.
    let flank_len = ((*parameter_mean_insert_int as f64 + 4.0 * *parameter_std_insert_int as f64)
        .ceil() as u64)
        .max(*parameter_length_reads_int as u64);

    let mut rng = RandomGenerator::new(*parameter_seed_int);

    let n_b = (n_total as f64 * parameter_reciprocal_rate_float).round() as u64;
    let n_a = n_total - n_b;
    println!(
        "fusion-simulate: producing {n_total} fusion read pair(s) ({n_a} primary->secondary \
         + {n_b} reciprocal, --parameter-reciprocal-rate-float {parameter_reciprocal_rate_float})"
    );

    let fusion_config = FusionConfig {
        length_reads: *parameter_length_reads_int,
        mean_insert_size: *parameter_mean_insert_int as f64,
        std_insert_size: *parameter_std_insert_int as f64,
        min_anchor: *parameter_minimum_anchor_int,
        profile_diversity,
        id_diversity: ref_primary.id_diversity.clone(),
        profile_sequencer,
    };
    let generator = FusionGenerator::new(fusion_config);

    let (mut forward, mut reverse, mut truth) = (Vec::new(), Vec::new(), Vec::new());
    let mut number = 0u64;
    // (5' partner, its reference), (3' partner, its reference), pairs.
    let orientations = [
        ((&primary, ref_primary), (&secondary, ref_secondary), n_a),
        ((&secondary, ref_secondary), (&primary, ref_primary), n_b),
    ];
    for ((five_bp, five_ref), (three_bp, three_ref), n) in orientations {
        // Each pair draws the strand of its two partners independently; the
        // pairs sharing a geometry are then generated from one junction.
        let mut per_geometry = [0u64; 4];
        for _ in 0..n {
            let five_reverse = rng.unit() >= *parameter_forward_strand_rate_float;
            let three_reverse = rng.unit() >= *parameter_forward_strand_rate_float;
            per_geometry[2 * five_reverse as usize + three_reverse as usize] += 1;
        }
        for (geometry, &count) in per_geometry.iter().enumerate() {
            if count == 0 {
                continue;
            }
            let strand = |reverse: bool| {
                if reverse {
                    fusion::Strand::Reverse
                } else {
                    fusion::Strand::Forward
                }
            };
            let five = fusion::Breakpoint { strand: strand(geometry >= 2), ..five_bp.clone() };
            let three = fusion::Breakpoint { strand: strand(geometry % 2 == 1), ..three_bp.clone() };
            let label = fusion::junction_label(&five, &three);
            let junction = fusion::build_junction(
                &fusion::JunctionPartner { faidx: &five_ref.faidx, bp: &five, variants: &[] },
                &fusion::JunctionPartner { faidx: &three_ref.faidx, bp: &three, variants: &[] },
                flank_len,
            );
            let (junction, junction_index) = match junction {
                Ok(v) => v,
                Err(e) => {
                    log::error!("failed to build the {label} junction sequence: {e}");
                    return 1;
                }
            };
            println!("fusion-simulate: {count} read pair(s) on junction {label}");
            let (f, r, t) = generator.generate(
                &junction,
                junction_index,
                count,
                &five,
                &three,
                &mut rng,
                number,
            );
            number += count;
            forward.extend(f);
            reverse.extend(r);
            truth.extend(t);
        }
    }

    let produced = forward.len() as u64;
    if produced < n_total {
        log::warn!(
            "{} fusion read pair(s) out of {n_total} requested were skipped",
            n_total - produced
        );
    }

    if let Err(e) = write_pool_fastq(output_forward_fastq, &forward) {
        log::error!("failed to write '{output_forward_fastq}': {e}");
        return 1;
    }
    if let Err(e) = write_pool_fastq(output_reverse_fastq, &reverse) {
        log::error!("failed to write '{output_reverse_fastq}': {e}");
        return 1;
    }

    let count_support =
        |support: fusion::Support| truth.iter().filter(|row| row.support == support).count();
    let metadata = vec![
        ("hmnrandomread", VERSION.to_string()),
        ("command", "fusion-simulate".to_string()),
        ("breakpoint_primary", primary.position()),
        ("breakpoint_secondary", secondary.position()),
        ("forward_strand_rate", parameter_forward_strand_rate_float.to_string()),
        ("reciprocal_rate", parameter_reciprocal_rate_float.to_string()),
        ("length_reads", parameter_length_reads_int.to_string()),
        ("mean_insert_size", parameter_mean_insert_int.to_string()),
        ("std_insert_size", parameter_std_insert_int.to_string()),
        ("minimum_anchor", parameter_minimum_anchor_int.to_string()),
        ("flank_length", flank_len.to_string()),
        ("seed", parameter_seed_int.to_string()),
        ("requested_pairs", n_total.to_string()),
        ("produced_pairs", produced.to_string()),
        ("split_read_pairs", count_support(fusion::Support::SplitRead).to_string()),
        ("short_anchor_pairs", count_support(fusion::Support::ShortAnchor).to_string()),
        ("discordant_pairs", count_support(fusion::Support::DiscordantPair).to_string()),
    ];
    if let Err(e) = truth::write_truth_tsv(output_truth_tsv, &metadata, &truth) {
        log::error!("failed to write '{output_truth_tsv}': {e}");
        return 1;
    }

    println!(
        "fusion-simulate: done ({produced} fusion read pair(s) written to \
         '{output_forward_fastq}', '{output_reverse_fastq}' and '{output_truth_tsv}')"
    );
    0
}

fn cmd_fusion_spike(command: &Commands) -> i32 {
    let Commands::FusionSpike {
        input_forward_fastq,
        input_reverse_fastq,
        input_bam,
        input_reference_fasta,
        input_fusion_forward_fastq,
        input_fusion_reverse_fastq,
        input_truth_tsv,
        parameter_allelic_fraction_float,
        output_forward_fastq,
        output_reverse_fastq,
        output_truth_tsv,
        parameter_seed_int,
    } = command
    else {
        unreachable!("cmd_fusion_spike is only called for Commands::FusionSpike");
    };

    println!("fusion-spike: starting");

    let fraction = *parameter_allelic_fraction_float;
    if !(fraction > 0.0 && fraction < 1.0) {
        log::error!(
            "--parameter-allelic-fraction-float must be within (0.0, 1.0), got {fraction}"
        );
        return 1;
    }
    for path in [
        input_forward_fastq,
        input_reverse_fastq,
        input_bam,
        input_fusion_forward_fastq,
        input_fusion_reverse_fastq,
        input_truth_tsv,
    ]
    .into_iter()
    .chain(input_reference_fasta)
    {
        if !Path::new(path).is_file() {
            log::error!("input file not found: {path}");
            return 1;
        }
    }

    let pool = match read_pool_metadata(input_truth_tsv) {
        Ok(pool) => pool,
        Err(e) => {
            log::error!("invalid --input-truth-tsv '{input_truth_tsv}': {e}");
            return 1;
        }
    };
    let (primary, secondary) = (&pool.primary, &pool.secondary);

    // Reference *fragments*, not pileup reads: one template molecule crossing
    // the breakpoint is one reference allele however many of its mates
    // overlap the position, and whether or not either of them does. `depth_at`
    // is reported alongside only so the ratio between the two stays visible.
    let mut reference = Vec::with_capacity(2);
    for bp in [primary, secondary] {
        match fusion::reference_fragments_at(input_bam, &bp.chrom, bp.pos) {
            Ok(fragments) => reference.push(fragments),
            Err(e) => {
                log::error!(
                    "failed to collect the reference fragments at breakpoint {}: {e}",
                    bp.position()
                );
                return 1;
            }
        }
    }
    let (reference_primary, reference_secondary) = (&reference[0], &reference[1]);
    let reads = match fusion::depth_at(input_bam, &primary.chrom, primary.pos) {
        Ok(d) => d,
        Err(e) => {
            log::error!("failed to compute depth at the primary breakpoint: {e}");
            return 1;
        }
    };
    let fragments = reference_primary.len() as u64;
    let fragments_secondary = reference_secondary.len() as u64;
    if fragments == 0 {
        log::error!(
            "no reference fragment crosses {} in '{input_bam}': nothing to define an allelic \
             fraction against",
            primary.position()
        );
        return 1;
    }
    println!(
        "fusion-spike: {fragments} reference fragment(s) ({reads} pileup read(s)) at primary \
         breakpoint {}, {fragments_secondary} at secondary breakpoint {}",
        primary.position(),
        secondary.position()
    );

    // n pairs injected, a share rho of them on the more represented
    // derivative, and as many reference fragments removed: the fraction
    // n / (F - rho n + n) is f for n = f F / (1 - f (1 - rho)). The pool's
    // own reciprocal rate stands in for rho, the sampled one is only known
    // once drawn and is what the realized fraction is computed from.
    let rho = pool.reciprocal_rate.max(1.0 - pool.reciprocal_rate);
    let n_target = (fraction * fragments as f64 / (1.0 - fraction * (1.0 - rho))).round() as u64;
    if n_target == 0 {
        log::error!(
            "an allelic fraction of {fraction} over {fragments} reference fragment(s) rounds to 0 \
             read pair(s): the breakpoint is not deep enough to carry that fraction"
        );
        return 1;
    }

    let mut rng = RandomGenerator::new(*parameter_seed_int);
    let (selected, pool_pairs) = match sample_truth_rows(input_truth_tsv, n_target, &mut rng) {
        Ok(v) => v,
        Err(e) => {
            log::error!("failed to sample '{input_truth_tsv}': {e}");
            return 1;
        }
    };
    if (selected.len() as u64) < n_target {
        log::error!(
            "the pool holds {pool_pairs} read pair(s), but {n_target} are needed for an allelic \
             fraction of {fraction} at this depth: generate a larger pool with fusion-simulate"
        );
        return 1;
    }

    // Whatever the strands drawn for its partners, a primary->secondary pair
    // has the primary breakpoint as its 5' partner.
    let injected_forward = selected
        .iter()
        .filter(|row| fusion::five_prime_position(&row.junction) == primary.position())
        .count() as u64;
    let injected_reciprocal = n_target - injected_forward;
    let to_remove = injected_forward.max(injected_reciprocal);
    if to_remove > fragments {
        log::error!(
            "{to_remove} reference fragment(s) would have to be removed at the primary breakpoint, \
             which only has {fragments}: lower --parameter-allelic-fraction-float"
        );
        return 1;
    }
    for (bp, reference) in [(primary, reference_primary), (secondary, reference_secondary)] {
        if reference.duplicate_pairs == 0 && reference.len() >= UNMARKED_DUPLICATES_DEPTH {
            log::warn!(
                "no read pair flagged as duplicate among the {} crossing {} in '{input_bam}': if \
                 the BAM wasn't duplicate-marked, PCR copies count as molecules and the realized \
                 fraction will be off once the output is deduplicated",
                reference.len(),
                bp.position()
            );
        }
    }
    let (removed, removed_molecules, removed_primary, removed_secondary) =
        pick_removed_fragments(reference_primary, reference_secondary, to_remove, &mut rng);
    let removed_duplicates = removed.len() as u64 - removed_molecules;
    if removed_secondary < to_remove {
        log::warn!(
            "only {removed_secondary} of the {to_remove} reference fragment(s) to remove cross the \
             secondary breakpoint {} in '{input_bam}' (is it outside the sequenced regions?)",
            secondary.position()
        );
    }
    println!(
        "fusion-spike: injecting {n_target} of the pool's {pool_pairs} fusion read pair(s) \
         ({injected_forward} primary->secondary + {injected_reciprocal} reciprocal) and removing \
         {} reference fragment(s) ({removed_primary} at the primary breakpoint, \
         {removed_secondary} at the secondary one, and {removed_duplicates} duplicate(s) of \
         them) for a target allelic fraction of {fraction}",
        removed_molecules
    );

    let wanted: HashSet<&str> = selected.iter().map(|row| row.qname.as_str()).collect();
    let (mut pool_forward, mut pool_reverse) = match collect_pool_records(
        input_fusion_forward_fastq,
        input_fusion_reverse_fastq,
        &wanted,
    ) {
        Ok(v) => v,
        Err(e) => {
            log::error!("{e}");
            return 1;
        }
    };

    // The fusion happened on one homolog of each partner: the injected
    // reads carry the sample's own alleles on it, called from its BAM.
    let mut references = Vec::with_capacity(input_reference_fasta.len());
    for path in input_reference_fasta {
        match FastaIndexedReader::open(path) {
            Ok(faidx) => references.push(faidx),
            Err(e) => {
                log::error!("failed to open '{path}': {e}");
                return 1;
            }
        }
    }
    let flank_len = pool.flank_length + HAPLOTYPE_FLANK_MARGIN;
    let mut haplotypes = Vec::with_capacity(2);
    for bp in [primary, secondary] {
        let variants = find_reference(&references, &bp.chrom).and_then(|faidx| {
            let (begin0, end0) = fusion::flank_region(faidx, bp, flank_len)?;
            let sequence = faidx.fetch(&bp.chrom, begin0, end0)?;
            haplotype::sample_haplotype(input_bam, &bp.chrom, begin0, end0, &sequence, &mut rng)
        });
        match variants {
            Ok(variants) => {
                println!(
                    "fusion-spike: {} sample variant(s) carried around breakpoint {}",
                    variants.len(),
                    bp.position()
                );
                haplotypes.push((bp.clone(), variants));
            }
            Err(e) => {
                log::error!(
                    "failed to call the sample's variants around breakpoint {}: {e}",
                    bp.position()
                );
                return 1;
            }
        }
    }
    let rewritten_pairs = match carry_sample_haplotype(
        &selected,
        &mut pool_forward,
        &mut pool_reverse,
        &references,
        &haplotypes,
        flank_len,
    ) {
        Ok(n) => n,
        Err(e) => {
            log::error!("failed to carry the sample's variants into the fusion reads: {e}");
            return 1;
        }
    };

    // Counted before anything is written: a mate-count mismatch, or removed
    // fragments the FASTQ doesn't hold, has to fail before it can leave a
    // pair of half-merged output files behind.
    let sample_pairs = match count_sample_pairs(input_forward_fastq, input_reverse_fastq, &removed)
    {
        Ok(n) => n,
        Err(e) => {
            log::error!("{e}");
            return 1;
        }
    };

    // One position list for both mates, so a pair stays a pair: the k-th
    // kept fusion record of R1 and of R2 land at the same output index.
    let output_pairs = sample_pairs - removed.len() as u64 + n_target;
    let positions = draw_positions(output_pairs, n_target, &mut rng);
    for (input, output, extra) in [
        (input_forward_fastq, output_forward_fastq, &pool_forward),
        (input_reverse_fastq, output_reverse_fastq, &pool_reverse),
    ] {
        if let Err(e) =
            write_merged_fastq(input, output, extra, &positions, sample_pairs, &removed)
        {
            log::error!("failed to write '{output}': {e}");
            return 1;
        }
    }

    let realized = n_target as f64 / (fragments - removed_primary + n_target) as f64;
    let realized_secondary = match fragments_secondary - removed_secondary + n_target {
        0 => 0.0,
        total => n_target as f64 / total as f64,
    };
    let metadata = vec![
        ("hmnrandomread", VERSION.to_string()),
        ("command", "fusion-spike".to_string()),
        ("breakpoint_primary", primary.position()),
        ("breakpoint_secondary", secondary.position()),
        ("reference_fragments", fragments.to_string()),
        ("reference_fragments_secondary", fragments_secondary.to_string()),
        ("reference_pileup_reads", reads.to_string()),
        ("target_allelic_fraction", fraction.to_string()),
        ("realized_allelic_fraction", format!("{realized:.9}")),
        ("realized_allelic_fraction_secondary", format!("{realized_secondary:.9}")),
        ("injected_pairs", n_target.to_string()),
        ("injected_primary_secondary_pairs", injected_forward.to_string()),
        ("injected_reciprocal_pairs", injected_reciprocal.to_string()),
        ("removed_reference_fragments", removed_molecules.to_string()),
        ("removed_reference_duplicates", removed_duplicates.to_string()),
        ("reference_duplicate_pairs", reference_primary.duplicate_pairs.to_string()),
        ("removed_reference_fragments_primary", removed_primary.to_string()),
        ("removed_reference_fragments_secondary", removed_secondary.to_string()),
        ("haplotype_variants_primary", variant_list(&haplotypes[0].1)),
        ("haplotype_variants_secondary", variant_list(&haplotypes[1].1)),
        ("haplotype_rewritten_pairs", rewritten_pairs.to_string()),
        ("pool_pairs", pool_pairs.to_string()),
        ("sample_pairs", sample_pairs.to_string()),
        ("output_pairs", output_pairs.to_string()),
        ("seed", parameter_seed_int.to_string()),
    ];
    if let Err(e) = truth::write_truth_tsv(output_truth_tsv, &metadata, &selected) {
        log::error!("failed to write '{output_truth_tsv}': {e}");
        return 1;
    }

    println!(
        "fusion-spike: done (realized allelic fraction {realized:.6}, written to \
         '{output_forward_fastq}', '{output_reverse_fastq}' and '{output_truth_tsv}')"
    );
    0
}

/// Comma-separated variants for a truth TSV metadata line, `.` for none.
fn variant_list(variants: &[haplotype::Variant]) -> String {
    if variants.is_empty() {
        return ".".to_string();
    }
    variants.iter().map(ToString::to_string).collect::<Vec<_>>().join(",")
}

/// Write a `fusion-simulate` pool: the produced reads and nothing else, so
/// the pool can be diluted into any sample at any fraction afterwards.
fn write_pool_fastq(output_path: &str, records: &[FastqRecord]) -> Result<(), Box<dyn Error>> {
    let mut writer = FastqWriter::create(output_path)?;
    for record in records {
        writer.write_record(record)?;
    }
    writer.finish()?;
    Ok(())
}

/// What `fusion-spike` needs to know about the pool it samples from.
struct PoolMetadata {
    primary: fusion::Breakpoint,
    secondary: fusion::Breakpoint,
    reciprocal_rate: f64,
    flank_length: u64,
}

/// The breakpoints and reciprocal rate a pool's truth TSV was generated
/// for.
fn read_pool_metadata(truth_path: &str) -> Result<PoolMetadata, Box<dyn Error>> {
    let reader = truth::TruthReader::open(truth_path)?;
    let primary = fusion::Breakpoint::parse(reader.require("breakpoint_primary")?)?;
    let secondary = fusion::Breakpoint::parse(reader.require("breakpoint_secondary")?)?;
    let reciprocal_rate: f64 = reader.require("reciprocal_rate")?.parse()?;
    let flank_length: u64 = reader.require("flank_length")?.parse()?;
    Ok(PoolMetadata { primary, secondary, reciprocal_rate, flank_length })
}

/// Extra bases beyond the pool's flanks that `fusion-spike` rebuilds each
/// junction with, so a deletion the sample carries in a flank still leaves
/// enough sequence for the farthest fragment.
const HAPLOTYPE_FLANK_MARGIN: u64 = 100;

/// The reference holding `chrom` among `references`.
fn find_reference<'a>(
    references: &'a [FastaIndexedReader],
    chrom: &str,
) -> Result<&'a FastaIndexedReader, Box<dyn Error>> {
    references
        .iter()
        .find(|faidx| faidx.seq_len(chrom).is_ok())
        .ok_or_else(|| format!("chromosome '{chrom}' not found in any --input-reference-fasta").into())
}

/// Make the sampled pool pairs carry the sample's variants, returning how
/// many pairs changed.
///
/// Each pair's junction is rebuilt twice from its truth row's label — from
/// the bare reference, as the pool was, and with the variants in
/// `haplotypes` (keyed by breakpoint) — and its fragment cut out of both at
/// the same offset from the junction: the reference version tells which
/// mate is the head read (drawn at random when the pool was built) and
/// where the pool added its own sequencing errors, the variant version
/// supplies the bases to write. Fragment length is kept, so an indel only
/// shifts the bases away from the junction. A pool built with a diversity
/// profile carrying indels loses the variants past such an indel in the
/// read, where it no longer lines up with the reference.
fn carry_sample_haplotype(
    rows: &[fusion::TruthRow],
    forward: &mut [FastqRawRecord],
    reverse: &mut [FastqRawRecord],
    references: &[FastaIndexedReader],
    haplotypes: &[(fusion::Breakpoint, Vec<haplotype::Variant>)],
    flank_len: u64,
) -> Result<u64, Box<dyn Error>> {
    type Junction = (String, usize);
    let by_name: HashMap<&str, &fusion::TruthRow> =
        rows.iter().map(|row| (row.qname.as_str(), row)).collect();
    let variants_of = |bp: &fusion::Breakpoint| {
        haplotypes
            .iter()
            .find(|(h, _)| h.position() == bp.position())
            .map_or(&[][..], |(_, v)| v.as_slice())
    };
    let mut junctions: HashMap<String, (Junction, Junction)> = HashMap::new();
    let mut rewritten = 0u64;

    for (r1, r2) in forward.iter_mut().zip(reverse.iter_mut()) {
        let row = by_name
            .get(r1.header.as_str())
            .ok_or_else(|| format!("'{}' has no truth row", r1.header))?;
        if !junctions.contains_key(&row.junction) {
            let (five, three) = fusion::parse_junction_label(&row.junction)?;
            let (five_faidx, three_faidx) =
                (find_reference(references, &five.chrom)?, find_reference(references, &three.chrom)?);
            let build = |carry: bool| {
                let pick = |bp| if carry { variants_of(bp) } else { &[][..] };
                fusion::build_junction(
                    &fusion::JunctionPartner { faidx: five_faidx, bp: &five, variants: pick(&five) },
                    &fusion::JunctionPartner {
                        faidx: three_faidx,
                        bp: &three,
                        variants: pick(&three),
                    },
                    flank_len,
                )
            };
            junctions.insert(row.junction.clone(), (build(false)?, build(true)?));
        }
        let ((clean, clean_index), (carried, carried_index)) = &junctions[&row.junction];

        let fragment_len = row.fragment_stop - row.fragment_start;
        let read_len = r1.sequence.len();
        let out_of_range = || format!("fragment of '{}' doesn't fit its junction", row.qname);
        let (clean_head, clean_tail) =
            fusion::fragment_reads(clean, *clean_index, row.junction_offset, fragment_len, read_len)
                .ok_or_else(out_of_range)?;
        let (carried_head, carried_tail) = fusion::fragment_reads(
            carried,
            *carried_index,
            row.junction_offset,
            fragment_len,
            read_len,
        )
        .ok_or_else(out_of_range)?;

        let r1_is_head = fusion::mismatches(&r1.sequence, &clean_head)
            + fusion::mismatches(&r2.sequence, &clean_tail)
            <= fusion::mismatches(&r1.sequence, &clean_tail)
                + fusion::mismatches(&r2.sequence, &clean_head);
        let ((c1, v1), (c2, v2)) = if r1_is_head {
            ((&clean_head, &carried_head), (&clean_tail, &carried_tail))
        } else {
            ((&clean_tail, &carried_tail), (&clean_head, &carried_head))
        };
        let s1 = fusion::carry_variants(&r1.sequence, c1, v1);
        let s2 = fusion::carry_variants(&r2.sequence, c2, v2);
        if s1 != r1.sequence || s2 != r2.sequence {
            rewritten += 1;
        }
        r1.sequence = s1;
        r2.sequence = s2;
    }
    Ok(rewritten)
}

/// Pick the reference fragments the fusion takes the place of: `n` crossing
/// the primary breakpoint and `n` crossing the secondary one (or as many as
/// there are), uniformly. Returns the read names to drop — each picked
/// molecule with every pair flagged as its duplicate, or a copy left in the
/// output would stand in for it once realigned and deduplicated — with the
/// number of molecules picked and how many of them cross each breakpoint.
///
/// A fragment crossing both breakpoints (two breakpoints closer than an
/// insert) counts for both, and the secondary ones are drawn among the
/// fragments not crossing the primary breakpoint, so the primary count is
/// exactly `n` and the fraction computed from it holds. Candidates are
/// sorted first so the draw depends on the seed only, not on hash order.
fn pick_removed_fragments(
    primary: &fusion::ReferenceFragments,
    secondary: &fusion::ReferenceFragments,
    n: u64,
    rng: &mut RandomGenerator,
) -> (HashSet<String>, u64, u64, u64) {
    let mut at_primary: Vec<&String> = primary.fragments.iter().collect();
    at_primary.sort();
    let picked_primary: Vec<&String> = draw_positions(at_primary.len() as u64, n, rng)
        .into_iter()
        .map(|i| at_primary[i as usize])
        .collect();

    let already = picked_primary
        .iter()
        .filter(|q| secondary.fragments.contains(**q))
        .count() as u64;
    let mut at_secondary: Vec<&String> = secondary
        .fragments
        .iter()
        .filter(|q| !primary.fragments.contains(*q))
        .collect();
    at_secondary.sort();
    let needed = n.saturating_sub(already).min(at_secondary.len() as u64);
    let picked_secondary: Vec<&String> = draw_positions(at_secondary.len() as u64, needed, rng)
        .into_iter()
        .map(|i| at_secondary[i as usize])
        .collect();

    let mut removed = HashSet::new();
    for qname in &picked_primary {
        removed.extend(primary.with_duplicates(qname).map(String::from));
        removed.extend(secondary.with_duplicates(qname).map(String::from));
    }
    for qname in &picked_secondary {
        removed.extend(secondary.with_duplicates(qname).map(String::from));
    }
    let molecules = (picked_primary.len() + picked_secondary.len()) as u64;
    (removed, molecules, n, already + needed)
}

/// Molecules at a breakpoint from which not a single duplicate-flagged pair
/// among them points to a BAM that was never duplicate-marked.
const UNMARKED_DUPLICATES_DEPTH: usize = 100;

/// The read name a FASTQ header shares with its BAM QNAME: the header up to
/// the first whitespace, without a trailing `/1` or `/2` mate suffix.
fn read_name(header: &str) -> &str {
    let name = header.split_whitespace().next().unwrap_or("");
    name.strip_suffix("/1")
        .or_else(|| name.strip_suffix("/2"))
        .unwrap_or(name)
}

/// Draw `n` truth rows uniformly without replacement from `truth_path`,
/// returning them alongside the pool's total size.
///
/// Reservoir sampling (Vitter's algorithm R) rather than a contiguous slice
/// of the pool: a pool split between the two derivative junctions holds all
/// of one orientation before any of the other, so a slice would be a biased
/// draw. Memory stays proportional to `n`, not to the pool.
fn sample_truth_rows(
    truth_path: &str,
    n: u64,
    rng: &mut RandomGenerator,
) -> Result<(Vec<fusion::TruthRow>, u64), Box<dyn Error>> {
    let mut reader = truth::TruthReader::open(truth_path)?;
    let mut reservoir: Vec<fusion::TruthRow> = Vec::new();
    let mut seen: u64 = 0;
    while let Some(row) = reader.next_row()? {
        if (reservoir.len() as u64) < n {
            reservoir.push(row);
        } else {
            let j = rng.range(0u64, seen);
            if j < n {
                reservoir[j as usize] = row;
            }
        }
        seen += 1;
    }
    reservoir.sort_by(|a, b| a.qname.cmp(&b.qname));
    Ok((reservoir, seen))
}

/// Pull the pool records named in `wanted` out of the pool's paired FASTQ,
/// in pool order and in lockstep so the k-th kept forward record is the mate
/// of the k-th kept reverse one.
fn collect_pool_records(
    forward_path: &str,
    reverse_path: &str,
    wanted: &HashSet<&str>,
) -> Result<(Vec<FastqRawRecord>, Vec<FastqRawRecord>), Box<dyn Error>> {
    let mut forward_reader = FastqReader::open(forward_path)?;
    let mut reverse_reader = FastqReader::open(reverse_path)?;
    let mut forward = Vec::with_capacity(wanted.len());
    let mut reverse = Vec::with_capacity(wanted.len());

    loop {
        let (r1, r2) = (forward_reader.next_record()?, reverse_reader.next_record()?);
        let (r1, r2) = match (r1, r2) {
            (None, None) => break,
            (Some(r1), Some(r2)) => (r1, r2),
            (Some(_), None) => {
                return Err(format!(
                    "'{forward_path}' holds more records than '{reverse_path}'"
                )
                .into())
            }
            (None, Some(_)) => {
                return Err(format!(
                    "'{reverse_path}' holds more records than '{forward_path}'"
                )
                .into())
            }
        };
        // fusion-simulate gives both mates the byte-identical QNAME, so a
        // mismatch here means the two files aren't two mates of one pool.
        if r1.header != r2.header {
            return Err(format!(
                "pool records out of step: '{}' in '{forward_path}' vs '{}' in '{reverse_path}'",
                r1.header, r2.header
            )
            .into());
        }
        if wanted.contains(r1.header.as_str()) {
            forward.push(r1);
            reverse.push(r2);
        }
    }
    if forward.len() != wanted.len() {
        return Err(format!(
            "{} of the {} sampled read pair(s) are missing from '{forward_path}': is it the pool \
             '{}' was written for?",
            wanted.len() - forward.len(),
            wanted.len(),
            forward_path
        )
        .into());
    }
    Ok((forward, reverse))
}

/// Number of records in `path`, and how many of them are named in `names`.
fn count_fastq_records(path: &str, names: &HashSet<String>) -> Result<(u64, u64), Box<dyn Error>> {
    let mut reader = FastqReader::open(path)?;
    let (mut n, mut named) = (0u64, 0u64);
    while let Some(record) = reader.next_record()? {
        n += 1;
        if names.contains(read_name(&record.header)) {
            named += 1;
        }
    }
    Ok((n, named))
}

/// Record count shared by the sample's two mate FASTQs, refusing a pair that
/// doesn't hold the same number of records — the merge splices the injected
/// pairs into both files at one shared list of positions, which only keeps
/// R1/R2 in step if the two files start out in step.
fn count_sample_pairs(
    forward_path: &str,
    reverse_path: &str,
    removed: &HashSet<String>,
) -> Result<u64, Box<dyn Error>> {
    let (forward, forward_removed) = count_fastq_records(forward_path, removed)
        .map_err(|e| format!("failed to read '{forward_path}': {e}"))?;
    let (reverse, reverse_removed) = count_fastq_records(reverse_path, removed)
        .map_err(|e| format!("failed to read '{reverse_path}': {e}"))?;
    for (path, found) in [(forward_path, forward_removed), (reverse_path, reverse_removed)] {
        if found != removed.len() as u64 {
            return Err(format!(
                "'{path}' holds {found} of the {} reference fragment(s) picked for removal from \
                 the BAM: was --input-bam aligned from these FASTQ?",
                removed.len()
            )
            .into());
        }
    }
    if forward != reverse {
        return Err(format!(
            "'{forward_path}' holds {forward} record(s) but '{reverse_path}' holds {reverse}: the \
             sample's two mate FASTQs must hold the same number of records"
        )
        .into());
    }
    Ok(forward)
}

/// Pick `n` distinct output slots out of `total_slots`, in increasing order.
///
/// Knuth's selection sampling: slot `i` is taken with probability
/// `needed / remaining`, which is exactly uniform over all `n`-subsets and
/// yields them already sorted, with no retry loop and no set to hold.
fn draw_positions(total_slots: u64, n: u64, rng: &mut RandomGenerator) -> Vec<u64> {
    let mut positions = Vec::with_capacity(n as usize);
    let mut needed = n;
    for i in 0..total_slots {
        if needed == 0 {
            break;
        }
        let remaining = total_slots - i;
        if rng.unit() * (remaining as f64) < needed as f64 {
            positions.push(i);
            needed -= 1;
        }
    }
    positions
}

/// Copy every record of `input_path` into `output_path` except those named
/// in `removed`, splicing `extra` in at `positions` (increasing output
/// indices) — the input file is never modified.
///
/// Scattering the injected pairs rather than appending them matters beyond
/// tidiness: appended reads all land in one aligner batch, and a caller or
/// duplicate-marker that sees them arrive together does not see what a real
/// library looks like.
fn write_merged_fastq(
    input_path: &str,
    output_path: &str,
    extra: &[FastqRawRecord],
    positions: &[u64],
    expected_records: u64,
    removed: &HashSet<String>,
) -> Result<(), Box<dyn Error>> {
    let mut reader = FastqReader::open(input_path)?;
    let mut writer = FastqWriter::create(output_path)?;
    let mut next_extra = 0usize;
    let mut out_index = 0u64;
    let mut copied = 0u64;
    let mut dropped = 0u64;

    loop {
        while next_extra < positions.len() && positions[next_extra] == out_index {
            writer.write_raw(&extra[next_extra])?;
            next_extra += 1;
            out_index += 1;
        }
        match reader.next_record()? {
            Some(record) if removed.contains(read_name(&record.header)) => dropped += 1,
            Some(record) => {
                writer.write_raw(&record)?;
                copied += 1;
                out_index += 1;
            }
            None => break,
        }
    }
    while next_extra < positions.len() {
        writer.write_raw(&extra[next_extra])?;
        next_extra += 1;
    }
    writer.finish()?;

    if copied + dropped != expected_records || dropped != removed.len() as u64 {
        return Err(format!(
            "'{input_path}' held {expected_records} record(s), {} of them to remove, when \
             counted but {} with {dropped} removed when merged: did it change under us?",
            removed.len(),
            copied + dropped
        )
        .into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use rust_htslib::bam::{self, Read as BamRead};

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
    fn statistics_bam_parses() {
        let cli = Cli::parse_from([
            APP_NAME,
            "statistics-bam",
            "--input-bam",
            "a.bam",
            "b.bam",
            "--output-statistics-csv",
            "stats.csv",
        ]);
        let Commands::StatisticsBam {
            input_bam,
            output_statistics_csv,
        } = &cli.command
        else {
            panic!("expected Commands::StatisticsBam");
        };
        assert_eq!(input_bam, &vec!["a.bam".to_string(), "b.bam".to_string()]);
        assert_eq!(output_statistics_csv.as_deref(), Some("stats.csv"));
    }

    #[test]
    fn statistics_bam_output_csv_is_optional() {
        let cli = Cli::parse_from([
            APP_NAME,
            "statistics-bam",
            "--input-bam",
            "a.bam",
        ]);
        let Commands::StatisticsBam {
            output_statistics_csv,
            ..
        } = &cli.command
        else {
            panic!("expected Commands::StatisticsBam");
        };
        assert!(output_statistics_csv.is_none());
    }

    #[test]
    fn statistics_bam_csv_reports_basename_not_full_path() {
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
            "statistics-bam",
            "--input-bam",
            bam_path.to_str().unwrap(),
            "--output-statistics-csv",
            csv_path.to_str().unwrap(),
        ]);
        assert_eq!(cmd_statistics_bam(&cli.command), 0);

        let contents = std::fs::read_to_string(&csv_path).unwrap();
        let row = contents.lines().nth(1).unwrap();
        assert_eq!(row, "reads.sam,200.00,0.00,1,chr1:101");
    }

    #[test]
    fn statistics_bam_requires_input_bam() {
        let result = Cli::try_parse_from([APP_NAME, "statistics-bam"]);
        assert!(result.is_err());
    }

    /// A contig of `len` bases with no N, wrapped at 60 columns so htslib
    /// can index it.
    fn write_reference(dir: &Path, name: &str, len: usize) -> String {
        let bases: String = "ACGTTGCAGGTCATCGAT".chars().cycle().take(len).collect();
        let mut fasta = format!(">{name}\n");
        for chunk in bases.as_bytes().chunks(60) {
            fasta.push_str(std::str::from_utf8(chunk).unwrap());
            fasta.push('\n');
        }
        let path = dir.join("ref.fa");
        std::fs::write(&path, fasta).unwrap();
        path.to_str().unwrap().to_string()
    }

    /// An indexed BAM holding `pairs` properly-paired fragments `p<i>`
    /// whose two mates both cover 1-based position 300 of `chr1` — so the
    /// pileup shows `2 * pairs` reads where only `pairs` molecules are
    /// present — and as many fragments `s<i>` over position 900.
    fn write_indexed_bam(dir: &Path, pairs: usize) -> String {
        write_indexed_bam_with(dir, pairs, false)
    }

    /// [`write_indexed_bam`], plus, with `duplicates`, a pair `pd<i>`
    /// flagged as a duplicate of each `p<i>` (same ends, flag 1024); the
    /// `p<i>` then start one base apart, so each molecule has its own
    /// duplicate key, as in a real library.
    fn write_indexed_bam_with(dir: &Path, pairs: usize, duplicates: bool) -> String {
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        let mut records: Vec<(usize, String)> = Vec::new();
        let mut pair = |name: String, first: usize, dup: u16| {
            let second = first + 11;
            records.push((
                first,
                format!("{name}\t{}\tchr1\t{first}\t60\t50M\t=\t{second}\t61\t{seq}\t{qual}\n", 99 + dup),
            ));
            records.push((
                second,
                format!("{name}\t{}\tchr1\t{second}\t60\t50M\t=\t{first}\t-61\t{seq}\t{qual}\n", 147 + dup),
            ));
        };
        for i in 0..pairs {
            let first = if duplicates { 255 + i } else { 260 };
            pair(format!("p{i}"), first, 0);
            if duplicates {
                pair(format!("pd{i}"), first, 1024);
            }
            pair(format!("s{i}"), 860, 0);
        }
        // Coordinate-sorted, as an index requires.
        records.sort_by_key(|(pos, _)| *pos);
        let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n");
        for (_, line) in records {
            sam.push_str(&line);
        }
        let sam_path = dir.join("sample.sam");
        std::fs::write(&sam_path, sam).unwrap();

        let bam_path = dir.join("sample.bam");
        {
            let template = bam::Reader::from_path(&sam_path).unwrap();
            let header = bam::Header::from_template(template.header());
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut reader = bam::Reader::from_path(&sam_path).unwrap();
            let mut record = bam::Record::new();
            while let Some(result) = reader.read(&mut record) {
                result.unwrap();
                writer.write(&record).unwrap();
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
        bam_path.to_str().unwrap().to_string()
    }

    /// A sample FASTQ holding the `pairs` fragments of [`write_indexed_bam`]
    /// at each breakpoint, then `filler` unrelated records, named as
    /// Illumina does (`<name> <mate>:N:0:1`).
    fn write_sample_fastq(dir: &Path, name: &str, pairs: usize, filler: usize, mate: u8) -> String {
        let names = (0..pairs)
            .map(|i| format!("p{i}"))
            .chain((0..pairs).map(|i| format!("s{i}")))
            .chain((0..filler).map(|i| format!("sample_{i}")));
        write_named_fastq(dir, name, names, mate)
    }

    /// A sample FASTQ holding one record per name, in order.
    fn write_named_fastq(
        dir: &Path,
        name: &str,
        names: impl Iterator<Item = String>,
        mate: u8,
    ) -> String {
        let mut fastq = String::new();
        for read in names {
            fastq.push_str(&format!(
                "@{read} {mate}:N:0:1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n+\n\
                 IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
            ));
        }
        let path = dir.join(name);
        std::fs::write(&path, fastq).unwrap();
        path.to_str().unwrap().to_string()
    }

    /// Every FASTQ header of a (gzip or plain) file, in file order.
    fn fastq_headers(path: &str) -> Vec<String> {
        let mut reader = FastqReader::open(path).unwrap();
        let mut headers = Vec::new();
        while let Some(record) = reader.next_record().unwrap() {
            headers.push(record.header);
        }
        headers
    }

    fn truth_rows(path: &str) -> (Vec<(String, String)>, Vec<hmnrandomread::TruthRow>) {
        let mut reader = truth::TruthReader::open(path).unwrap();
        let metadata = [
            "breakpoint_primary",
            "injected_pairs",
            "produced_pairs",
            "removed_reference_fragments_primary",
            "removed_reference_fragments_secondary",
            "realized_allelic_fraction",
            "haplotype_variants_primary",
            "haplotype_rewritten_pairs",
            "forward_strand_rate",
        ]
        .iter()
            .filter_map(|k| reader.get(k).map(|v| (k.to_string(), v.to_string())))
            .collect();
        let mut rows = Vec::new();
        while let Some(row) = reader.next_row().unwrap() {
            rows.push(row);
        }
        (metadata, rows)
    }

    /// Build a fusion pool of `pairs` pairs on `chr1:300`, returning
    /// (forward, reverse, truth) paths.
    fn run_fusion_simulate(dir: &Path, pairs: u64) -> (String, String, String) {
        run_fusion_simulate_with(dir, pairs, &[])
    }

    /// [`run_fusion_simulate`], with `extra` flags appended.
    fn run_fusion_simulate_with(dir: &Path, pairs: u64, extra: &[&str]) -> (String, String, String) {
        let reference = write_reference(dir, "chr1", 2000);
        let r1 = dir.join("pool_R1.fastq.gz").to_str().unwrap().to_string();
        let r2 = dir.join("pool_R2.fastq.gz").to_str().unwrap().to_string();
        let tsv = dir.join("pool.tsv").to_str().unwrap().to_string();
        let pool_spec = format!("{reference},{pairs}");

        let mut args = vec![
            APP_NAME,
            "fusion-simulate",
            "--input-reference-fasta",
            &pool_spec,
            "--parameter-breakpoint-primary-roi",
            "chr1:300",
            "--parameter-breakpoint-secondary-roi",
            "chr1:900",
            "--parameter-reciprocal-rate-float",
            "0.0",
            "--output-forward-fastq",
            &r1,
            "--output-reverse-fastq",
            &r2,
            "--output-truth-tsv",
            &tsv,
            "--parameter-length-reads-int",
            "50",
            "--parameter-mean-insert-int",
            "150",
            "--parameter-std-insert-int",
            "10",
            "--parameter-minimum-anchor-int",
            "10",
            "--parameter-seed-int",
            "42",
        ];
        args.extend_from_slice(extra);
        let cli = Cli::parse_from(args);
        assert_eq!(cmd_fusion_simulate(&cli.command), 0);
        (r1, r2, tsv)
    }

    #[test]
    fn parses_fusion_simulate_flags() {
        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-simulate",
            "--input-reference-fasta",
            "genome.fa,1000,human",
            "--parameter-breakpoint-primary-roi",
            "chr9:130854064",
            "--parameter-breakpoint-secondary-roi",
            "chr22:23632600",
            "--output-forward-fastq",
            "pool_r1.fastq.gz",
            "--output-reverse-fastq",
            "pool_r2.fastq.gz",
            "--output-truth-tsv",
            "pool.tsv",
        ]);
        let Commands::FusionSimulate {
            input_reference_fasta,
            parameter_breakpoint_primary_roi,
            parameter_forward_strand_rate_float,
            parameter_reciprocal_rate_float,
            output_truth_tsv,
            parameter_length_reads_int,
            parameter_minimum_anchor_int,
            ..
        } = &cli.command
        else {
            panic!("expected Commands::FusionSimulate");
        };
        assert_eq!(input_reference_fasta, &vec!["genome.fa,1000,human".to_string()]);
        assert_eq!(parameter_breakpoint_primary_roi, "chr9:130854064");
        assert_eq!(*parameter_forward_strand_rate_float, 1.0);
        assert_eq!(*parameter_reciprocal_rate_float, 0.5);
        assert_eq!(output_truth_tsv, "pool.tsv");
        assert_eq!(*parameter_length_reads_int, 150);
        assert_eq!(*parameter_minimum_anchor_int, 20);
    }

    #[test]
    fn fusion_simulate_rejects_length_reads_too_short_for_minimum_anchor() {
        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-simulate",
            "--input-reference-fasta",
            "genome.fa,10",
            "--parameter-breakpoint-primary-roi",
            "chr9:130854064",
            "--parameter-breakpoint-secondary-roi",
            "chr22:23632600",
            "--parameter-length-reads-int",
            "30",
            "--parameter-minimum-anchor-int",
            "20",
            "--output-forward-fastq",
            "out_r1.fastq.gz",
            "--output-reverse-fastq",
            "out_r2.fastq.gz",
            "--output-truth-tsv",
            "out.tsv",
        ]);
        assert_eq!(cmd_fusion_simulate(&cli.command), 1);
    }

    #[test]
    fn fusion_simulate_rejects_reciprocal_rate_out_of_range() {
        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-simulate",
            "--input-reference-fasta",
            "genome.fa,10",
            "--parameter-breakpoint-primary-roi",
            "chr9:130854064",
            "--parameter-breakpoint-secondary-roi",
            "chr22:23632600",
            "--parameter-reciprocal-rate-float",
            "1.5",
            "--output-forward-fastq",
            "out_r1.fastq.gz",
            "--output-reverse-fastq",
            "out_r2.fastq.gz",
            "--output-truth-tsv",
            "out.tsv",
        ]);
        assert_eq!(cmd_fusion_simulate(&cli.command), 1);
    }

    #[test]
    fn fusion_simulate_rejects_a_reference_with_no_read_count() {
        let dir = tempfile::tempdir().unwrap();
        let reference = write_reference(dir.path(), "chr1", 2000);
        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-simulate",
            "--input-reference-fasta",
            &format!("{reference},,human"),
            "--parameter-breakpoint-primary-roi",
            "chr1:300",
            "--parameter-breakpoint-secondary-roi",
            "chr1:900",
            "--output-forward-fastq",
            "out_r1.fastq.gz",
            "--output-reverse-fastq",
            "out_r2.fastq.gz",
            "--output-truth-tsv",
            "out.tsv",
        ]);
        assert_eq!(cmd_fusion_simulate(&cli.command), 1);
    }

    #[test]
    fn fusion_simulate_writes_a_pool_whose_qnames_match_its_truth_file() {
        let dir = tempfile::tempdir().unwrap();
        let (r1, r2, tsv) = run_fusion_simulate(dir.path(), 40);

        let (metadata, rows) = truth_rows(&tsv);
        assert!(metadata.contains(&("breakpoint_primary".to_string(), "chr1:300".to_string())));
        assert_eq!(rows.len(), 40);

        let forward = fastq_headers(&r1);
        let reverse = fastq_headers(&r2);
        // Both mates carry the same whitespace-free, self-identifying name,
        // and the truth file lists exactly those names.
        assert_eq!(forward, reverse);
        assert_eq!(
            forward,
            rows.iter().map(|r| r.qname.clone()).collect::<Vec<_>>()
        );
        for header in &forward {
            assert!(header.starts_with(fusion::QNAME_PREFIX), "{header}");
            assert!(!header.contains(char::is_whitespace), "{header}");
        }
        // Only the primary->secondary junction, as asked.
        assert!(rows.iter().all(|r| r.junction == "chr1:300(+)>chr1:900(+)"));
    }

    #[test]
    fn fusion_simulate_keeps_every_kind_of_supporting_pair() {
        let dir = tempfile::tempdir().unwrap();
        let (_, _, tsv) = run_fusion_simulate(dir.path(), 400);
        let (_, rows) = truth_rows(&tsv);
        // 150bp inserts, 50bp reads: the junction often falls between the
        // two reads, and those pairs are kept as discordant pairs.
        for support in [
            fusion::Support::SplitRead,
            fusion::Support::ShortAnchor,
            fusion::Support::DiscordantPair,
        ] {
            assert!(rows.iter().any(|r| r.support == support), "no {support:?}");
        }
    }

    #[test]
    fn fusion_simulate_mixes_strands_at_the_requested_rate() {
        let dir = tempfile::tempdir().unwrap();
        let (r1, _, tsv) = run_fusion_simulate_with(
            dir.path(),
            400,
            &["--parameter-forward-strand-rate-float", "0.5"],
        );
        let (_, rows) = truth_rows(&tsv);
        assert_eq!(rows.len(), 400);
        // Four geometries at ~25% each.
        for label in [
            "chr1:300(+)>chr1:900(+)",
            "chr1:300(+)>chr1:900(-)",
            "chr1:300(-)>chr1:900(+)",
            "chr1:300(-)>chr1:900(-)",
        ] {
            let n = rows.iter().filter(|r| r.junction == label).count();
            assert!((60..=140).contains(&n), "{label}: {n}");
        }
        // Numbered once across geometries: every QNAME is unique.
        let names: HashSet<String> = fastq_headers(&r1).into_iter().collect();
        assert_eq!(names.len(), 400);
    }

    #[test]
    fn fusion_simulate_takes_every_partner_on_the_minus_strand_at_rate_zero() {
        let dir = tempfile::tempdir().unwrap();
        let (_, _, tsv) = run_fusion_simulate_with(
            dir.path(),
            20,
            &["--parameter-forward-strand-rate-float", "0.0"],
        );
        let (metadata, rows) = truth_rows(&tsv);
        assert!(metadata.contains(&("forward_strand_rate".to_string(), "0".to_string())));
        assert!(rows.iter().all(|r| r.junction == "chr1:300(-)>chr1:900(-)"));
    }

    #[test]
    fn fusion_simulate_rejects_a_strand_rate_out_of_range() {
        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-simulate",
            "--input-reference-fasta",
            "genome.fa,10",
            "--parameter-breakpoint-primary-roi",
            "chr9:130854064",
            "--parameter-breakpoint-secondary-roi",
            "chr22:23632600",
            "--parameter-forward-strand-rate-float",
            "1.5",
            "--output-forward-fastq",
            "out_r1.fastq.gz",
            "--output-reverse-fastq",
            "out_r2.fastq.gz",
            "--output-truth-tsv",
            "out.tsv",
        ]);
        assert_eq!(cmd_fusion_simulate(&cli.command), 1);
    }

    #[test]
    fn fusion_spike_injects_the_pairs_the_target_fraction_calls_for() {
        let dir = tempfile::tempdir().unwrap();
        let (pool_r1, pool_r2, pool_tsv) = run_fusion_simulate(dir.path(), 40);
        let reference = dir.path().join("ref.fa").to_str().unwrap().to_string();
        // 20 fragments, both mates over chr1:300 — 40 pileup reads — and
        // 20 more over chr1:900.
        let bam = write_indexed_bam(dir.path(), 20);
        let sample_r1 = write_sample_fastq(dir.path(), "sample_R1.fastq", 20, 60, 1);
        let sample_r2 = write_sample_fastq(dir.path(), "sample_R2.fastq", 20, 60, 2);
        let out_r1 = dir.path().join("out_R1.fastq.gz").to_str().unwrap().to_string();
        let out_r2 = dir.path().join("out_R2.fastq.gz").to_str().unwrap().to_string();
        let out_tsv = dir.path().join("out.tsv").to_str().unwrap().to_string();

        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-spike",
            "--input-forward-fastq",
            &sample_r1,
            "--input-reverse-fastq",
            &sample_r2,
            "--input-bam",
            &bam,
            "--input-reference-fasta",
            &reference,
            "--input-fusion-forward-fastq",
            &pool_r1,
            "--input-fusion-reverse-fastq",
            &pool_r2,
            "--input-truth-tsv",
            &pool_tsv,
            "--parameter-allelic-fraction-float",
            "0.2",
            "--output-forward-fastq",
            &out_r1,
            "--output-reverse-fastq",
            &out_r2,
            "--output-truth-tsv",
            &out_tsv,
            "--parameter-seed-int",
            "1",
        ]);
        assert_eq!(cmd_fusion_spike(&cli.command), 0);

        // 0.2 * 20 fragments = 4 pairs — derived from the 20 molecules, not
        // from the 40 alignments they contribute to the pileup — replacing
        // 4 reference fragments at each breakpoint.
        let (metadata, rows) = truth_rows(&out_tsv);
        let has = |k: &str, v: &str| metadata.contains(&(k.to_string(), v.to_string()));
        assert!(has("injected_pairs", "4"));
        assert!(has("removed_reference_fragments_primary", "4"));
        assert!(has("removed_reference_fragments_secondary", "4"));
        assert!(has("realized_allelic_fraction", "0.200000000"));
        assert_eq!(rows.len(), 4);
        // The BAM's reads (ACGT repeats) disagree with the reference
        // (ACGTTGCAGGTCATCGAT repeats) around chr1:300, where only one mate
        // covers a position: the sample's homozygous variants, carried by
        // the injected reads. Where the two mates overlap they disagree
        // with each other, and each such fragment counts for nothing.
        let variants = metadata
            .iter()
            .find(|(k, _)| k == "haplotype_variants_primary")
            .map(|(_, v)| v.clone())
            .unwrap();
        assert!(variants.contains("chr1:261:G>C:hom"), "{variants}");
        assert!(!variants.contains(":het"), "{variants}");
        let rewritten = metadata
            .iter()
            .find(|(k, _)| k == "haplotype_rewritten_pairs")
            .map(|(_, v)| v.parse::<u64>().unwrap())
            .unwrap();
        assert!(rewritten > 0);

        let forward = fastq_headers(&out_r1);
        let reverse = fastq_headers(&out_r2);
        assert_eq!(forward.len(), 100 - 8 + 4);
        let names = |headers: &[String]| -> Vec<String> {
            headers.iter().map(|h| read_name(h).to_string()).collect()
        };
        assert_eq!(names(&forward), names(&reverse));

        let names = names(&forward);
        let kept = |prefix: &str| {
            (0..20).filter(|i| names.contains(&format!("{prefix}{i}"))).count()
        };
        assert_eq!(kept("p"), 16);
        assert_eq!(kept("s"), 16);
        assert_eq!((0..60).filter(|i| names.contains(&format!("sample_{i}"))).count(), 60);

        let injected: Vec<usize> = forward
            .iter()
            .enumerate()
            .filter(|(_, h)| h.starts_with(fusion::QNAME_PREFIX))
            .map(|(i, _)| i)
            .collect();
        assert_eq!(injected.len(), 4);
        // Scattered through the output, not appended to the end of it.
        assert_ne!(injected, vec![92, 93, 94, 95]);
        // ...and every injected name is one the truth file claims.
        let claimed: HashSet<&str> = rows.iter().map(|r| r.qname.as_str()).collect();
        assert!(injected.iter().all(|&i| claimed.contains(forward[i].as_str())));
    }

    #[test]
    fn fusion_spike_rejects_a_fraction_the_pool_cannot_cover() {
        let dir = tempfile::tempdir().unwrap();
        let (pool_r1, pool_r2, pool_tsv) = run_fusion_simulate(dir.path(), 2);
        let reference = dir.path().join("ref.fa").to_str().unwrap().to_string();
        let bam = write_indexed_bam(dir.path(), 20);
        let sample_r1 = write_sample_fastq(dir.path(), "sample_R1.fastq", 20, 0, 1);
        let sample_r2 = write_sample_fastq(dir.path(), "sample_R2.fastq", 20, 0, 2);

        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-spike",
            "--input-forward-fastq",
            &sample_r1,
            "--input-reverse-fastq",
            &sample_r2,
            "--input-bam",
            &bam,
            "--input-reference-fasta",
            &reference,
            "--input-fusion-forward-fastq",
            &pool_r1,
            "--input-fusion-reverse-fastq",
            &pool_r2,
            "--input-truth-tsv",
            &pool_tsv,
            "--parameter-allelic-fraction-float",
            "0.2",
            "--output-forward-fastq",
            "out_r1.fastq.gz",
            "--output-reverse-fastq",
            "out_r2.fastq.gz",
            "--output-truth-tsv",
            "out.tsv",
            "--parameter-seed-int",
            "1",
        ]);
        // 4 pairs needed, 2 in the pool.
        assert_eq!(cmd_fusion_spike(&cli.command), 1);
    }

    #[test]
    fn fusion_spike_removes_the_duplicates_of_the_removed_fragments() {
        let dir = tempfile::tempdir().unwrap();
        let (pool_r1, pool_r2, pool_tsv) = run_fusion_simulate(dir.path(), 40);
        let reference = dir.path().join("ref.fa").to_str().unwrap().to_string();
        let bam = write_indexed_bam_with(dir.path(), 20, true);
        let names = || {
            (0..20)
                .map(|i| format!("p{i}"))
                .chain((0..20).map(|i| format!("pd{i}")))
                .chain((0..20).map(|i| format!("s{i}")))
        };
        let sample_r1 = write_named_fastq(dir.path(), "sample_R1.fastq", names(), 1);
        let sample_r2 = write_named_fastq(dir.path(), "sample_R2.fastq", names(), 2);
        let out_r1 = dir.path().join("out_R1.fastq.gz").to_str().unwrap().to_string();
        let out_r2 = dir.path().join("out_R2.fastq.gz").to_str().unwrap().to_string();
        let out_tsv = dir.path().join("out.tsv").to_str().unwrap().to_string();

        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-spike",
            "--input-forward-fastq",
            &sample_r1,
            "--input-reverse-fastq",
            &sample_r2,
            "--input-bam",
            &bam,
            "--input-reference-fasta",
            &reference,
            "--input-fusion-forward-fastq",
            &pool_r1,
            "--input-fusion-reverse-fastq",
            &pool_r2,
            "--input-truth-tsv",
            &pool_tsv,
            "--parameter-allelic-fraction-float",
            "0.2",
            "--output-forward-fastq",
            &out_r1,
            "--output-reverse-fastq",
            &out_r2,
            "--output-truth-tsv",
            &out_tsv,
            "--parameter-seed-int",
            "1",
        ]);
        assert_eq!(cmd_fusion_spike(&cli.command), 0);

        // 20 molecules at the primary breakpoint, duplicates not counted:
        // 4 pairs injected, 4 molecules removed there — each with its copy.
        let names: HashSet<String> = fastq_headers(&out_r1)
            .iter()
            .map(|h| read_name(h).to_string())
            .collect();
        let removed: Vec<usize> = (0..20).filter(|i| !names.contains(&format!("p{i}"))).collect();
        assert_eq!(removed.len(), 4);
        for i in 0..20 {
            assert_eq!(
                names.contains(&format!("p{i}")),
                names.contains(&format!("pd{i}")),
                "p{i} and its duplicate pd{i} must go together"
            );
        }
        let (metadata, _) = truth_rows(&out_tsv);
        assert!(metadata.contains(&("injected_pairs".to_string(), "4".to_string())));
    }

    #[test]
    fn fusion_spike_rejects_a_fastq_missing_the_fragments_to_remove() {
        let dir = tempfile::tempdir().unwrap();
        let (pool_r1, pool_r2, pool_tsv) = run_fusion_simulate(dir.path(), 40);
        let reference = dir.path().join("ref.fa").to_str().unwrap().to_string();
        let bam = write_indexed_bam(dir.path(), 20);
        // Not the FASTQ the BAM was aligned from: none of its fragments.
        let sample_r1 = write_sample_fastq(dir.path(), "sample_R1.fastq", 0, 50, 1);
        let sample_r2 = write_sample_fastq(dir.path(), "sample_R2.fastq", 0, 50, 2);
        let out_r1 = dir.path().join("out_R1.fastq.gz").to_str().unwrap().to_string();
        let out_r2 = dir.path().join("out_R2.fastq.gz").to_str().unwrap().to_string();
        let out_tsv = dir.path().join("out.tsv").to_str().unwrap().to_string();

        let cli = Cli::parse_from([
            APP_NAME,
            "fusion-spike",
            "--input-forward-fastq",
            &sample_r1,
            "--input-reverse-fastq",
            &sample_r2,
            "--input-bam",
            &bam,
            "--input-reference-fasta",
            &reference,
            "--input-fusion-forward-fastq",
            &pool_r1,
            "--input-fusion-reverse-fastq",
            &pool_r2,
            "--input-truth-tsv",
            &pool_tsv,
            "--parameter-allelic-fraction-float",
            "0.2",
            "--output-forward-fastq",
            &out_r1,
            "--output-reverse-fastq",
            &out_r2,
            "--output-truth-tsv",
            &out_tsv,
        ]);
        assert_eq!(cmd_fusion_spike(&cli.command), 1);
        assert!(!Path::new(&out_r1).exists());
    }

    #[test]
    fn read_name_matches_the_bam_qname() {
        assert_eq!(read_name("SRR1.7 1:N:0:ACGT"), "SRR1.7");
        assert_eq!(read_name("frag_3/1"), "frag_3");
        assert_eq!(read_name("frag_3/2 extra"), "frag_3");
        assert_eq!(read_name("frag_3"), "frag_3");
    }

    #[test]
    fn fusion_spike_rejects_an_allelic_fraction_out_of_range() {
        for fraction in ["0.0", "1.0"] {
            let cli = Cli::parse_from([
                APP_NAME,
                "fusion-spike",
                "--input-forward-fastq",
                "sample_r1.fastq.gz",
                "--input-reverse-fastq",
                "sample_r2.fastq.gz",
                "--input-bam",
                "sample.bam",
                "--input-reference-fasta",
                "genome.fa",
                "--input-fusion-forward-fastq",
                "pool_r1.fastq.gz",
                "--input-fusion-reverse-fastq",
                "pool_r2.fastq.gz",
                "--input-truth-tsv",
                "pool.tsv",
                "--parameter-allelic-fraction-float",
                fraction,
                "--output-forward-fastq",
                "out_r1.fastq.gz",
                "--output-reverse-fastq",
                "out_r2.fastq.gz",
                "--output-truth-tsv",
                "out.tsv",
            ]);
            assert_eq!(cmd_fusion_spike(&cli.command), 1);
        }
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
