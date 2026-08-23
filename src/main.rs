use std::path::Path;

use clap::{Parser, Subcommand};

use hmnrandomread::{BuiltProfile, Config, Generator, ProfileDiversity, ProfileSequencer, Reference};

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
        Commands::Version => println!("{VERSION}"),
    }
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
    let mut references = Vec::with_capacity(input_reference_fasta.len());
    for spec in input_reference_fasta {
        let (path, nb_reads, id_diversity) = match Reference::parse_spec(spec) {
            Ok(v) => v,
            Err(e) => {
                log::error!("invalid --input-reference-fasta '{spec}': {e}");
                return 1;
            }
        };
        if !Path::new(&path).is_file() {
            log::error!("reference file not found: {path}");
            return 1;
        }
        match Reference::open(path.clone(), nb_reads, id_diversity, min_scaffold_len) {
            Ok(reference) => references.push(reference),
            Err(e) => {
                log::error!("failed to open reference '{path}': {e}");
                return 1;
            }
        }
    }
    let total_reads: u64 = references.iter().map(|r| r.nb_reads).sum();
    println!(
        "simulate: loaded {} reference(s), {total_reads} read pair(s) requested",
        references.len()
    );

    let profile_diversity = match input_profile_diversity_csv {
        None => None,
        Some(path) => {
            if !Path::new(path).is_file() {
                log::error!("profile diversity file not found: {path}");
                return 1;
            }
            match ProfileDiversity::parse_csv(path) {
                Ok(profile) => {
                    println!("simulate: loaded diversity profile from '{path}'");
                    Some(profile)
                }
                Err(e) => {
                    log::error!("failed to parse profile diversity '{path}': {e}");
                    return 1;
                }
            }
        }
    };

    let profile_sequencer = match input_profile_sequencer_csv {
        None => None,
        Some(path) => {
            let Some(id) = parameter_profile_sequencer_id_str else {
                log::error!(
                    "--parameter-profile-sequencer-id-str is required when \
                     --input-profile-sequencer-csv is set"
                );
                return 1;
            };
            if !Path::new(path).is_file() {
                log::error!("profile sequencer file not found: {path}");
                return 1;
            }
            match ProfileSequencer::parse_csv(path, id, true) {
                Ok(profile) => {
                    println!("simulate: loaded sequencer profile from '{path}' (id={id})");
                    Some(profile)
                }
                Err(e) => {
                    log::error!("failed to parse profile sequencer '{path}': {e}");
                    return 1;
                }
            }
        }
    };

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
