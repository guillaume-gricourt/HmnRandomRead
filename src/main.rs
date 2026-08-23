use std::path::Path;

use clap::{Parser, Subcommand};

use hmnrandomread::{Config, Generator, ProfileDiversity, ProfileError, Reference};

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
        input_profile_error_csv: Option<String>,

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

        /// Identifier to select within --input-profile-error-csv. Required
        /// if --input-profile-error-csv is set.
        #[arg(long)]
        parameter_profile_error_id_str: Option<String>,

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
        Commands::Version => println!("{VERSION}"),
    }
}

fn cmd_simulate(command: &Commands) -> i32 {
    let Commands::Simulate {
        input_reference_fasta,
        input_profile_diversity_csv,
        input_profile_error_csv,
        output_forward_fastq,
        output_reverse_fastq,
        parameter_length_reads_int,
        parameter_mean_insert_int,
        parameter_std_insert_int,
        parameter_profile_error_id_str,
        parameter_seed_int,
    } = command
    else {
        unreachable!("cmd_simulate is only called for Commands::Simulate");
    };

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

    let profile_diversity = match input_profile_diversity_csv {
        None => None,
        Some(path) => {
            if !Path::new(path).is_file() {
                log::error!("profile diversity file not found: {path}");
                return 1;
            }
            match ProfileDiversity::parse_csv(path) {
                Ok(profile) => Some(profile),
                Err(e) => {
                    log::error!("failed to parse profile diversity '{path}': {e}");
                    return 1;
                }
            }
        }
    };

    let profile_error = match input_profile_error_csv {
        None => None,
        Some(path) => {
            let Some(id) = parameter_profile_error_id_str else {
                log::error!(
                    "--parameter-profile-error-id-str is required when \
                     --input-profile-error-csv is set"
                );
                return 1;
            };
            if !Path::new(path).is_file() {
                log::error!("profile error file not found: {path}");
                return 1;
            }
            match ProfileError::parse_csv(path, id, true) {
                Ok(profile) => Some(profile),
                Err(e) => {
                    log::error!("failed to parse profile error '{path}': {e}");
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
        profile_error,
    };

    let generator = match Generator::new(references, config) {
        Ok(generator) => generator,
        Err(e) => {
            log::error!("{e}");
            return 1;
        }
    };

    match generator.run(output_forward_fastq, output_reverse_fastq) {
        Ok(()) => 0,
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
}
