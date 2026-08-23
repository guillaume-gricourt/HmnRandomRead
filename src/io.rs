//! Indexed FASTA reading and gzip FASTQ writing.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::faidx;

use crate::fastq::FastqRecord;

/// Indexed access to a (optionally bgzipped) reference FASTA file.
///
/// Wraps `rust_htslib::faidx::Reader`, which builds/loads a standard
/// samtools-compatible `.fai` sidecar via htslib's own `fai_load`.
pub struct FastaIndexedReader {
    path: PathBuf,
    inner: faidx::Reader,
}

impl FastaIndexedReader {
    /// Open `path`, building its `.fai` index if it doesn't already exist.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let inner = faidx::Reader::from_path(&path).map_err(htslib_err)?;
        Ok(FastaIndexedReader { path, inner })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Every contig name in the index, in file order.
    pub fn seq_names(&self) -> io::Result<Vec<String>> {
        self.inner.seq_names().map_err(htslib_err)
    }

    /// Length of a contig. Only call this with names obtained from
    /// [`FastaIndexedReader::seq_names`] — see the note on [`fetch`].
    pub fn seq_len(&self, name: &str) -> io::Result<u64> {
        // rust-htslib casts htslib's signed "not found" sentinel (-1) into an
        // unsigned value here without checking it, so a huge result means the
        // name isn't actually in the index rather than a real length.
        let len = self.inner.fetch_seq_len(name);
        if len == u64::MAX {
            Err(not_found(name))
        } else {
            Ok(len)
        }
    }

    /// Fetch the whole (uppercased) sequence of a contig.
    pub fn fetch_all(&self, name: &str) -> io::Result<String> {
        let len = self.seq_len(name)?;
        if len == 0 {
            return Ok(String::new());
        }
        self.fetch(name, 0, len - 1)
    }

    /// Fetch `[begin, end]` (0-based, inclusive) of a contig, uppercased.
    ///
    /// `name` must come from [`FastaIndexedReader::seq_names`]: rust-htslib
    /// 0.50's `fetch_seq` does not validate htslib's returned pointer/length
    /// before building a `Vec` from them, so fetching an unknown contig name
    /// is undefined behavior rather than a clean error. Calling
    /// [`seq_len`](Self::seq_len) first (as this method does) keeps that
    /// invariant enforced in one place.
    pub fn fetch(&self, name: &str, begin: u64, end: u64) -> io::Result<String> {
        self.seq_len(name)?;
        let seq = self
            .inner
            .fetch_seq_string(name, begin as usize, end as usize)
            .map_err(htslib_err)?;
        Ok(seq.to_uppercase())
    }
}

fn not_found(name: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::NotFound,
        format!("sequence '{name}' not found in reference index"),
    )
}

fn htslib_err(e: impl std::fmt::Display) -> io::Error {
    io::Error::other(e.to_string())
}

/// Gzip-compressed FASTQ writer.
pub struct FastqWriter {
    encoder: GzEncoder<BufWriter<File>>,
}

impl FastqWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(FastqWriter {
            encoder: GzEncoder::new(BufWriter::new(file), Compression::fast()),
        })
    }

    pub fn write_record(&mut self, record: &FastqRecord) -> io::Result<()> {
        write!(self.encoder, "{record}")
    }

    pub fn finish(self) -> io::Result<()> {
        self.encoder.finish()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_fasta(dir: &Path, contents: &str) -> PathBuf {
        let path = dir.join("ref.fa");
        let mut f = File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    #[test]
    fn open_builds_index_and_fetches_sequence() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGTACGTNN\n>chr2\nTTTT\n");

        let reader = FastaIndexedReader::open(&path).unwrap();
        assert!(path.with_extension("fa.fai").exists() || dir.path().join("ref.fa.fai").exists());
        assert_eq!(reader.seq_names().unwrap(), vec!["chr1", "chr2"]);
        assert_eq!(reader.seq_len("chr1").unwrap(), 10);
        assert_eq!(reader.fetch_all("chr1").unwrap(), "ACGTACGTNN");
        assert_eq!(reader.fetch("chr1", 0, 3).unwrap(), "ACGT");
    }

    #[test]
    fn seq_len_reports_not_found_for_unknown_contig() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">chr1\nACGT\n");
        let reader = FastaIndexedReader::open(&path).unwrap();
        let err = reader.seq_len("chrX").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
    }
}
