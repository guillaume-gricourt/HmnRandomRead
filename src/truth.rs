//! The truth TSV shared by `fusion-simulate` and `fusion-spike`.
//!
//! Every spiked read pair is recorded here at the moment it is created, so a
//! quantification is read off what was *injected* rather than reconstructed
//! from the alignment it is supposed to evaluate — recovering the numerator
//! by keeping only the reads an aligner placed as expected conditions the
//! measurement on the very thing being measured.
//!
//! The format is a plain TSV preceded by a `#key<TAB>value` metadata block,
//! so `awk`/`pandas`/`cut` read it with no special handling while
//! `fusion-spike` can still pick the breakpoint and pool size back out of a
//! pool written by `fusion-simulate`.

use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::fusion::{Support, TruthRow};

const COLUMNS: [&str; 6] = [
    "qname",
    "junction",
    "fragment_start",
    "fragment_stop",
    "junction_offset",
    "support",
];

fn malformed(msg: impl Into<String>) -> Box<dyn Error> {
    io::Error::new(io::ErrorKind::InvalidData, msg.into()).into()
}

/// Write the metadata block, the column header, then one row per read pair.
pub fn write_truth_tsv<P: AsRef<Path>>(
    path: P,
    metadata: &[(&str, String)],
    rows: &[TruthRow],
) -> io::Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    for (key, value) in metadata {
        writeln!(out, "#{key}\t{value}")?;
    }
    writeln!(out, "{}", COLUMNS.join("\t"))?;
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            row.qname,
            row.junction,
            row.fragment_start,
            row.fragment_stop,
            row.junction_offset,
            row.support.as_str()
        )?;
    }
    out.flush()
}

/// Streaming reader over a truth TSV: the metadata block is parsed up front
/// by [`TruthReader::open`], rows are then pulled one at a time so a pool of
/// any size can be sampled in memory proportional to the sample, not to the
/// pool.
pub struct TruthReader {
    inner: Box<dyn BufRead>,
    metadata: Vec<(String, String)>,
}

impl TruthReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let mut inner: Box<dyn BufRead> = Box::new(BufReader::new(File::open(path)?));
        let mut metadata = Vec::new();

        loop {
            let mut line = String::new();
            if inner.read_line(&mut line)? == 0 {
                return Err(malformed("truth file ends before its column header"));
            }
            let line = line.trim_end_matches(['\n', '\r']);
            let Some(entry) = line.strip_prefix('#') else {
                if line.split('\t').ne(COLUMNS) {
                    return Err(malformed(format!(
                        "truth file column header is '{line}', expected '{}'",
                        COLUMNS.join("\t")
                    )));
                }
                break;
            };
            let (key, value) = entry
                .split_once('\t')
                .ok_or_else(|| malformed(format!("malformed truth metadata line '#{entry}'")))?;
            metadata.push((key.to_string(), value.to_string()));
        }

        Ok(TruthReader { inner, metadata })
    }

    /// Value of a `#key<TAB>value` metadata entry, if present.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.metadata
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// Value of a metadata entry that must be there for the file to be usable.
    pub fn require(&self, key: &str) -> Result<&str, Box<dyn Error>> {
        self.get(key)
            .ok_or_else(|| malformed(format!("truth file has no '#{key}' metadata line")))
    }

    /// The next data row, or `None` at end of file.
    pub fn next_row(&mut self) -> Result<Option<TruthRow>, Box<dyn Error>> {
        let mut line = String::new();
        if self.inner.read_line(&mut line)? == 0 {
            return Ok(None);
        }
        let line = line.trim_end_matches(['\n', '\r']);
        let fields: Vec<&str> = line.split('\t').collect();
        let [qname, junction, fragment_start, fragment_stop, junction_offset, support] = fields[..]
        else {
            return Err(malformed(format!(
                "truth row has {} field(s), expected {}: '{line}'",
                fields.len(),
                COLUMNS.len()
            )));
        };
        Ok(Some(TruthRow {
            qname: qname.to_string(),
            junction: junction.to_string(),
            fragment_start: fragment_start.parse()?,
            fragment_stop: fragment_stop.parse()?,
            junction_offset: junction_offset.parse()?,
            support: Support::parse(support)?,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(qname: &str) -> TruthRow {
        TruthRow {
            qname: qname.to_string(),
            junction: "chr9:130854064(+)>chr22:23632600(+)".to_string(),
            fragment_start: 412,
            fragment_stop: 630,
            junction_offset: 138,
            support: Support::DiscordantPair,
        }
    }

    #[test]
    fn round_trips_metadata_and_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("truth.tsv");
        let rows = vec![row("hmnrr_fusion_a_0000000000"), row("hmnrr_fusion_a_0000000001")];
        write_truth_tsv(
            &path,
            &[
                ("breakpoint_primary", "chr9:130854064".to_string()),
                ("produced_pairs", "2".to_string()),
            ],
            &rows,
        )
        .unwrap();

        let mut reader = TruthReader::open(&path).unwrap();
        assert_eq!(reader.get("breakpoint_primary"), Some("chr9:130854064"));
        assert_eq!(reader.require("produced_pairs").unwrap(), "2");
        assert!(reader.get("absent").is_none());
        assert!(reader.require("absent").is_err());

        let mut read_back = Vec::new();
        while let Some(row) = reader.next_row().unwrap() {
            read_back.push(row);
        }
        assert_eq!(read_back, rows);
    }

    #[test]
    fn writes_no_rows_for_an_empty_pool() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("truth.tsv");
        write_truth_tsv(&path, &[("produced_pairs", "0".to_string())], &[]).unwrap();
        let mut reader = TruthReader::open(&path).unwrap();
        assert!(reader.next_row().unwrap().is_none());
    }

    #[test]
    fn rejects_a_file_whose_header_is_not_the_expected_columns() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("truth.tsv");
        std::fs::write(&path, "#key\tvalue\nqname\tjunction\n").unwrap();
        assert!(TruthReader::open(&path).is_err());
    }

    #[test]
    fn rejects_a_file_with_no_header_at_all() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("truth.tsv");
        std::fs::write(&path, "#key\tvalue\n").unwrap();
        assert!(TruthReader::open(&path).is_err());
    }

    #[test]
    fn rejects_a_row_with_a_missing_field() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("truth.tsv");
        std::fs::write(&path, format!("{}\nonly\ttwo\n", COLUMNS.join("\t"))).unwrap();
        let mut reader = TruthReader::open(&path).unwrap();
        assert!(reader.next_row().is_err());
    }
}
