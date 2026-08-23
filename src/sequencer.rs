//! Guess a sequencer model from a FASTQ read header, a BAM query name, or a
//! BAM header's `@RG` platform tags.
//!
//! Illumina instrument serial numbers carry a prefix tied to the instrument
//! model (e.g. `NB552023` is a NextSeq 550, `M00947` is a MiSeq). This is the
//! same heuristic long used by QC tools to report "what sequencer produced
//! this data" from the read names alone.

/// Ordered by prefix specificity: longer/more specific prefixes are checked
/// before the shorter ones they could otherwise be mistaken for.
const INSTRUMENT_PREFIXES: &[(&str, &str)] = &[
    ("HWI-EAS", "genome analyzer"),
    ("HWUSI", "genome analyzer iix"),
    ("HWI-ST", "hiseq 2000"),
    ("LH", "novaseq x"),
    ("VH", "nextseq 2000"),
    ("NB", "nextseq 550"),
    ("NS", "nextseq 500"),
    ("K0", "hiseq 4000"),
    ("D0", "hiseq 2500"),
    ("C0", "hiseq 1500"),
    ("A0", "novaseq 6000"),
];

/// The name reported when no known instrument prefix or platform tag
/// matches.
pub const UNKNOWN: &str = "unknown";

/// Guess the sequencer model from a FASTQ header (with or without its
/// leading `@`) or a BAM `QNAME`, both of which share Illumina's
/// `<instrument>:<run>:<flowcell>:...` naming convention.
pub fn from_header(header: &str) -> String {
    let header = header.strip_prefix('@').unwrap_or(header);
    let Some(instrument_id) = header.split(':').next().filter(|s| !s.is_empty()) else {
        return UNKNOWN.to_string();
    };
    from_instrument_id(instrument_id)
}

/// Guess the sequencer model from an instrument serial number/id.
pub fn from_instrument_id(instrument_id: &str) -> String {
    let upper = instrument_id.to_ascii_uppercase();
    for (prefix, name) in INSTRUMENT_PREFIXES {
        if upper.starts_with(prefix) {
            return name.to_string();
        }
    }
    if matches!(upper.as_bytes(), [b'M', d, ..] if d.is_ascii_digit()) {
        return "miseq".to_string();
    }
    UNKNOWN.to_string()
}

/// Look for a `PM:` (platform model) tag on an `@RG` line of a BAM header's
/// text, returning it lowercased. `None` if no `@RG`/`PM` tag is present, in
/// which case the caller should fall back to [`from_header`] on a read's
/// `QNAME`.
pub fn from_bam_header(header_text: &str) -> Option<String> {
    header_text.lines().find_map(|line| {
        if !line.starts_with("@RG") {
            return None;
        }
        line.split('\t')
            .find_map(|field| field.strip_prefix("PM:"))
            .map(|pm| pm.to_ascii_lowercase())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recognizes_nextseq_550_from_fastq_header() {
        assert_eq!(
            from_header("@NB552023:15:H2N7YAFX2:1:11101:5000:1000 1:N:0:1"),
            "nextseq 550"
        );
    }

    #[test]
    fn recognizes_miseq_from_bare_qname() {
        assert_eq!(from_header("M00947:34:000000000-A1B2C:1:1101:1:2"), "miseq");
    }

    #[test]
    fn recognizes_novaseq_6000() {
        assert_eq!(from_instrument_id("A00111"), "novaseq 6000");
    }

    #[test]
    fn unknown_prefix_falls_back() {
        assert_eq!(from_instrument_id("ZZZ123"), UNKNOWN);
    }

    #[test]
    fn empty_header_is_unknown() {
        assert_eq!(from_header(""), UNKNOWN);
    }

    #[test]
    fn bam_header_reads_platform_model_tag() {
        let header = "@HD\tVN:1.6\n@RG\tID:1\tPL:ILLUMINA\tPM:NextSeq550\n";
        assert_eq!(from_bam_header(header).as_deref(), Some("nextseq550"));
    }

    #[test]
    fn bam_header_without_rg_returns_none() {
        assert_eq!(from_bam_header("@HD\tVN:1.6\n"), None);
    }
}
