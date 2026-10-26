//! Non-destructive chromatogram edits. Raw channels and imported call families
//! remain unchanged; each manual call keeps an explicit original sample position.
use crate::{Document, Error, Molecule, Segment, Trace, TraceCallSet};

fn error(message: &str) -> Error {
    Error::Operation(message.into())
}
pub fn validate(document: &Document) -> Result<(), Error> {
    let Some(trace) = &document.trace else {
        return Ok(());
    };
    let samples = trace.channels.first().map_or(0, Vec::len);
    let valid_calls = |sequence: &str, qualities: &[u8], peaks: &[usize]| {
        !sequence.is_empty()
            && Molecule::Dna
                .normalize(sequence)
                .is_ok_and(|s| s == sequence)
            && qualities.len() == sequence.len()
            && peaks.len() == sequence.len()
            && peaks.iter().all(|&p| p < samples)
    };
    if document.molecule != Molecule::Dna
        || document.placeholder
        || trace.channels.len() != 4
        || samples == 0
        || samples > 2_000_000
        || trace.channels.iter().any(|c| c.len() != samples)
        || trace.channel_order.len() != 4
        || !b"ACGT"
            .iter()
            .all(|b| trace.channel_order.as_bytes().contains(b))
        || !valid_calls(&document.sequence, &trace.qualities, &trace.peaks)
    {
        return Err(error("Invalid trace signal or active call arrays"));
    }
    if (!trace.manual_calls.is_empty() && trace.manual_calls.len() != document.sequence.len())
        || trace.original_calls.len() > 16
        || trace
            .original_calls
            .iter()
            .any(|calls| !valid_calls(&calls.sequence, &calls.qualities, &calls.peaks))
    {
        return Err(error("Invalid preserved trace call family"));
    }
    if trace.sample_range.as_ref().is_some_and(|r| {
        r.start >= r.end
            || r.end > samples
            || trace.peaks.iter().any(|&p| p < r.start || p >= r.end)
    }) {
        return Err(error(
            "Trace sample crop excludes active calls or exceeds the signal",
        ));
    }
    Ok(())
}
fn prepare(document: &Document) -> Result<(Document, Trace), Error> {
    validate(document)?;
    let mut result = document.clone();
    let mut trace = result
        .trace
        .take()
        .ok_or_else(|| error("This document has no chromatogram"))?;
    if trace.original_calls.is_empty() {
        trace.original_calls.push(TraceCallSet {
            label: format!("Original calls {}", trace.calls_version),
            version: trace.calls_version,
            sequence: document.sequence.to_string(),
            qualities: trace.qualities.clone(),
            qualities_known: trace.qualities_known,
            peaks: trace.peaks.clone(),
        });
    }
    if trace.manual_calls.is_empty() {
        trace.manual_calls = vec![false; document.sequence.len()];
    }
    Ok((result, trace))
}
fn finish(document: &mut Document, mut result: Document, mut trace: Trace) -> Result<(), Error> {
    trace.edited = true;
    result.trace = Some(trace);
    validate(&result)?;
    result.validate()?;
    *document = result;
    Ok(())
}
/// Replace one called symbol at its existing peak. Quality zero is recorded and
/// manual provenance is explicit; no confidence is invented for the new call.
pub fn replace_call(document: &mut Document, index: usize, base: &str) -> Result<(), Error> {
    let base = Molecule::Dna.normalize(base)?;
    if base.len() != 1 || index >= document.sequence.len() {
        return Err(error("Choose one existing call and one nucleotide symbol"));
    }
    let (mut result, mut trace) = prepare(document)?;
    result.replace(index, index + 1, &base)?;
    trace.qualities[index] = 0;
    trace.manual_calls[index] = true;
    finish(document, result, trace)
}
/// Insert a manually specified base at a real sample between its adjacent calls.
pub fn insert_call(
    document: &mut Document,
    index: usize,
    base: &str,
    sample: usize,
) -> Result<(), Error> {
    let base = Molecule::Dna.normalize(base)?;
    if base.len() != 1 || index > document.sequence.len() {
        return Err(error("Choose a call boundary and one nucleotide symbol"));
    }
    let (mut result, mut trace) = prepare(document)?;
    let range = trace.sample_range.clone().unwrap_or(Segment {
        start: 0,
        end: trace.channels[0].len(),
    });
    if sample < range.start
        || sample >= range.end
        || index > 0 && sample <= trace.peaks[index - 1]
        || index < trace.peaks.len() && sample >= trace.peaks[index]
    {
        return Err(error(
            "Insertion sample must be within the visible signal and strictly between adjacent call peaks",
        ));
    }
    result.replace(index, index, &base)?;
    trace.peaks.insert(index, sample);
    trace.qualities.insert(index, 0);
    trace.manual_calls.insert(index, true);
    finish(document, result, trace)
}
pub fn delete_calls(document: &mut Document, range: Segment) -> Result<(), Error> {
    if range.start >= range.end
        || range.end > document.sequence.len()
        || range.end - range.start == document.sequence.len()
    {
        return Err(error(
            "Choose an existing call interval and retain at least one call",
        ));
    }
    let (mut result, mut trace) = prepare(document)?;
    result.replace(range.start, range.end, "")?;
    trace.peaks.drain(range.start..range.end);
    trace.qualities.drain(range.start..range.end);
    trace.manual_calls.drain(range.start..range.end);
    finish(document, result, trace)
}
/// Keep the selected calls and crop the display at neighboring peak midpoints.
/// This is a retained range in original sample coordinates; all signal data survives.
pub fn trim_calls(document: &mut Document, keep: Segment) -> Result<(), Error> {
    if keep.start >= keep.end || keep.end > document.sequence.len() {
        return Err(error("Choose a non-empty call interval to retain"));
    }
    let (mut result, mut trace) = prepare(document)?;
    let old = trace.sample_range.clone().unwrap_or(Segment {
        start: 0,
        end: trace.channels[0].len(),
    });
    let start = if keep.start == 0 {
        old.start
    } else {
        (trace.peaks[keep.start - 1] + trace.peaks[keep.start]).div_ceil(2)
    };
    let end = if keep.end == trace.peaks.len() {
        old.end
    } else {
        (trace.peaks[keep.end - 1] + trace.peaks[keep.end])
            .div_ceil(2)
            .max(trace.peaks[keep.end - 1] + 1)
    };
    if keep.end < result.sequence.len() {
        result.replace(keep.end, result.sequence.len(), "")?;
    }
    if keep.start > 0 {
        result.replace(0, keep.start, "")?;
    }
    trace.peaks = trace.peaks[keep.start..keep.end].to_vec();
    trace.qualities = trace.qualities[keep.start..keep.end].to_vec();
    trace.manual_calls = trace.manual_calls[keep.start..keep.end].to_vec();
    // Edited chromatograms can place adjacent calls out of sample order.
    let start = start.min(*trace.peaks.iter().min().unwrap());
    let end = end.max(trace.peaks.iter().max().unwrap() + 1);
    trace.sample_range = Some(Segment { start, end });
    finish(document, result, trace)
}
/// Restore a preserved call family against the original raw channels. Existing
/// annotations prevent changing coordinate systems without an explicit remapping.
pub fn choose_call_set(document: &mut Document, index: usize) -> Result<(), Error> {
    let (mut result, mut trace) = prepare(document)?;
    let calls = trace
        .original_calls
        .get(index)
        .ok_or_else(|| error("Call family is unavailable"))?
        .clone();
    if !result.features.is_empty() && result.sequence.as_ref() != calls.sequence {
        return Err(error(
            "Remove or export annotations before switching to a call family with different sequence coordinates",
        ));
    }
    if result.sequence.as_ref() != calls.sequence {
        result.invalidate_methylation();
    }
    result.sequence = calls.sequence.into();
    result.sequence_qualities = None;
    trace.calls_version = calls.version;
    trace.qualities = calls.qualities;
    trace.qualities_known = calls.qualities_known;
    trace.peaks = calls.peaks;
    trace.manual_calls = vec![false; result.sequence.len()];
    trace.sample_range = None;
    finish(document, result, trace)
}
/// Longest contiguous interval whose sliding-window mean meets a chosen Phred
/// threshold. This is a user-controlled quality crop, not a new base-calling model.
pub fn quality_trim_range(
    document: &Document,
    minimum: u8,
    window: usize,
) -> Result<Segment, Error> {
    validate(document)?;
    let trace = document
        .trace
        .as_ref()
        .ok_or_else(|| error("This document has no chromatogram"))?;
    if !trace.qualities_known {
        return Err(error("This trace has no recorded quality estimates"));
    }
    if window == 0 || window > trace.qualities.len() {
        return Err(error(
            "Quality window must contain between one and all active calls",
        ));
    }
    let mut sum: usize = trace.qualities[..window].iter().map(|&v| v as usize).sum();
    let threshold = minimum as usize * window;
    let mut best = Segment { start: 0, end: 0 };
    let mut start = None;
    for index in 0..=trace.qualities.len() - window {
        if index > 0 {
            sum = sum + trace.qualities[index + window - 1] as usize
                - trace.qualities[index - 1] as usize;
        }
        if sum >= threshold {
            let begin = *start.get_or_insert(index);
            let end = index + window;
            if end - begin > best.end - best.start {
                best = Segment { start: begin, end };
            }
        } else {
            start = None;
        }
    }
    if best.start == best.end {
        return Err(error(
            "No call interval meets this quality threshold and window",
        ));
    }
    Ok(best)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn doc() -> Document {
        let mut d = Document::new("trace", "ACGTAC").unwrap();
        d.trace = Some(Trace {
            calls_version: 2,
            qualities: vec![5, 30, 35, 35, 30, 5],
            qualities_known: true,
            peaks: vec![10, 20, 30, 40, 50, 60],
            channel_order: "ACGT".into(),
            channels: vec![vec![1; 80]; 4],
            original_calls: vec![],
            edited: false,
            sample_range: None,
            manual_calls: vec![],
        });
        d
    }
    #[test]
    fn edits_retain_original_calls_signal_positions_and_honest_quality() {
        let mut d = doc();
        let original = d.clone();
        replace_call(&mut d, 1, "N").unwrap();
        insert_call(&mut d, 2, "A", 25).unwrap();
        delete_calls(&mut d, Segment { start: 3, end: 4 }).unwrap();
        assert_eq!(d.sequence.as_ref(), "ANATAC");
        let trace = d.trace.as_ref().unwrap();
        assert_eq!(trace.peaks, vec![10, 20, 25, 40, 50, 60]);
        assert_eq!(trace.qualities, vec![5, 0, 0, 35, 30, 5]);
        assert_eq!(
            trace.manual_calls,
            vec![false, true, true, false, false, false]
        );
        assert_eq!(trace.channels, original.trace.unwrap().channels);
        assert_eq!(trace.original_calls[0].sequence, "ACGTAC");
        choose_call_set(&mut d, 0).unwrap();
        assert_eq!(d.sequence.as_ref(), "ACGTAC");
        assert_eq!(d.trace.unwrap().peaks, vec![10, 20, 30, 40, 50, 60]);
    }
    #[test]
    fn trim_retains_raw_samples_and_uses_original_sample_coordinates() {
        let mut d = doc();
        assert_eq!(
            quality_trim_range(&d, 30, 2).unwrap(),
            Segment { start: 1, end: 5 }
        );
        trim_calls(&mut d, Segment { start: 1, end: 5 }).unwrap();
        assert_eq!(d.sequence.as_ref(), "CGTA");
        let trace = d.trace.as_ref().unwrap();
        assert_eq!(trace.peaks, vec![20, 30, 40, 50]);
        assert_eq!(trace.sample_range, Some(Segment { start: 15, end: 55 }));
        assert_eq!(trace.channels[0].len(), 80);
        assert_eq!(trace.original_calls[0].sequence, "ACGTAC");
        trim_calls(&mut d, Segment { start: 1, end: 3 }).unwrap();
        assert_eq!(
            d.trace.unwrap().sample_range,
            Some(Segment { start: 25, end: 45 })
        );
    }
    #[test]
    fn invalid_mutations_are_atomic_and_no_synthetic_peak_is_invented() {
        let mut d = doc();
        let original = d.clone();
        assert!(insert_call(&mut d, 2, "A", 20).is_err());
        assert_eq!(d, original);
        assert!(delete_calls(&mut d, Segment { start: 0, end: 6 }).is_err());
        assert_eq!(d, original);
        assert!(replace_call(&mut d, 0, "ZZ").is_err());
        assert_eq!(d, original);
        d.trace.as_mut().unwrap().sample_range = Some(Segment { start: 20, end: 80 });
        assert!(validate(&d).is_err());
    }
    #[test]
    fn valid_out_of_order_peak_calls_retain_sequence_order_when_trimmed() {
        let mut d = doc();
        d.trace.as_mut().unwrap().peaks = vec![10, 35, 25, 40, 50, 60];
        validate(&d).unwrap();
        trim_calls(&mut d, Segment { start: 1, end: 3 }).unwrap();
        assert_eq!(d.sequence.as_ref(), "CG");
        let trace = d.trace.unwrap();
        assert_eq!(trace.peaks, vec![35, 25]);
        assert!(trace.sample_range.unwrap().end > 35);
    }
}
