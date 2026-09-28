//! Minimal single-page PDF writer, used only to wrap an already-
//! rendered RGB raster (a `seqlogo_png::Canvas`) as `--oformat pdf`
//! output for `sc_seqLogo.py` (DIV-0016). This is deliberately NOT a
//! general vector-graphics PDF writer -- the page contains exactly one
//! full-page image XObject, nothing else. Reusing the PNG renderer's
//! own pixel buffer this way is a much smaller lift than a real
//! vector-graphics or font-embedding PDF backend would be, while still
//! producing a genuine, standards-conforming, viewer-openable PDF file
//! (not a renamed PNG or another format masquerading as one).
//!
//! The image is embedded RAW (uncompressed 8-bit DeviceRGB, no
//! `/Filter`) rather than DCT/Flate-compressed -- simpler and lower
//! risk to get byte-exact than implementing PDF stream compression,
//! at the cost of a larger file than a "real" PDF writer would
//! produce. For sequence-logo-sized images this is a few hundred KB
//! at most, not a practical problem.
//!
//! PDF object/xref-table mechanics are simple enough that this doesn't
//! need a PDF-writing crate dependency: a handful of indirect objects
//! (Catalog, Pages, Page, an Image XObject, a Content stream), an xref
//! table with the byte OFFSET of each object (tracked while writing),
//! and a trailer. Written directly against the PDF 1.4 spec's own
//! object syntax, not derived from any implementation.

use crate::seqlogo::StackOrder;
use crate::seqlogo_png::{build_frequency_canvas, build_mean_centered_canvas, Canvas};
use std::fmt::Write as _;

/// Renders the plain (non-centered) sequence logo as a PDF -- same
/// data as `seqlogo::render_frequency_logo_svg`/`seqlogo_png::
/// render_frequency_logo_png`, wrapped as a single full-page raster
/// image (see this module's own doc comment for why).
pub fn render_frequency_logo_pdf(
    bases: &[char],
    rows: &[Vec<i64>],
    stack_order: StackOrder,
    highlight: Option<(i64, i64)>,
) -> Vec<u8> {
    canvas_to_pdf(&build_frequency_canvas(bases, rows, stack_order, highlight))
}

/// Renders the mean-centered sequence logo as a PDF -- same data as
/// `seqlogo::render_mean_centered_logo_svg`/`seqlogo_png::
/// render_mean_centered_logo_png`.
pub fn render_mean_centered_logo_pdf(
    bases: &[char],
    rows: &[Vec<i64>],
    stack_order: StackOrder,
    highlight: Option<(i64, i64)>,
) -> Vec<u8> {
    canvas_to_pdf(&build_mean_centered_canvas(
        bases,
        rows,
        stack_order,
        highlight,
    ))
}

/// Wraps `canvas`'s RGB pixel buffer as a single-page PDF, with the
/// image scaled to fill the whole page (one PDF point per pixel).
pub(crate) fn canvas_to_pdf(canvas: &Canvas) -> Vec<u8> {
    let width = canvas.width();
    let height = canvas.height();
    let pixels = canvas.pixels();

    let content_stream = format!("q {width} 0 0 {height} 0 0 cm /Im0 Do Q");

    let mut objects: Vec<Vec<u8>> = Vec::new();
    // 1: Catalog
    objects.push(b"<< /Type /Catalog /Pages 2 0 R >>".to_vec());
    // 2: Pages
    objects.push(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>".to_vec());
    // 3: Page
    objects.push(
        format!("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width} {height}] /Resources << /XObject << /Im0 4 0 R >> >> /Contents 5 0 R >>")
            .into_bytes(),
    );
    // 4: Image XObject (raw, uncompressed DeviceRGB) -- object body is
    // binary (the pixel buffer), so it's built separately below rather
    // than going through the `objects` list of dictionary-only bodies.
    // 5: Content stream
    let content_obj = format!(
        "<< /Length {} >>\nstream\n{content_stream}\nendstream",
        content_stream.len()
    );

    let mut out: Vec<u8> = Vec::new();
    let mut offsets: Vec<usize> = Vec::with_capacity(5);

    out.extend_from_slice(b"%PDF-1.4\n");

    for (i, body) in objects.iter().enumerate() {
        let obj_num = i + 1;
        offsets.push(out.len());
        out.extend_from_slice(format!("{obj_num} 0 obj\n").as_bytes());
        out.extend_from_slice(body);
        out.extend_from_slice(b"\nendobj\n");
    }

    // Object 4: image, written manually to interleave binary pixel data.
    offsets.push(out.len());
    out.extend_from_slice(
        format!("4 0 obj\n<< /Type /XObject /Subtype /Image /Width {width} /Height {height} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Length {} >>\nstream\n", pixels.len())
            .as_bytes(),
    );
    out.extend_from_slice(pixels);
    out.extend_from_slice(b"\nendstream\nendobj\n");

    // Object 5: content stream.
    offsets.push(out.len());
    out.extend_from_slice(b"5 0 obj\n");
    out.extend_from_slice(content_obj.as_bytes());
    out.extend_from_slice(b"\nendobj\n");

    let xref_offset = out.len();
    let object_count = offsets.len() + 1; // +1 for the free-list head entry (object 0)
    let mut xref = String::new();
    let _ = writeln!(xref, "xref\n0 {object_count}\n0000000000 65535 f ");
    for offset in &offsets {
        let _ = writeln!(xref, "{offset:010} 00000 n ");
    }
    out.extend_from_slice(xref.as_bytes());

    out.extend_from_slice(
        format!("trailer\n<< /Size {object_count} /Root 1 0 R >>\nstartxref\n{xref_offset}\n%%EOF")
            .as_bytes(),
    );

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqlogo::StackOrder;
    use crate::seqlogo_png::build_frequency_canvas;

    #[test]
    fn produces_a_well_formed_pdf_header_and_trailer() {
        let bases = vec!['A', 'C'];
        let rows = vec![vec![3, 1], vec![0, 4]];
        let canvas = build_frequency_canvas(&bases, &rows, StackOrder::BigOnTop, None);
        let pdf = canvas_to_pdf(&canvas);

        assert!(pdf.starts_with(b"%PDF-1.4\n"));
        let tail = String::from_utf8_lossy(&pdf[pdf.len().saturating_sub(200)..]);
        assert!(tail.contains("%%EOF"));
        assert!(tail.contains("trailer"));
        assert!(tail.contains("/Root 1 0 R"));
    }

    #[test]
    fn xref_offsets_point_at_real_object_headers() {
        let bases = vec!['A'];
        let rows = vec![vec![5]];
        let canvas = build_frequency_canvas(&bases, &rows, StackOrder::Fixed, None);
        let pdf = canvas_to_pdf(&canvas);
        let text = String::from_utf8_lossy(&pdf);

        // Find the xref table and confirm each non-free-list offset
        // really does point at "<n> 0 obj" in the file.
        let xref_pos = text.find("\nxref\n").unwrap();
        let xref_section = &text[xref_pos + 1..];
        let lines: Vec<&str> = xref_section.lines().skip(3).take(5).collect(); // skip "xref" + count line + the free-list head entry
        for (i, line) in lines.iter().enumerate() {
            let offset: usize = line[..10].parse().unwrap();
            let obj_num = i + 1;
            let expected_prefix = format!("{obj_num} 0 obj");
            assert!(
                pdf[offset..].starts_with(expected_prefix.as_bytes()),
                "object {obj_num} offset {offset} doesn't point at its header"
            );
        }
    }

    #[test]
    fn embeds_the_exact_pixel_buffer_bytes() {
        let bases = vec!['A'];
        let rows = vec![vec![5]];
        let canvas = build_frequency_canvas(&bases, &rows, StackOrder::Fixed, None);
        let pixels = canvas.pixels().to_vec();
        let pdf = canvas_to_pdf(&canvas);
        // The raw pixel bytes should appear verbatim in the image
        // object's (uncompressed) stream.
        assert!(pdf.windows(pixels.len()).any(|w| w == pixels.as_slice()));
    }
}
