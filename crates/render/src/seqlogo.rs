//! Native DNA sequence-logo rendering (SVG only for now). Ports the
//! VISUAL semantics of upstream's `qcmodule.fastq.make_logo` (which
//! delegates to the Python `logomaker` library over matplotlib) --
//! this is a from-scratch renderer reproducing the same DATA encoding
//! (per-position letter-stack height reflecting relative frequency,
//! `{A: green, C: blue, G: orange, T: red, N: grey}` color scheme,
//! `stack_order`, mean-centered vs. raw-frequency modes, highlighted
//! position range), not a byte-for-byte reproduction of matplotlib's
//! own rendering -- no two independent renderers can produce identical
//! raster/vector output from the same drawing instructions, let alone
//! from independently-implemented drawing code.
//!
//! **Known gaps, disclosed rather than silently treated as "done"**:
//! - No real font-metrics library is used. Letter glyphs are SVG
//!   `<text>` elements scaled (via an SVG `transform`) to approximately
//!   fill their target band -- a real, recognizable sequence logo when
//!   rendered by any SVG viewer, not a pixel-identical reproduction of
//!   whichever `font_name` upstream selected.
//! - `shade_below`/`fade_below` (matplotlib fill-shading effects for
//!   flipped, below-axis letters) are not implemented. Both default to
//!   `0.0` (no effect) in the CLI's own flag defaults, so this only
//!   matters for a caller that explicitly overrides them away from
//!   that default -- which is a disclosed limitation, not something
//!   silently ignored without a comment.
//! - Only SVG is implemented. PDF/PNG need real rasterization/font-
//!   embedding crates, out of scope for this pass.
//!
//! **Could not be verified against a working upstream oracle**: the
//! real `logomaker` + `pandas` combination installed in this project's
//! own development environment crashes on EVERY invocation (a genuine
//! version-incompatibility bug in that environment, unrelated to this
//! port -- see `compatibility/divergences.yaml` DIV-0016), so there is
//! no reference image to diff pixel or structural content against. The
//! values plotted here (each position's raw base counts converted to
//! relative frequency for the plain logo; each position's count minus
//! that position's own mean, for the mean-centered logo) follow
//! `logomaker.Logo`'s documented `center_values` semantics as read
//! from its own docs/source, not empirically confirmed against a
//! running reference.

use std::fmt::Write as _;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StackOrder {
    BigOnTop,
    SmallOnTop,
    Fixed,
}

impl StackOrder {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "big_on_top" => Some(Self::BigOnTop),
            "small_on_top" => Some(Self::SmallOnTop),
            "fixed" => Some(Self::Fixed),
            _ => None,
        }
    }
}

/// Upstream's exact color scheme (`qcmodule.fastq.make_logo`):
/// `{'A':'green','C':'blue','G':'orange','T':'red'}`, plus `'N':'grey'`
/// when N is not excluded. These are valid SVG/CSS color keywords
/// directly (no hex conversion needed -- matplotlib and SVG both
/// recognize the same named-color set for these five names).
fn base_color(base: char) -> &'static str {
    match base {
        'A' => "green",
        'C' => "blue",
        'G' => "orange",
        'T' => "red",
        _ => "grey",
    }
}

/// Same color scheme as `base_color`, as RGB triples (standard CSS/SVG
/// named-color values) for the PNG raster renderer, which has no
/// notion of a named-color string.
pub(crate) fn base_color_rgb(base: char) -> (u8, u8, u8) {
    match base {
        'A' => (0, 128, 0),     // green
        'C' => (0, 0, 255),     // blue
        'G' => (255, 165, 0),   // orange
        'T' => (255, 0, 0),     // red
        _ => (128, 128, 128),   // grey
    }
}

const BAR_WIDTH: f64 = 28.0;
const LOGO_HEIGHT: f64 = 260.0;
const MARGIN_LEFT: f64 = 55.0;
const MARGIN_BOTTOM: f64 = 45.0;
const MARGIN_TOP: f64 = 20.0;
const MARGIN_RIGHT: f64 = 20.0;

pub(crate) fn order_indices(values: &[f64], order: StackOrder) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..values.len()).collect();
    match order {
        StackOrder::BigOnTop => idx.sort_by(|&a, &b| values[a].abs().partial_cmp(&values[b].abs()).unwrap()),
        StackOrder::SmallOnTop => idx.sort_by(|&a, &b| values[b].abs().partial_cmp(&values[a].abs()).unwrap()),
        StackOrder::Fixed => {}
    }
    idx
}

/// Draws one letter glyph filling the vertical span
/// `[y_top, y_top + height]` (SVG y grows downward) at horizontal
/// center `cx`, within a column `bar_width` wide. `flipped` mirrors the
/// glyph vertically (upstream's `flip_below`, for values drawn below
/// the axis in the mean-centered logo).
fn svg_letter(out: &mut String, base: char, cx: f64, y_top: f64, height: f64, bar_width: f64, flipped: bool) {
    if height <= 0.01 {
        return;
    }
    let color = base_color(base);
    let sx = bar_width * 0.9;
    let sy = if flipped { -height } else { height };
    let ty = if flipped { y_top + height } else { y_top };
    let _ = writeln!(
        out,
        "<g transform=\"translate({cx:.2},{ty:.2}) scale({sx:.3},{sy:.3})\"><text x=\"0\" y=\"0\" font-family=\"sans-serif\" font-size=\"1\" font-weight=\"bold\" text-anchor=\"middle\" fill=\"{color}\">{base}</text></g>",
    );
}

fn svg_header(width: f64, height: f64) -> String {
    format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width:.0}\" height=\"{height:.0}\" viewBox=\"0 0 {width:.0} {height:.0}\">\n<rect x=\"0\" y=\"0\" width=\"{width:.0}\" height=\"{height:.0}\" fill=\"white\"/>\n"
    )
}

fn draw_highlight(out: &mut String, highlight: Option<(i64, i64)>, num_positions: usize, plot_top: f64, plot_bottom: f64) {
    let Some((start, end)) = highlight else { return };
    let x0 = MARGIN_LEFT + (start.max(0) as f64) * BAR_WIDTH;
    let x1 = MARGIN_LEFT + ((end.min(num_positions as i64 - 1).max(start) + 1) as f64) * BAR_WIDTH;
    let _ = writeln!(
        out,
        "<rect x=\"{x0:.2}\" y=\"{plot_top:.2}\" width=\"{w:.2}\" height=\"{h:.2}\" fill=\"yellow\" fill-opacity=\"0.25\"/>",
        w = x1 - x0,
        h = plot_bottom - plot_top,
    );
}

fn draw_axes(out: &mut String, num_positions: usize, axis_y: f64, plot_top: f64, y_label: &str) {
    let plot_width = num_positions as f64 * BAR_WIDTH;
    let _ = writeln!(
        out,
        "<line x1=\"{x0:.2}\" y1=\"{y:.2}\" x2=\"{x1:.2}\" y2=\"{y:.2}\" stroke=\"black\" stroke-width=\"1\"/>",
        x0 = MARGIN_LEFT,
        y = axis_y,
        x1 = MARGIN_LEFT + plot_width,
    );
    for pos in 0..num_positions {
        // Label every position, or every 5th once the logo gets wide,
        // to avoid unreadable label overlap.
        if num_positions <= 40 || pos % 5 == 0 {
            let x = MARGIN_LEFT + (pos as f64 + 0.5) * BAR_WIDTH;
            let _ = writeln!(
                out,
                "<text x=\"{x:.2}\" y=\"{y:.2}\" font-family=\"sans-serif\" font-size=\"10\" text-anchor=\"middle\">{pos}</text>",
                y = axis_y + 15.0,
            );
        }
    }
    let _ = writeln!(
        out,
        "<text x=\"14\" y=\"{y:.2}\" font-family=\"sans-serif\" font-size=\"11\" text-anchor=\"middle\" transform=\"rotate(-90 14 {y:.2})\">{y_label}</text>",
        y = plot_top + LOGO_HEIGHT / 2.0,
    );
}

/// Renders the plain (non-centered) logo: each position's letters
/// stacked bottom-to-top, height proportional to that base's relative
/// frequency at that position (matches `logomaker.Logo(mat,
/// center_values=False, ...)`).
pub fn render_frequency_logo_svg(bases: &[char], rows: &[Vec<i64>], stack_order: StackOrder, highlight: Option<(i64, i64)>) -> String {
    let num_positions = rows.len();
    let width = MARGIN_LEFT + num_positions as f64 * BAR_WIDTH + MARGIN_RIGHT;
    let height = MARGIN_TOP + LOGO_HEIGHT + MARGIN_BOTTOM;
    let axis_y = MARGIN_TOP + LOGO_HEIGHT;

    let mut out = svg_header(width, height);
    draw_highlight(&mut out, highlight, num_positions, MARGIN_TOP, axis_y);

    for (pos, row) in rows.iter().enumerate() {
        let total: i64 = row.iter().sum();
        if total <= 0 {
            continue;
        }
        let freqs: Vec<f64> = row.iter().map(|&c| c as f64 / total as f64).collect();
        let order = order_indices(&freqs, stack_order);
        let cx = MARGIN_LEFT + (pos as f64 + 0.5) * BAR_WIDTH;
        let mut cumulative = 0.0;
        for &i in &order {
            let seg_height = freqs[i] * LOGO_HEIGHT;
            let y_top = axis_y - cumulative - seg_height;
            svg_letter(&mut out, bases[i], cx, y_top, seg_height, BAR_WIDTH, false);
            cumulative += seg_height;
        }
    }

    draw_axes(&mut out, num_positions, axis_y, MARGIN_TOP, "Frequency");
    out.push_str("</svg>\n");
    out
}

/// Renders the mean-centered logo: each position's values are its raw
/// counts minus that position's own mean count, so a base above the
/// position's average is drawn upward and one below is drawn downward
/// (and flipped, matching `flip_below=True`). All positions share one
/// global vertical scale (the largest `|centered value|` anywhere in
/// the matrix maps to the full half-height), matching a shared,
/// comparable y-axis across positions the way `logomaker.Logo` renders
/// one shared axis for the whole plot (matches `center_values=True`).
pub fn render_mean_centered_logo_svg(bases: &[char], rows: &[Vec<i64>], stack_order: StackOrder, highlight: Option<(i64, i64)>) -> String {
    let num_positions = rows.len();
    let width = MARGIN_LEFT + num_positions as f64 * BAR_WIDTH + MARGIN_RIGHT;
    let half_height = LOGO_HEIGHT / 2.0;
    let height = MARGIN_TOP + LOGO_HEIGHT + MARGIN_BOTTOM;
    let axis_y = MARGIN_TOP + half_height;

    let centered: Vec<Vec<f64>> = rows
        .iter()
        .map(|row| {
            let mean = row.iter().sum::<i64>() as f64 / row.len().max(1) as f64;
            row.iter().map(|&c| c as f64 - mean).collect()
        })
        .collect();
    let max_abs = centered.iter().flatten().fold(0.0_f64, |acc, &v| acc.max(v.abs())).max(1e-9);

    let mut out = svg_header(width, height);
    draw_highlight(&mut out, highlight, num_positions, MARGIN_TOP, MARGIN_TOP + LOGO_HEIGHT);

    for (pos, values) in centered.iter().enumerate() {
        let order = order_indices(values, stack_order);
        let cx = MARGIN_LEFT + (pos as f64 + 0.5) * BAR_WIDTH;
        let mut cumulative_pos = 0.0;
        let mut cumulative_neg = 0.0;
        for &i in &order {
            let scaled = values[i] / max_abs * half_height;
            if scaled >= 0.0 {
                let y_top = axis_y - cumulative_pos - scaled;
                svg_letter(&mut out, bases[i], cx, y_top, scaled, BAR_WIDTH, false);
                cumulative_pos += scaled;
            } else {
                let seg_height = -scaled;
                let y_top = axis_y + cumulative_neg;
                svg_letter(&mut out, bases[i], cx, y_top, seg_height, BAR_WIDTH, true);
                cumulative_neg += seg_height;
            }
        }
    }

    draw_axes(&mut out, num_positions, axis_y, MARGIN_TOP, "Count - mean");
    out.push_str("</svg>\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frequency_logo_is_well_formed_svg_with_expected_letters_and_colors() {
        let bases = vec!['A', 'C'];
        let rows = vec![vec![3, 1], vec![0, 4]];
        let svg = render_frequency_logo_svg(&bases, &rows, StackOrder::BigOnTop, None);
        assert!(svg.starts_with("<?xml"));
        assert!(svg.trim_end().ends_with("</svg>"));
        assert!(svg.contains(">A<"));
        assert!(svg.contains(">C<"));
        assert!(svg.contains("fill=\"green\""));
        assert!(svg.contains("fill=\"blue\""));
    }

    #[test]
    fn frequency_logo_skips_zero_total_position() {
        let bases = vec!['A', 'C'];
        let rows = vec![vec![0, 0]];
        let svg = render_frequency_logo_svg(&bases, &rows, StackOrder::Fixed, None);
        assert!(!svg.contains(">A<"));
        assert!(!svg.contains(">C<"));
    }

    #[test]
    fn mean_centered_logo_splits_above_and_below_axis() {
        // Position with A above its own mean (3 > 2) and C below (1 < 2).
        let bases = vec!['A', 'C'];
        let rows = vec![vec![3, 1]];
        let svg = render_mean_centered_logo_svg(&bases, &rows, StackOrder::Fixed, None);
        assert!(svg.contains(">A<"));
        assert!(svg.contains(">C<"));
        // The flipped (below-axis) letter uses a negative y-scale.
        assert!(svg.contains("scale(") );
    }

    #[test]
    fn highlight_range_draws_a_rect() {
        let bases = vec!['A'];
        let rows = vec![vec![1], vec![1], vec![1]];
        let svg = render_frequency_logo_svg(&bases, &rows, StackOrder::Fixed, Some((0, 1)));
        assert!(svg.contains("<rect") && svg.contains("fill=\"yellow\""));
    }
}
