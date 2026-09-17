//! Native PNG sequence-logo rendering (DIV-0016). Shares the stacking/
//! ordering computation (`order_indices`) and color scheme
//! (`base_color_rgb`) with `seqlogo`'s SVG renderer, so both backends
//! agree on WHAT to draw; this module is only responsible for turning
//! that into actual pixels.
//!
//! Letters are drawn with a tiny embedded 5x7 bitmap font
//! (`crate::bitmap_font`), nearest-neighbor-scaled to fill each
//! segment's pixel box -- a real, recognizable glyph shape, not a
//! plain colored rectangle, but deliberately not attempting real
//! font-hinting/anti-aliasing (out of scope; see this crate's other
//! doc comments on what "correct" means for an image artifact here).

use crate::bitmap_font::{GLYPH_HEIGHT, GLYPH_WIDTH, glyph_rows};
use crate::seqlogo::{StackOrder, base_color_rgb, order_indices};

const BAR_WIDTH: u32 = 24;
const LOGO_HEIGHT: u32 = 200;
const MARGIN_LEFT: u32 = 30;
const MARGIN_BOTTOM: u32 = 20;
const MARGIN_TOP: u32 = 10;
const MARGIN_RIGHT: u32 = 10;
const WHITE: (u8, u8, u8) = (255, 255, 255);
const HIGHLIGHT: (u8, u8, u8) = (255, 250, 205); // pale yellow, matches the SVG renderer's highlight color at low opacity

#[derive(Clone, Copy)]
struct Rect {
    x0: u32,
    y0: u32,
    w: u32,
    h: u32,
}

struct Canvas {
    width: u32,
    height: u32,
    pixels: Vec<u8>, // RGB8, row-major
}

impl Canvas {
    fn new(width: u32, height: u32, fill: (u8, u8, u8)) -> Self {
        let mut pixels = Vec::with_capacity((width * height * 3) as usize);
        for _ in 0..(width * height) {
            pixels.push(fill.0);
            pixels.push(fill.1);
            pixels.push(fill.2);
        }
        Self { width, height, pixels }
    }

    fn set(&mut self, x: u32, y: u32, color: (u8, u8, u8)) {
        if x >= self.width || y >= self.height {
            return;
        }
        let i = ((y * self.width + x) * 3) as usize;
        self.pixels[i] = color.0;
        self.pixels[i + 1] = color.1;
        self.pixels[i + 2] = color.2;
    }

    fn fill_rect(&mut self, x0: u32, y0: u32, w: u32, h: u32, color: (u8, u8, u8)) {
        for y in y0..(y0 + h).min(self.height) {
            for x in x0..(x0 + w).min(self.width) {
                self.set(x, y, color);
            }
        }
    }

    /// Draws `c` nearest-neighbor-scaled to fill `rect`. `flipped`
    /// reverses the glyph vertically (the mean-centered logo's
    /// below-axis letters).
    fn draw_glyph(&mut self, c: char, rect: Rect, color: (u8, u8, u8), flipped: bool) {
        if rect.w == 0 || rect.h == 0 {
            return;
        }
        let rows = glyph_rows(c);
        for py in 0..rect.h {
            let gy = (py * GLYPH_HEIGHT as u32 / rect.h).min(GLYPH_HEIGHT as u32 - 1);
            let row = rows[if flipped { GLYPH_HEIGHT - 1 - gy as usize } else { gy as usize }];
            for px in 0..rect.w {
                let gx = (px * GLYPH_WIDTH as u32 / rect.w).min(GLYPH_WIDTH as u32 - 1);
                let bit = (row >> (GLYPH_WIDTH - 1 - gx as usize)) & 1;
                if bit == 1 {
                    self.set(rect.x0 + px, rect.y0 + py, color);
                }
            }
        }
    }

    fn encode_png(self) -> io::Result<Vec<u8>> {
        let mut buf = Vec::new();
        {
            let mut encoder = png::Encoder::new(&mut buf, self.width, self.height);
            encoder.set_color(png::ColorType::Rgb);
            encoder.set_depth(png::BitDepth::Eight);
            let mut writer = encoder.write_header().map_err(|e| io::Error::other(e.to_string()))?;
            writer.write_image_data(&self.pixels).map_err(|e| io::Error::other(e.to_string()))?;
        }
        Ok(buf)
    }
}

use std::io;

fn canvas_size(num_positions: usize) -> (u32, u32) {
    let width = MARGIN_LEFT + num_positions as u32 * BAR_WIDTH + MARGIN_RIGHT;
    let height = MARGIN_TOP + LOGO_HEIGHT + MARGIN_BOTTOM;
    (width, height)
}

fn draw_highlight(canvas: &mut Canvas, highlight: Option<(i64, i64)>, num_positions: usize, top: u32, plot_height: u32) {
    let Some((start, end)) = highlight else { return };
    let x0 = MARGIN_LEFT + (start.max(0) as u32) * BAR_WIDTH;
    let x1 = MARGIN_LEFT + ((end.min(num_positions as i64 - 1).max(start) + 1) as u32) * BAR_WIDTH;
    canvas.fill_rect(x0, top, x1.saturating_sub(x0), plot_height, HIGHLIGHT);
}

fn draw_position_ticks(canvas: &mut Canvas, num_positions: usize, axis_y: u32) {
    for pos in 0..num_positions {
        if num_positions <= 40 || pos % 5 == 0 {
            let label = pos.to_string();
            let mut x = MARGIN_LEFT + pos as u32 * BAR_WIDTH + 2;
            for ch in label.chars() {
                canvas.draw_glyph(ch, Rect { x0: x, y0: axis_y + 4, w: 6, h: 8 }, (0, 0, 0), false);
                x += 7;
            }
        }
    }
}

/// Renders the plain (non-centered) sequence logo as a PNG, matching
/// `render_frequency_logo_svg`'s data (see that function's own doc
/// comment for the semantics: letter-stack height by relative
/// frequency per position).
pub fn render_frequency_logo_png(bases: &[char], rows: &[Vec<i64>], stack_order: StackOrder, highlight: Option<(i64, i64)>) -> io::Result<Vec<u8>> {
    let num_positions = rows.len();
    let (width, height) = canvas_size(num_positions);
    let axis_y = MARGIN_TOP + LOGO_HEIGHT;
    let mut canvas = Canvas::new(width, height, WHITE);
    draw_highlight(&mut canvas, highlight, num_positions, MARGIN_TOP, LOGO_HEIGHT);

    for (pos, row) in rows.iter().enumerate() {
        let total: i64 = row.iter().sum();
        if total <= 0 {
            continue;
        }
        let freqs: Vec<f64> = row.iter().map(|&c| c as f64 / total as f64).collect();
        let order = order_indices(&freqs, stack_order);
        let x0 = MARGIN_LEFT + pos as u32 * BAR_WIDTH;
        let mut cumulative_px = 0u32;
        for &i in &order {
            let seg_height_px = (freqs[i] * LOGO_HEIGHT as f64).round() as u32;
            if seg_height_px == 0 {
                continue;
            }
            let y0 = axis_y.saturating_sub(cumulative_px + seg_height_px);
            canvas.draw_glyph(bases[i], Rect { x0: x0 + 1, y0, w: BAR_WIDTH.saturating_sub(2), h: seg_height_px }, base_color_rgb(bases[i]), false);
            cumulative_px += seg_height_px;
        }
    }

    draw_position_ticks(&mut canvas, num_positions, axis_y);
    canvas.encode_png()
}

/// Renders the mean-centered sequence logo as a PNG, matching
/// `render_mean_centered_logo_svg`'s data (positions' counts minus
/// their own mean; below-axis letters flipped).
pub fn render_mean_centered_logo_png(bases: &[char], rows: &[Vec<i64>], stack_order: StackOrder, highlight: Option<(i64, i64)>) -> io::Result<Vec<u8>> {
    let num_positions = rows.len();
    let (width, height) = canvas_size(num_positions);
    let half_height = LOGO_HEIGHT / 2;
    let axis_y = MARGIN_TOP + half_height;
    let mut canvas = Canvas::new(width, height, WHITE);
    draw_highlight(&mut canvas, highlight, num_positions, MARGIN_TOP, LOGO_HEIGHT);

    let centered: Vec<Vec<f64>> = rows
        .iter()
        .map(|row| {
            let mean = row.iter().sum::<i64>() as f64 / row.len().max(1) as f64;
            row.iter().map(|&c| c as f64 - mean).collect()
        })
        .collect();
    let max_abs = centered.iter().flatten().fold(0.0_f64, |acc, &v| acc.max(v.abs())).max(1e-9);

    for (pos, values) in centered.iter().enumerate() {
        let order = order_indices(values, stack_order);
        let x0 = MARGIN_LEFT + pos as u32 * BAR_WIDTH;
        let mut cumulative_pos_px = 0i64;
        let mut cumulative_neg_px = 0i64;
        for &i in &order {
            let scaled_px = (values[i] / max_abs * half_height as f64).round() as i64;
            if scaled_px == 0 {
                continue;
            }
            if scaled_px > 0 {
                let seg = scaled_px as u32;
                let y0 = axis_y as i64 - cumulative_pos_px - scaled_px;
                canvas.draw_glyph(bases[i], Rect { x0: x0 + 1, y0: y0.max(0) as u32, w: BAR_WIDTH.saturating_sub(2), h: seg }, base_color_rgb(bases[i]), false);
                cumulative_pos_px += scaled_px;
            } else {
                let seg = (-scaled_px) as u32;
                let y0 = axis_y as i64 + cumulative_neg_px;
                canvas.draw_glyph(bases[i], Rect { x0: x0 + 1, y0: y0.max(0) as u32, w: BAR_WIDTH.saturating_sub(2), h: seg }, base_color_rgb(bases[i]), true);
                cumulative_neg_px += seg as i64;
            }
        }
    }

    draw_position_ticks(&mut canvas, num_positions, MARGIN_TOP + LOGO_HEIGHT);
    canvas.encode_png()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn decode_png(bytes: &[u8]) -> (u32, u32, Vec<u8>) {
        let decoder = png::Decoder::new(std::io::Cursor::new(bytes));
        let mut reader = decoder.read_info().unwrap();
        let mut buf = vec![0; reader.output_buffer_size().unwrap()];
        let info = reader.next_frame(&mut buf).unwrap();
        (info.width, info.height, buf[..info.buffer_size()].to_vec())
    }

    #[test]
    fn frequency_logo_png_is_well_formed_and_has_colored_pixels() {
        let bases = vec!['A', 'C'];
        let rows = vec![vec![3, 1], vec![0, 4]];
        let png_bytes = render_frequency_logo_png(&bases, &rows, StackOrder::BigOnTop, None).unwrap();
        assert_eq!(&png_bytes[..8], &[0x89, b'P', b'N', b'G', b'\r', b'\n', 0x1a, b'\n']);
        let (w, h, pixels) = decode_png(&png_bytes);
        assert!(w > 0 && h > 0);
        // At least one non-white pixel should exist (a drawn glyph).
        assert!(pixels.chunks(3).any(|p| p != [255, 255, 255]));
    }

    #[test]
    fn mean_centered_logo_png_is_well_formed() {
        let bases = vec!['A', 'C'];
        let rows = vec![vec![3, 1]];
        let png_bytes = render_mean_centered_logo_png(&bases, &rows, StackOrder::Fixed, None).unwrap();
        let (w, h, pixels) = decode_png(&png_bytes);
        assert!(w > 0 && h > 0);
        assert!(pixels.chunks(3).any(|p| p != [255, 255, 255]));
    }

    #[test]
    fn zero_total_position_draws_no_letter_glyph() {
        // The position-0 axis tick label ("0", drawn in black) is
        // still expected -- only the colored LETTER glyph should be
        // absent for a position with no data.
        let bases = vec!['A'];
        let rows = vec![vec![0]];
        let png_bytes = render_frequency_logo_png(&bases, &rows, StackOrder::Fixed, None).unwrap();
        let (_, _, pixels) = decode_png(&png_bytes);
        assert!(!pixels.chunks(3).any(|p| p == [base_color_rgb('A').0, base_color_rgb('A').1, base_color_rgb('A').2]));
    }
}
