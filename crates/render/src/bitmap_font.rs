//! A tiny embedded 5x7 pixel bitmap font, covering only the characters
//! `seqlogo`'s PNG renderer needs: the six nucleotide-code letters
//! (`A C G T N X`) and digits `0-9` for axis position labels. This is
//! a classic, public-domain-style dot-matrix letterform (the same
//! basic 5x7 glyph shapes used by countless embedded/LED-matrix
//! projects) hand-encoded here, not derived from or copying any
//! specific font file -- there is no licensing concern with
//! reproducing standard letter/digit SHAPES at this resolution.
//!
//! Each glyph is 7 rows of 5 bits (MSB-first, bit 4 = leftmost pixel).

pub const GLYPH_WIDTH: usize = 5;
pub const GLYPH_HEIGHT: usize = 7;

pub fn glyph_rows(c: char) -> [u8; GLYPH_HEIGHT] {
    match c {
        'A' => [0b01110, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001],
        'C' => [0b01111, 0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b01111],
        'G' => [0b01110, 0b10001, 0b10000, 0b10111, 0b10001, 0b10001, 0b01110],
        'T' => [0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100],
        'N' => [0b10001, 0b11001, 0b10101, 0b10101, 0b10011, 0b10001, 0b10001],
        'X' => [0b10001, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001, 0b10001],
        '0' => [0b01110, 0b10011, 0b10101, 0b10101, 0b10101, 0b11001, 0b01110],
        '1' => [0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110],
        '2' => [0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111],
        '3' => [0b11110, 0b00001, 0b00001, 0b01110, 0b00001, 0b00001, 0b11110],
        '4' => [0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010],
        '5' => [0b11111, 0b10000, 0b11110, 0b00001, 0b00001, 0b10001, 0b01110],
        '6' => [0b00110, 0b01000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110],
        '7' => [0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b01000, 0b01000],
        '8' => [0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110],
        '9' => [0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00010, 0b01100],
        _ => [0; GLYPH_HEIGHT],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_supported_glyph_has_at_least_one_lit_pixel() {
        for c in "ACGTNX0123456789".chars() {
            let rows = glyph_rows(c);
            assert!(rows.iter().any(|&r| r != 0), "glyph for {c:?} is blank");
        }
    }

    #[test]
    fn unsupported_character_is_blank_not_a_panic() {
        assert_eq!(glyph_rows('?'), [0; GLYPH_HEIGHT]);
    }
}
