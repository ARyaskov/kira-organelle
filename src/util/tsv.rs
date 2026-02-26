use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

pub struct TsvReader {
    reader: BufReader<File>,
    line_buf: String,
}

impl TsvReader {
    pub fn open(path: &Path) -> io::Result<Self> {
        let file = File::open(path)?;
        Ok(Self {
            reader: BufReader::new(file),
            line_buf: String::with_capacity(16 * 1024),
        })
    }

    pub fn read_record(&mut self, fields: &mut Vec<Range<usize>>) -> io::Result<bool> {
        self.line_buf.clear();
        fields.clear();

        let n = self.reader.read_line(&mut self.line_buf)?;
        if n == 0 {
            return Ok(false);
        }

        let line = self.line_buf.trim_end_matches(['\n', '\r']);
        split_ranges(line, fields);
        Ok(true)
    }

    pub fn field<'a>(&'a self, fields: &'a [Range<usize>], idx: usize) -> Option<&'a str> {
        let range = fields.get(idx)?;
        let line = self.line_buf.trim_end_matches(['\n', '\r']);
        line.get(range.clone())
    }

    pub fn fields_len(fields: &[Range<usize>]) -> usize {
        fields.len()
    }
}

pub fn split_ranges(line: &str, out: &mut Vec<Range<usize>>) {
    out.clear();
    let mut start = 0usize;
    for (idx, b) in line.as_bytes().iter().enumerate() {
        if *b == b'\t' {
            out.push(start..idx);
            start = idx + 1;
        }
    }
    out.push(start..line.len());
}
