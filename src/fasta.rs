// A (still partial) re-implementation of
// https://github.com/galaxyproject/sequence_utils/blob/master/galaxy_utils/sequence/fasta.py
// in Rust with using code adapted from 
// https://github.com/rust-bio/rust-bio/blob/master/src/io/fasta.rs

use std::fmt;
use std::io;

/// Trait for FASTA readers.
pub trait FastaRead {
    fn read(&mut self, record: &mut FastaSequence) -> io::Result<()>;
}

/// A FASTA reader.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct FastaReader<B> {
    reader: B,
    error_has_occured: bool,
    line_cache: String, // cache last (header) line obtained from reader
}

impl<B> FastaReader<B>
where
    B: io::BufRead,
{
    /// Create a new Fasta reader with an object that implements `io::BufRead`.
    pub fn new(bufreader: B) -> Self {
        FastaReader {
            reader: bufreader,
            error_has_occured: false,
            line_cache: String::new(),
        }
    }
}

impl<B> Iterator for FastaReader<B>
where
    B: io::BufRead,
{
    type Item = io::Result<FastaSequence>;

    fn next(&mut self) -> Option<io::Result<FastaSequence>> {
        if self.error_has_occured {
            None
        } else {
            let mut record = FastaSequence::new();
            match self.read(&mut record) {
                Ok(()) if record.is_empty() => None,
                Ok(()) => Some(Ok(record)),
                Err(err) => {
                    self.error_has_occured = true;
                    Some(Err(err))
                }
            }
        }
    }
}

impl<B> FastaRead for FastaReader<B>
where
    B: io::BufRead,
{
    /// Read the next FASTA record.
    /// An Ok, but empty result indicates that there are no more records in
    /// the input.
    fn read(&mut self, record: &mut FastaSequence) -> io::Result<()> {
        record.clear();
        if self.line_cache.is_empty() {
            self.reader.read_line(&mut self.line_cache)?;
            if self.line_cache.is_empty() {
                return Ok(());
            }
        }

        if !self.line_cache.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at record start.",
            ));
        }
        let mut header_fields = self.line_cache[1..].trim_end().splitn(2, char::is_whitespace);
        record.id = header_fields.next().map(|s| s.to_owned()).unwrap();
        record.desc = header_fields.next().map(|s| s.to_owned());
        loop {
            self.line_cache.clear();
            self.reader.read_line(&mut self.line_cache)?;
            if self.line_cache.is_empty() || self.line_cache.starts_with('>') {
                break;
            }
            record.seq.push_str(self.line_cache.trim_end());
        }

        Ok(())
    }
}


/// A FASTA record.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct FastaSequence {
    pub id: String,
    pub desc: Option<String>,
    pub seq: String,
}

impl FastaSequence {
    /// Create a new instance.
    pub fn new() -> Self {
        FastaSequence {
            id: String::new(),
            desc: None,
            seq: String::new(),
        }
    }

    /// Create a new `FastaSequence` from given attributes.
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: &str) -> Self {
        let desc = desc.map(|desc| desc.to_owned());
        FastaSequence {
            id: id.to_owned(),
            desc,
            seq: seq.to_owned(),
        }
    }
    
    /// Get the length of the sequence in bases.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Check if record is empty.
    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_none() && self.seq.is_empty()
    }

    /// Check validity of Fasta record.
    pub fn check(&self) -> Result<(), &str> {
        if self.id.is_empty() {
            return Err("Expecting id for Fasta record.");
        }
        if !self.seq.is_ascii() {
            return Err("Non-ascii character found in sequence.");
        }

        Ok(())
    }

    /// Clear the record.
    fn clear(&mut self) {
        self.id.clear();
        self.desc = None;
        self.seq.clear();
    }
}

impl fmt::Display for FastaSequence {
    /// Allows for using `FastaSequence` in a given formatter `f`. In general this is for
    /// creating a `String` representation of a `FastaSequence` and, optionally, writing it to
    /// a file.
    ///
    /// # Errors
    /// Returns [`std::fmt::Error`](https://doc.rust-lang.org/std/fmt/struct.Error.html)
    /// if there is an issue formatting to the stream.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let header = match self.desc.to_owned() {
            Some(d) => format!("{} {}", self.id.to_owned(), d),
            None => self.id.to_owned(),
        };
        write!(
            f,
            ">{}\n{}\n",
            header,
            self.seq.to_owned(),
        )
    }
}

