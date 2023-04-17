use gx_sequence_utils_rs::fasta::FastaReader;

const FASTA_FILE: &[u8] = b">id desc
ACCGTAGGCTGA
CCGTAGGCTGAA
CGTAGGCTGAAA
GTAGGCTGAAAA
CCCC
>id2
ATTGTTGTTTTA
ATTGTTGTTTTA
ATTGTTGTTTTA
GGGG
";

fn main() {
    let reader = FastaReader::new(FASTA_FILE);
    for record in reader {
        let faseq = record.unwrap();
        println!("{}length: {}", faseq, faseq.len());
    }
}
