use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();

    let input_filename = &args[1];
    let output_filename = &args[2];

    let num_reads = 0;
    let fastq_read: Option<String> = None;

    if num_reads == 0 {
        println!("No valid FASTQ reads could be processed from {input_filename}");
    } else {
        println!("{num_reads} FASTQ reads were converted to FASTA.");
    }
}
