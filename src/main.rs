use std::io::prelude::*;
use std::io::{BufWriter,
              BufReader,
              BufRead};
use std::path::PathBuf;
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::error::Error;
use itertools::Itertools;
use rayon::prelude::*;
use crossbeam_channel::unbounded;
use std::thread;
use gzp::{par::compress::{ParCompress, ParCompressBuilder},
         deflate::Bgzf, ZWriter};
use clap::{Arg, Command};

fn is_fastq(v: &str) -> Result<(), String> {
    if v.ends_with("fastq.gz") || v.ends_with("fq.gz") {
        return Ok(());
    }
    Err(format!("{} is not a fastq.gz file", v))
 }

fn cli() -> Command<'static> {
    Command::new("turbotrim")
        .args(&[
            Arg::new("read1")
                .long("read1")
                .short('a')
                .help("first fastq.gz file")
                .takes_value(true)
                .value_name("FILE")
                .validator(is_fastq)
                .required(true),
            Arg::new("read2")
                .long("read2")
                .short('b')
                .help("second fastq.gz file")
                .takes_value(true)
                .value_name("FILE")
                .validator(is_fastq)
                .required(true),
            Arg::new("trim-len")
                .long("trim-len")
                .short('l')
                .help("number of bases to trim")
                .default_value("10")
                .value_name("NUMB")
                .validator(|v| v.parse::<usize>()),
            Arg::new("out1")
                .long("out1")
                .short('x')
                .help("first output fastq.gz file")
                .takes_value(true)
                .value_name("FILE")
                .validator(is_fastq)
                .required(true),
            Arg::new("out2")
                .long("out2")
                .short('y')
                .help("second output fastq.gz file")
                .takes_value(true)
                .value_name("FILE")
                .validator(is_fastq)
                .required(true), 
            ])

}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = cli().get_matches();
    let read1_path = matches.value_of("read1")
                            .expect("error parsing read1");
    let read2_path = matches.value_of("read2")
                            .expect("error parsing read2");

    let out1_path = matches.value_of("out1")
                            .expect("error parsing out1")
                            .to_string();
    let out2_path = matches.value_of("out2")
                            .expect("error parsing out2")
                            .to_string();

    let trim_len = matches.value_of("trim-len")
                            .expect("error parsing trim-len")
                            .parse::<usize>()
                            .unwrap();

    let fq1_f = BufReader::new(File::open(PathBuf::from(read1_path))?);
    let fq1_gz =  MultiGzDecoder::new(fq1_f);
    let fq1 = BufReader::with_capacity(100 * 1024_usize.pow(2), fq1_gz);
    let fq2_f = BufReader::new(File::open(PathBuf::from(read2_path))?);
    let fq2_gz = MultiGzDecoder::new(fq2_f);
    let fq2 = BufReader::with_capacity(100 * 1024_usize.pow(2), fq2_gz);


    let (s, r) = unbounded();
    let (s_chunk, r_chunk) = unbounded();

    
    let trim_hdl = thread::Builder::new()
                            .name("trimer".to_string())
                            .spawn(move || {
        r_chunk.into_iter().par_bridge().for_each(|(mut read1, mut read2): (Vec<String>, Vec<String>)| {
            let drain1 = read1[1].len().saturating_sub(trim_len);
            let drain2 = read2[1].len().saturating_sub(trim_len);
            
            read1[1].truncate(drain1);
            read1[3].truncate(drain1);
            let read1 = read1.join("\n") + "\n";

            read2[1].truncate(drain2);
            read2[3].truncate(drain2);
            let read2 = read2.join("\n") + "\n";

            s.send((read1, read2)).unwrap();
        })}).unwrap();

    let write_hdl = thread::Builder::new()
                            .name("writer".to_string())
                            .spawn(move || {

                                let out_f1 = File::create(PathBuf::from(out1_path)).unwrap();
                                let out_buf1 = BufWriter::with_capacity(100 * 1024 * 1024, out_f1);
                                let mut out_gz1: ParCompress<Bgzf> = ParCompressBuilder::new()
                                                        .num_threads(4).unwrap()
                                                        .from_writer(out_buf1);

                                let out_f2 = File::create(PathBuf::from(out2_path)).unwrap();
                                let out_buf2 = BufWriter::with_capacity(100 * 1024 * 1024, out_f2);
                                let mut out_gz2: ParCompress<Bgzf> = ParCompressBuilder::new()
                                                        .num_threads(4).unwrap()
                                                        .from_writer(out_buf2);

                                for (read1, read2) in r {
                                    out_gz1.write_all(read1.as_bytes())
                                            .unwrap();
                                    out_gz2.write_all(read2.as_bytes())
                                            .unwrap();
                                }

                                out_gz1.finish()
                                        .unwrap();
                                out_gz2.finish()
                                        .unwrap();
                                
                            }).unwrap();

    fq1.lines()
        .zip(fq2.lines())
        .map(|(l1, l2)| (l1.unwrap(), l2.unwrap()))
        .chunks(4)
        .into_iter()
        .for_each(|chunk| {
            let (read1, read2): (Vec<String>, Vec<String>) = chunk.multiunzip();
            s_chunk.send((read1, read2)).unwrap();
        });

    drop(s_chunk);
    let trim_status = trim_hdl.join();
    let write_status = write_hdl.join();
    
    assert!(write_status.is_ok());
    assert!(trim_status.is_ok());
    
    Ok(())
}
