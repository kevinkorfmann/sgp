#!/usr/bin/env python


from simple_genomic_pipeline import *
from clint.textui import puts, colored, indent
import os
import argparse

parser = argparse.ArgumentParser(description='Simple genomic processor arguments.')
parser.add_argument('--srr',  required=True)
parser.add_argument('--reference_genome', required=True)
parser.add_argument('--download_directory')
parser.add_argument('--data_directory')
parser.add_argument('--sample_name')

parser.add_argument('--suffix_fastq_1')
parser.add_argument('--suffix_fastq_2')
parser.add_argument('--bwa_threads')

args = parser.parse_args()

if args.download_directory is None:
    args.download_directory = "./pysradb_downloads"
    
if args.data_directory is None:
    args.data_directory = "./"

if args.sample_name is None:
    args.sample_name = args.srr

if args.suffix_fastq_1 is not None:
    suffix_fastq_1 = args.suffix_fastq_1


if args.suffix_fastq_2 is not None:
    suffix_fastq_2 = args.suffix_fastq_2

if args.bwa_threads is None:
    bwa_threads = 1
else:
    bwa_threads = args.bwa_threads

print(""" 
  A_A
 (-.-)    ~+ Simple Genomic Processor, meow +~
  |-|
 /   \
|     |   __
|  || |  |  \__
 \_||_/_/

""")

print("Make sure all necessary dependencies are accessible in your $PATH")
print("Dependencies: ")
print("fastq-dump bwa samtools picard gatk")

kwargs = {
    "srr"                 :   args.srr,
    "download_directory"  :   args.download_directory,
    "data_directory"      :   args.data_directory,
    "sample_name"         :   args.sample_name,
    "reference_genome"    :   args.reference_genome
}

print("")
print('\033[1mCommand line parameters: \033[0;0m')

for k,v in zip(kwargs.keys(), kwargs.values()):
    print("\t>>>", end= " ")
    print(colored.blue(k) + " : " + str(v))

print("")
#print(kwargs)
"""
kwargs = {
    "srr"                 :   "ERR3997272",
    "download_directory"  :   "/data/home/users/k.korfmann/tests/pysradb_downloads",
    "data_directory"      :   "/data/home/users/k.korfmann/tests",
    "sample_name"         :   "ERR3997272",
    "reference_genome"    :   "/data/home/users/k.korfmann/reference/p_vivax_ref.fasta"
}
"""


def srr2vcf(suffix_fastq_1="_pass_1.fastq.gz", suffix_fastq_2="_pass_2.fastq.gz", bwa_threads=1, **kwargs):
    sgp = Pipeline(
        [
            SRR_Downloader(**kwargs),
            StaticFileGenerator(**kwargs),
            BWA_Mapper(threads=bwa_threads, 
                suffix_fastq_1 = suffix_fastq_1,
                suffix_fastq_2 = suffix_fastq_1, **kwargs),
            MarkDuplicates(**kwargs),
            AddOrReplaceReadGroups(**kwargs),
            FixMateInformation(**kwargs),
            HaplotypeCaller(**kwargs)
        ]
    ).run()

if args.suffix_fastq_1 is not None and args.suffix_fastq_2 is not None:
    srr2vcf(suffix_fastq_1, suffix_fastq_2, bwa_threads=bwa_threads, **kwargs)
else:
    srr2vcf(bwa_threads=bwa_threads, **kwargs)

