import re
import os
import sys
import subprocess
from pathlib import Path

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from pysradb.sraweb import SRAweb
db = SRAweb()

from contextlib import contextmanager
from time import time
import datetime
from clint.textui import puts, colored, indent


@contextmanager
def timing():
    start = time()
    yield
    ellapsed_time = int(time() - start)
    ellapsed_time = str(datetime.timedelta(seconds=ellapsed_time))
    print(f" - Step finished in {ellapsed_time}.")


class Pipeline():
    def __init__(self, steps=[]):
        self.steps = steps
        
    def run(self):
        print('\033[1mRunning \'Simple Genomic Processor\' with parameters: \033[0;0m')
        for step in self.steps:
            print("\t>>>", end= " ")
            with timing():
                print(colored.yellow(step.__class__.__name__), end=" ")
                step() 
        
    def show_commands(self):
        for step in self.steps:
            print(step.command, end="\n\n")

class Step():
    def __init__(self, **kwargs):
        if 'data_directory' in kwargs and 'sample_name' in kwargs and 'reference_genome' in kwargs:
            self.data_directory = Path(kwargs['data_directory'])
            self.sample_name = Path(kwargs['sample_name'])
            self.reference_genome = Path(kwargs['reference_genome'])
            self.verbose = kwargs['verbose']
    
    
    def __call__(self, stdout=False, stderr=False):
        p = subprocess.Popen(self.command,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)

        if self.verbose:
            stdout = True
            stderr = True

        if stdout:
            try:
                for line in p.stdout.readlines():
                    print(line.decode("ascii"))
            except AttributeError:
                pass

        if stderr:
            try:
                for line in p.stderr.readlines():
                    print(line.decode("ascii"))
            except AttributeError:
                pass

        retval = p.wait()      

class StaticFileGenerator(Step):
    def __init__(self, bwa=True, picard=True, **kwargs):
        super().__init__(**kwargs)
        self.bwa = bwa
        self.picard = picard
        
        self.command = [
            "samtools faidx",
            str(self.reference_genome),
            ";"
        ]

        if self.bwa and not self.picard:
            self.command += ["bwa index",
                             str(self.reference_genome),
                            ]
            
        if "fasta" in str(self.reference_genome):
            dict_file = re.sub("fasta", "dict", str(self.reference_genome))
        elif "fa" in str(self.reference_genome):
            dict_file = re.sub("fa", "dict", str(self.reference_genome))


        if self.picard and not self.bwa:           
            self.command += [
                "picard CreateSequenceDictionary",
                "REFERENCE=" + str(self.reference_genome),
                "OUTPUT=" + dict_file
             ]
            
        if self.bwa and self.picard:
            self.command += ["bwa index",
                             str(self.reference_genome),
                             ";"]
            dict_file = re.sub("fasta", "dict", str(self.reference_genome))
            self.command += [
                "picard CreateSequenceDictionary",
                "REFERENCE=" + str(self.reference_genome),
                "OUTPUT=" + dict_file
             ]
            
        self.command = " ".join(self.command)  


class BWA_Mapper(Step):
    def __init__(self, threads=1, 
                 suffix_fastq_1 = "_pass_1.fastq.gz", 
                 suffix_fastq_2 = "_pass_2.fastq.gz",
                 suffix_output = ".sorted.bam",
                 **kwargs
                 ):
        
        
        super().__init__(**kwargs)
        self.threads = threads
        self.command = ("bwa mem -t", 
                        str(self.threads),
                        str(self.reference_genome),
                        str(self.data_directory/self.sample_name) + suffix_fastq_1,
                        str(self.data_directory/self.sample_name) + suffix_fastq_2,
                        "|",
                        "samtools view -@",
                        str(self.threads),
                        "-b -o",
                        str(self.data_directory/self.sample_name) + ".mapped.bam",
                        ";",
                        "samtools sort -@",
                        str(self.threads),
                        "-T",
                        str(self.data_directory/self.sample_name/"temp"),
                        "-o",
                        str(self.data_directory/self.sample_name) + suffix_output,
                        "-O bam",
                        str(self.data_directory/self.sample_name) + ".mapped.bam",
                        "; rm -f",
                        str(self.data_directory/self.sample_name) + ".bam"
                       )
        self.command = " ".join(self.command)

# remove duplicates also?
class MarkDuplicates(Step):
    def __init__(self, suffix_input=".sorted.bam", suffix_output=".markdup.bam",**kwargs):
        super().__init__(**kwargs)
        self.command = ("picard MarkDuplicates", 
                        "I=" + str(self.data_directory/self.sample_name) + suffix_input,
                        "O=" + str(self.data_directory/self.sample_name) + suffix_output,
                        "M=" + str(self.data_directory/self.sample_name) + ".metrics.txt",
                        "CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 MAX_RECORDS_IN_RAM=10000000 VALIDATION_STRINGENCY=SILENT >>",
                        str(self.data_directory/self.sample_name) + "postmap.log.txt"
                       )
        self.command = " ".join(self.command)

class AddOrReplaceReadGroups(Step):
    def __init__(self, suffix_input=".markdup.bam", suffix_output=".readgroups.bam",**kwargs):
        super().__init__(**kwargs)
        self.command = ("picard AddOrReplaceReadGroups", 
                        "I=" + str(self.data_directory/self.sample_name) + suffix_input,
                        "O=" + str(self.data_directory/self.sample_name) + suffix_output,
                        "RGID=" + str(self.sample_name),
                        "RGLB=" + str(self.sample_name),
                        "RGPL=" + "ILLUMINA",
                        "RGPU=" + str(self.sample_name),
                        "RGSM=" + str(self.sample_name),
                        ">>",
                        str(self.data_directory/self.sample_name) + "postmap.log.txt"
                       )
        self.command = " ".join(self.command)

class FixMateInformation(Step):
    def __init__(self, suffix_input=".readgroups.bam", suffix_output=".fixm.bam",**kwargs):
        super().__init__(**kwargs)
        self.command = ("picard FixMateInformation", 
                        "I=" + str(self.data_directory/self.sample_name) + suffix_input,
                        "O=" + str(self.data_directory/self.sample_name) + suffix_output,
                        "TMP_DIR=" + str(self.data_directory/self.sample_name/"temp"),
                        ">>",
                        str(self.data_directory/self.sample_name) + "postmap.log.txt"
                       )
        self.command = " ".join(self.command)

class HaplotypeCaller(Step):
    def __init__(self, suffix_input=".fixm.bam", suffix_output= ".g.vcf.gz", java_memory="-Xmx15g",interval="",remove_genomic_from_gvcf=False,**kwargs):
        super().__init__(**kwargs)
        self.remove_genomic_from_gvcf=remove_genomic_from_gvcf,
        self.interval=interval
        self.command = [
            "samtools index",
            str(self.data_directory/self.sample_name) + suffix_input,
            ";"
        ]

        if interval:
            L = interval
        else:
            L = ""

        self.command += ["gatk",
                        "--java-options",
                        java_memory,
                        "HaplotypeCaller",
                        "--reference",
                        str(self.reference_genome),
                        "--input",
                        str(self.data_directory/self.sample_name) + suffix_input,
                        "--output",
                        str(self.data_directory/self.sample_name) + suffix_output,
                        "--emit-ref-confidence GVCF",
                        L,
                        ">>",
                        str(self.data_directory/self.sample_name) + "_map_call.log.txt"
                        ]

        if remove_genomic_from_gvcf:
            self.command += ["; gatk",
                    "--java-options",
                    java_memory,
                    "GenotypeGVCFs",
                    "-R",
                    str(self.reference_genome),
                    "-V",
                    str(self.data_directory/self.sample_name) + suffix_output,
                    "-O", 
                    str(self.data_directory/self.sample_name) + ".vcf.gz"
                ]
        self.command = " ".join(self.command)

#class ENA_Downloader():
#   def __init__(self, **kwargs):


class SRR_Downloader():
    def __init__(self, delete_download_directory=True,**kwargs):
        self.srr = kwargs['srr']
        self.download_directory = kwargs['download_directory']
        self.data_directory = kwargs['data_directory']
        self.command = "pysradb"
        self.delete_download_directory = delete_download_directory
    
    def __call__(self):
        
        
        df = db.sra_metadata(self.srr, detailed=True)
        download_directory = os.path.expanduser(self.download_directory)
        
        # redirecting stdout for a second
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = open(os.devnull, "w")
        sys.stderr = open(os.devnull, "w")
        try:
            db.download(df=df, skip_confirmation=True, out_dir=str(download_directory), protocol="ftp")
        #finally:
        #    sys.stdout.close()
        #    sys.stdout = old_stdout
        # reverting stdout 

            for currentpath, folders, files in os.walk(download_directory):
                for file in files:
                    srr_path = os.path.join(currentpath, file)
            #print("Unpacking...")
            cmd_unpack = "fastq-dump "+ srr_path  +" --outdir " + self.data_directory + " --gzip --split-files --skip-technical --readids --read-filter pass --clip > /dev/null" 
            #cmd_unpack = "fasterq-dump -O " + self.data_directory + " --split-files " + srr_path
            #print(cmd_unpack)
            if self.delete_download_directory:
                os.system("rm -rf " + str(self.download_directory))
            os.system(cmd_unpack)
        finally:
            sys.stdout.close()
            sys.stdout = old_stdout
            sys.stderr.close()
            sys.stderr = old_stderr
            # reverting stdout

