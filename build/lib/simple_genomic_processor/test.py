#!/data/home/users/k.korfmann/software/anaconda3/bin/python

from simple_genomic_pipeline import *

kwargs = {
    "srr"                 :   "ERR3997272",
    "download_directory"  :   "/data/home/users/k.korfmann/tests/pysradb_downloads",
    "data_directory"      :   "/data/home/users/k.korfmann/tests",
    "sample_name"         :   "ERR3997272",
    "reference_genome"    :   "/data/home/users/k.korfmann/reference/p_vivax_ref.fasta"
}



"""sgp = Pipeline(
    [
        #SRR_Downloader(srr, download_directory, data_directory),
        StaticFileGenerator(**kwargs),
        BWA_Mapper(**kwargs),
        MarkDuplicates(**kwargs),
        AddOrReplaceReadGroups(**kwargs),
        FixMateInformation(**kwargs),
        HaplotypeCaller(**kwargs)
   ]
)
"""

#sgp.show_commands()


""" Testing """
def test_download(**kwargs):
    sgp = Pipeline(
        [
            SRR_Downloader(**kwargs)
            ]
    )
    sgp.run()


def test_static_file_generator(**kwargs):
    sgp = Pipeline(
        [
            StaticFileGenerator(**kwargs)
        ]
    ).run()

def test_mapper(**kwargs):
    sgp = Pipeline(
        [
            BWA_Mapper(**kwargs)
        ]
    ).run()

def test_mark_duplicates(**kwargs):
    sgp = Pipeline(
        [
            MarkDuplicates(**kwargs)
        ]
    ).run()

def test_add_or_replace_read_groups(**kwargs):
    sgp = Pipeline(
        [
            AddOrReplaceReadGroups(**kwargs)
        ]
    ).run()

def test_fix_mate_information(**kwargs):
    sgp = Pipeline(
        [
            FixMateInformation(**kwargs)
        ]
    ).run()

def test_haplotype_caller(**kwargs):
    sgp = Pipeline(
        [
            HaplotypeCaller(**kwargs)
        ]
    ).run()


#test_download(**kwargs)
#test_static_file_generator(**kwargs)
#test_mapper(**kwargs)
#test_mark_duplicates(**kwargs)
#test_add_or_replace_read_groups(**kwargs)
#test_fix_mate_information(**kwargs)
#test_haplotype_caller(**kwargs)


def show_command(step):
    sgp = Pipeline([
            step
        ]
    ).show_commands()

#sgp = Pipeline(
#    [
#        BWA_Mapper(**kwargs)
#    ]
#)
#sgp.show_commands()



