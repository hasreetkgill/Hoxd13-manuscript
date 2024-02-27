import os
import glob

TRIMMED_READS_PATH = "outFiles/"

CONDITIONS = {
    'S1':'HindA',
    'S2':'HindB',
    'S3':'HindC',
    'S4':'HindD',
    'S5':'HoxA',
    'S6':'HoxB',
    'S7':'HoxC',
    'S8':'HoxD',
    'S9':'MidA',
    'S10':'MidB',
    'S11':'MidC',
    'S12':'MidD',
    }

trimmed_reads = sorted(
    glob.glob(TRIMMED_READS_PATH + "*_R1_trimmed.fastq.gz")
)

for cond in list(CONDITIONS.keys()):

    cond_files = [file for file in trimmed_reads if cond in file]
    LS_cond = CONDITIONS.get(cond)
    
    catstring = "cat %s > %s_%s_Merged_R1.fastq.gz" % (
        ' '.join(cond_files),
        LS_cond,
        cond,
    )
    print(catstring)
    os.system(catstring)
