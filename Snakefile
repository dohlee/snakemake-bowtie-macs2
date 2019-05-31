from pathlib import Path
import pandas as pd

configfile: 'config.yaml'

include: 'rules/bowtie.smk'
include: 'rules/macs2.smk'

manifest = pd.read_csv(config['manifest'])
RESULT_DIR = Path(config['result_dir'])

SAMPLES = manifest.name.values
SE_SAMPLES = manifest[manifest.library_layout == 'single'].name.values
PE_SAMPLES = manifest[manifest.library_layout == 'paired'].name.values

ALIGNED_BAM_SE = expand(str(RESULT_DIR / '01_bowtie' / 'se' / '{sample}.sorted.bam'), sample=SE_SAMPLES)
ALIGNED_BAM_PE = expand(str(RESULT_DIR / '01_bowtie' / 'pe' / '{sample}.sorted.bam'), sample=PE_SAMPLES)

FILTERDUPED = expand(str(RESULT_DIR / '02_macs2_filterdup' / '{sample}.sorted.filterdup.bed'), sample=SAMPLES)

RESULT_FILES = []
RESULT_FILES.append(ALIGNED_BAM_SE)
RESULT_FILES.append(ALIGNED_BAM_PE)
RESULT_FILES.append(FILTERDUPED)

rule all:
    input: RESULT_FILES
