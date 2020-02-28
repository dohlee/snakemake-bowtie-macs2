ASCP_BIN = config['ascp_bin']
ASCP_KEY = config['ascp_key']

rule prefetch_accession:
    output:
        temp('{run}.sra')
    resources:
        network = 1
    shell:
        f'prefetch --ascp-path "{ASCP_BIN}|{ASCP_KEY}" -v {{wildcards.run}} && mv {{wildcards.run}}/{{wildcards.run}}.sra . && rm -r {{wildcards.run}}'

c = config['parallel_fastq_dump']
rule parallel_fastq_dump_single:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        '{run}.sra'
    output:
        # Required output.
        temp(DATA_DIR / '{run}.fastq.gz')
    params:
        extra = c['extra']
    threads: config['threads']['parallel_fastq_dump']
    wrapper:
        'http://dohlee-bio.info:9193/parallel-fastq-dump'

rule parallel_fastq_dump_paired:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        DATA_DIR / '{run}.sra'
    output:
        # Required output.
        temp(DATA_DIR / '{run}.read1.fastq.gz'),
        temp(DATA_DIR / '{run}.read2.fastq.gz'),
    params:
        # Optional parameters. Omit if unused.
        extra = c['extra']
    threads: config['threads']['parallel_fastq_dump']
    wrapper:
        'http://dohlee-bio.info:9193/parallel-fastq-dump'

