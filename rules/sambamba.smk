rule sambamba_index:
    input:
        RESULT_DIR / '02_bowtie' / '{run}.sorted.bam'
    output:
        RESULT_DIR / '02_bowtie' / '{run}.sorted.bam.bai'
    threads: 1
    log: 'logs/sambamba_index/{run}.log'
    benchmark: 'benchmarks/sambamba_index/{run}.benchmark'
    wrapper:
        'http://dohlee-bio.info:9193/sambamba/index'

