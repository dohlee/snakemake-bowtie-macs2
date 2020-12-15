rule to_tdf:
    input:
        RESULT_DIR / '05_macs2_bdgcmp' / '{name}_ppois.bdg'
    output:
        RESULT_DIR / '05_macs2_bdgcmp' / '{name}_ppois.bdg.tdf'
    params:
        genome_version = 'hg38'
    threads: 1
    log: 'logs/igvtools/to_tdf/{name}.log'
    benchmark: repeat('benchmarks/igvtools/to_tdf/{name}.benchmark', 1)
    wrapper:
        'http://dohlee-bio.info:9193/igvtools/totdf'

