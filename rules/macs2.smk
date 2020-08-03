c = config['macs2_filterdup']
rule macs2_filterdup:
    input:
        RESULT_DIR / '02_bowtie' / '{sample}.sorted.bam'
    output:
        RESULT_DIR / '03_macs2_filterdup' / '{sample}.sorted.filterdup.bed'
    params:
        # Extra options.
        extra = c['extra'],
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = c['gsize'],
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = c['tsize'],
        # Pvalue cutoff for binomial distribution test.
        # Default: 1e-5
        pvalue = c['pvalue'],
        # MACS2 filterdup's behavior towards duplicate tags/pairs at the exact
        # same location.
        # 'auto': Calculate the maximum tags at the exact same location based on
        # binomial distribution.
        # integer value: Keep at most that much reads.
        # 'all': Keep all duplicates.
        keep_dup = c['keep_dup'],
        # Set verbose level.
        # 0: only show critical message.
        # 1: show additional warning message.
        # 2: show process information.
        # 3: show debug messages.
        # If you want to know where are the duplicate reads, use 3.
        # Default: 2
        verbose = c['verbose'],
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_filterdup/{sample}.tsv', 1)
    log: 'logs/macs2_filterdup/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/filterdup'

def macs2_callpeak_narrow_input(wildcards):
    run_t, run_c = name2pair[wildcards.name]
    ret = {}
    ret['treatment'] = str(RESULT_DIR / '03_macs2_filterdup' / f'{run_t}.sorted.filterdup.bed')
    ret['control'] = str(RESULT_DIR / '03_macs2_filterdup' / f'{run_c}.sorted.filterdup.bed')

    return ret

c = config['macs2_callpeak_narrow']
checkpoint macs2_callpeak_narrow:
    input: unpack(macs2_callpeak_narrow_input)
    output:
        peak = RESULT_DIR / '04_macs2_callpeak' / 'narrow' / '{name}_peaks.narrowPeak',
        excel = RESULT_DIR / '04_macs2_callpeak' / 'narrow' / '{name}_peaks.xls',
        summits = RESULT_DIR / '04_macs2_callpeak' / 'narrow' / '{name}_summits.bed',
        model_script = RESULT_DIR / '04_macs2_callpeak' / 'narrow' / '{name}_model.r',
        control_bdg = RESULT_DIR / '04_macs2_callpeak' / 'narrow' / '{name}_control_lambda.bdg',
        treat_bdg = RESULT_DIR / '04_macs2_callpeak' / 'narrow' / '{name}_treat_pileup.bdg',
    params:
        # Extra options.
        extra = c['extra'],
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = c['gsize'],
        # It controls the MACS behavior towards duplicate tags at the
        # exact same location -- the same coordination and the same
        # strand. The 'auto' option makes MACS calculate the maximum tags
        # at the exact same location based on binomial distribution using
        # 1e-5 as pvalue cutoff; and the 'all' option keeps every tags. If
        # an integer is given, at most this number of tags will be kept at
        # the same location. Note, if you've used samtools or picard to flag reads
        # as 'PCR/Optical duplicate' in bit 1024, MACS2 will still read them
        # although the reads may be decided by MACS2 as duplicate later.
        # The default is to keep one tag at the same location.
        # Default: 1
        keep_dup = c['keep_dup'],
        # Buffer size for incrementally increasing internal array size to store
        # reads alignment information. In most cases, you don't have to change this
        # parameter. However, if there are large number of chromosomes/contigs/scaffolds
        # in your alignment, it's recommended to specify a smaller buffer size in order
        # to decrease memory useage (but it will take longer time to read alignment files).
        # Minimum memory requested for reading an alignment file is about # of CHROMOSOME *
        # BUFFER_SIZE * 2 Bytes.
        # Default: 100000
        buffer_size = c['buffer_size'],
        # Experiment name, which will be used to generate output file names.
        # Default: "NA"
        name = c['name'],
        # Whether or not to save extended fragment pileup, and local lambda tracks (two files)
        # at every bp into a bedGraph file.
        # Default: False
        bdg = c['bdg'],
        # Tells MACS to include trackline with bedGraph files. To include this trackline while
        # displaying bedGraph at UCSC genome browser, can show name and description of the
        # file as well. However my suggestion is to convert bedGraph to bigWig, then show the
        # smaller and faster binary bigWig file at UCSC genome browser, as well as downstream
        # analysis. Require -B to be set.
        # Default: Not include trackline.
        trackline = c['trackline'],
        # If True, MACS will save signal per million reads for fragment pileup profiles.
        # Require -B to be set.
        # Default: False
        SPMR = c['SPMR'],
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = c['tsize'],
        # Band width for picking regions to compute fragment size. This value is only used
        # while building the shifting model.
        # Default: 300
        bw = c['bw'],
        # Select the regions within MFOLD range of high-confidence enrichment ratio against
        # background to build model. Fold-enrichment in regions must be lower than upper limit,
        # and higher than the lower limit. Use as "-m 10 30"
        # Default: 5 50
        mfold = c['mfold'],
        # Whether turn on the auto pair model process. If se, when MACS failed to build paired model,
        # it will use the noodel settings, the --exsize parameter to extend each tags towards 3'
        # direction. Not to use this automate fixation is a default behavior now.
        # Default: False
        fix_bimodal = c['fix_bimodal'],
        # Wheter or not to build the shifting model. If True, MACS will not build model.
        # By default it means shifting size = 100, try to set extsize to change it.
        # Default: False
        nomodel = c['nomodel'],
        # The arbitrary shift in bp. Use discretion while setting it other than default value. 
        # When NOMODEL is set, MACS will user this value to move cutting ends (5') towards 5'->3'
        # direction then apply EXTSIZE to extend them to fragments. When this value is negative,
        # ends will be moved toward 3'->5' direction. Recommended to keep it as default 0 for
        # ChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for detecting
        # enriched cuttin gloci such as certain DNAseI-Seq datasets. Note, you can't set values other
        # than 0 if format is BAMPE or BEDPE for paired-end data.
        # Default: 0
        shift = c['shift'],
        # The arbitrary extension size in bp. When nomodel is tru, MACS will use this value as
        # fragment size to extend each read towards 3' end, then pile them up.
        # It's exactly twice then umber of obsolete SHIFTSIZE. In previous language,
        # each read is moved 5'->3' direction to middle of fragment by 1/2 d, then extended to
        # both direction with 1/2 d. This is equivalent to say each read is extended towards 5'->3'
        # into a d size fragment. EXTSIZE and SHIFT can be combined when necessary. Check SHIFT option.
        # Default: 200
        extsize = c['extsize'],
        # Minium FDR (q-value) cutoff for peak detection. -q and -p are mutually exclusive.
        # Default; 0.05
        qvalue = c['qvalue'],
        # Pvalue cutoff for peak detection. -q and -p are mutually exclusive. If pvalue cutoff is
        # set, qvalue will not be calculated and reported as -1 in the final .xls file.
        # Default: not set.
        pvalue = c['pvalue'],
        # When set, scale the small sample up to the bigger sample. By default, the bigger dataset
        # will be scaled down towards the smaller dataset, which will lead to smaller p/qvalues and
        # more specific results. Keep in mind that scaling down will bring down background noise more.
        # Default: False
        to_large = c['to_large'],
        # When set, use a custome scaling ratio of ChIP/control (e.g. calculated using NCIS) for 
        # linear scaling.
        # Default: ignore
        ratio = c['ratio'],
        # When set, random sampling method will scale down the bigger sample. By default, MACS
        # uses linear scaling.
        # Default: False
        down_sample = c['down_sample'],
        # Set the random seed while down sampling data. Must be a non-negative integer.
        # Default: 0
        seed = c['seed'],
        # Optional directory to store temp files.
        # Default: tmp
        tempdir = c['tempdir'],
        # If True, MACS will use fixed background lambda as local lambda for every peak region.
        # Normally, MACS calculates a dynamic local lambda to reflect the local bias due to 
        # potential chromatin structure.
        nolambda = c['nolambda'],
        # The small nearby region in basepairs to calculate dynamic lambda. This is used to
        # cap ture the bias near the peak summit region. Invalid if there is no control data.
        # If you set this to 0, MACS will skip slocal lambda calculation. *Note* that MACS will
        # always perform a d-size local lambda calculation. The final local bias should be the
        # maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 1000
        slocal = c['slocal'],
        # The large nearby region in basepairs to calculate dynamic lambda. This is used to 
        # capture the surround bias. If you set this to 0, MACS will skip llocal lambda calculation.
        # *Note* that MACS will always perform a d-size local lambda calculation. The final local
        # bias should be the maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 10000
        llocal = c['llocal'],
        # If set, MACS will try to call broad peaks by linking nearby highly enriched regions.
        # The linking region is controlled by another cutoff through --linking-cutoff.
        # The maximum linking region length is 4 times of d from MACS.
        # Default: False
        broad = c['broad'],
        # Cutoff for broad region. This option is not available unless --broad is set.
        # If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff.
        # Default: 0.1
        broad_cutoff = c['broad_cutoff'],
        # While set, MACS2 will analyze number or total length of peaks that can be called
        # by different p-value cutoff then output a summary table to help user decide
        # a better cutoff. The table will be saved in NAME_cutoff_analysis.txt file.
        # Note, minlen and maxgap may affect the results.
        # WARNING: May take ~30 folds longer time to finish.
        # Default: False
        cutoff_analysis = c['cutoff_analysis'],
        # If set, MACS will use a more sophisticated signal processing approach to find subpeak
        # summits in each enriched peak region.
        # Default: False
        call_summits = c['call_summits'],
        # When set, the value will be used to filter out peaks with low fold-enrichment.
        # Note, MACS2 use 1.0 as pseudocount while calculating fold-enrichment.
        # Default: 1.0
        fe_cutoff = c['fe_cutoff'],
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_callpeak/{name}.tsv', 1)
    log: 'logs/macs2_callpeak/{name}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/callpeak'

def macs2_callpeak_broad_input(wildcards):
    run_t, run_c = name2pair[wildcards.name]
    ret = {}
    ret['treatment'] = str(RESULT_DIR / '03_macs2_filterdup' / f'{run_t}.sorted.filterdup.bed')
    ret['control'] = str(RESULT_DIR / '03_macs2_filterdup' / f'{run_c}.sorted.filterdup.bed')
    return ret

c = config['macs2_callpeak_broad']
checkpoint macs2_callpeak_broad:
    input: unpack(macs2_callpeak_broad_input)
    output:
        peak = RESULT_DIR / '04_macs2_callpeak' / 'broad' / '{name}_peaks.broadPeak',
        excel = RESULT_DIR / '04_macs2_callpeak' / 'broad' / '{name}_peaks.xls',
        summits = RESULT_DIR / '04_macs2_callpeak' / 'broad' / '{name}_summits.bed',
        model_script = RESULT_DIR / '04_macs2_callpeak' / 'broad' / '{name}_model.r',
        control_bdg = RESULT_DIR / '04_macs2_callpeak' / 'broad' / '{name}_control_lambda.bdg',
        treat_bdg = RESULT_DIR / '04_macs2_callpeak' / 'broad' / '{name}_treat_pileup.bdg',
    params:
        # Extra options.
        extra = c['extra'],
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = c['gsize'],
        # It controls the MACS behavior towards duplicate tags at the
        # exact same location -- the same coordination and the same
        # strand. The 'auto' option makes MACS calculate the maximum tags
        # at the exact same location based on binomial distribution using
        # 1e-5 as pvalue cutoff; and the 'all' option keeps every tags. If
        # an integer is given, at most this number of tags will be kept at
        # the same location. Note, if you've used samtools or picard to flag reads
        # as 'PCR/Optical duplicate' in bit 1024, MACS2 will still read them
        # although the reads may be decided by MACS2 as duplicate later.
        # The default is to keep one tag at the same location.
        # Default: 1
        keep_dup = c['keep_dup'],
        # Buffer size for incrementally increasing internal array size to store
        # reads alignment information. In most cases, you don't have to change this
        # parameter. However, if there are large number of chromosomes/contigs/scaffolds
        # in your alignment, it's recommended to specify a smaller buffer size in order
        # to decrease memory useage (but it will take longer time to read alignment files).
        # Minimum memory requested for reading an alignment file is about # of CHROMOSOME *
        # BUFFER_SIZE * 2 Bytes.
        # Default: 100000
        buffer_size = c['buffer_size'],
        # Experiment name, which will be used to generate output file names.
        # Default: "NA"
        name = c['name'],
        # Whether or not to save extended fragment pileup, and local lambda tracks (two files)
        # at every bp into a bedGraph file.
        # Default: False
        bdg = c['bdg'],
        # Tells MACS to include trackline with bedGraph files. To include this trackline while
        # displaying bedGraph at UCSC genome browser, can show name and description of the
        # file as well. However my suggestion is to convert bedGraph to bigWig, then show the
        # smaller and faster binary bigWig file at UCSC genome browser, as well as downstream
        # analysis. Require -B to be set.
        # Default: Not include trackline.
        trackline = c['trackline'],
        # If True, MACS will save signal per million reads for fragment pileup profiles.
        # Require -B to be set.
        # Default: False
        SPMR = c['SPMR'],
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = c['tsize'],
        # Band width for picking regions to compute fragment size. This value is only used
        # while building the shifting model.
        # Default: 300
        bw = c['bw'],
        # Select the regions within MFOLD range of high-confidence enrichment ratio against
        # background to build model. Fold-enrichment in regions must be lower than upper limit,
        # and higher than the lower limit. Use as "-m 10 30"
        # Default: 5 50
        mfold = c['mfold'],
        # Whether turn on the auto pair model process. If se, when MACS failed to build paired model,
        # it will use the noodel settings, the --exsize parameter to extend each tags towards 3'
        # direction. Not to use this automate fixation is a default behavior now.
        # Default: False
        fix_bimodal = c['fix_bimodal'],
        # Wheter or not to build the shifting model. If True, MACS will not build model.
        # By default it means shifting size = 100, try to set extsize to change it.
        # Default: False
        nomodel = c['nomodel'],
        # The arbitrary shift in bp. Use discretion while setting it other than default value. 
        # When NOMODEL is set, MACS will user this value to move cutting ends (5') towards 5'->3'
        # direction then apply EXTSIZE to extend them to fragments. When this value is negative,
        # ends will be moved toward 3'->5' direction. Recommended to keep it as default 0 for
        # ChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for detecting
        # enriched cuttin gloci such as certain DNAseI-Seq datasets. Note, you can't set values other
        # than 0 if format is BAMPE or BEDPE for paired-end data.
        # Default: 0
        shift = c['shift'],
        # The arbitrary extension size in bp. When nomodel is tru, MACS will use this value as
        # fragment size to extend each read towards 3' end, then pile them up.
        # It's exactly twice then umber of obsolete SHIFTSIZE. In previous language,
        # each read is moved 5'->3' direction to middle of fragment by 1/2 d, then extended to
        # both direction with 1/2 d. This is equivalent to say each read is extended towards 5'->3'
        # into a d size fragment. EXTSIZE and SHIFT can be combined when necessary. Check SHIFT option.
        # Default: 200
        extsize = c['extsize'],
        # Minium FDR (q-value) cutoff for peak detection. -q and -p are mutually exclusive.
        # Default; 0.05
        qvalue = c['qvalue'],
        # Pvalue cutoff for peak detection. -q and -p are mutually exclusive. If pvalue cutoff is
        # set, qvalue will not be calculated and reported as -1 in the final .xls file.
        # Default: not set.
        pvalue = c['pvalue'],
        # When set, scale the small sample up to the bigger sample. By default, the bigger dataset
        # will be scaled down towards the smaller dataset, which will lead to smaller p/qvalues and
        # more specific results. Keep in mind that scaling down will bring down background noise more.
        # Default: False
        to_large = c['to_large'],
        # When set, use a custome scaling ratio of ChIP/control (e.g. calculated using NCIS) for 
        # linear scaling.
        # Default: ignore
        ratio = c['ratio'],
        # When set, random sampling method will scale down the bigger sample. By default, MACS
        # uses linear scaling.
        # Default: False
        down_sample = c['down_sample'],
        # Set the random seed while down sampling data. Must be a non-negative integer.
        # Default: 0
        seed = c['seed'],
        # Optional directory to store temp files.
        # Default: tmp
        tempdir = c['tempdir'],
        # If True, MACS will use fixed background lambda as local lambda for every peak region.
        # Normally, MACS calculates a dynamic local lambda to reflect the local bias due to 
        # potential chromatin structure.
        nolambda = c['nolambda'],
        # The small nearby region in basepairs to calculate dynamic lambda. This is used to
        # cap ture the bias near the peak summit region. Invalid if there is no control data.
        # If you set this to 0, MACS will skip slocal lambda calculation. *Note* that MACS will
        # always perform a d-size local lambda calculation. The final local bias should be the
        # maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 1000
        slocal = c['slocal'],
        # The large nearby region in basepairs to calculate dynamic lambda. This is used to 
        # capture the surround bias. If you set this to 0, MACS will skip llocal lambda calculation.
        # *Note* that MACS will always perform a d-size local lambda calculation. The final local
        # bias should be the maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 10000
        llocal = c['llocal'],
        # If set, MACS will try to call broad peaks by linking nearby highly enriched regions.
        # The linking region is controlled by another cutoff through --linking-cutoff.
        # The maximum linking region length is 4 times of d from MACS.
        # Default: False
        broad = c['broad'],
        # Cutoff for broad region. This option is not available unless --broad is set.
        # If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff.
        # Default: 0.1
        broad_cutoff = c['broad_cutoff'],
        # While set, MACS2 will analyze number or total length of peaks that can be called
        # by different p-value cutoff then output a summary table to help user decide
        # a better cutoff. The table will be saved in NAME_cutoff_analysis.txt file.
        # Note, minlen and maxgap may affect the results.
        # WARNING: May take ~30 folds longer time to finish.
        # Default: False
        cutoff_analysis = c['cutoff_analysis'],
        # If set, MACS will use a more sophisticated signal processing approach to find subpeak
        # summits in each enriched peak region.
        # Default: False
        call_summits = c['call_summits'],
        # When set, the value will be used to filter out peaks with low fold-enrichment.
        # Note, MACS2 use 1.0 as pseudocount while calculating fold-enrichment.
        # Default: 1.0
        fe_cutoff = c['fe_cutoff'],
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_callpeak/{name}.tsv', 1)
    log: 'logs/macs2_callpeak/{name}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/callpeak'

def macs2_bdgcmp_input(wildcards):
    peak_type = name2peaktype[wildcards.name]

    if peak_type.upper().startswith('NARROW'):
        return {
            'treatment': str(RESULT_DIR / '04_macs2_callpeak' / 'narrow' / '{name}_treat_pileup.bdg'),
            'control': str(RESULT_DIR / '04_macs2_callpeak' / 'narrow' / '{name}_control_lambda.bdg'),
        }
    else:
        return {
            'treatment': str(RESULT_DIR / '04_macs2_callpeak' / 'broad' / '{name}_treat_pileup.bdg'),
            'control': str(RESULT_DIR / '04_macs2_callpeak' / 'broad' / '{name}_control_lambda.bdg'),
        }

c = config['macs2_bdgcmp']
rule macs2_bdgcmp:
    input:
        unpack(macs2_bdgcmp_input)
    output:
        # Output suffixes automatically determines the '-m' parameters.
        # '_ppois.bdg': use '-m ppois', Poisson P-value
        # '_qpois.bdg': use '-m qpois', Poisson Q-value
        # '_subtract.bdg': use '-m subtract', subtraction from treatment
        # '_FE.bdg': use '-m FE', linear scale fold enrichment
        # '_logFE.bdg': use '-m logFE', log10 fold enrichment (need to set pseudocount)
        # '_logLR.bdg': use '-m logLR', log10 likelihood between ChIP-enriched model
        # and open chromatin model (need to set pseudocount)
        # '_slogLR.bdg': use '-m slogLR', symmetric log10 likelihood between two ChIP-enrichment
        # models, or maximum value between the two tracks.
        [
            temp(RESULT_DIR / '05_macs2_bdgcmp' / '{name}_ppois.bdg'),
            temp(RESULT_DIR / '05_macs2_bdgcmp' / '{name}_qpois.bdg'),
            temp(RESULT_DIR / '05_macs2_bdgcmp' / '{name}_FE.bdg'),
        ]
    params:
        extra = c['extra'],
        # Scaling factor for treatment and control track. Keep it as
        # 1.0 or default in most cases. Set it ONLY while you have SPMR
        # output from MACS2 callpeak, and plan to calculate scores as MACS2
        # callpeak module. If you want to simulate 'callpeak' w/o '--to-large',
        # calculate effective smaller sample size after filtering redundant
        # reads in million (e.g., put 31.415926 if effective reads are
        # 31,415,926) and input if for '-S'; for 'callpeak --to-large',
        # calculate effective reads in larger sample.
        # Default: 1.0
        scaling_factor = c['scaling_factor'],
        # The pseudocount used for calculating logLR, logFE or FE.
        # The count will be applied after normalization of sequencing
        # depth.
        # Default: False
        pseudocount = c['pseudocount'],
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_bdgcmp/{name}.tsv', 1)
    log: 'logs/macs2_callpeak/{name}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/bdgcmp'

