from pathlib import Path
import pandas as pd

DATA_DIR = Path(config['data_dir'])
RESULT_DIR = Path(config['result_dir'])
m = pd.read_csv('manifest.csv')

def get_library_layout_of_sample(s):
    return m[m.name == s].library_layout.values[0]

def macs2_filterdup_input(wildcards):
    l = get_library_layout_of_sample(wildcards.sample)
    if l == 'single':
        return str(RESULT_DIR / '01_bowtie' / 'se' / '%s.sorted.bam') % wildcards.sample
    else:
        return str(RESULT_DIR / '01_bowtie' / 'pe' / '%s.sorted.bam') % wildcards.sample

# Configurations for macs2_filterdup.
MF = config['macs2_filterdup']
rule macs2_filterdup:
    input:
        macs2_filterdup_input
    output:
        RESULT_DIR / '02_macs2_filterdup' / '{sample}.sorted.filterdup.bed'
    params:
        # Extra options.
        extra = MF['extra'],
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = MF['gsize'],
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = MF['tsize'],
        # Pvalue cutoff for binomial distribution test.
        # Default: 1e-5
        pvalue = MF['pvalue'],
        # MACS2 filterdup's behavior towards duplicate tags/pairs at the exact
        # same location.
        # 'auto': Calculate the maximum tags at the exact same location based on
        # binomial distribution.
        # integer value: Keep at most that much reads.
        # 'all': Keep all duplicates.
        keep_dup = MF['keep_dup'],
        # Set verbose level.
        # 0: only show critical message.
        # 1: show additional warning message.
        # 2: show process information.
        # 3: show debug messages.
        # If you want to know where are the duplicate reads, use 3.
        # Default: 2
        verbose = MF['verbose'],
    threads: config['threads']['macs2_filterdup']  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_filterdup/{sample}.tsv', 1)
    log: 'logs/macs2_filterdup/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/filterdup'

MCN = config['macs2_callpeak_narrow']
rule macs2_callpeak_narrow:
    input:
        # Required input.
        treatment = RESULT_DIR / '02_macs2_filterdup' / '{treatment}.sorted.filterdup.bed',
        # Optional input.
        control = RESULT_DIR / '02_macs2_filterdup' / '{control}.sorted.filterdup.bed',
    output:
        peak = RESULT_DIR /  '03_macs2_callpeak' / '{treatment}_vs_{control}_peaks.narrowPeak',
        excel = RESULT_DIR / '03_macs2_callpeak' / '{treatment}_vs_{control}_peaks.xls',
        summits = RESULT_DIR /  '03_macs2_callpeak' / '{treatment}_vs_{control}_summits.bed',
        model_script = RESULT_DIR / '03_macs2_callpeak' / '{treatment}_vs_{control}_model.r',
    params:
        # Extra options.
        extra = MCN['extra'],
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = MCN['gsize'],
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
        keep_dup = MCN['keep_dup'],
        # Buffer size for incrementally increasing internal array size to store
        # reads alignment information. In most cases, you don't have to change this
        # parameter. However, if there are large number of chromosomes/contigs/scaffolds
        # in your alignment, it's recommended to specify a smaller buffer size in order
        # to decrease memory useage (but it will take longer time to read alignment files).
        # Minimum memory requested for reading an alignment file is about # of CHROMOSOME *
        # BUFFER_SIZE * 2 Bytes.
        # Default: 100000
        buffer_size = MCN['buffer_size'],
        # Experiment name, which will be used to generate output file names.
        # Default: "NA"
        name = MCN['name'],
        # Whether or not to save extended fragment pileup, and local lambda tracks (two files)
        # at every bp into a bedGraph file.
        # Default: False
        bdg = MCN['bdg'],
        # Tells MACS to include trackline with bedGraph files. To include this trackline while
        # displaying bedGraph at UCSC genome browser, can show name and description of the
        # file as well. However my suggestion is to convert bedGraph to bigWig, then show the
        # smaller and faster binary bigWig file at UCSC genome browser, as well as downstream
        # analysis. Require -B to be set.
        # Default: Not include trackline.
        trackline = MCN['trackline'],
        # If True, MACS will save signal per million reads for fragment pileup profiles.
        # Require -B to be set.
        # Default: False
        SPMR = MCN['SPMR'],
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = MCN['tsize'],
        # Band width for picking regions to compute fragment size. This value is only used
        # while building the shifting model.
        # Default: 300
        bw = MCN['bw'],
        # Select the regions within MFOLD range of high-confidence enrichment ratio against
        # background to build model. Fold-enrichment in regions must be lower than upper limit,
        # and higher than the lower limit. Use as "-m 10 30"
        # Default: 5 50
        mfold = MCN['mfold'],
        # Whether turn on the auto pair model process. If se, when MACS failed to build paired model,
        # it will use the noodel settings, the --exsize parameter to extend each tags towards 3'
        # direction. Not to use this automate fixation is a default behavior now.
        # Default: False
        fix_bimodal = MCN['fix_bimodal'],
        # Wheter or not to build the shifting model. If True, MACS will not build model.
        # By default it means shifting size = 100, try to set extsize to change it.
        # Default: False
        nomodel = MCN['nomodel'],
        # The arbitrary shift in bp. Use discretion while setting it other than default value. 
        # When NOMODEL is set, MACS will user this value to move cutting ends (5') towards 5'->3'
        # direction then apply EXTSIZE to extend them to fragments. When this value is negative,
        # ends will be moved toward 3'->5' direction. Recommended to keep it as default 0 for
        # ChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for detecting
        # enriched cuttin gloci such as certain DNAseI-Seq datasets. Note, you can't set values other
        # than 0 if format is BAMPE or BEDPE for paired-end data.
        # Default: 0
        shift = MCN['shift'],
        # The arbitrary extension size in bp. When nomodel is tru, MACS will use this value as
        # fragment size to extend each read towards 3' end, then pile them up.
        # It's exactly twice then umber of obsolete SHIFTSIZE. In previous language,
        # each read is moved 5'->3' direction to middle of fragment by 1/2 d, then extended to
        # both direction with 1/2 d. This is equivalent to say each read is extended towards 5'->3'
        # into a d size fragment. EXTSIZE and SHIFT can be combined when necessary. Check SHIFT option.
        # Default: 200
        extsize = MCN['extsize'],
        # Minium FDR (q-value) cutoff for peak detection. -q and -p are mutually exclusive.
        # Default; 0.05
        qvalue = MCN['qvalue'],
        # Pvalue cutoff for peak detection. -q and -p are mutually exclusive. If pvalue cutoff is
        # set, qvalue will not be calculated and reported as -1 in the final .xls file.
        # Default: not set.
        pvalue = MCN['pvalue'],
        # When set, scale the small sample up to the bigger sample. By default, the bigger dataset
        # will be scaled down towards the smaller dataset, which will lead to smaller p/qvalues and
        # more specific results. Keep in mind that scaling down will bring down background noise more.
        # Default: False
        to_large = MCN['to_large'],
        # When set, use a custome scaling ratio of ChIP/control (e.g. calculated using NCIS) for 
        # linear scaling.
        # Default: ignore
        ratio = MCN['ratio'],
        # When set, random sampling method will scale down the bigger sample. By default, MACS
        # uses linear scaling.
        # Default: False
        down_sample = MCN['down_sample'],
        # Set the random seed while down sampling data. Must be a non-negative integer.
        # Default: 0
        seed = MCN['seed'],
        # Optional directory to store temp files.
        # Default: tmp
        tempdir = MCN['tempdir'],
        # If True, MACS will use fixed background lambda as local lambda for every peak region.
        # Normally, MACS calculates a dynamic local lambda to reflect the local bias due to 
        # potential chromatin structure.
        nolambda = MCN['nolambda'],
        # The small nearby region in basepairs to calculate dynamic lambda. This is used to
        # cap ture the bias near the peak summit region. Invalid if there is no control data.
        # If you set this to 0, MACS will skip slocal lambda calculation. *Note* that MACS will
        # always perform a d-size local lambda calculation. The final local bias should be the
        # maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 1000
        slocal = MCN['slocal'],
        # The large nearby region in basepairs to calculate dynamic lambda. This is used to 
        # capture the surround bias. If you set this to 0, MACS will skip llocal lambda calculation.
        # *Note* that MACS will always perform a d-size local lambda calculation. The final local
        # bias should be the maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 10000
        llocal = MCN['llocal'],
        # If set, MACS will try to call broad peaks by linking nearby highly enriched regions.
        # The linking region is controlled by another cutoff through --linking-cutoff.
        # The maximum linking region length is 4 times of d from MACS.
        # Default: False
        broad = MCN['broad'],
        # Cutoff for broad region. This option is not available unless --broad is set.
        # If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff.
        # Default: 0.1
        broad_cutoff = MCN['broad_cutoff'],
        # While set, MACS2 will analyze number or total length of peaks that can be called
        # by different p-value cutoff then output a summary table to help user decide
        # a better cutoff. The table will be saved in NAME_cutoff_analysis.txt file.
        # Note, minlen and maxgap may affect the results.
        # WARNING: May take ~30 folds longer time to finish.
        # Default: False
        cutoff_analysis = MCN['cutoff_analysis'],
        # If set, MACS will use a more sophisticated signal processing approach to find subpeak
        # summits in each enriched peak region.
        # Default: False
        call_summits = MCN['call_summits'],
        # When set, the value will be used to filter out peaks with low fold-enrichment.
        # Note, MACS2 use 1.0 as pseudocount while calculating fold-enrichment.
        # Default: 1.0
        fe_cutoff = MCN['fe_cutoff'],
    threads: config['threads']['macs2_callpeak_narrow']  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_callpeak/{treatment}_vs_{control}.tsv', 1)
    log: 'logs/macs2_callpeak/{treatment}_vs_{control}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/callpeak'

# Configurations for macs2_callpeak_broad
MCB = config['macs2_callpeak_broad']
rule macs2_callpeak_broad:
    input:
        # Required input.
        treatment = 'YOUR_TREATMENT_SAMPLE',
        # Optional input.
        control = 'YOUR_CONTROL_SAMPLE',
    output:
        peak = RESULT_DIR / '03_macs2_callpeak' / '{treatment}_vs_{control}_peaks.broadPeak',
        excel = RESULT_DIR / '03_macs2_callpeak' / '{treatment}_vs_{control}_peaks.xls',
        summits = RESULT_DIR / '03_macs2_callpeak' / '{treatment}_vs_{control}_summits.bed',
        model_script = RESULT_DIR / '03_macs2_callpeak' / '{treatment}_vs_{control}_model.r',
    params:
        # Extra options.
        extra = MCB['extra'],
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = MCB['gsize'],
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
        keep_dup = MCB['keep_dup'],
        # Buffer size for incrementally increasing internal array size to store
        # reads alignment information. In most cases, you don't have to change this
        # parameter. However, if there are large number of chromosomes/contigs/scaffolds
        # in your alignment, it's recommended to specify a smaller buffer size in order
        # to decrease memory useage (but it will take longer time to read alignment files).
        # Minimum memory requested for reading an alignment file is about # of CHROMOSOME *
        # BUFFER_SIZE * 2 Bytes.
        # Default: 100000
        buffer_size = MCB['buffer_size'],
        # Experiment name, which will be used to generate output file names.
        # Default: "NA"
        name = MCB['name'],
        # Whether or not to save extended fragment pileup, and local lambda tracks (two files)
        # at every bp into a bedGraph file.
        # Default: False
        bdg = MCB['bdg'],
        # Tells MACS to include trackline with bedGraph files. To include this trackline while
        # displaying bedGraph at UCSC genome browser, can show name and description of the
        # file as well. However my suggestion is to convert bedGraph to bigWig, then show the
        # smaller and faster binary bigWig file at UCSC genome browser, as well as downstream
        # analysis. Require -B to be set.
        # Default: Not include trackline.
        trackline = MCB['trackline'],
        # If True, MACS will save signal per million reads for fragment pileup profiles.
        # Require -B to be set.
        # Default: False
        SPMR = MCB['SPMR'],
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = MCB['tsize'],
        # Band width for picking regions to compute fragment size. This value is only used
        # while building the shifting model.
        # Default: 300
        bw = MCB['bw'],
        # Select the regions within MFOLD range of high-confidence enrichment ratio against
        # background to build model. Fold-enrichment in regions must be lower than upper limit,
        # and higher than the lower limit. Use as "-m 10 30"
        # Default: 5 50
        mfold = MCB['mfold'],
        # Whether turn on the auto pair model process. If se, when MACS failed to build paired model,
        # it will use the noodel settings, the --exsize parameter to extend each tags towards 3'
        # direction. Not to use this automate fixation is a default behavior now.
        # Default: False
        fix_bimodal = MCB['fix_bimodal'],
        # Wheter or not to build the shifting model. If True, MACS will not build model.
        # By default it means shifting size = 100, try to set extsize to change it.
        # Default: False
        nomodel = MCB['nomodel'],
        # The arbitrary shift in bp. Use discretion while setting it other than default value. 
        # When NOMODEL is set, MACS will user this value to move cutting ends (5') towards 5'->3'
        # direction then apply EXTSIZE to extend them to fragments. When this value is negative,
        # ends will be moved toward 3'->5' direction. Recommended to keep it as default 0 for
        # ChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for detecting
        # enriched cuttin gloci such as certain DNAseI-Seq datasets. Note, you can't set values other
        # than 0 if format is BAMPE or BEDPE for paired-end data.
        # Default: 0
        shift = MCB['shift'],
        # The arbitrary extension size in bp. When nomodel is tru, MACS will use this value as
        # fragment size to extend each read towards 3' end, then pile them up.
        # It's exactly twice then umber of obsolete SHIFTSIZE. In previous language,
        # each read is moved 5'->3' direction to middle of fragment by 1/2 d, then extended to
        # both direction with 1/2 d. This is equivalent to say each read is extended towards 5'->3'
        # into a d size fragment. EXTSIZE and SHIFT can be combined when necessary. Check SHIFT option.
        # Default: 200
        extsize = MCB['extsize'],
        # Minium FDR (q-value) cutoff for peak detection. -q and -p are mutually exclusive.
        # Default; 0.05
        qvalue = MCB['qvalue'],
        # Pvalue cutoff for peak detection. -q and -p are mutually exclusive. If pvalue cutoff is
        # set, qvalue will not be calculated and reported as -1 in the final .xls file.
        # Default: not set.
        pvalue = MCB['pvalue'],
        # When set, scale the small sample up to the bigger sample. By default, the bigger dataset
        # will be scaled down towards the smaller dataset, which will lead to smaller p/qvalues and
        # more specific results. Keep in mind that scaling down will bring down background noise more.
        # Default: False
        to_large = MCB['to_large'],
        # When set, use a custome scaling ratio of ChIP/control (e.g. calculated using NCIS) for 
        # linear scaling.
        # Default: ignore
        ratio = MCB['ratio'],
        # When set, random sampling method will scale down the bigger sample. By default, MACS
        # uses linear scaling.
        # Default: False
        down_sample = MCB['down_sample'],
        # Set the random seed while down sampling data. Must be a non-negative integer.
        # Default: 0
        seed = MCB['seed'],
        # Optional directory to store temp files.
        # Default: tmp
        tempdir = MCB['tempdir'],
        # If True, MACS will use fixed background lambda as local lambda for every peak region.
        # Normally, MACS calculates a dynamic local lambda to reflect the local bias due to 
        # potential chromatin structure.
        nolambda = MCB['nolambda'],
        # The small nearby region in basepairs to calculate dynamic lambda. This is used to
        # cap ture the bias near the peak summit region. Invalid if there is no control data.
        # If you set this to 0, MACS will skip slocal lambda calculation. *Note* that MACS will
        # always perform a d-size local lambda calculation. The final local bias should be the
        # maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 1000
        slocal = MCB['slocal'],
        # The large nearby region in basepairs to calculate dynamic lambda. This is used to 
        # capture the surround bias. If you set this to 0, MACS will skip llocal lambda calculation.
        # *Note* that MACS will always perform a d-size local lambda calculation. The final local
        # bias should be the maximum of the lambda value from d, slocal, and llocal size windows.
        # Default: 10000
        llocal = MCB['llocal'],
        # If set, MACS will try to call broad peaks by linking nearby highly enriched regions.
        # The linking region is controlled by another cutoff through --linking-cutoff.
        # The maximum linking region length is 4 times of d from MACS.
        # Default: False
        broad = MCB['broad'],
        # Cutoff for broad region. This option is not available unless --broad is set.
        # If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff.
        # Default: 0.1
        broad_cutoff = MCB['broad_cutoff'],
        # While set, MACS2 will analyze number or total length of peaks that can be called
        # by different p-value cutoff then output a summary table to help user decide
        # a better cutoff. The table will be saved in NAME_cutoff_analysis.txt file.
        # Note, minlen and maxgap may affect the results.
        # WARNING: May take ~30 folds longer time to finish.
        # Default: False
        cutoff_analysis = MCB['cutoff_analysis'],
        # If set, MACS will use a more sophisticated signal processing approach to find subpeak
        # summits in each enriched peak region.
        # Default: False
        call_summits = MCB['call_summits'],
        # When set, the value will be used to filter out peaks with low fold-enrichment.
        # Note, MACS2 use 1.0 as pseudocount while calculating fold-enrichment.
        # Default: 1.0
        fe_cutoff = MCB['fe_cutoff'],
    threads: config['threads']['macs2_callpeak_broad']  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_callpeak/{treatment}_vs_{control}.tsv', 1)
    log: 'logs/macs2_callpeak/{treatment}_vs_{control}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/callpeak'

