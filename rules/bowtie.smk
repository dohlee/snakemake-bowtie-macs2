from os.path import join

c = config['bowtie_build']
rule bowtie_build:
    input:
        # Required input.
        # Reference genome fasta.
        # e.g.
        # reference = config['reference']['genome']['fasta']
        config['reference']['fasta'],
    output:
        # e.g.
        # reference_indices = config['reference']['genome']['bowtie_index_dir']
        index_dir = directory(config['bowtie_index_dir']),
    params:
        # Additional parameters go here.
        extra = c['extra'],
        # Build a colorspace index.
        color = c['color'],
        # Disable automatic -p/--bmax/--dcv memory-fitting.
        noauto = c['noauto'],
        # Use packed strings internally; slower, uses less mem.
        packed = c['packed'], 
        # Max bucket size for blockwise suffix-array builder.
        bmax = c['bmax'],
        # Max bucket size as divisor of ref len.
        # Default: 4
        bmaxdivn = c['bmaxdivn'],
        # diff-cover period for blockwise
        # Default: 1024
        dcv = c['dcv'],
        # disable diff-cover (algorithm becomes quadratic)
        nodc = c['nodc'],
        # Don't build .3/.4.ebwt (packed reference) portion.
        noref = c['noref'],
        # Just build .3/.4.ebwt (packed reference) portion.
        justref = c['justref'],
        # SA is sampled every 2^offRate BWT chars.
        # Default: 5
        offrate = c['offrate'],
        # Number of chars consumed in initial lookup.
        # Default: 10
        ftabchars = c['ftabchars'],
        # convert Ns in reference to As.
        ntoa = c['ntoa'],
        # Seed for random number generator.
        seed = c['seed'],
        # Verbose output (for debugging).
        quiet = c['quiet'],
    threads: config['threads']['bowtie_build']
    benchmark:
        repeat('benchmarks/bowtie_build/mm10.tsv', 1)
    log: 'logs/bowtie_build/mm10.log'
    wrapper:
        'http://dohlee-bio.info:9193/bowtie/build'

def bowtie_input(wildcards):
    lib = run2lib[wildcards.run]

    ret = {}
    ret['index_dir'] = directory(config['bowtie_index_dir'])
    if lib.upper().startswith('SINGLE'):
        ret['reads'] = [
            str(RESULT_DIR / '01_trim_galore' / f'{wildcards.run}.trimmed.fastq.gz'),
        ]
    else:
        ret['reads'] = [
            str(RESULT_DIR / '01_trim_galore' / f'{wildcards.run}.read1.fastq.gz'),
            str(RESULT_DIR / '01_trim_galore' / f'{wildcards.run}.read2.fastq.gz'),
        ]
    
    return ret

c = config['bowtie']
rule bowtie:
    input: unpack(bowtie_input)
    output:
        # It automatically sorts the output bam file if its file name ends with '.sorted.bam',
        # e.g.
        # '{sample}.sorted.bam'
        RESULT_DIR / '02_bowtie' / '{run}.sorted.bam'
    params:
        # Additional parameters go here.
        extra = c['extra'],
        # Seed for random number generator.
        seed = c['seed'],
        # Skip the first <int> reads/pairs in the input.
        skip = c['skip'],
        # Stop after first <int> reads/pairs, excluding skipped reads.
        qupto = c['qupto'],
        # Trim <int> bases from 5' (left) end of reads
        trim5 = c['trim5'],
        # Trim <int> bases from 3' (right) end of reads
        trim3 = c['trim3'],
        # Input quals are Phred+33 (default)
        phred33_quals = c['phred33_quals'],
        # Input quals are Phred+64 (same as --solexa1.3-quals)
        phred64_quals = c['phred64_quals'],
        # Input quals are from GA Pipeline ver. < 1.3
        solexa_quals = c['solexa_quals'],
        # Input quals are from GA Pipeline ver. >= 1.3
        solexa13_quals = c['solexa13_quals'],
        # Qualities are given as space-separated integers (not ASCII)
        integer_quals = c['integer_quals'],
        # Force usage of a 'arge' index, even if a small one is present.
        large_index = c['large_index'],
        # Report end-to-end hits with <=v mismatches; ignore qualities.
        v = c['v'],
        # Max mismatches in seed (can be 0-3)
        # Deafult: 2
        seedmms = c['seedmms'],
        # Max sum of mismatch quals across alignment for -n
        # Default: 70
        maqerr = c['maqerr'],
        # Seed length for -n
        # Default: 28
        seedlen = c['seedlen'],
        # Disable Maq-like quality rounding for -n (nearest 10 <= 30)
        nomaqround = c['nomaqround'],
        # Minimum insert size for paired-end alignment
        # Default: 0
        minins = c['minins'],
        # Maximum insert size for paired-end alignment
        # Default: 250
        maxins = c['maxins'],
        # -1, -2 mates align fw/rev, rev/fw, fw/fw.
        # Default: fr
        fr = c['fr'],
        rf = c['rf'],
        ff = c['ff'],
        # Do not align to forward/reverse-complement reference strand.
        nofw = c['nofw'],
        norc = c['norc'],
        # Max # backtracks for -n 2/3
        # Default: 125, 800 for --best
        maxbts = c['maxbts'],
        # Max # attempts to find mate for anchor hit.
        # Default: 100
        pairtries = c['pairtries'],
        # Try hard to find valid alignments, at the expense of speed.
        tryhard = c['tryhard'],
        # Max megabytes of RAM for best-first search frames.
        # Default: 64
        chunkmbs = c['chunkmbs'],
        # # of reads to read from input file at once.
        # Default: 16
        reads_per_batch = c['reads_per_batch'],
        # Report up to <int> good alignments per read.
        # Default: 1
        k = c['k'],
        # Report all alignment sper read (much slower than low -k)
        all = c['all'],
        # Suppress all alignments if > <int> exist.
        # Default: no limit
        m = c['m'],
        # Like -m, but reports 1 random hit (MAPQ=0); requires --best
        M = c['M'],
        # Hist guaranteed best stratu; ties broken by quality.
        best = c['best'],
        # Hits in sub-optimal strata arent' reported (requires --best)
        strata = c['strata'],
        # Print wall-clock time taken by search phases.
        time = c['time'],
        # Leftmost ref offset = <int> in bowtie output.
        # Default: 0
        offbase = c['offbase'],
        # Print nothing but the alignments.
        quiet = c['quiet'],
        # Refer to ref. seqs by 0-based index rather than name.
        refidx = c['refidx'],
        # Write aligned reads/pairs to file(s) <fname>
        al = c['al'],
        # Write unaligned reads/pairs to file(s) <fname>
        un = c['un'],
        # Suppress SAM records for unaligned reads.
        no_unal = c['no_unal'],
        # Write reads/pairs over -m limit to file(s) <fname>
        max = c['max'],
        # Suppresses given columns (comma-delim'ed) in default output.
        suppress = c['suppress'],
        # Write entire ref name.
        # Default: Only up to 1st space
        fullref = c['fullref'],
        # Phred penalty for SNP when decoding colorspace.
        # Default: 30
        snpphred = c['snpphred'],
        # Approximate fraction of SNP bases (e.g. 0.001); sets --snpphred
        snpfrac = c['snpfrac'],
        # Print aligned colorspace seqs as colors, not decoded bases.
        col_cseq = c['col_cseq'],
        # Print original colorspace quals, not decoded quals.
        col_cqual = c['col_cqual'],
        # Keep nucleotides at extreme ends of decoded alignment.
        col_keepends = c['col_keepends'],
        # Write hits in SAM format.
        # Default: True (for this wrapper)
        sam = c['sam'],
        # Default mapping quality (MAPQ) to print for SAM alignments.
        mapq = c['mapq'],
        # Suppress header lines (starting with @) for SAM output.
        sam_nohead = c['sam_nohead'],
        # Suppress @SQ header lines for SAM output.
        sam_nosq = c['sam_nosq'],
        # ADD <text> (usually "lab=value") to @RG line of SAM header.
        sam_rg = c['sam_rg'],
        # Override offrate of index; must be >= index's offrate.
        offrate = c['offrate'],
        # Use memory-mapped I/O for index; many 'bowtie's can share.
        mm = c['mm'],
        # Use shared mem for index; many 'bowtie's can share.
        shmem = c['shmem'],
        # Discard mapped reads having mapping quality (MAPQ) below this value.
        # NOTE: This will be done after the alignment, with `samtools view -bS -q <int>` command.
        mapq_cutoff = c['mapq_cutoff'],
    threads: config['threads']['bowtie']
    benchmark:
        repeat('benchmarks/bowtie/{run}.tsv', 1)
    log: 'logs/bowtie/{run}.log'
    wrapper:
        'http://dohlee-bio.info:9193/bowtie'


