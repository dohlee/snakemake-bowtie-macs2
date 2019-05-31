from pathlib import Path
from os.path import join

DATA_DIR = Path(config['data_dir'])
RESULT_DIR = Path(config['result_dir'])

# Configurations for bowtie_build
BB = config['bowtie_build']
rule bowtie_build:
    input:
        # Required input.
        # Reference genome fasta.
        # e.g.
        # reference = config['reference']['genome']['fasta']
        reference = config['reference']['fasta']
    output:
        # e.g.
        # reference_indices = config['reference']['genome']['bowtie_index_dir']
        index_dir = directory(config['bowtie_index_dir'])
    params:
        # Additional parameters go here.
        extra = BB['extra'],
        # Build a colorspace index.
        color = BB['color'],
        # Disable automatic -p/--bmax/--dcv memory-fitting.
        noauto = BB['noauto'],
        # Use packed strings internally; slower, uses less mem.
        packed = BB['packed'], 
        # Max bucket size for blockwise suffix-array builder.
        bmax = BB['bmax'],
        # Max bucket size as divisor of ref len.
        # Default: 4
        bmaxdivn = BB['bmaxdivn'],
        # diff-cover period for blockwise
        # Default: 1024
        dcv = BB['dcv'],
        # disable diff-cover (algorithm becomes quadratic)
        nodc = BB['nodc'],
        # Don't build .3/.4.ebwt (packed reference) portion.
        noref = BB['noref'],
        # Just build .3/.4.ebwt (packed reference) portion.
        justref = BB['justref'],
        # SA is sampled every 2^offRate BWT chars.
        # Default: 5
        offrate = BB['offrate'],
        # Number of chars consumed in initial lookup.
        # Default: 10
        ftabchars = BB['ftabchars'],
        # convert Ns in reference to As.
        ntoa = BB['ntoa'],
        # Seed for random number generator.
        seed = BB['seed'],
        # Verbose output (for debugging).
        quiet = BB['quiet'],
    threads: config['threads']['bowtie_build']
    benchmark:
        repeat('benchmarks/bowtie_build/%s.tsv' % config['reference']['name'], 1)
    log: 'logs/bowtie_build/%s.log' % config['reference']['name']
    wrapper:
        'http://dohlee-bio.info:9193/bowtie/build'

# Configurations for bowtie_single.
BS = config['bowtie_single']
rule bowtie_single:
    input:
        index_dir = config['bowtie_index_dir'],
        reads = [
            DATA_DIR / '{sample}.fastq.gz',
        ],
    output:
        # It automatically sorts the output bam file if its file name ends with '.sorted.bam',
        # e.g.
        # '{sample}.sorted.bam'
        bam = RESULT_DIR / '01_bowtie' / 'se' / '{sample}.sorted.bam'
    params:
        # Additional parameters go here.
        extra = BS['extra'],
        # Seed for random number generator.
        seed = BS['seed'],
        # Skip the first <int> reads/pairs in the input.
        skip = BS['skip'],
        # Stop after first <int> reads/pairs, excluding skipped reads.
        qupto = BS['qupto'],
        # Trim <int> bases from 5' (left) end of reads
        trim5 = BS['trim5'],
        # Trim <int> bases from 3' (right) end of reads
        trim3 = BS['trim3'],
        # Input quals are Phred+33 (default)
        phred33_quals = BS['phred33_quals'],
        # Input quals are Phred+64 (same as --solexa1.3-quals)
        phred64_quals = BS['phred64_quals'],
        # Input quals are from GA Pipeline ver. < 1.3
        solexa_quals = BS['solexa_quals'],
        # Input quals are from GA Pipeline ver. >= 1.3
        solexa13_quals = BS['solexa13_quals'],
        # Qualities are given as space-separated integers (not ASCII)
        integer_quals = BS['integer_quals'],
        # Force usage of a 'arge' index, even if a small one is present.
        large_index = BS['large_index'],
        # Report end-to-end hits with <=v mismatches; ignore qualities.
        v = BS['v'],
        # Max mismatches in seed (can be 0-3)
        # Deafult: 2
        seedmms = BS['seedmms'],
        # Max sum of mismatch quals across alignment for -n
        # Default: 70
        maqerr = BS['maqerr'],
        # Seed length for -n
        # Default: 28
        seedlen = BS['seedlen'],
        # Disable Maq-like quality rounding for -n (nearest 10 <= 30)
        nomaqround = BS['nomaqround'],
        # Minimum insert size for paired-end alignment
        # Default: 0
        minins = BS['minins'],
        # Maximum insert size for paired-end alignment
        # Default: 250
        maxins = BS['maxins'],
        # -1, -2 mates align fw/rev, rev/fw, fw/fw.
        # Default: fr
        fr = BS['fr'],
        rf = BS['rf'],
        ff = BS['ff'],
        # Do not align to forward/reverse-complement reference strand.
        nofw = BS['nofw'],
        norc = BS['norc'],
        # Max # backtracks for -n 2/3
        # Default: 125, 800 for --best
        maxbts = BS['maxbts'],
        # Max # attempts to find mate for anchor hit.
        # Default: 100
        pairtries = BS['pairtries'],
        # Try hard to find valid alignments, at the expense of speed.
        tryhard = BS['tryhard'],
        # Max megabytes of RAM for best-first search frames.
        # Default: 64
        chunkmbs = BS['chunkmbs'],
        # # of reads to read from input file at once.
        # Default: 16
        reads_per_batch = BS['reads_per_batch'],
        # Report up to <int> good alignments per read.
        # Default: 1
        k = BS['k'],
        # Report all alignment sper read (much slower than low -k)
        all = BS['all'],
        # Suppress all alignments if > <int> exist.
        # Default: no limit
        m = BS['m'],
        # Like -m, but reports 1 random hit (MAPQ=0); requires --best
        M = BS['M'],
        # Hist guaranteed best stratu; ties broken by quality.
        best = BS['best'],
        # Hits in sub-optimal strata arent' reported (requires --best)
        strata = BS['strata'],
        # Print wall-clock time taken by search phases.
        time = BS['time'],
        # Leftmost ref offset = <int> in bowtie output.
        # Default: 0
        offbase = BS['offbase'],
        # Print nothing but the alignments.
        quiet = BS['quiet'],
        # Refer to ref. seqs by 0-based index rather than name.
        refidx = BS['refidx'],
        # Write aligned reads/pairs to file(s) <fname>
        al = BS['al'],
        # Write unaligned reads/pairs to file(s) <fname>
        un = BS['un'],
        # Suppress SAM records for unaligned reads.
        no_unal = BS['no_unal'],
        # Write reads/pairs over -m limit to file(s) <fname>
        max = BS['max'],
        # Suppresses given columns (comma-delim'ed) in default output.
        suppress = BS['suppress'],
        # Write entire ref name.
        # Default: Only up to 1st space
        fullref = BS['fullref'],
        # Phred penalty for SNP when decoding colorspace.
        # Default: 30
        snpphred = BS['snpphred'],
        # Approximate fraction of SNP bases (e.g. 0.001); sets --snpphred
        snpfrac = BS['snpfrac'],
        # Print aligned colorspace seqs as colors, not decoded bases.
        col_cseq = BS['col_cseq'],
        # Print original colorspace quals, not decoded quals.
        col_cqual = BS['col_cqual'],
        # Keep nucleotides at extreme ends of decoded alignment.
        col_keepends = BS['col_keepends'],
        # Write hits in SAM format.
        # Default: True (for this wrapper)
        sam = BS['sam'],
        # Default mapping quality (MAPQ) to print for SAM alignments.
        mapq = BS['mapq'],
        # Suppress header lines (starting with @) for SAM output.
        sam_nohead = BS['sam_nohead'],
        # Suppress @SQ header lines for SAM output.
        sam_nosq = BS['sam_nosq'],
        # ADD <text> (usually "lab=value") to @RG line of SAM header.
        sam_rg = BS['sam_rg'],
        # Override offrate of index; must be >= index's offrate.
        offrate = BS['offrate'],
        # Use memory-mapped I/O for index; many 'bowtie's can share.
        mm = BS['mm'],
        # Use shared mem for index; many 'bowtie's can share.
        shmem = BS['shmem'],
        # Discard mapped reads having mapping quality (MAPQ) below this value.
        # NOTE: This will be done after the alignment, with `samtools view -bS -q <int>` command.
        mapq_cutoff = BS['mapq_cutoff'],
    threads: config['threads']['bowtie_single']
    benchmark:
        repeat('benchmarks/bowtie/{sample}.tsv', 1)
    log: 'logs/bowtie/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/bowtie'

# Configurations for bowtie_paired.
BP = config['bowtie_paired']
rule bowtie_paired:
    input:
        index_dir = config['bowtie_index_dir'],
        reads = [
            DATA_DIR / '{sample}.read1.fastq.gz',
            DATA_DIR / '{sample}.read2.fastq.gz',
        ],
    output:
        # It automatically sorts the output bam file if its file name ends with '.sorted.bam',
        # e.g.
        # '{sample}.sorted.bam'
        bam = RESULT_DIR / '01_bowtie' / 'pe' / '{sample}.sorted.bam'
    params:
        # Additional parameters go here.
        extra = BP['extra'],
        # Seed for random number generator.
        seed = BP['seed'],
        # Skip the first <int> reads/pairs in the input.
        skip = BP['skip'],
        # Stop after first <int> reads/pairs, excluding skipped reads.
        qupto = BP['qupto'],
        # Trim <int> bases from 5' (left) end of reads
        trim5 = BP['trim5'],
        # Trim <int> bases from 3' (right) end of reads
        trim3 = BP['trim3'],
        # Input quals are Phred+33 (default)
        phred33_quals = BP['phred33_quals'],
        # Input quals are Phred+64 (same as --solexa1.3-quals)
        phred64_quals = BP['phred64_quals'],
        # Input quals are from GA Pipeline ver. < 1.3
        solexa_quals = BP['solexa_quals'],
        # Input quals are from GA Pipeline ver. >= 1.3
        solexa13_quals = BP['solexa13_quals'],
        # Qualities are given as space-separated integers (not ASCII)
        integer_quals = BP['integer_quals'],
        # Force usage of a 'arge' index, even if a small one is present.
        large_index = BP['large_index'],
        # Report end-to-end hits with <=v mismatches; ignore qualities.
        v = BP['v'],
        # Max mismatches in seed (can be 0-3)
        # Deafult: 2
        seedmms = BP['seedmms'],
        # Max sum of mismatch quals across alignment for -n
        # Default: 70
        maqerr = BP['maqerr'],
        # Seed length for -n
        # Default: 28
        seedlen = BP['seedlen'],
        # Disable Maq-like quality rounding for -n (nearest 10 <= 30)
        nomaqround = BP['nomaqround'],
        # Minimum insert size for paired-end alignment
        # Default: 0
        minins = BP['minins'],
        # Maximum insert size for paired-end alignment
        # Default: 250
        maxins = BP['maxins'],
        # -1, -2 mates align fw/rev, rev/fw, fw/fw.
        # Default: fr
        fr = BP['fr'],
        rf = BP['rf'],
        ff = BP['ff'],
        # Do not align to forward/reverse-complement reference strand.
        nofw = BP['nofw'],
        norc = BP['norc'],
        # Max # backtracks for -n 2/3
        # Default: 125, 800 for --best
        maxbts = BP['maxbts'],
        # Max # attempts to find mate for anchor hit.
        # Default: 100
        pairtries = BP['pairtries'],
        # Try hard to find valid alignments, at the expense of speed.
        tryhard = BP['tryhard'],
        # Max megabytes of RAM for best-first search frames.
        # Default: 64
        chunkmbs = BP['chunkmbs'],
        # # of reads to read from input file at once.
        # Default: 16
        reads_per_batch = BP['reads_per_batch'],
        # Report up to <int> good alignments per read.
        # Default: 1
        k = BP['k'],
        # Report all alignment sper read (much slower than low -k)
        all = BP['all'],
        # Suppress all alignments if > <int> exist.
        # Default: no limit
        m = BP['m'],
        # Like -m, but reports 1 random hit (MAPQ=0); requires --best
        M = BP['M'],
        # Hist guaranteed best stratu; ties broken by quality.
        best = BP['best'],
        # Hits in sub-optimal strata arent' reported (requires --best)
        strata = BP['strata'],
        # Print wall-clock time taken by search phases.
        time = BP['time'],
        # Leftmost ref offset = <int> in bowtie output.
        # Default: 0
        offbase = BP['offbase'],
        # Print nothing but the alignments.
        quiet = BP['quiet'],
        # Refer to ref. seqs by 0-based index rather than name.
        refidx = BP['refidx'],
        # Write aligned reads/pairs to file(s) <fname>
        al = BP['al'],
        # Write unaligned reads/pairs to file(s) <fname>
        un = BP['un'],
        # Suppress SAM records for unaligned reads.
        no_unal = BP['no_unal'],
        # Write reads/pairs over -m limit to file(s) <fname>
        max = BP['max'],
        # Suppresses given columns (comma-delim'ed) in default output.
        suppress = BP['suppress'],
        # Write entire ref name.
        # Default: Only up to 1st space
        fullref = BP['fullref'],
        # Phred penalty for SNP when decoding colorspace.
        # Default: 30
        snpphred = BP['snpphred'],
        # Approximate fraction of SNP bases (e.g. 0.001); sets --snpphred
        snpfrac = BP['snpfrac'],
        # Print aligned colorspace seqs as colors, not decoded bases.
        col_cseq = BP['col_cseq'],
        # Print original colorspace quals, not decoded quals.
        col_cqual = BP['col_cqual'],
        # Keep nucleotides at extreme ends of decoded alignment.
        col_keepends = BP['col_keepends'],
        # Write hits in SAM format.
        # Default: True (for this wrapper)
        sam = BP['sam'],
        # Default mapping quality (MAPQ) to print for SAM alignments.
        mapq = BP['mapq'],
        # Suppress header lines (starting with @) for SAM output.
        sam_nohead = BP['sam_nohead'],
        # Suppress @SQ header lines for SAM output.
        sam_nosq = BP['sam_nosq'],
        # ADD <text> (usually "lab=value") to @RG line of SAM header.
        sam_rg = BP['sam_rg'],
        # Override offrate of index; must be >= index's offrate.
        offrate = BP['offrate'],
        # Use memory-mapped I/O for index; many 'bowtie's can share.
        mm = BP['mm'],
        # Use shared mem for index; many 'bowtie's can share.
        shmem = BP['shmem'],
        # Discard mapped reads having mapping quality (MAPQ) below this value.
        # NOTE: This will be done after the alignment, with `samtools view -bS -q <int>` command.
        mapq_cutoff = BP['mapq_cutoff'],
    threads: config['threads']['bowtie_paired']
    benchmark:
        repeat('benchmarks/bowtie/{sample}.tsv', 1)
    log: 'logs/bowtie/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/bowtie'
