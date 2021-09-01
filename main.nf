nextflow.enable.dsl=2

params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------

    STRling-nf
    ==========

    Documentation and issues can be found at:
    https://github.com/brwnj/...

    Required arguments:
    -------------------
    --crams               Aligned sequences in .bam and/or .cram format.
                          Indexes (.bai/.crai) must be present.
    --reference           Reference FASTA. Index (.fai) must exist in same
                          directory.

    Optional:
    -------------------
    --joint               Perform joint calling and find outliers.
                          Default: false
    --proportion_repeat   Proportion of read that is repetitive to be considered as STR.
                          Default: 0.8

    // extract and call
    --min_mapq            Minimum mapping quality (does not apply to STR reads).
                          Default: 40
    --min_support         Minimum number of supporting reads for a locus to be reported
                          Default: 5
    --min_clip            Minimum number of supporting clipped reads for each side of a locus
                          Default: 0
    --min_clip_total      Minimum total number of supporting clipped reads for a locus
                          Default: 0

    // merge
    --window              Number of bp within which to search for reads supporting the
                          other side of a bound. Estimated from the insert size distribution
                          by default.
                          Default: -1

    // outliers
    --control             Input file for median and standard deviation estimates at
                          each locus from a set of control samples. This file can be
                          produced by this script using the emit option. If this
                          option is not set, all samples in the current batch will
                          be used as controls by default.
                          Default: false
    --slop                Merge loci that are within this many bp of each other and
                          have the same repeat unit.
                          Default: 50
    --min_clips           In the individual sample files, only report loci with at
                          least many soft-clipped reads in that sample.
                          Default: 0
    --min_size            In the individual sample files, only report loci with at
                          least this allele2_est size in that sample.
                          Default: 0

    -----------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

params.crams = false
params.reference = false

if(!params.crams) {
    exit 1, "--crams argument like '/path/to/*.cram' is required"
}
if(!params.reference) {
    exit 1, "--reference argument is required"
}

process strling_index {
    input:
    path(reference)
    path(fai)

    output:
    path("${reference}.str"), emit: str

    script:
    """
    strling index $reference -g ${reference}.str
    """
}

process strling_extract {
    input:
    tuple val(sample), path(cram), path(crai)
    path(reference)
    path(fai)
    path(str)
    val(proportion_repeat)
    val(min_mapq)

    output:
    tuple val("${sample}"), path(cram), path(crai), path("${sample}.bin"), emit: bin
    path("${sample}.bin"), emit: bin_only

    script:
    """
    strling extract -f $reference -g $str -p $proportion_repeat -q $min_mapq $cram ${sample}.bin
    """
}

process strling_merge {
    input:
    path(bin)
    path(reference)
    path(fai)
    val(window)
    val(min_support)
    val(min_clip)
    val(min_clip_total)
    val(min_mapq)

    output:
    path("joint-bounds.txt"), emit: bounds

    script:
    """
    strling merge -f $reference -w $window -m $min_support -c $min_clip \
        -t $min_clip_total -q $min_mapq -o joint $bin
    """
}

process strling_call {
    input:
    tuple val(sample), path(cram), path(crai), path(bin)
    path(reference)
    path(fai)
    path(bounds)
    val(min_mapq)
    val(min_support)
    val(min_clip)
    val(min_clip_total)

    output:
    path("${sample}-bounds.txt"), emit: bounds
    path("${sample}-genotype.txt"), emit: genotypes
    path("${sample}-unplaced.txt"), emit: unplaced

    script:
    b = bounds ? "-b $bounds" : ""
    """
    strling call -o $sample $b -m $min_support -c $min_clip -t $min_clip_total \
        -q $min_mapq -v -f $reference $cram $bin
    """
}

process strling_outliers {
    publishDir "${params.outdir}/outliers", mode: 'symlink'

    input:
    path(genotypes)
    path(unplaced)
    path(control)
    val(slop)
    val(min_clips)
    val(min_size)

    output:
    path("*STRs.tsv")
    path("control.tsv")
    path("depths.tsv")
    path("unplaced.tsv")

    script:
    c = control ? "--control $control" : ""
    """
    strling-outliers.py --genotypes $genotypes --unplaced $unplaced \
        --emit control.tsv --slop $slop --min_clips $min_clips \
        --min_size $min_size $c
    """
}

workflow {

    crams = Channel.fromPath(params.crams, checkIfExists: true)
        .map { file -> tuple(file.simpleName, file, file + ("${file}".endsWith('.cram') ? '.crai' : '.bai')) }
    fai = "${params.reference}.fai"

    strling_index(params.reference, fai)
    strling_extract(
        crams,
        params.reference,
        fai,
        strling_index.out.str,
        params.proportion_repeat,
        params.min_mapq
    )
    if (params.joint) {
        strling_merge(
            strling_extract.out.bin_only.collect(),
            params.reference,
            fai,
            params.window,
            params.min_support,
            params.min_clip,
            params.min_clip_total,
            params.min_mapq
        )
        strling_call(
            strling_extract.out.bin,
            params.reference,
            fai,
            strling_merge.out.bounds,
            params.min_mapq,
            params.min_support,
            params.min_clip,
            params.min_clip_total
        )
    } else {
        strling_call(
            strling_extract.out.bin,
            params.reference,
            fai,
            [],
            params.min_mapq,
            params.min_support,
            params.min_clip,
            params.min_clip_total
        )
    }
    strling_outliers(
        strling_call.out.genotypes.collect(),
        strling_call.out.unplaced.collect(),
        params.control,
        params.slop,
        params.min_clips,
        params.min_size
    )
}
