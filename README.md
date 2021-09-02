A Nextflow implementation of the STRling data flow:

https://github.com/quinlan-lab/STRling

STRling: a method to detect novel (and reference) STR
expandsions from short-read data.

# Usage

To run joint calls across more than one sample on version 1.0 of the workflow:

```
nextflow run quinlan-lab/STRling-nf -r 1.0 \
    --joint --crams 'preprocessing/*.bam' --reference GRCh38.fasta \
```

The workflow expects the alignments to be indexed (contain a .bai/.crai)
in the same directory as the .bam/.cram. It also expects the reference
.fasta to be indexed (.fai) and be in the same directory as the .fasta.

Single sample, or non-joint calls, omit `--joint`:

```
nextflow run quinlan-lab/STRling-nf -r 1.0 \
    --crams 'preprocessing/*.bam' --reference GRCh38.fasta \
```

Complete `--help` docs available using:

```
nextflow run quinlan-lab/STRling-nf -r 1.0 --help
```

# Output

Results of the workflow are written to `<outdir>/outliers` and
consist of files defined in the STRling docs:

https://strling.readthedocs.io/en/latest/outputs.html
