# Usage

To run joint calls across more than one sample on version 1.0 of the workflow:

```
nextflow run quinlan-lab/STRling-nf -r 1.0 \
    --joint --crams 'preprocessing/*.bam' --reference GRCh38.fasta \
```

Single sample, or non-joint calls, omit `--joint`:

```
nextflow run quinlan-lab/STRling-nf -r 1.0 \
    --crams 'preprocessing/*.bam' --reference GRCh38.fasta \
```

Complete `--help` docs available using:

```
nextflow run quinlan-lab/STRling-nf -r 1.0 --help
```
