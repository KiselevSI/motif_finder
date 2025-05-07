# motif_finder

## A C-based tool for detecting sequence motifs in bacteriophage genomes

`motif_finder` is a high-performance C program designed for the identification and analysis of specific DNA sequence motifs in bacteriophage genomes, with particular emphasis on finding Start Associated Sequences (SAS) and Extended Start Associated Sequences (ESAS) in mycobacteriophages.

## Overview

This tool was developed to analyze regulatory sequence features found in mycobacteriophage genomes such as those in Cluster K (TM4, ZoeJ, etc.). It efficiently scans genomic sequences to identify important regulatory motifs that control gene expression during phage life cycles.

## Features

- Fast and memory-efficient scanning of DNA sequences using optimized C code
- Support for IUPAC ambiguity codes for flexible motif definition
- Detection of SAS motifs (13bp sequences typically found upstream of gene start codons)
- Identification of ESAS motifs (paired inverted repeats separated by variable spacers)
- Calculation of Hamming distance to allow for imperfect matches
- Support for gapped motifs
- FASTA format support for sequence input
- Comprehensive output reporting with sequence positions and match statistics

## Compilation

The program can be easily compiled using the provided Makefile:

```
make
```

This will produce the executable `motif_finder`.

## Usage

```
./motif_finder -i input.fa -o output.tsv --left MOTIF1 --right MOTIF2 --gap-min N --gap-max M [options]
```

### Required arguments:
- `-i, --input`: Input FASTA file(s) containing sequences to search
- `-o, --output`: Output file (TSV format)
- `--left`: Left motif sequence to search for (supports IUPAC codes)
- `--right`: Right motif sequence to search for (supports IUPAC codes)
- `--gap-min`: Minimum gap size between motifs
- `--gap-max`: Maximum gap size between motifs

### Optional arguments:
- `--err-left`: Maximum number of mismatches allowed in left motif (default: 0)
- `--err-right`: Maximum number of mismatches allowed in right motif (default: 0)
- `-a, --after`: Number of bases to include after the match in output (default: 0)

## Example

To search for SAS-like sequences (consensus: GGGATAGGAGCCC) allowing up to 2 mismatches:

```
./motif_finder -i mycobacteriophage.fa -o sas_hits.tsv --left GGGATAGGAGCCC --right "" --gap-min 0 --gap-max 0 --err-left 2
```

To search for ESAS motifs with variable spacers:

```
./motif_finder -i phage_genome.fa -o esas_results.tsv --left TGTTGACNNNTCAACA --right TGTTGANNNGTCAACA --gap-min 4 --gap-max 13 --err-left 1 --err-right 1
```

## Applications

Beyond identifying SAS and ESAS sequences in mycobacteriophages, this tool can be useful for:

1. Detecting promoter regions in bacterial genomes
2. Finding transcription factor binding sites
3. Identifying conserved regulatory elements across multiple phage genomes
4. Characterizing the structure and distribution of repetitive elements
5. Analyzing evolutionary relationships between phage regulatory systems
6. Predicting gene expression patterns based on upstream sequence features
7. Supporting phage genome annotation by identifying regulatory motifs
8. Studying the mechanisms of lysogeny/lytic cycle control

## Scientific Background

Recent research on Cluster K mycobacteriophages has identified important regulatory sequences:

- **SAS (Start Associated Sequences)**: 13 bp asymmetric sequences positioned 3-8 bp upstream of translation start sites
- **ESAS (Extended Start Associated Sequences)**: Complex motifs with paired inverted repeats separated by variable spacers

These sequences appear to play crucial roles in phage gene regulation during both lytic and lysogenic growth cycles. Studies of phages like ZoeJ suggest that ESAS sites may promote transcription during lysogeny while being repressed during lytic growth.

## Performance

The program is optimized for speed and memory efficiency, implementing:
- Efficient bit-masking for IUPAC code handling
- Linear-time scanning algorithms
- Memory-conscious data structures

## Requirements

- C compiler (gcc recommended)
- Standard C library
- POSIX-compliant system

## License

[MIT](LICENSE)
