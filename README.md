# Motif Finder in FASTA Genomes

This Python script searches for occurrences of a given DNA motif in a FASTA genome file and outputs the results in BED format.

## Features

* Search motif occurrences on the forward strand by default.
* Optionally search the reverse complement strand.
* Supports DNA motifs with IUPAC ambiguity codes (e.g. N, V, H, D, B, M, R, W, S, Y, K).
* Outputs BED file with motif locations and matched sequences.
* Supports output compression with `.gz` extension.
* Progress bar shows processing progress based on total bases for smooth tracking.
* Prints summary of total matches and elapsed runtime.
* Minimal dependencies (`biopython` required).

## Installation

```bash
# Python 3.x
pip install biopython tqdm
```

## Usage

```bash
python find_motif.py --fasta genome.fa --out motifs.bed --motif GATC
```

### Options

| Parameter   | Description                                                         |
| ----------- | ------------------------------------------------------------------- |
| `--fasta`   | Input genome FASTA file (required)                                  |
| `--out`     | Output BED file path (required). Use `.gz` for gzip output          |
| `--motif`   | Motif sequence with IUPAC codes (required). Examples: `GATC`, `NVH` |
| `--reverse` | If set, also searches reverse complement strand                     |

## IUPAC Ambiguity Codes Supported

| Code  | Meaning          | Bases Included |
| ----- | ---------------- | -------------- |
| A     | Adenine          | A              |
| C     | Cytosine         | C              |
| G     | Guanine          | G              |
| T (U) | Thymine (Uracil) | T              |
| N     | Any base         | A, C, G, T     |
| V     | Not T            | A, C, G        |
| H     | Not G            | A, C, T        |
| D     | Not C            | A, G, T        |
| B     | Not A            | C, G, T        |
| M     | Amino            | A, C           |
| R     | Purine           | A, G           |
| W     | Weak             | A, T           |
| S     | Strong           | C, G           |
| Y     | Pyrimidine       | C, T           |
| K     | Keto             | G, T           |

## Output BED Format

The output file contains tab-separated columns:

```
chrom  start  end  name  score  strand  matched_sequence
```

* `start` is 0-based coordinate.
* `end` is non-inclusive.
* `name` includes motif, chromosome, and match index.
* `score` is always 0.
* `strand` is `+` or `-`.
* `matched_sequence` is the actual matched motif sequence from the genome.

## Example

Search for motif "GATC" on both strands and save compressed output:

```bash
python find_motif.py --fasta dm6.fa --out dm6_GATC_motifs.bed.gz --motif GATC
```

Search for a motif with ambiguity codes:

```bash
python find_motif.py --fasta genome.fa --out motifs.bed --motif NVH --reverse
```

## License

MIT License

## Contact

For suggestions or issues, please submit a GitHub issue.

