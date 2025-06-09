import argparse
from Bio import SeqIO
import re
import gzip
import time

# IUPAC nucleotide ambiguity codes mapping to regex classes
IUPAC_CODES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",  # treat U as T
    "N": "[ACGT]",
    "V": "[ACG]",
    "H": "[ACT]",
    "D": "[AGT]",
    "B": "[CGT]",
    "M": "[AC]",
    "R": "[AG]",
    "W": "[AT]",
    "S": "[CG]",
    "Y": "[CT]",
    "K": "[GT]"
}

def iupac_to_regex(motif):
    regex_parts = []
    for letter in motif.upper():
        if letter in IUPAC_CODES:
            regex_parts.append(IUPAC_CODES[letter])
        else:
            raise ValueError(f"Invalid IUPAC code in motif: {letter}")
    return "".join(regex_parts)

def reverse_complement(seq):
    complement = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(complement)[::-1]

def find_motifs(fasta_file, output_bed, motif, search_reverse=False):
    regex_motif = iupac_to_regex(motif)
    pattern = re.compile(regex_motif)
    total_matches = 0

    open_func = gzip.open if output_bed.endswith(".gz") else open

    start_time = time.time()
    with open_func(output_bed, "wt") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            chrom = record.id
            seq = str(record.seq).upper()

            # forward strand
            for match in pattern.finditer(seq):
                start = match.start()
                end = match.end()
                matched_seq = seq[start:end]
                total_matches += 1
                name = f"{motif}_{chrom}_{total_matches}"
                score = 0
                strand = "+"
                out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{matched_seq}\n")

            # reverse strand
            if search_reverse:
                rev_seq = reverse_complement(seq)
                for match in pattern.finditer(rev_seq):
                    start = len(seq) - match.end()
                    end = len(seq) - match.start()
                    matched_seq = seq[start:end]
                    total_matches += 1
                    name = f"{motif}_rev_{chrom}_{total_matches}"
                    score = 0
                    strand = "-"
                    out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{matched_seq}\n")

    elapsed = time.time() - start_time
    print(f"Finished searching motif '{motif}' in {fasta_file}")
    print(f"Total matches found: {total_matches}")
    print(f"Output saved to: {output_bed}")
    print(f"Elapsed time: {elapsed:.2f} seconds")

def main():
    parser = argparse.ArgumentParser(description="Find motif occurrences in a FASTA genome and output as BED.")
    parser.add_argument("--fasta", required=True, help="Input genome FASTA file")
    parser.add_argument("--out", required=True, help="Output BED file (supports .gz compression)")
    parser.add_argument("--motif", required=True, help="Motif sequence (IUPAC codes allowed, e.g. GATC, N, V, H)")
    parser.add_argument("--reverse", action="store_true", help="Also search motif on reverse strand")

    args = parser.parse_args()
    find_motifs(args.fasta, args.out, args.motif, args.reverse)

if __name__ == "__main__":
    main()
