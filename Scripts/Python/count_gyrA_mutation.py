# Valencia, 27-04-2026
# This script was written by Matilde Ibba
# This script calcultes the exact number of mutations in an
# alignment (.faa) previously seen in AliView

from Bio import SeqIO
from collections import Counter

alignment = "gyrA_aligned.faa"

# Keyword of the NCBI reference
ref_keyword = "WP_002515815.1"

# Amino acid position in the reference without counting gaps
target_ref_pos = 72

records = list(SeqIO.parse(alignment, "fasta"))

# --------------------------------------------------
# Find reference sequence
# --------------------------------------------------
ref = None
for r in records:
    if ref_keyword in r.id or ref_keyword in r.description:
        ref = r
        break
        
if ref is None:
    raise ValueError("Reference not found")

# --------------------------------------------------
# Find the column corresponding to the aa position
# --------------------------------------------------
ungapped_pos = 0
target_col = None

for i, aa in enumerate(str(ref.seq)):
    if aa != "-":
        ungapped_pos += 1
    if ungapped_pos == target_ref_pos:
        target_col = i
        break

if target_col is None:
    raise ValueError("Target position not found in the alignment")

print("Column of alignment:", target_col +1)
print("Amino acid of reference:", ref.seq[target_col])

# --------------------------------------------------
# Count amino acid in that column
# --------------------------------------------------
counts = Counter()

for r in records:
    if r.id == ref.id:
        continue
    aa = str(r.seq)[target_col]
    counts[aa] += 1
    
print("\nAmino acid counts at position:")
print(counts)

# --------------------------------------------------
# Listing samples with mutation (D > G)
# --------------------------------------------------
print("\nSamples with D > G:")
for r in records:
    if r.id == ref.id:
        continue
    aa = str(r.seq)[target_col]
    if aa == "G":
        print(r.id)
        
# --------------------------------------------------
# Calculte frecuency of the mutation
# --------------------------------------------------
total = sum(counts.values())
mut = counts.get("G", 0)

print(f"\nD > G frecuency: {mut}/{total} ({(mut/total)*100:.2f}%)")