from Bio import SeqIO
import os
import re

# Referans CDS dosyasını oku
record = SeqIO.read("../data/processed/BRCA1_reference_CDS.fasta", "fasta")
cds_sequence = str(record.seq)

# ClinVar varyant dosyasını oku
with open("../data/raw/clinvar.txt") as f:
    variants = [line.strip() for line in f if line.strip()]

# İşlenecek diziler için klasör oluştur
os.makedirs("../data/processed/variants", exist_ok=True)

for var in variants:
    # Örn: "c.*4056C>A"
    parts = var.split("*")
    if len(parts) != 2:
        print(f"Atlandı: {var}")
        continue

    pos_change = parts[1]  # örn: 4056C>A
    match = re.match(r"(\d+)([ACGT])>([ACGT])", pos_change)
    if not match:
        print(f"Format hatası: {var}")
        continue

    pos, ref_nt, alt_nt = int(match.group(1)), match.group(2), match.group(3)
    cds_index = pos - 1  # Python 0-index

    # CDS sınır kontrolü
    if cds_index < 0 or cds_index >= len(cds_sequence):
        print(f"ATLANDI: CDS dışında: {var}")
        continue

    # Referans ile uyuşma kontrolü
    if cds_sequence[cds_index] != ref_nt:
        print(f"UYARI: Referans ile uyuşmayan nükleotid: {var} ({cds_sequence[cds_index]} != {ref_nt})")

    # Varyantı uygula
    var_seq = cds_sequence[:cds_index] + alt_nt + cds_sequence[cds_index+1:]

    # FASTA olarak kaydet
    filename = f"../data/processed/variants/{var.replace('*','_')}.fasta"
    with open(filename, "w") as f_out:
        f_out.write(f">{var}\n")
        f_out.write(var_seq)

print("Tüm varyantlar için CDS dizileri oluşturuldu.")