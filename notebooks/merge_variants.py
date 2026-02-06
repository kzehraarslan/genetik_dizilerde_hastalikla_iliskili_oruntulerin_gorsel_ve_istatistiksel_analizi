from Bio import SeqIO
import os

# Klasör yolları
variants_folder = "../data/processed/variants/"
reference_cds_file = "../data/processed/BRCA1_reference_CDS.fasta"
output_file = "../data/processed/BRCA1_variants_combined.fasta"

# 1️⃣ Tüm varyant FASTA dosyalarını oku
variant_files = [f for f in os.listdir(variants_folder) if f.endswith(".fasta")]

# 2️⃣ Tek bir listeye FASTA kayıtlarını topla
records = []

# Referans CDS'i ekle
records.append(SeqIO.read(reference_cds_file, "fasta"))

# Varyantları ekle
for vf in variant_files:
    var_record = SeqIO.read(os.path.join(variants_folder, vf), "fasta")
    records.append(var_record)

# 3️⃣ Tek FASTA dosyası olarak kaydet
SeqIO.write(records, output_file, "fasta")

print(f"Tüm varyantlar ve referans CDS tek dosyada birleştirildi: {output_file}")
print(f"Toplam {len(records)} kayıt eklendi (1 referans + {len(records)-1} varyant).")