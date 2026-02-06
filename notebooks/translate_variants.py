from Bio import SeqIO
import csv
import os

# Klasör kontrolü
os.makedirs("data/processed/variants", exist_ok=True)

# Girdi FASTA (combined varyant CDS dizileri)
input_fasta = "data/processed/variants/BRCA1_variants_combined.fasta"

# Çıktılar
output_fasta = "data/processed/variants/BRCA1_variants_proteins.fasta"
output_csv = "data/processed/variants/BRCA1_variants_proteins.csv"

# FASTA ve CSV aç
with open(output_fasta, "w") as fasta_out, open(output_csv, "w", newline="") as csv_out:
    writer = csv.writer(csv_out)
    writer.writerow(["Variant", "Protein_sequence"])
    
    # Her varyant için
    for record in SeqIO.parse(input_fasta, "fasta"):
        protein_seq = record.seq.translate(to_stop=True)  # Stop kodonu ile dur
        # FASTA olarak yaz
        fasta_out.write(f">{record.id}\n{protein_seq}\n")
        # CSV olarak yaz
        writer.writerow([record.id, protein_seq])

print("Tüm varyantların protein dizileri oluşturuldu ve kaydedildi.")