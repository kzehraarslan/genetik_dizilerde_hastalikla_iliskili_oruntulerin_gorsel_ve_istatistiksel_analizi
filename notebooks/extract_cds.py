from Bio import SeqIO

# FASTA dosyasını oku
record = SeqIO.read("data/raw/BRCA1_reference_mRNA.fasta", "fasta")

sequence = str(record.seq)

# CDS sınırları (NCBI'dan alındı.)
cds_start = 114 - 1  # Python 0-index
cds_end = 5705

cds_sequence = sequence[cds_start:cds_end]

# CDS'i yeni FASTA olarak kaydet
with open("data/processed/BRCA1_reference_CDS.fasta", "w") as f:
    f.write(">BRCA1_reference_CDS\n")
    f.write(cds_sequence)

print("CDS extraction tamamlandı.")
print("CDS uzunluğu:", len(cds_sequence))