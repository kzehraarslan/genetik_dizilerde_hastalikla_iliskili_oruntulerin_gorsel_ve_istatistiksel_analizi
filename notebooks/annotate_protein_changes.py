import pandas as pd

# CSV dosyasını oku
input_csv = "../data/processed/BRCA1_variants_proteins.csv"
df = pd.read_csv(input_csv)

# Referans protein dizisini al
ref_seq = df.loc[df['Variant'] == "BRCA1_reference_CDS", "Protein_sequence"].values[0]

# Fonksiyon: değişiklikleri tespit et
def detect_changes(ref, var):
    changes = []
    for i, (r, v) in enumerate(zip(ref, var), start=1):  # i = pozisyon (1-index)
        if r != v:
            changes.append(f"{r}{i}{v}")  # Örn: K123N
    return ";".join(changes) if changes else "-"

# Her varyant için protein değişikliklerini hesapla
protein_changes = []
for idx, row in df.iterrows():
    if row['Variant'] == "BRCA1_reference_CDS":
        protein_changes.append("-")
    else:
        var_seq = row['Protein_sequence']
        protein_changes.append(detect_changes(ref_seq, var_seq))

# Yeni sütunu ekle
df['Protein_changes'] = protein_changes

# CSV olarak kaydet
output_csv = "../data/processed/BRCA1_variants_proteins_annotated.csv"
df.to_csv(output_csv, index=False)

print(f"Tüm protein değişiklikleri tespit edildi ve '{output_csv}' olarak kaydedildi.")