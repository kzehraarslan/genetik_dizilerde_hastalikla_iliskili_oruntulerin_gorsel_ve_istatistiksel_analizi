import pandas as pd

# proteins.csv dosyasını oku
df = pd.read_csv("data/processed/variants/BRCA1_variants_proteins.csv")

# Referans proteini al
ref_protein = df.loc[df['Variant'] == "BRCA1_reference_CDS", "Protein_sequence"].values[0]

# Yeni bir sütun oluştur: değişim açıklaması
variant_changes = []

for idx, row in df.iterrows():
    if row['Variant'] == "BRCA1_reference_CDS":
        variant_changes.append("-")  # referansın kendisi değişim yok
        continue
    
    var_protein = row['Protein_sequence']
    changes = []
    
    # Kısa proteinlerde karşılaştırmak için zip kullanıyoruz
    for pos, (r, v) in enumerate(zip(ref_protein, var_protein), start=1):
        if r != v:
            changes.append(f"{r}{pos}{v}")  # Örn: M1V
    
    if not changes:
        variant_changes.append("No change")
    else:
        variant_changes.append(", ".join(changes))

# Yeni sütunu ekle
df['Protein_changes'] = variant_changes

# Sonucu kaydet
df.to_csv("data/processed/proteins_with_changes.csv", index=False)

print("Protein değişim analizi tamamlandı. Sonuç: data/processed/proteins_with_changes.csv")