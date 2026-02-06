import pandas as pd
import matplotlib.pyplot as plt

# 1️⃣ Veri setini oku

input_csv = "data/processed/variants/proteins_with_changes.csv"
df = pd.read_csv(input_csv)

# Referansı çıkar
df = df[df["Variant"] != "BRCA1_reference_CDS"]


# 2️⃣ Değişim sayısını hesapla

def count_changes(change_str):
    if change_str == "-" or change_str.lower() == "no change":
        return 0
    return len(change_str.split(","))

df["Change_Count"] = df["Protein_changes"].apply(count_changes)


# 3️⃣ Sırala (etkisi yüksek olanlar üstte)

df_sorted = df.sort_values("Change_Count", ascending=False)


# 4️⃣ Grafik

plt.figure(figsize=(14, 6))
plt.bar(
    df_sorted["Variant"],
    df_sorted["Change_Count"]
)

plt.xticks(rotation=90, fontsize=7)
plt.xlabel("BRCA1 Varyantları")
plt.ylabel("Amino Asit Değişim Sayısı")
plt.title("BRCA1 Varyantlarının Protein Üzerindeki Etkisi")

plt.tight_layout()


# 5️⃣ Kaydet ve göster

plt.savefig("data/processed/BRCA1_variant_protein_impact.png", dpi=300)
plt.show()

print("Grafik oluşturuldu ve kaydedildi:")
print("data/processed/BRCA1_variant_protein_impact.png")