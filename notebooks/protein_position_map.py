import pandas as pd
import matplotlib.pyplot as plt
import re
from collections import Counter

# ===============================
# 1️⃣ CSV dosyasını oku
# ===============================
input_csv = "data/processed/variants/proteins_with_changes.csv"
df = pd.read_csv(input_csv)

# Referansı çıkar
df = df[df["Variant"] != "BRCA1_reference_CDS"]

# ===============================
# 2️⃣ Tüm amino asit pozisyonlarını topla
# ===============================
positions = []

for change_str in df["Protein_changes"]:
    if change_str == "-" or change_str.lower() == "no change":
        continue
    
    # Örnek: K123N → 123
    changes = change_str.split(",")
    for ch in changes:
        match = re.search(r"\d+", ch)
        if match:
            positions.append(int(match.group()))

# ===============================
# 3️⃣ Pozisyon frekanslarını hesapla
# ===============================
position_counts = Counter(positions)

# DataFrame'e çevir
pos_df = pd.DataFrame(
    position_counts.items(),
    columns=["Protein_Position", "Mutation_Count"]
)

pos_df = pos_df.sort_values("Protein_Position")

# ===============================
# 4️⃣ Grafik
# ===============================
plt.figure(figsize=(16, 6))
plt.bar(
    pos_df["Protein_Position"],
    pos_df["Mutation_Count"]
)

plt.xlabel("Protein Pozisyonu (Amino Asit)")
plt.ylabel("Varyant Sayısı")
plt.title("BRCA1 Proteininde Amino Asit Değişimlerinin Pozisyon Dağılımı")

plt.tight_layout()

# ===============================
# 5️⃣ Kaydet ve göster
# ===============================
plt.savefig(
    "data/processed/BRCA1_protein_position_hotspots.png",
    dpi=300
)

plt.show()

print("Protein pozisyon haritası oluşturuldu:")
print("data/processed/BRCA1_protein_position_hotspots.png")