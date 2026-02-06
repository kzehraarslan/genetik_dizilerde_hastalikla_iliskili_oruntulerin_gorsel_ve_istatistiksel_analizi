import pandas as pd
import re

# ===============================
# 1️⃣ Veriyi oku
# ===============================
df = pd.read_csv("data/processed/variants/proteins_with_changes.csv")
df = df[df["Variant"] != "BRCA1_reference_CDS"]

# ===============================
# 2️⃣ Domain tanımları
# ===============================
DOMAINS = {
    "RING": (24, 65),
    "BRCT_1": (1642, 1736),
    "BRCT_2": (1760, 1855)
}

# ===============================
# 3️⃣ Amino asit grupları
# ===============================
AA_GROUPS = {
    "hydrophobic": set("AILMFWV"),
    "polar": set("STNQY"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "special": set("CGP")
}

def aa_group(aa):
    for group, aas in AA_GROUPS.items():
        if aa in aas:
            return group
    return "unknown"

# ===============================
# 4️⃣ Tek varyant için risk hesabı
# ===============================
def calculate_risk(change_str):
    if change_str == "-" or change_str.lower() == "no change":
        return 0, 0, 0

    changes = change_str.split(",")
    domain_hits = 0
    radical_changes = 0

    for ch in changes:
        match = re.match(r"([A-Z])(\d+)([A-Z])", ch.strip())
        if not match:
            continue

        ref_aa, pos, var_aa = match.groups()
        pos = int(pos)

        # Domain kontrolü
        for start, end in DOMAINS.values():
            if start <= pos <= end:
                domain_hits += 1

        # Amino asit tipi değişimi
        if aa_group(ref_aa) != aa_group(var_aa):
            radical_changes += 1

    return len(changes), domain_hits, radical_changes

# ===============================
# 5️⃣ Tüm varyantlar için uygula
# ===============================
results = df["Protein_changes"].apply(calculate_risk)
df[["Total_Changes", "Domain_Hits", "Radical_Changes"]] = pd.DataFrame(results.tolist(), index=df.index)

# ===============================
# 6️⃣ Final Risk Skoru
# ===============================
df["Risk_Score"] = (
    df["Total_Changes"]
    + 2 * df["Domain_Hits"]
    + df["Radical_Changes"]
)

# ===============================
# 7️⃣ Kaydet
# ===============================
output = "data/processed/BRCA1_variant_risk_scores.csv"
df.sort_values("Risk_Score", ascending=False).to_csv(output, index=False)

print("Final varyant risk analizi tamamlandı:")
print(output)