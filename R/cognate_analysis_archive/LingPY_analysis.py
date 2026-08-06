import pandas as pd

# Define orthography-to-IPA mapping
def convert_to_ipa(word):
    if pd.isna(word) or word == "":
        return ""
    return (word
        .replace("ng", "ŋ")
        .replace("7", "ʔ")
        .replace("N", "ŋ")
        .replace("3", "ə")
        .replace("sh", "ʃ")
        .replace("ch", "tʃ")
        .replace("j", "dʒ")
        .replace("ñ", "ɲ")
        .replace("y", "j")
        .replace("x", "x") 
        .replace("q", "ʔ")
        .replace("%", "")
        .replace("~", "")
    )

# Load your CSV
df = pd.read_csv(r"C:\Users\ramma\Box\Ram_Ximena_Nicole\Indp Research Phillipine Languages\Cognate Analysis\WORDLISTdf_forlingpy_unflitered.csv")

# Melt to long format
df_long = df.melt(id_vars=["concept"], var_name="Language", value_name="form")

# Replace NA with blank strings
df_long["form"] = df_long["form"].fillna("")

# Apply IPA conversion
df_long["form"] = df_long["form"].apply(convert_to_ipa)

# Add unique ID column
df_long.insert(0, "ID", range(1, len(df_long) + 1))

# Rename 'concept' to 'Concept' for LingPy compatibility
df_long = df_long.rename(columns={"concept": "Concept"})

# Save as LingPy-ready TSV
df_long.to_csv("lingpy_input_unfiltered.tsv", sep="\t", index=False)

print("LingPy-ready file with IPA conversion saved as 'lingpy_input_unfiltered.tsv'")