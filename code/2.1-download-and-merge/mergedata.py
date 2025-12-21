import os
import re

import pandas as pd

# --- CONFIGURATION ---
INPUT_FOLDER = "data"
OUTPUT_FILE = "merged_data.xlsx"

# Map 3-letter codes to Sheet Names
SPECIES_MAP = {
    "hsa": "Human",
    "mmu": "House mouse",
    "dme": "Fruit fly",
    "cel": "Roundworm",
}

# The strict column list
FINAL_COLUMNS = [
    "MirGeneDB ID",
    "MiRBase ID",
    "Family",
    "Seed",
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Precursor sequence",
    "Mature sequence",
    "Star sequence",
    "Mature location",
    "Precursor length",
    "Mature length",
    "Star length",
]


def clean_csv_id(raw_id):
    s = str(raw_id).strip()
    if s.endswith(" V"):
        s = s[:-2]
    return s


def clean_fasta_header(header):
    # Remove '>', '*', and version suffixes like '_v1'
    raw = header.strip().replace(">", "").replace("*", "")
    clean = re.sub(r"[-_]v\d+", "", raw)

    if clean.endswith("_5p"):
        return clean[:-3], "5p"
    elif clean.endswith("_3p"):
        return clean[:-3], "3p"
    elif clean.endswith("_pre"):
        return clean[:-4], None
    else:
        return clean, None


def parse_fasta(file_path):
    if not os.path.exists(file_path):
        print(f"  Warning: Missing {file_path}")
        return {}

    seqs = {}
    current_key = None
    current_seq = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_key:
                    seqs[current_key] = "".join(current_seq)
                base_id, arm = clean_fasta_header(line)
                current_key = f"{base_id}_{arm}" if arm else base_id
                current_seq = []
            else:
                current_seq.append(line)
        if current_key:
            seqs[current_key] = "".join(current_seq)
    return seqs


def resolve_mature(mir_id, seed_seq, fasta_dict):
    key_5p, key_3p = f"{mir_id}_5p", f"{mir_id}_3p"
    seq_5p, seq_3p = fasta_dict.get(key_5p), fasta_dict.get(key_3p)

    if seq_5p and not seq_3p:
        return seq_5p, "5p"
    if seq_3p and not seq_5p:
        return seq_3p, "3p"

    if seq_5p and seq_3p:
        seed = str(seed_seq).strip() if pd.notna(seed_seq) else ""
        if seed in seq_5p and seed not in seq_3p:
            return seq_5p, "5p"
        if seed in seq_3p and seed not in seq_5p:
            return seq_3p, "3p"
        return seq_5p, "5p"  # Fallback

    return "N/A", "N/A"


# --- MAIN EXECUTION ---
print(f"Starting merge... Saving to '{OUTPUT_FILE}'")

# Use xlsxwriter for formatting
with pd.ExcelWriter(OUTPUT_FILE, engine="xlsxwriter") as writer:
    for code, sheet_name in SPECIES_MAP.items():
        print(f"Processing {sheet_name} ({code})...")

        # 1. Load CSV
        csv_path = os.path.join(INPUT_FOLDER, f"{code}.csv")
        if not os.path.exists(csv_path):
            print(f"  Skipping {sheet_name} (CSV not found)")
            continue

        try:
            # Try header 1 (common), fallback to 0
            df = pd.read_csv(csv_path, header=1)
            id_col = next((c for c in df.columns if "MirGeneDB ID" in str(c)), None)
            if not id_col:
                df = pd.read_csv(csv_path, header=0)
                id_col = next((c for c in df.columns if "MirGeneDB ID" in str(c)), None)

            if not id_col:
                print("  ID column missing")
                continue
            df.rename(columns={id_col: "MirGeneDB ID"}, inplace=True)
        except Exception as e:
            print(f"  Error reading CSV: {e}")
            continue

        # 2. Load FASTAs
        pre_dict = parse_fasta(os.path.join(INPUT_FOLDER, f"{code}_pre.fas"))
        mat_dict = parse_fasta(os.path.join(INPUT_FOLDER, f"{code}_mature.fas"))
        star_dict = parse_fasta(os.path.join(INPUT_FOLDER, f"{code}_star.fas"))

        # 3. Build Data Lists
        pre_seqs, mat_seqs, star_seqs, mat_locs = [], [], [], []

        for index, row in df.iterrows():
            mir_id = clean_csv_id(row["MirGeneDB ID"])
            seed = row.get("Seed", "")

            # Precursor
            pre_seqs.append(pre_dict.get(mir_id, "N/A"))

            # Mature
            m_seq, m_loc = resolve_mature(mir_id, seed, mat_dict)
            mat_seqs.append(m_seq)
            mat_locs.append(m_loc)

            # Star Logic (Rescue Integrated)
            target_star = "3p" if m_loc == "5p" else "5p"
            star_key = f"{mir_id}_{target_star}"

            # --- THE FIX: Check Mature dict first for the 'rejected' candidate ---
            s_seq = mat_dict.get(star_key)

            # If not found in Mature dict, check Star dict
            if not s_seq:
                s_seq = star_dict.get(star_key, "")

            # Star Fallback (Original logic)
            if not s_seq and f"{mir_id}_{m_loc}" in star_dict:
                s_seq = star_dict[f"{mir_id}_{m_loc}"]

            star_seqs.append(s_seq)

        # 4. Assign & Calc Lengths
        df["Precursor sequence"] = pre_seqs
        df["Mature sequence"] = mat_seqs
        df["Star sequence"] = star_seqs
        df["Mature location"] = mat_locs

        df["Precursor length"] = df["Precursor sequence"].apply(
            lambda x: len(x) if x != "N/A" and x is not None else 0
        )
        df["Mature length"] = df["Mature sequence"].apply(
            lambda x: len(x) if x != "N/A" and x is not None else 0
        )
        df["Star length"] = df["Star sequence"].apply(
            lambda x: len(x) if x != "N/A" and x is not None else 0
        )

        # 5. Filter Columns
        available_cols = [c for c in FINAL_COLUMNS if c in df.columns]
        final_df = df[available_cols]

        # 6. Save & Format
        final_df.to_excel(writer, sheet_name=sheet_name, index=False)

        worksheet = writer.sheets[sheet_name]
        for idx, col in enumerate(final_df.columns):
            # Calculate width (with safety check for non-string data)
            col_len = len(str(col))
            data_len = (
                final_df[col].astype(str).map(len).max()
                if not final_df[col].empty
                else 0
            )
            max_len = max(data_len, col_len)
            worksheet.set_column(idx, idx, min(max_len + 2, 50))

print(f"Done! Saved to {os.path.abspath(OUTPUT_FILE)}")
