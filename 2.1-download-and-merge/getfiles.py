import os
import urllib.request
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import requests

url = "https://mirgenedb.org/download"
base_url = "https://mirgenedb.org"

# Load the table
df = pd.read_html(url, extract_links="body")[0]

# Wanted dict: common name -> latin name
wanted = {
    "Human": "Homo sapiens",
    "House mouse": "Mus musculus",
    "Fruit fly": "Drosophila melanogaster",
    "Roundworm": "Caenorhabditis elegans",
}

# Only match on latin names
wanted_vals = set(wanted.values())

cols_needed = [
    "Precursor sequences",
    "Mature sequences",
    "Star sequences",
    "Genomic coordinates",
]

# Get relevant rows
rows = []

for _, row in df.iterrows():
    first_col = row.iloc[0]  # may be ('Human (Homo sapiens)', None)

    first_col = str(first_col)

    for latin in wanted_vals:
        if latin in first_col:
            row.iloc[0] = row.iloc[0][0]
            rows.append(row)
            break

sub = pd.DataFrame(rows, columns=df.columns)

alllinks = {}

# Get relevant links
for _, row in sub.iterrows():
    links = {}
    species_name = row.iloc[0]

    for col in cols_needed:
        cell = row[col]  # tuple: (text, href)
        href = cell[1] if isinstance(cell, tuple) else None

        if href:
            # make absolute URL
            if href.startswith("http"):
                links[col] = href
            else:
                links[col] = base_url + href
        else:
            links[col] = None
    alllinks[species_name] = links

# Create folder if needed
out_dir = Path("data")
os.makedirs(out_dir, exist_ok=True)

nametocode = {
    "Human": "hsa",
    "House mouse": "mmu",
    "Fruit fly": "dme",
    "Roundworm": "cel",
}

typetosuffix = {
    "Precursor sequences": "pre",
    "Mature sequences": "mature",
    "Star sequences": "star",
    "Genomic coordinates": "gff",
}

# Download files from links
for org_label, link_dict in alllinks.items():
    # find prefix code (hsa/mmu/...)
    matches = [code for name, code in nametocode.items() if name in org_label]
    if not matches:
        raise ValueError(f"no prefix found for {org_label}")
    prefix = matches[0]

    for kind, url in link_dict.items():
        suffix = typetosuffix[kind]

        if suffix == "gff":
            filename = out_dir / f"{prefix}_{suffix}.gff"
        else:
            filename = out_dir / f"{prefix}_{suffix}.fas"

        print(f"downloading {org_label} - {kind} -> {filename}")

        resp = requests.get(url)
        resp.raise_for_status()
        filename.write_bytes(resp.content)
