from pathlib import Path

import pandas as pd

out_dir = Path("data")

nametocode = {
    "Human": "hsa",
    "House mouse": "mmu",
    "Fruit fly": "dme",
    "Roundworm": "cel",
}
# Download browse tabs as csv
for org_name, code in nametocode.items():
    url = f"https://mirgenedb.org/browse/{code}"
    print(f"processing {org_name} from {url}")

    # read all tables in the page
    tables = pd.read_html(url, header=None)

    if not tables:
        print(f"no tables found for {org_name} at {url}")
        continue

    df = tables[0]  # first table on the page

    out_name = out_dir / f"{code}.csv"
    df.to_csv(out_name, index=False)
