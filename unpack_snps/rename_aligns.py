# python3 scripts/unpack_snps/rename_aligns.py --root svinoserver --meta meta/meta_6_may.tsv --temp temp

import re
import os
import shutil
import argparse
import requests
import pandas as pd
from pathlib import Path

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
PATTERN = re.compile(r"^(ALIGNS\d{6,7})_(.+)\.vcf\.gz$")

def fetch_uniprot_name(uniprot_id, cache: dict) -> str | None:
    """Вернуть «короткое» имя фактора (до первого _), кэшируя запросы."""
    if uniprot_id in cache:
        return cache[uniprot_id]

    r = requests.get(UNIPROT_URL.format(uniprot_id), timeout=15)
    if r.ok:
        name = r.json().get("uniProtkbId")
        cache[uniprot_id] = name
        return name
    print(f"[WARN] UniProt {uniprot_id}: HTTP {r.status_code}")
    cache[uniprot_id] = None
    return None

def main(root_dir: Path, meta_path: Path, temp_dir: Path) -> None:
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    meta = meta.set_index("algn_id")

    uniprot_cache: dict[str, str | None] = {}

    for path in root_dir.glob("ALIGNS*_*.vcf.gz"):
        m = PATTERN.match(path.name)
        if not m:
            continue

        algn_id, cell_type = m.groups()
        if algn_id not in meta.index:
            print(f"[SKIP] {path.name}: {algn_id} нет в метаданных")
            continue

        tf_uniprot_id = meta.at[algn_id, "tf_uniprot_id"]
        factor = fetch_uniprot_name(tf_uniprot_id, uniprot_cache)
        if not factor:
            print(f"[SKIP] {path.name}: не смогли получить UniProt-имя для {tf_uniprot_id}")
            continue

        new_name = f"{factor.split('_')[0]}_{cell_type}_{algn_id}.vcf.gz"
        dst = temp_dir / new_name
        shutil.copy2(path, dst) 
        print(f"✔ {path.name}\t→\ttemp/{new_name}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Копирует файлы ALIGNS.._*.vcf.gz в temp с фактором-префиксом.")
    p.add_argument("--root", default="svinoserver",
                   help="Каталог, где лежат файлы (по умолчанию svinoserver)")
    p.add_argument("--meta", default="meta/meta_6_may.tsv",
                   help="Путь к meta-файлу (TSV)")
    p.add_argument("--temp", default="temp",
                   help="Куда класть результат (по умолчанию ./temp)")
    args = p.parse_args()

    main(Path(args.root), Path(args.meta), Path(args.temp))
