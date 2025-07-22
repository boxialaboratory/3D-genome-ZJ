import os
import re

def extract_microc_number(filename):
    match = re.search(r'_MicroC(\d+)_', filename)
    return int(match.group(1)) if match else -1

def parse_stat_file(filepath):
    values = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                key = parts[0]
                try:
                    value = int(parts[1])
                except ValueError:
                    value = float(parts[1])
                values[key] = value
    return values

# collect and sort
data = []
for filename in os.listdir("."):
    if filename.endswith(".stat"):
        microc_num = extract_microc_number(filename)
        stats = parse_stat_file(filename)
        total = stats.get("total", 0)
        total_dups = stats.get("total_dups", 0)
        total_nodups = stats.get("total_nodups", 0)
        cis = stats.get("cis", 0)
        trans = stats.get("trans", 0)
        nodup_ratio = total_nodups / total if total else 0
        cis_trans_ratio = cis / trans if trans else float("inf")
        data.append({
            "filename": filename,
            "microc_num": microc_num,
            "total": int(total),
            "total_dups": int(total_dups),
            "total_nodups": int(total_nodups),
            "nodup_ratio": f"{nodup_ratio:.4f}",
            "cis": int(cis),
            "trans": int(trans),
            "cis_trans_ratio": f"{cis_trans_ratio:.4f}",
        })

# sort
data.sort(key=lambda x: x["microc_num"])

# print by columns
columns = [
    "filename", "total", "total_dups", "total_nodups",
    "nodup_ratio", "cis", "trans", "cis_trans_ratio"
]

for col in columns:
    print(f"\n### {col}")
    for row in data:
        print(row[col])
