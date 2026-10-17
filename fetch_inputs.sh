#!/usr/bin/env bash
# fetch_inputs.sh -- download the public input files of manifests/EXTERNAL_INPUTS.tsv to their listed paths and
# check each against its sha256. Needs curl, gunzip, unzip, sha256sum and python3.
# Empty manifest fields are read as "-". Files already present with the right checksum are skipped. Exits non-zero if any file is missing or differs.
# Live sources (KEGG REST, GEO text pages) can change after the listed fetch date; a checksum mismatch on those
# is reported, not hidden.
set -u
cd "$(dirname "$0")"
MAN=manifests/EXTERNAL_INPUTS.tsv
fail=0; manual=()

check() { [ -f "$1" ] && [ "$(sha256sum "$1" | cut -c1-64)" = "$2" ]; }

pincus() {   # first six columns of sheet 1 of s002.xlsx, header spaces -> underscores
  local dest=$1 tmp; tmp=$(mktemp -d)
  curl -fsSL --retry 3 -o "$tmp/supp.zip" "$2" && unzip -q -o "$tmp/supp.zip" mbc-29-3168-s002.xlsx -d "$tmp" &&
  python3 - "$tmp/mbc-29-3168-s002.xlsx" "$dest" <<'PY'
import sys, zipfile, xml.etree.ElementTree as ET
z = zipfile.ZipFile(sys.argv[1]); a = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
ss = ["".join(t.itertext()) for t in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(a + "si")]
out = []
for r in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).iter(a + "row"):
    row = []
    for c in r.findall(a + "c"):
        v = c.find(a + "v"); x = v.text if v is not None else ""
        row.append(ss[int(x)] if c.get("t") == "s" else x)
    out.append(row[:6])
out[0] = [h.replace(" ", "_") for h in out[0]]
def num(x):
    try: return repr(float(x)) if "." in x else x
    except ValueError: return x
with open(sys.argv[2], "w", newline="") as fh:
    for i, row in enumerate(out):
        fh.write("\t".join(row if i == 0 else row[:2] + [num(x) for x in row[2:]]) + "\r\n")
PY
  local rc=$?; rm -rf "$tmp"; return $rc
}

while IFS=$'\t' read -r file url fetched bytes sha note; do
  [ "$file" = "file" ] && continue
  if check "$file" "$sha"; then echo "  ok        $file"; continue; fi
  [ "$url" = "-" ] && url=""
  if [ -z "$url" ]; then manual+=("$file  ($note)"); continue; fi
  mkdir -p "$(dirname "$file")"
  case "$file" in
    *pincus2018_nacseq_counts.tsv) pincus "$file" "$url" ;;
    *.gtf) curl -fsSL --retry 3 -o "$file.gz" "$url" && gunzip -f "$file.gz" ;;
    *) curl -fsSL --retry 3 -o "$file.part" "$url" && mv "$file.part" "$file" ;;
  esac
  if check "$file" "$sha"; then echo "  fetched   $file"
  elif [ -f "$file" ]; then echo "  MISMATCH  $file (sha256 differs from the manifest; the source may have changed since $fetched)"; fail=1
  else echo "  FAILED    $file ($url)"; fail=1; fi
done < <(sed -e 's/\t\t/\t-\t/g; s/\t\t/\t-\t/g' "$MAN")

if [ ${#manual[@]} -gt 0 ]; then
  echo "Not downloadable automatically; place these at the listed path:"
  printf '  %s\n' "${manual[@]}"; fail=1
fi
exit $fail
