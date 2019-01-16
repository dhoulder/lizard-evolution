#!/bin/sh
cd "$(dirname "$0")" || exit 1
for d in */; do mkdir -p "$d/out/"; done
rm -f */out/*.{yml,csv,nc}
