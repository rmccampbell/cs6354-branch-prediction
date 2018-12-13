#!/bin/bash

if [[ "$#" -eq 0 || "$1" =~ ^(-h|--help)$ ]]; then
    echo "usage: run_traces.sh <out_dir> [-u num-uops-to-simulate]" >&2
    exit 1
fi

make

dir="$(dirname "${BASH_SOURCE[0]}")"
outdir="$1"

mkdir -p "$outdir"
i=1
for f in "$dir/traces/"*.bz2; do
    echo "$i: $f"
    ((i++))
    ./cbp3 -t "$f" "${@:2}" | tee "$outdir/$(basename $f .bz2).out"
    echo
    echo
done
