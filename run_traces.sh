#!/bin/bash

if [ $# -eq 0 ]; then
	echo "usage: run_traces.sh <out_dir>" >&2
	exit 1
fi

dir="$(dirname "${BASH_SOURCE[0]}")"
outdir="$1"

mkdir -p "$outdir"
for f in "$dir/traces/"*.bz2; do
    echo "$f"
	./cbp3 -t "$f" "${@:2}" | tee "$outdir/$(basename $f .bz2).out"
	echo
	echo
done
