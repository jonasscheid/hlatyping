#!/usr/bin/env bash
# Regenerate data/references/{alleles.h5, hla_reference_dna.fasta, hla_reference_rna.fasta}
# from a fresh IPD-IMGT/HLA release. Self-contained — clones IMGTHLA, builds a
# pinned conda env, and writes outputs in place.
#
# Usage:  bash data/references/regenerate.sh [IMGT_TAG]
#         (default IMGT_TAG: v3.63.1-alpha)
#
# Env vars:
#   WORKDIR   scratch dir (default: mktemp)
#   ENV_NAME  conda env name (default: optitype-refgen-1.3.5)
#
# Requires: git, mamba (or conda), ~6 GB free disk for IMGTHLA.

set -euo pipefail

IMGT_TAG="${1:-v3.63.1-alpha}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$(mktemp -d -t optitype-refgen.XXXX)}"
ENV_NAME="${ENV_NAME:-optitype-refgen-1.3.5}"
CONDA="${CONDA:-mamba}"

echo "[regenerate] IMGT tag : $IMGT_TAG"
echo "[regenerate] workdir  : $WORKDIR"
echo "[regenerate] env name : $ENV_NAME"
echo "[regenerate] output   : $HERE"

# 1) Conda env — pinned to optitype 1.3.5 (matches the runtime container).
#    Bioconda's optitype=1.3.5 recipe lands on Python 2.7 + biopython 1.70 +
#    pandas 0.24 + pytables 3.5; hlatyper.py was written against this stack.
if ! $CONDA env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "[regenerate] creating conda env $ENV_NAME ..."
    $CONDA create -n "$ENV_NAME" -c bioconda -c conda-forge \
        optitype=1.3.5 'python<3.11' biopython 'pytables<3.10' -y
fi

# 2) IMGTHLA clone (skip if already present in WORKDIR).
if [ ! -d "$WORKDIR/IMGTHLA" ]; then
    echo "[regenerate] cloning ANHIG/IMGTHLA @ $IMGT_TAG ..."
    git clone --depth 1 --branch "$IMGT_TAG" \
        https://github.com/ANHIG/IMGTHLA.git "$WORKDIR/IMGTHLA"
fi
(
    cd "$WORKDIR/IMGTHLA"
    [ -f hla.dat       ] || unzip -o hla.dat.zip       > /dev/null
    [ -f hla_gen.fasta ] || unzip -o hla_gen.fasta.zip > /dev/null
)

# 3) Build outputs via the Python helper (see regenerate.py).
echo "[regenerate] running regenerate.py ..."
conda run --no-capture-output -n "$ENV_NAME" python "$HERE/regenerate.py" \
    --imgt-dir   "$WORKDIR/IMGTHLA" \
    --output-dir "$HERE"

echo
echo "[regenerate] done. Outputs:"
ls -lh "$HERE"/alleles.h5 "$HERE"/hla_reference_dna.fasta "$HERE"/hla_reference_rna.fasta
echo
echo "[regenerate] MD5s (paste these into assets/software_meta.json):"
md5sum "$HERE"/alleles.h5 "$HERE"/hla_reference_dna.fasta "$HERE"/hla_reference_rna.fasta
