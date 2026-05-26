#!/bin/bash
# run-cm-fda.sh — Launch cardiomyocyte FDA drug virtual screening pipeline
set -e

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$PROJECT_ROOT"

PYTHON="$HOME/miniconda3/envs/xpert/bin/python"
DEVICE="${1:-cuda:0}"

echo "=== Cardiomyocyte FDA Virtual Screening Pipeline ==="
echo "Project root: $PROJECT_ROOT"
echo "Device: $DEVICE"

# --- Step 1: FDA drug identification ---
if [ ! -f "vcms/fda/fda_drugs.csv" ]; then
    echo ""
    echo "=== Step 1: Identifying FDA-approved drugs ==="
    $PYTHON vcms/src/fda_drugs.py
else
    echo ""
    echo "=== Step 1: vcms/fda/fda_drugs.csv exists, skipping ==="
fi

# --- Step 2: Control transcriptome ---
if [ ! -f "vcms/Xcm.csv" ]; then
    echo ""
    echo "=== Step 2: Building control cardiomyocyte transcriptome ==="
    $PYTHON vcms/src/Xcm.py
else
    echo ""
    echo "=== Step 2: vcms/Xcm.csv exists, skipping ==="
fi

# --- Step 3: Virtual screening ---
echo ""
echo "=== Step 3: Running ensemble virtual screening ==="
$PYTHON vcms/src/vcms_fda_screen.py --device "$DEVICE" "${@:2}"

echo ""
echo "=== Pipeline complete ==="
echo "Output files:"
echo "  vcms/fda/fda_drugs.csv       — FDA-approved drug list"
echo "  vcms/Xcm.csv                 — Unified control transcriptome"
echo "  vcms/fda/cm_fda_pert_deg.csv — Predicted perturbation deg matrix"
