"""Only annotation defaults; callers may supply an external data root explicitly."""
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
HGNC_PATH=ROOT/'pert_data/Annotation/source/hgnc_complete_set.txt'
COMPOUND_ANNOTATION_PATH=ROOT/'pert_data/Annotation/local/cm_compound_annotation.csv'
