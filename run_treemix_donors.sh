#!/usr/bin/env bash

set -euo pipefail

mkdir -p "results/treemix"

# ==================== Spanish ====================

treemix -i "data/treemix/spanishtreemix_phonemes96.gz" -tf "data/treemix/spanishtreemix_fixed_tree.nwk" -root "Spanish" -o "results/treemix/spanish_m0" -m 0 -k 1 -noss

treemix -i "data/treemix/spanishtreemix_phonemes96.gz" -tf "data/treemix/spanishtreemix_fixed_tree.nwk" -root "Spanish" -o "results/treemix/spanish_m1" -m 1 -k 1 -noss


# ==================== English ====================

treemix -i "data/treemix/englishtreemix_phonemes96.gz" -tf "data/treemix/englishtreemix_fixed_tree.nwk" -root "English" -o "results/treemix/english_m0" -m 0 -k 1 -noss

treemix -i "data/treemix/englishtreemix_phonemes96.gz" -tf "data/treemix/englishtreemix_fixed_tree.nwk" -root "English" -o "results/treemix/english_m1" -m 1 -k 1 -noss


# ==================== Japanese ====================

treemix -i "data/treemix/japanesetreemix_phonemes96.gz" -tf "data/treemix/japanesetreemix_fixed_tree.nwk" -root "Japanese" -o "results/treemix/japanese_m0" -m 0 -k 1 -noss

treemix -i "data/treemix/japanesetreemix_phonemes96.gz" -tf "data/treemix/japanesetreemix_fixed_tree.nwk" -root "Japanese" -o "results/treemix/japanese_m1" -m 1 -k 1 -noss


