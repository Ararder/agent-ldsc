"""Semantic comparison of two ldsR output directories (not byte identity)."""
import json
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

a, b = map(Path, sys.argv[1:3])
for f in ("ld.parquet", "annot.parquet", "annot_ref.parquet"):
    pd.testing.assert_frame_equal(pq.read_table(a / f).to_pandas(), pq.read_table(b / f).to_pandas(),
                                  check_exact=True)
ma, mb = (json.loads((d / "manifest.json").read_text()) for d in (a, b))
assert ma["scientific_key"] == mb["scientific_key"], "scientific keys differ"
print("parity ok:", ma["scientific_key"])
