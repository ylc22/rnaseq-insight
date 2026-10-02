import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
import pandas as pd
from rnaseq_insight.pipeline import normalize_log_cpm, differential_expression


def test_normalization_shape():
    x = pd.DataFrame({"s1":[10,20],"s2":[30,40]}, index=["g1","g2"])
    assert normalize_log_cpm(x).shape == x.shape


def test_de_ranks_strong_signal():
    x = pd.DataFrame({"a":[1,5],"b":[1.1,5],"c":[8,5],"d":[9,5]}, index=["g1","g2"])
    y = pd.Series({"a":"ctrl","b":"ctrl","c":"case","d":"case"})
    out = differential_expression(x, y, "case", "ctrl")
    assert out.iloc[0].gene == "g1"
