import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


def filter_low_expression(counts: pd.DataFrame, min_cpm: float = 1.0, min_samples: int = 3) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= min_samples
    return counts.loc[keep]


def normalize_log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1)


def _bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def differential_expression(log_expr: pd.DataFrame, labels: pd.Series, case: str, control: str) -> pd.DataFrame:
    case_cols = labels.index[labels == case]
    ctrl_cols = labels.index[labels == control]
    rows = []
    for gene, row in log_expr.iterrows():
        a, b = row[case_cols].values, row[ctrl_cols].values
        stat, p = ttest_ind(a, b, equal_var=False)
        lfc = np.nanmean(a) - np.nanmean(b)
        rows.append((gene, lfc, p))
    out = pd.DataFrame(rows, columns=["gene", "log2_fc", "p_value"])
    out["fdr"] = _bh(out["p_value"].fillna(1.0))
    out["abs_log2_fc"] = out["log2_fc"].abs()
    return out.sort_values(["fdr", "abs_log2_fc"], ascending=[True, False]).reset_index(drop=True)
