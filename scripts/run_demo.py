from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from rnaseq_insight.pipeline import filter_low_expression, normalize_log_cpm, differential_expression

rng = np.random.default_rng(42)
ng, ns = 1200, 16
genes = [f"GENE_{i:04d}" for i in range(ng)]
samples = [f"S{i:02d}" for i in range(ns)]
labels = pd.Series(["control"]*8 + ["treated"]*8, index=samples, name="condition")
base = rng.gamma(2.2, 85, size=ng)
counts = np.zeros((ng, ns), int)
for j in range(ns):
    mu = base.copy()
    if labels.iloc[j] == "treated":
        mu[:60] *= 2.8
        mu[60:110] *= 0.35
    counts[:, j] = rng.poisson(mu)
counts = pd.DataFrame(counts, index=genes, columns=samples)

(ROOT/"data").mkdir(exist_ok=True)
outdir = ROOT/"outputs"; outdir.mkdir(exist_ok=True)
counts.to_csv(ROOT/"data"/"demo_counts.csv")
labels.to_csv(ROOT/"data"/"demo_metadata.csv", header=True)

filtered = filter_low_expression(counts)
logcpm = normalize_log_cpm(filtered)
de = differential_expression(logcpm, labels, "treated", "control")
de.to_csv(outdir/"differential_expression.csv", index=False)

pca = PCA(n_components=2).fit_transform(logcpm.T)
fig, ax = plt.subplots(figsize=(7,5))
for group in labels.unique():
    m = labels.values == group
    ax.scatter(pca[m,0], pca[m,1], label=group, s=55)
ax.set(title="RNA-seq sample separation", xlabel="PC1", ylabel="PC2"); ax.legend(); fig.tight_layout(); fig.savefig(outdir/"pca.png", dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(7,5))
y = -np.log10(de.p_value.clip(lower=1e-300))
sig = (de.fdr < .05) & (de.abs_log2_fc > 1)
ax.scatter(de.log2_fc[~sig], y[~sig], s=10, alpha=.45)
ax.scatter(de.log2_fc[sig], y[sig], s=14, alpha=.75)
ax.axvline(-1, ls="--", lw=1); ax.axvline(1, ls="--", lw=1)
ax.set(title="Differential expression", xlabel="log2 fold-change", ylabel="-log10(p-value)"); fig.tight_layout(); fig.savefig(outdir/"volcano.png", dpi=160); plt.close(fig)
print(de.head(15).to_string(index=False))
print(f"\nSignificant genes: {sig.sum()} / {len(de)}")
