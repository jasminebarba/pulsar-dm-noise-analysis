from scipy.stats import pearsonr, spearmanr
import numpy as np

rms, reddening = [], []
with open('dust_reddening_all76.txt') as f:
    next(f)
    for line in f:
        parts = line.split()
        rms.append(float(parts[1]))
        reddening.append(float(parts[4]))

rms = np.array(rms)
reddening = np.array(reddening)
r_p, p_p = pearsonr(reddening, rms)
r_s, p_s = spearmanr(reddening, rms)
print(f'n={len(rms)}')
print(f'Pearson: r={r_p:.4f}, p={p_p:.4f}')
print(f'Spearman: r={r_s:.4f}, p={p_s:.4f}')
