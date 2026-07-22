import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

pulsars, rms, dist, reddening = [], [], [], []
with open("dust_reddening_all76.txt") as f:
    next(f)
    for line in f:
        parts = line.split()
        pulsars.append(parts[0])
        rms.append(float(parts[1]))
        dist.append(float(parts[2]))
        reddening.append(float(parts[4]))

pulsars = np.array(pulsars)
rms_log = np.log10(np.array(rms))
dist = np.array(dist)
reddening = np.array(reddening)

def residualize(y, x):
    slope, intercept = np.polyfit(x, y, 1)
    return y - (slope * x + intercept)

rms_resid = residualize(rms_log, dist)
reddening_resid = residualize(reddening, dist)

r_partial, p_partial = pearsonr(reddening_resid, rms_resid)

plt.figure(figsize=(9, 7))
plt.scatter(reddening_resid, rms_resid, s=60, alpha=0.7, color='teal')
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.axvline(0, color='gray', linestyle='--', linewidth=1)

plt.xlabel("Dust reddening residual (distance removed)")
plt.ylabel("log(DM RMS) residual (distance removed)")
plt.title(f"Dust vs DM RMS, controlling for distance\nPartial r={r_partial:.3f}, p={p_partial:.4f}")

textstr = f"Partial Pearson r={r_partial:.3f}, p={p_partial:.4f}\n(not significant)"
plt.gca().text(0.02, 0.98, textstr, transform=plt.gca().transAxes,
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

plt.tight_layout()
plt.savefig("dust_dmrms_distance_controlled.png", dpi=150)
plt.show()

print(f"Partial correlation: r={r_partial:.4f}, p={p_partial:.4f}")