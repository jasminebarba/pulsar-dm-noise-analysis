import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

pulsars, rms, reddening = [], [], []
with open("dust_reddening_all76.txt") as f:
    next(f)
    for line in f:
        parts = line.split()
        pulsars.append(parts[0])
        rms.append(float(parts[1]))
        reddening.append(float(parts[4]))

rms = np.array(rms)
reddening = np.array(reddening)
pulsars = np.array(pulsars)

r_p, p_p = pearsonr(reddening, rms)
r_s, p_s = spearmanr(reddening, rms)

# Pulsars to label: anything with E(g-r) > 1, plus your three named ones from group discussion
always_label = {"J1903+0327", "J1713+0747", "B1953+29"}
label_mask = (reddening > 1.0) | np.isin(pulsars, list(always_label))

plt.figure(figsize=(11, 8))
plt.scatter(reddening[~label_mask], rms[~label_mask], s=60, alpha=0.6, color='darkred', label='Other pulsars')
plt.scatter(reddening[label_mask], rms[label_mask], s=90, alpha=0.9, color='gold', edgecolor='black', zorder=5, label='Labeled pulsars')

for name, x, y in zip(pulsars[label_mask], reddening[label_mask], rms[label_mask]):
    plt.annotate(name, (x, y), textcoords="offset points", xytext=(6, 6), fontsize=9)

plt.yscale('log')
plt.xlabel("Dust reddening E(g-r) (mag)")
plt.ylabel("Pulsar DM RMS (pc/cm$^3$)")
plt.title(f"DM RMS vs Dust Reddening — {len(pulsars)} pulsars")

# Correlation coefficients printed directly on the plot
textstr = f"Pearson r={r_p:.3f}, p={p_p:.4f}\nSpearman r={r_s:.3f}, p={p_s:.4f}"
plt.gca().text(0.02, 0.98, textstr, transform=plt.gca().transAxes,
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("dust_vs_dmrms_labeled.png", dpi=150)
plt.show()

print(f"Labeled pulsars: {list(pulsars[label_mask])}")
print(f"Pearson r={r_p:.4f}, p={p_p:.4f}")
print(f"Spearman r={r_s:.4f}, p={p_s:.4f}")