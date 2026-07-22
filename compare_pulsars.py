import numpy as np
import matplotlib.pyplot as plt
import csv, glob

def get_ratios(pattern):
    ratios = []
    for f in sorted(glob.glob(pattern)):
        mag, magerr, catflags = [], [], []
        try:
            with open(f) as fh:
                for row in csv.DictReader(fh):
                    mag.append(float(row["mag"]))
                    magerr.append(float(row["magerr"]))
                    catflags.append(int(row["catflags"]))
        except Exception:
            continue
        mag = np.array(mag); magerr = np.array(magerr); catflags = np.array(catflags)
        good = catflags == 0
        if good.sum() < 2:
            continue
        ratios.append(np.std(mag[good]) / np.median(magerr[good]))
    return ratios

pulsars = {
    "J1713": ("star*_lightcurve.csv", 1.30e-4, 'blue'),
    "J1903": ("j1903_star*_lightcurve.csv", 6.44e-3, 'red'),
    "J1946": ("J1946_stars/*.csv", 2.5e-3, 'green'),
    "B1953+29": ("B1953_stars/*.csv", 9.59e-4, 'purple'),
    "J2043": ("J2043_stars/*.csv", 8.32e-5, 'orange'),
    "J1640": ("J1640_stars/*.csv", 1.77e-4, 'brown'),
    "J1911": ("J1911_stars/*.csv", 3.04e-4, 'teal'),
    "J1918": ("J1918_stars/*.csv", 3.98e-4, 'magenta'),
    "J2017": ("J2017_stars/*.csv", 3.99e-4, 'gold'),
    "B1855+09": ("B1855_stars/*.csv", 5.16e-4, 'cyan'),
    "J1910": ("J1910_stars/*.csv", 7.62e-4, 'darkgreen'),
    "J1832": ("J1832_stars/*.csv", 1.42e-3, 'crimson'),
    "J1738": ("J1738_stars/*.csv", 3.98e-4, 'navy'),
    "J1745": ("J1745_stars/*.csv", 1.44e-3, 'olive'),
    "J0030": ("J0030_stars/*.csv", 1.34e-4, 'darkorange'),
    "J2214+3000": ("J2214_stars/*.csv", 9.23e-4, 'slateblue'),
    "J0613-0200": ("J0613_stars/*.csv", 2.45e-4, 'deeppink'),
    "J0645+5158": ("J0645_stars/*.csv", 1.09e-4, 'sienna'),
}

plt.figure(figsize=(13, 7.5))

for name, (pattern, rms, color) in pulsars.items():
    ratios = get_ratios(pattern)
    if not ratios:
        print(f"WARNING: {name} has no data, skipping")
        continue
    plt.scatter([rms]*len(ratios), ratios, s=50, color=color, alpha=0.6, label=f'{name}')
    plt.scatter([rms], [np.median(ratios)], s=200, color=color, marker='_', linewidths=3)

annotations = [
    (2.5e-3, 4.83, "J1946: pulsar\ncontamination"),
    (5.16e-4, 4.59, "B1855+09: eclipsing\nbinary"),
    (3.98e-4, 2.99, "J1918: open lead\n(unexplained)"),
]
for x, y, note in annotations:
    plt.annotate(note, (x, y), xytext=(30, 15), textcoords="offset points",
                 fontsize=8, style='italic',
                 arrowprops=dict(arrowstyle='->', color='black', lw=0.8))

plt.axhline(2, color='gray', linestyle='--', label='variability threshold')
plt.xscale('log')
plt.xlabel("Pulsar DM RMS (pc/cm$^3$)")
plt.ylabel("Star variability ratio (scatter / error)")
plt.title("Stellar variability vs pulsar DM variability — 18 pulsars")
plt.legend(fontsize=8, ncol=2, loc='upper left', bbox_to_anchor=(1.01, 1.0))
plt.tight_layout()
plt.savefig("pulsar_comparison_full.png", dpi=120, bbox_inches='tight')
plt.show()

print("\n--- Summary ---")
total_stars = 0
for name, (pattern, rms, color) in pulsars.items():
    ratios = get_ratios(pattern)
    if ratios:
        n_variable = sum(1 for r in ratios if r > 2)
        total_stars += len(ratios)
        print(f"{name} (DM RMS {rms:.2e}): {len(ratios)} stars, median {np.median(ratios):.2f}, {n_variable} above threshold")
print(f"\nTotal: {len(pulsars)} pulsars, {total_stars} stars")