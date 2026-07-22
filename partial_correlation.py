import numpy as np
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

rms = np.log10(np.array(rms))  # log since DM RMS spans orders of magnitude
dist = np.array(dist)
reddening = np.array(reddening)

# Partial correlation: regress out distance from both variables, then correlate the residuals
def residualize(y, x):
    slope, intercept = np.polyfit(x, y, 1)
    return y - (slope * x + intercept)

rms_resid = residualize(rms, dist)
reddening_resid = residualize(reddening, dist)

r_partial, p_partial = pearsonr(reddening_resid, rms_resid)
print(f"Partial correlation (dust vs DM RMS, controlling for distance): r={r_partial:.4f}, p={p_partial:.4f}")