import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, get_sun
from astropy.time import Time
import astropy.units as u

# 1. Read the raw DM data
data_file = "wideband/par/J1909-3744gbt_PINT_20230131.wb.par"
dm_vals, mjd_vals = [], []
base_dm = None
with open(data_file, "r") as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("DM ") or line.startswith("DM\t"):
            try:
                base_dm = float(line.split()[1])
            except (ValueError, IndexError):
                pass
        if line.startswith("DMX_") and not line.startswith("DMXR"):
            parts = line.split()
            dm_vals.append(float(parts[1]))
            mjd_vals.append((float(lines[i+1].split()[1]) + float(lines[i+2].split()[1])) / 2)

mjd_vals = np.array(mjd_vals)
dm_vals = np.array(dm_vals)

# 2. Sun angle at each MJD
ra = "19:09:47.4"
dec = "-37:44:14.5"
pulsar_coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
times = Time(mjd_vals, format='mjd')
sun = get_sun(times)
sun_icrs = SkyCoord(sun.ra, sun.dec, frame="icrs")
solar_angles = pulsar_coord.separation(sun_icrs).deg

# 3. Select points more than 30 degrees from the Sun
far = np.where(solar_angles > 30)
near = np.where(solar_angles < 30)
selected_mjd = mjd_vals[far]
selected_dm = dm_vals[far]

# 4. Linear fit and detrend
slope, intercept = np.polyfit(selected_mjd, selected_dm, 1)
fit_line = slope * selected_mjd + intercept
detrended_dm = selected_dm - fit_line

# Find the MJDs where the sun angle is at a local minimum (closest approach)
order = np.argsort(mjd_vals)
m_sorted = mjd_vals[order]
a_sorted = solar_angles[order]
min_mjds = []
for k in range(1, len(a_sorted) - 1):
    if a_sorted[k] < a_sorted[k-1] and a_sorted[k] < a_sorted[k+1] and a_sorted[k] < 40:
        min_mjds.append(m_sorted[k])

# 5. Four-panel figure
fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
title = "J1909-3744"
if base_dm is not None:
    title += f"   (DM $\\approx$ {base_dm:.2f} pc/cm$^3$)"
axs[0].set_title(title)

axs[0].plot(mjd_vals, solar_angles, 'b.', markersize=3)
axs[0].set_ylabel("Sun angle (deg)")

axs[1].plot(mjd_vals[far], dm_vals[far], 'k.', markersize=3, label="Measured $\\Delta$DM")
axs[1].plot(mjd_vals[near], dm_vals[near], 'r.', markersize=3, label=r"$\theta_{sun}$ < 30°")
axs[1].set_ylabel("$\\Delta$DM (pc/cm$^3$)")
axs[1].legend(loc="upper right", fontsize=8)

axs[2].plot(selected_mjd, selected_dm, 'k.', markersize=3, label="Selected $\\Delta$DM")
axs[2].plot(selected_mjd, fit_line, 'r-', lw=1.5, label="Linear fit")
axs[2].set_ylabel("Selected $\\Delta$DM")
axs[2].legend(loc="upper right", fontsize=8)

axs[3].plot(selected_mjd, detrended_dm, 'k.', markersize=3)
axs[3].set_ylabel("Detrended $\\Delta$DM")
axs[3].set_xlabel("MJD")

for mj in min_mjds:
    for ax in axs:
        ax.axvline(mj, color="gray", lw=0.6, alpha=0.5, zorder=0)

plt.tight_layout()
plt.savefig("five_panel.png", dpi=120)

# 6. Save data
dm_selected_full = np.full(len(mjd_vals), np.nan)
dm_detrended_full = np.full(len(mjd_vals), np.nan)
dm_selected_full[far] = selected_dm
dm_detrended_full[far] = detrended_dm
np.savetxt("j1909_data.txt",
           np.column_stack([mjd_vals, solar_angles, dm_vals,
                            dm_selected_full, dm_detrended_full]),
           header="MJD  Sun_angle  dDM  dDM_selected  dDM_detrended", fmt="%.6f")

plt.show()