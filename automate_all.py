import matplotlib.pyplot as plt
import numpy as np
import glob, os
from astropy.coordinates import SkyCoord, get_sun, BarycentricTrueEcliptic
from astropy.time import Time
import astropy.units as u

# make a folder to hold all the output plots and data
os.makedirs("output_plots", exist_ok=True)
os.makedirs("output_data", exist_ok=True)

# helper: read DMX values and ecliptic coords out of one par file
def read_par(par_path):
    dm_vals, mjd_vals = [], []
    elong = elat = None
    with open(par_path, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("ELONG"):
            elong = float(line.split()[1])
        elif line.startswith("ELAT"):
            elat = float(line.split()[1])
        elif line.startswith("DMX_") and not line.startswith("DMXR"):
            parts = line.split()
            try:
                dm = float(parts[1])
                r1 = float(lines[i+1].split()[1])
                r2 = float(lines[i+2].split()[1])
                dm_vals.append(dm)
                mjd_vals.append((r1 + r2) / 2)
            except (ValueError, IndexError):
                continue
    return np.array(mjd_vals), np.array(dm_vals), elong, elat

# process one pulsar
def process(par_path):
    name = os.path.basename(par_path).split("_")[0]
    mjd, dm, elong, elat = read_par(par_path)

    # need coords and enough points to be useful
    if elong is None or elat is None or len(mjd) < 10:
        print(f"  skipping {name} (missing coords or too few points)")
        return

    # convert ecliptic coords from par file to equatorial
    ecl = SkyCoord(lon=elong*u.deg, lat=elat*u.deg, frame=BarycentricTrueEcliptic)
    pulsar_coord = ecl.icrs

    # sun angle at each MJD
    times = Time(mjd, format="mjd")
    sun = get_sun(times)
    sun_icrs = SkyCoord(sun.ra, sun.dec, frame="icrs")
    sun_angle = pulsar_coord.separation(sun_icrs).deg

    # select points more than 30 deg from the Sun
    far = np.where(sun_angle > 30)
    near = np.where(sun_angle < 30)
    sel_mjd = mjd[far]
    sel_dm = dm[far]

    # detrend the selected DM
    if len(sel_mjd) < 2:
        print(f"  skipping {name} (not enough far points to detrend)")
        return
    slope, intercept = np.polyfit(sel_mjd, sel_dm, 1)
    detr = sel_dm - (slope * sel_mjd + intercept)

    # five-panel figure
    fig, axs = plt.subplots(5, 1, figsize=(10, 14), sharex=True)
    axs[0].plot(mjd, dm, 'k.', markersize=3); axs[0].set_ylabel("DM"); axs[0].set_title(name)
    axs[1].plot(mjd, sun_angle, 'b.', markersize=3); axs[1].set_ylabel("Sun angle (deg)")
    axs[2].plot(mjd[far], dm[far], 'k.', markersize=3)
    axs[2].plot(mjd[near], dm[near], 'r.', markersize=3); axs[2].set_ylabel("DM (red < 30)")
    axs[3].plot(sel_mjd, sel_dm, 'k.', markersize=3); axs[3].set_ylabel("Selected DM")
    axs[4].plot(sel_mjd, detr, 'k.', markersize=3); axs[4].axhline(0, color="gray", lw=0.7)
    axs[4].set_ylabel("Detrended DM"); axs[4].set_xlabel("MJD")
    plt.tight_layout()
    plt.savefig(f"output_plots/{name}_five_panel.png")
    plt.close()

    # save data with all five columns
    sel_full = np.full(len(mjd), np.nan)
    detr_full = np.full(len(mjd), np.nan)
    sel_full[far] = sel_dm
    detr_full[far] = detr
    np.savetxt(f"output_data/{name}_data.txt",
               np.column_stack([mjd, sun_angle, dm, sel_full, detr_full]),
               header="MJD  Sun_angle  DM  DM_selected  DM_detrended", fmt="%.6f")
    print(f"  done: {name}  ({len(mjd)} points)")

# run on every par file in the wideband folder
par_files = sorted(glob.glob("wideband/par/*.par"))
print(f"found {len(par_files)} par files")
for p in par_files:
    try:
        process(p)
    except Exception as e:
        print(f"  ERROR on {p}: {e}")

print("all done")