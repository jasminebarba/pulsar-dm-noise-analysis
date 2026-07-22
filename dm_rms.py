import numpy as np
import glob, os

results = []
for path in sorted(glob.glob("output_data/*_data.txt")):
    name = os.path.basename(path).replace("_data.txt", "")
    # columns: MJD, Sun_angle, DM, DM_selected, DM_detrended
    data = np.loadtxt(path)
    if data.ndim < 2 or data.shape[1] < 5:
        continue
    detrended = data[:, 4]
    detrended = detrended[~np.isnan(detrended)]   # drop the nan (near-Sun) rows
    if len(detrended) < 2:
        continue
    rms = np.std(detrended)
    results.append((name, rms, len(detrended)))

# sort by RMS, most variable first
results.sort(key=lambda x: x[1], reverse=True)

print(f"{'Pulsar':<20} {'DM_RMS':<12} {'N_points'}")
print("-" * 45)
for name, rms, n in results:
    print(f"{name:<20} {rms:<12.6e} {n}")

# save to a file
with open("dm_rms.txt", "w") as f:
    f.write("Pulsar  DM_RMS  N_points\n")
    for name, rms, n in results:
        f.write(f"{name}  {rms:.6e}  {n}\n")

print(f"\nSaved {len(results)} pulsars to dm_rms.txt")