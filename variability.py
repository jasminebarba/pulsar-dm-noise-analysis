import numpy as np
import csv, sys

csv_file = sys.argv[1]
mjd, mag, magerr, catflags = [], [], [], []
oid = None
with open(csv_file, "r") as f:
    for row in csv.DictReader(f):
        oid = row["oid"]
        mjd.append(float(row["mjd"])); mag.append(float(row["mag"]))
        magerr.append(float(row["magerr"])); catflags.append(int(row["catflags"]))

mag = np.array(mag); magerr = np.array(magerr); catflags = np.array(catflags)
good = catflags == 0
mag, magerr = mag[good], magerr[good]

scatter = np.std(mag)                 # how much the brightness actually varies
typical_error = np.median(magerr)     # typical measurement uncertainty
ratio = scatter / typical_error       # >1 means real variation beyond noise

print(f"star {oid}")
print(f"scatter (std of mag): {scatter:.4f}")
print(f"typical error:        {typical_error:.4f}")
print(f"ratio (scatter/error): {ratio:.2f}")
print("-> likely VARIABLE" if ratio > 2 else "-> consistent with noise")