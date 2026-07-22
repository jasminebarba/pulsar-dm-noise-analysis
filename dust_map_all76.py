import subprocess, re, glob
from astropy.coordinates import SkyCoord, BarycentricTrueEcliptic
import astropy.units as u
from dustmaps.bayestar import BayestarQuery

# Distances in kpc, from Agazie et al. 2023 (NANOGrav 15yr, ApJL 951, L50), Table 2.
# Key = base pulsar name (no ao/gbt suffix). Value = (distance_kpc, method)
DISTANCES = {
    "B1855+09": (1.18, "PX"), "B1937+21": (3.10, "PX"), "B1953+29": (4.64, "DM"),
    "J0023+0923": (1.02, "PX"), "J0030+0451": (0.3296, "PX"), "J0340+4130": (1.71, "DM"),
    "J0406+3039": (1.72, "DM"), "J0437-4715": (0.1549, "PX"), "J0509+0856": (1.45, "DM"),
    "J0557+1551": (2.87, "DM"), "J0605+3757": (0.70, "DM"), "J0610-2100": (1.38, "PX"),
    "J0613-0200": (1.07, "PX"), "J0636+5128": (1.27, "PX"), "J0645+5158": (1.37, "PX"),
    "J0709+0458": (1.80, "DM"), "J0740+6620": (0.94, "PX"), "J0931-1902": (1.88, "DM"),
    "J1012+5307": (0.862, "PX"), "J1012-4235": (2.50, "DM"), "J1022+1001": (0.706, "PX"),
    "J1024-0719": (1.080, "PX"), "J1125+7819": (0.63, "DM"), "J1312+0051": (0.84, "DM"),
    "J1453+1902": (1.15, "DM"), "J1455-3330": (1.01, "PX"), "J1600-3053": (1.84, "PX"),
    "J1614-2230": (0.699, "PX"), "J1630+3734": (0.089, "PX"), "J1640+2224": (1.404, "PX"),
    "J1643-1224": (0.835, "PX"), "J1705-1903": (1.62, "DM"), "J1713+0747": (1.138, "PX"),
    "J1719-1438": (1.21, "DM"), "J1730-2304": (0.529, "PX"), "J1738+0333": (1.64, "PX"),
    "J1741+1351": (2.36, "PX"), "J1744-1134": (0.4141, "PX"), "J1745+1017": (1.27, "DM"),
    "J1747-4036": (3.50, "DM"), "J1751-2857": (1.11, "DM"), "J1802-2124": (2.96, "DM"),
    "J1811-2405": (1.79, "DM"), "J1832-0836": (2.00, "PX"), "J1843-1113": (1.72, "DM"),
    "J1853+1303": (1.91, "PX"), "J1903+0327": (6.49, "DM"), "J1909-3744": (1.159, "PX"),
    "J1910+1256": (3.52, "PX"), "J1911+1347": (2.08, "DM"), "J1918-0642": (1.44, "PX"),
    "J1923+2515": (0.94, "PX"), "J1944+0907": (1.38, "PX"), "J1946+3417": (5.12, "DM"),
    "J2010-1323": (1.94, "PX"), "J2017+0603": (1.57, "DM"), "J2033+1734": (1.99, "DM"),
    "J2043+1711": (1.58, "PX"), "J2124-3358": (0.413, "PX"), "J2145-0750": (0.624, "PX"),
    "J2214+3000": (1.54, "DM"), "J2229+2643": (1.43, "DM"), "J2234+0611": (1.23, "PX"),
    "J2234+0944": (1.00, "DM"), "J2302+4442": (1.18, "DM"), "J2317+1439": (1.57, "PX"),
    "J2322+2057": (1.00, "PX"),
}

def base_name(name):
    return re.sub(r'(ao|gbt)$', '', name)

def get_radec(base):
    matches = glob.glob(f"wideband/par/{base}_PINT*.par")
    if not matches:
        return None
    with open(matches[0]) as f:
        elong, elat = None, None
        for line in f:
            if line.startswith("ELONG"):
                elong = float(line.split()[1])
            elif line.startswith("ELAT"):
                elat = float(line.split()[1])
    if elong is None or elat is None:
        return None
    c = SkyCoord(lon=elong*u.deg, lat=elat*u.deg, frame=BarycentricTrueEcliptic).icrs
    return c.ra.deg, c.dec.deg

pulsars = []
with open("dm_rms.txt") as f:
    next(f)
    for line in f:
        parts = line.split()
        name, rms = parts[0], float(parts[1])
        pulsars.append((name, rms))

print(f"Loaded {len(pulsars)} entries from dm_rms.txt")

bayestar = BayestarQuery(map_fname='./dustmaps_data/bayestar/bayestar2019.h5')

results = []
skipped = []
seen_base = set()

for name, rms in pulsars:
    base = base_name(name)
    if base in seen_base:
        continue
    seen_base.add(base)

    if base not in DISTANCES:
        skipped.append((name, "no distance in table"))
        continue

    radec = get_radec(base)
    if radec is None:
        skipped.append((name, "no par file / coords found"))
        continue
    ra, dec = radec

    if dec < -30:
        skipped.append((name, f"dec={dec:.1f} outside Bayestar footprint"))
        continue

    dist_kpc, method = DISTANCES[base]
    coord = SkyCoord(ra*u.deg, dec*u.deg, distance=dist_kpc*u.kpc, frame='icrs')
    try:
        reddening = float(bayestar(coord, mode='median'))
    except Exception as e:
        skipped.append((name, f"bayestar query failed: {e}"))
        continue

    results.append((base, rms, dist_kpc, method, reddening))
    print(f"{base}: DM_RMS={rms:.2e}, dist={dist_kpc} kpc ({method}), E(g-r)={reddening:.4f}")

with open("dust_reddening_all76.txt", "w") as f:
    f.write("# Pulsar  DM_RMS  Distance_kpc  DistMethod  E(g-r)\n")
    for base, rms, dist_kpc, method, reddening in results:
        f.write(f"{base} {rms:.6e} {dist_kpc} {method} {reddening:.4f}\n")

with open("dust_skipped.txt", "w") as f:
    f.write("# Pulsar  Reason\n")
    for name, reason in skipped:
        f.write(f"{name}\t{reason}\n")

print(f"\nDone. {len(results)} pulsars with dust values -> dust_reddening_all76.txt")
print(f"{len(skipped)} skipped -> dust_skipped.txt (see reasons)")