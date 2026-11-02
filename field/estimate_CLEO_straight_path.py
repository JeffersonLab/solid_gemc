"""Run: python3 estimate_CLEO_straight_path.py [9 10 11] --condition pvdis|pvdis_baffle|sidis|jpsi_full|jpsi_lite|jpsi_la

Read maps beside this script and write CLEOvN_CONDITION_straight_path.csv.
Defaults: v9, pvdis. Multiple maps also write CLEO_CONDITION_straight_path_compare.csv.
PVDIS: target Z=0 cm, last plane Z=400 cm, angles 22/35 deg.
PVDIS baffle: target Z=0 cm, last plane Z=200 cm, angles 22/35 deg.
SIDIS: target Z=-350 cm, last plane Z=200 cm, angles 8/13 deg.
J/psi full (CLEOv9): target Z=-315 cm, last plane Z=90 cm, angles 8/13/14/24 deg.
J/psi lite (CLEOv10): target Z=-315 cm, last plane Z=190 cm, angles 8/13/14/24 deg.
J/psi LA: target Z=-315 cm, last plane Z=-70 cm, angles 8/13/14/24 deg.
Both use target R=0, momenta 2/5 GeV/c, and angles measured from +Z.
Signed bending uses electrons (q=-1), with azimuthal deflection q*(t cross B).
Ignore material, energy loss, and changes to the reference trajectory.
"""

import argparse
import csv
import math
from pathlib import Path

import numpy as np


def interpolate(field, r, z, rgrid, zgrid):
    if np.any((r < rgrid[0]) | (r > rgrid[-1]) | (z < zgrid[0]) | (z > zgrid[-1])):
        raise ValueError("Reference path leaves the field map")
    i = np.clip(np.searchsorted(rgrid, r, side="right") - 1, 0, len(rgrid) - 2)
    j = np.clip(np.searchsorted(zgrid, z, side="right") - 1, 0, len(zgrid) - 2)
    u = (r - rgrid[i]) / (rgrid[i + 1] - rgrid[i])
    v = (z - zgrid[j]) / (zgrid[j + 1] - zgrid[j])
    return ((1-u)*(1-v)*field[i,j] + u*(1-v)*field[i+1,j]
            + (1-u)*v*field[i,j+1] + u*v*field[i+1,j+1])


def integrate(angle_deg, step_cm, rgrid, zgrid, br, bz, target_z_cm, endpoint_z_cm):
    angle = math.radians(angle_deg)
    z = np.linspace(target_z_cm, endpoint_z_cm, round((endpoint_z_cm-target_z_cm) / step_cm) + 1)
    r = (z-target_z_cm) * math.tan(angle)
    s_m = (z-target_z_cm) / (100 * math.cos(angle))
    # Defined even at the target by the specified initial track direction.
    bperp_t = (interpolate(br,r,z,rgrid,zgrid) * math.cos(angle)
               - interpolate(bz,r,z,rgrid,zgrid) * math.sin(angle)) * 1e-4
    positive = np.trapz(np.maximum(bperp_t,0), s_m)
    negative = np.trapz(np.minimum(bperp_t,0), s_m)
    return {
        "angle_deg": angle_deg,
        "target_z_cm": target_z_cm,
        "target_r_cm": 0,
        "endpoint_z_cm": endpoint_z_cm,
        "endpoint_r_cm": r[-1],
        "path_length_m": s_m[-1],
        "signed_integral_Tm": positive + negative,
        "absolute_integral_Tm": positive - negative,
        "positive_integral_Tm": positive,
        "negative_integral_Tm": negative,
        "weighted_integral_Tm2": np.trapz((s_m[-1]-s_m)*bperp_t, s_m),
    }


def main(version, condition):
    target_z_cm, endpoint_z_cm, angles = {
        "pvdis": (0,400,(22,35)),
        "pvdis_baffle": (0,200,(22,35)),
        "sidis": (-350,200,(8,13)),
        "jpsi_full": (-315,90,(8,13,14,24)),
        "jpsi_lite": (-315,190,(8,13,14,24)),
        "jpsi_la": (-315,-70,(8,13,14,24)),
    }[condition]
    source = Path(__file__).with_name(f"solenoid_CLEOv{version}.dat")
    with source.open() as stream:
        for header_lines, line in enumerate(stream, 1):
            if "</mfield>" in line:
                break
        else:
            raise ValueError("Missing GEMC header terminator")
    data = np.loadtxt(source, skiprows=header_lines)
    rgrid, zgrid = np.unique(data[:,0]), np.unique(data[:,1])
    assert len(data) == len(rgrid) * len(zgrid)
    assert np.array_equal(data[:,0], np.repeat(rgrid,len(zgrid)))
    assert np.array_equal(data[:,1], np.tile(zgrid,len(rgrid)))
    br = data[:,2].reshape(len(rgrid),len(zgrid))
    bz = data[:,3].reshape(len(rgrid),len(zgrid))
    # Verify interpolation at map nodes, including both boundary corners.
    indices = np.linspace(0,len(data)-1,1001,dtype=int)
    for field, column in ((br,2),(bz,3)):
        assert np.allclose(interpolate(field,data[indices,0],data[indices,1],rgrid,zgrid),
                           data[indices,column], rtol=0, atol=1e-10)

    charge = -1
    rows = []
    for angle in angles:
        coarse = integrate(angle,0.25,rgrid,zgrid,br,bz,target_z_cm,endpoint_z_cm)
        fine = integrate(angle,0.125,rgrid,zgrid,br,bz,target_z_cm,endpoint_z_cm)
        for key in ("signed_integral_Tm","absolute_integral_Tm","weighted_integral_Tm2"):
            assert abs(fine[key]-coarse[key]) <= 1e-4 * max(abs(fine[key]),1e-6), key
        print(f"{condition}, CLEOv{version}, {angle} deg: signed={fine['signed_integral_Tm']:.6f} T m, "
              f"absolute={fine['absolute_integral_Tm']:.6f} T m; refinement verified")
        for momentum in (2,5):
            row = {"condition": condition, "map_version": version, **fine, "momentum_GeV_c": momentum, "charge_e": charge,
                   "deflection_deg": math.degrees(0.299792458*charge*fine['signed_integral_Tm']/momentum),
                   "azimuthal_displacement_cm": 100*0.299792458*charge*fine['weighted_integral_Tm2']/momentum}
            rows.append(row)
            print(f"  {momentum} GeV/c: angle={row['deflection_deg']:.4f} deg, "
                  f"displacement={row['azimuthal_displacement_cm']:.4f} cm")
    output = source.with_name(f"CLEOv{version}_{condition}_straight_path.csv")
    with output.open("w",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {output.name}")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("versions",nargs="*",type=int,choices=(9,10,11),default=[9])
    parser.add_argument("--condition",choices=("pvdis","pvdis_baffle","sidis","jpsi_full","jpsi_lite","jpsi_la"),default="pvdis")
    args = parser.parse_args()
    versions = args.versions
    results = [row for version in versions for row in main(version,args.condition)]
    if len(versions)>1:
        output = Path(__file__).with_name(f"CLEO_{args.condition}_straight_path_compare.csv")
        with output.open("w",newline="") as stream:
            writer = csv.DictWriter(stream,fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"Saved {output.name}")
