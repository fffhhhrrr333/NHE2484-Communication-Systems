"""Publish and audit the optimized complete E3 antenna.

Run after refine_full.py. Writes E3.nec with full 361-point 4NEC2 pattern
requests, frequency data, complete-model director scans, and results.json.
"""

from __future__ import annotations

import csv
import argparse
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

import optimize_full as full


CODE = Path(__file__).resolve().parent
ROOT = CODE.parent if (CODE.parent / "E3.nec").exists() else CODE / "output"
DEST = ROOT / "Experiment3"
OUTPUT = ROOT / "E3.nec"


def delivered_deck(x: np.ndarray) -> str:
    raw = full.deck(x).splitlines()
    assert raw[-2] == "RP 0 1 2 1000 90 0 0 180"
    raw[:3] = [
        "CM NHE2484 objective 3: optimized complete 14-wire Yagi/corner antenna",
        "CM Driven dipole plus nine corner-reflector rods and four directors.",
        "CM Free space; all rods 8 mm diameter; 300 MHz; 50 ohm reference.",
        "CE",
    ]
    raw[-2:-1] = ["RP 0 361 1 1000 -180 0 1 0",
                  "RP 0 1 361 1000 90 -180 0 1"]
    return "\r\n".join(raw) + "\r\n"


def run_exact_deck(deck: str) -> str:
    with tempfile.TemporaryDirectory(prefix="nec_delivery_audit_") as tmp:
        inp, out = Path(tmp)/"verify.inp", Path(tmp)/"verify.out"
        inp.write_bytes(deck.encode("ascii"))
        proc = subprocess.run([str(full.ENGINE)], input=f"{inp}\n{out}\n",
                              text=True, capture_output=True, cwd=full.ENGINE.parent, timeout=90)
        if proc.returncode or not out.exists():
            raise RuntimeError(f"NEC2 delivery verification failed: {proc.returncode} {proc.stdout[-300:]}")
        text = out.read_text(encoding="latin1")
    if "ERROR" in text.upper():
        raise RuntimeError("NEC2 delivery output contains an ERROR")
    return text


def independent_metrics(output: str) -> dict:
    """Separate parser for the actual delivered deck's 361-point patterns."""
    import re
    lines = output.splitlines()
    z = None
    elevation = []
    azimuth = []
    mode = None
    for line in lines:
        if "ANTENNA INPUT PARAMETERS" in line:
            mode = "input"
        elif "RADIATION PATTERNS" in line:
            mode = "pattern"
            if elevation:
                mode = "azimuth"
        elif mode == "input" and re.match(r"^\s*1\s+\d+\s", line):
            nums = re.findall(r"[+-]?\d+\.\d+E[+-]\d+", line)
            if len(nums) >= 6:
                z = complex(float(nums[4]), float(nums[5]))
                mode = None
        elif mode in ("pattern", "azimuth"):
            parts = line.split()
            if len(parts) >= 5:
                try:
                    theta, phi, gain = float(parts[0]), float(parts[1]), float(parts[4])
                except ValueError:
                    continue
                if mode == "pattern" and phi == 0.0 and -180 <= theta <= 180:
                    elevation.append((theta, gain))
                elif mode == "azimuth" and theta == 90.0 and -180 <= phi <= 180:
                    azimuth.append((phi, gain))
    if z is None or len(elevation) != 361 or len(azimuth) != 361:
        raise RuntimeError(f"Delivery parser incomplete: z={z}, elevation={len(elevation)}, azimuth={len(azimuth)}")
    elev = dict(elevation)
    az = dict(azimuth)
    forward, backward = elev[90.0], az[180.0]
    gamma = abs((z-50)/(z+50))
    return {"r_ohm": z.real, "x_ohm": z.imag,
            "swr": (1+gamma)/(1-gamma), "s11_db": 20*math.log10(gamma),
            "gain_dbi": forward, "realized_gain_dbi": forward + 10*math.log10(1-gamma*gamma),
            "fb_db": forward-backward, "elevation": elevation, "azimuth": azimuth}


def hpbw(pattern: list[tuple[float, float]]) -> float:
    a = dict(pattern)
    peak = a[90.0]
    target = peak - 3.0
    left = 90
    while a[float(left)] >= target and left > -180:
        left -= 1
    right = 90
    while a[float(right)] >= target and right < 180:
        right += 1
    def cross(i, j):
        yi, yj = a[float(i)], a[float(j)]
        return i + (target-yi)*(j-i)/(yj-yi)
    return cross(right-1,right) - cross(left,left+1)


def main() -> None:
    global ROOT, DEST, OUTPUT
    parser = argparse.ArgumentParser(description="Publish and independently verify a full-model E3 candidate")
    parser.add_argument("--source", default="refinement_result.json",
                        help="Result JSON within code/ (default: refinement_result.json)")
    parser.add_argument("--output-root", default=None,
                        help="Directory for E3.nec and Experiment3/ (default: adjacent assignment directory, or code/output)")
    args = parser.parse_args()
    source = (CODE / args.source).resolve()
    if source.parent != CODE or not source.is_file():
        parser.error("--source must name an existing JSON file in code/")
    if args.output_root:
        ROOT = Path(args.output_root).resolve()
        DEST = ROOT / "Experiment3"
        OUTPUT = ROOT / "E3.nec"
    DEST.mkdir(parents=True, exist_ok=True)
    result = json.loads(source.read_text(encoding="utf-8"))
    x = np.array([result["variables"][name] for name in full.NAMES])
    model = delivered_deck(x)
    metrics = independent_metrics(run_exact_deck(model))
    direct = full.evaluate(x)
    for key in ("swr", "realized_gain_dbi", "fb_db"):
        if abs(direct[key] - metrics[key]) > 0.02:
            raise RuntimeError(f"Delivery differs from optimization: {key} {direct[key]} vs {metrics[key]}")
    if metrics["realized_gain_dbi"] <= 12.5331 or metrics["swr"] > 2 or metrics["fb_db"] < 8:
        raise RuntimeError("Optimized delivery does not improve the constrained baseline")
    if OUTPUT.exists() and not (CODE/"E3_previous.nec").exists():
        shutil.copy2(OUTPUT, CODE/"E3_previous.nec")
    OUTPUT.write_bytes(model.encode("ascii"))
    # A frequency point is a full 14-wire NEC2 simulation, including mismatch.
    freqs = np.arange(240., 361., 1.)
    sweep = [(float(f), full.evaluate(x, float(f))) for f in freqs]
    with (DEST/"frequency_sweep.csv").open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(["Frequency (MHz)", "Rin (ohm)", "Xin (ohm)", "SWR (50 ohm)",
                         "S11 (dB)", "Forward gain (dBi)", "Realized gain (dBi)", "Front-to-back (dB)"])
        for f, r in sweep:
            writer.writerow([f, r["r_ohm"], r["x_ohm"], r["swr"], r["s11_db"],
                             r["gain_dbi"], r["realized_gain_dbi"], r["fb_db"]])
    def cross(lo, hi):
        return brentq(lambda f: full.evaluate(x, f)["swr"]-2., lo, hi, xtol=1e-4)
    central = next(j for j,(f,_) in enumerate(sweep) if f == 300.)
    start = central
    while start > 0 and sweep[start-1][1]["swr"] <= 2:
        start -= 1
    end = central
    while end < len(sweep)-1 and sweep[end+1][1]["swr"] <= 2:
        end += 1
    low = cross(sweep[start-1][0], sweep[start][0]) if start > 0 else sweep[0][0]
    high = cross(sweep[end][0], sweep[end+1][0]) if end < len(sweep)-1 else sweep[-1][0]
    # Preserve individual director differences while changing the common length
    # or each of the four adjacent gaps by the same amount.
    for kind, indices, offsets, path in [
        ("length", range(2,6), np.linspace(-.020,.020,41), DEST/"director_length_sweep.csv"),
        ("spacing", range(6,10), np.linspace(-.040,.040,41), DEST/"director_spacing_sweep.csv"),
    ]:
        with path.open("w", newline="", encoding="utf-8") as fp:
            writer = csv.writer(fp)
            writer.writerow(["offset_m", "realized_gain_dbi", "swr", "fb_db"])
            for offset in offsets:
                trial = x.copy()
                trial[list(indices)] += offset
                r = full.evaluate(trial)
                writer.writerow([float(offset), r["realized_gain_dbi"], r["swr"], r["fb_db"]])
    pos = np.cumsum(x[6:10])
    with (DEST/"Director_Parameters.csv").open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(["Director", "Length (m)", "Gap from previous element (m)", "X position (m)"])
        for j in range(4):
            writer.writerow([j+1, f"{x[j+2]:.6f}", f"{x[j+6]:.6f}", f"{pos[j]:.6f}"])
    prior_results = DEST / "results.json"
    out = json.loads(prior_results.read_text(encoding="utf-8")) if prior_results.exists() else {}
    best_swr = min(sweep, key=lambda row: row[1]["swr"])
    best_gain = max(sweep, key=lambda row: row[1]["realized_gain_dbi"])
    out.update({
        "swr_min": best_swr[1]["swr"], "swr_min_mhz": best_swr[0],
        "band_contains_300_mhz": bool(metrics["swr"] <= 2),
        "bw_low_mhz": low, "bw_high_mhz": high, "bw_mhz": high-low,
        "fractional_bw_percent": 100*(high-low)/((high+low)/2),
        "gain_peak_dbi": best_gain[1]["realized_gain_dbi"], "gain_peak_mhz": best_gain[0],
        "at_300": {"z": [metrics["r_ohm"], metrics["x_ohm"]],
                   "swr": metrics["swr"], "s11": metrics["s11_db"],
                   "gain": metrics["gain_dbi"], "realized_gain": metrics["realized_gain_dbi"],
                   "fb": metrics["fb_db"]},
        "distance_m": float(x[0]), "reflector_length_m": float(x[1]),
        "director_lengths_m": x[2:6].tolist(), "director_gaps_m": x[6:10].tolist(),
        "director_positions_m": pos.tolist(),
        "hpbw_elevation_deg": hpbw(metrics["elevation"]),
        "optimization": {"method": result["method"], "evaluations": result["evaluations"],
                         "full_model": "14 wires; 9 reflector rods; 4 directors; 31 segments per wire"},
    })
    (DEST/"results.json").write_text(json.dumps(out, indent=2)+"\n", encoding="utf-8")
    audit = {k:v for k,v in metrics.items() if k not in ("elevation","azimuth")}
    (CODE/"delivery_audit.json").write_text(json.dumps(audit,indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"audit":audit, "band_mhz":[low,high], "hpbw_deg":out["hpbw_elevation_deg"],
                      "director_lengths_m":x[2:6].tolist(), "director_gaps_m":x[6:10].tolist(),
                      "director_positions_m":pos.tolist()}, indent=2), flush=True)


if __name__ == "__main__":
    main()
