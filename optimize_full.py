"""Optimize the complete 14-wire 300 MHz antenna with 4NEC2's NEC2 engine.

Each trial contains the driven element, all nine corner-reflector rods, and all
four directors. Only physical dimensions change; no reduced surrogate model is
used. Run with ``python optimize_full.py`` from any directory.
"""

from __future__ import annotations

import csv
import math
import subprocess
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution


ROOT = Path(__file__).resolve().parent
ENGINE = Path(r"C:\4nec2\exe\nec2dxs1K5.exe")
OUT = ROOT / "optimization.csv"

# Reflector apex distance and length; four director lengths and four gaps.
BASE = np.array([
    0.161456, 0.446935,
    0.4350, 0.4215, 0.4185, 0.4245,
    0.1525, 0.3150, 0.3400, 0.3650,
])
BOUNDS = [
    (0.12, 0.24), (0.42, 0.50),
    (0.38, 0.46), (0.38, 0.46), (0.38, 0.46), (0.38, 0.46),
    (0.12, 0.38), (0.18, 0.42), (0.18, 0.42), (0.18, 0.42),
]
NAMES = ["reflector_distance", "reflector_length",
         "director1_length", "director2_length", "director3_length", "director4_length",
         "gap1", "gap2", "gap3", "gap4"]
LOCK = threading.Lock()
CACHE: dict[tuple[float, ...], dict[str, float]] = {}
BEST: dict[str, object] = {"score": -math.inf, "x": BASE.copy()}


def deck(x: np.ndarray, freq_mhz: float = 300.0) -> str:
    d, lr = float(x[0]), float(x[1])
    ls = list(map(float, x[2:6]))
    gaps = list(map(float, x[6:10]))
    lines = [
        "CM NHE2484 full 14-wire optimization candidate",
        "CM 8 mm diameter, free space, 300 MHz, complete 120 degree corner",
        "CE",
        "GW 1 31 0 0 -0.231300 0 0 0.231300 0.004",
    ]
    for tag, k in enumerate(range(-4, 5), 2):
        xx = -d - abs(k) * 0.075
        yy = k * 0.15 * math.sin(math.pi / 3)
        lines.append(
            f"GW {tag} 31 {xx:.6f} {yy:.6f} {-lr/2:.6f} "
            f"{xx:.6f} {yy:.6f} {lr/2:.6f} 0.004"
        )
    pos = 0.0
    for tag, (length, gap) in enumerate(zip(ls, gaps), 11):
        pos += gap
        lines.append(
            f"GW {tag} 31 {pos:.6f} 0 {-length/2:.6f} "
            f"{pos:.6f} 0 {length/2:.6f} 0.004"
        )
    lines += [
        "GE 0", "EK 0", "EX 0 1 16 0 1 0",
        f"FR 0 1 0 0 {freq_mhz:.6f} 0", "RP 0 1 2 1000 90 0 0 180", "EN",
    ]
    return "\r\n".join(lines) + "\r\n"


def evaluate(x: np.ndarray, freq_mhz: float = 300.0) -> dict[str, float]:
    with tempfile.TemporaryDirectory(prefix="nec_full_") as folder:
        inp = Path(folder) / "trial.inp"
        out = Path(folder) / "trial.out"
        inp.write_bytes(deck(x, freq_mhz).encode("ascii"))
        proc = subprocess.run(
            [str(ENGINE)], input=f"{inp}\n{out}\n", text=True,
            capture_output=True, cwd=ENGINE.parent, timeout=60,
        )
        if proc.returncode or not out.is_file():
            raise RuntimeError(f"NEC failed: {proc.returncode} {proc.stdout[-300:]}")
        text = out.read_text(encoding="latin1")
    if "ERROR" in text.upper():
        raise RuntimeError("NEC output contains an error")
    import re
    z = None
    gains = []
    in_input = in_pattern = False
    for line in text.splitlines():
        if "ANTENNA INPUT PARAMETERS" in line:
            in_input, in_pattern = True, False
            continue
        if in_input and re.match(r"^\s*1\s+\d+\s", line):
            vals = re.findall(r"[+-]?\d+\.\d+E[+-]\d+", line)
            if len(vals) >= 6:
                z = complex(float(vals[4]), float(vals[5]))
                in_input = False
        if "RADIATION PATTERNS" in line:
            in_pattern = True
            continue
        if in_pattern:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    theta, phi, gain = float(parts[0]), float(parts[1]), float(parts[4])
                except ValueError:
                    continue
                if theta == 90 and phi in (0, 180):
                    gains.append(gain)
    if z is None or len(gains) != 2:
        raise RuntimeError("Incomplete NEC impedance or pattern output")
    gamma = abs((z - 50) / (z + 50))
    swr = (1 + gamma) / (1 - gamma)
    mismatch = 10 * math.log10(1 - gamma * gamma)
    return {
        "r_ohm": z.real, "x_ohm": z.imag,
        "swr": swr, "s11_db": 20 * math.log10(gamma),
        "gain_dbi": gains[0], "realized_gain_dbi": gains[0] + mismatch,
        "fb_db": gains[0] - gains[1],
    }


def objective(x: np.ndarray) -> float:
    key = tuple(round(float(v), 6) for v in x)
    with LOCK:
        if key in CACHE:
            return -CACHE[key]["score"]
    try:
        result = evaluate(np.asarray(key))
        # Realized gain already includes mismatch. Preserve usable matching and
        # front-to-back performance while searching for a higher forward gain.
        penalty = 3.0 * max(0.0, result["swr"] - 2.0)
        penalty += 0.5 * max(0.0, 8.0 - result["fb_db"])
        score = result["realized_gain_dbi"] - penalty
    except Exception:
        result = {k: float("nan") for k in
                  ("r_ohm", "x_ohm", "swr", "s11_db", "gain_dbi", "realized_gain_dbi", "fb_db")}
        score = -1000.0
    result["score"] = score
    with LOCK:
        CACHE[key] = result
        with OUT.open("a", newline="", encoding="utf-8") as handle:
            csv.writer(handle).writerow([len(CACHE), *key, *[result[k] for k in
                ("r_ohm", "x_ohm", "swr", "s11_db", "gain_dbi", "realized_gain_dbi", "fb_db", "score")]])
        if score > BEST["score"]:
            BEST.update(score=score, x=np.asarray(key), result=result)
            print(f"best eval={len(CACHE)} score={score:.4f} realized={result['realized_gain_dbi']:.4f} "
                  f"SWR={result['swr']:.3f} FB={result['fb_db']:.2f}", flush=True)
    return -score


def main() -> None:
    if not ENGINE.is_file():
        raise FileNotFoundError(ENGINE)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerow(["eval", *NAMES, "r_ohm", "x_ohm", "swr", "s11_db",
                                      "gain_dbi", "realized_gain_dbi", "fb_db", "score"])
    base = evaluate(BASE)
    print(f"baseline realized={base['realized_gain_dbi']:.4f} SWR={base['swr']:.3f} FB={base['fb_db']:.2f}", flush=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        result = differential_evolution(
            objective, BOUNDS, x0=BASE, seed=20261024,
            maxiter=14, popsize=8, mutation=(0.5, 1.0),
            recombination=0.8, workers=pool.map, updating="deferred",
            polish=True, tol=0.002,
            callback=lambda x, convergence: print(
                f"generation best={-objective(x):.4f} convergence={convergence:.5f}", flush=True),
        )
    final_x = np.asarray(BEST["x"])
    final = evaluate(final_x)
    print("evaluations", len(CACHE), "optimizer", result.message, flush=True)
    print("final x", [round(float(v), 6) for v in final_x], flush=True)
    print("final", final, flush=True)
    import json
    (ROOT / "optimization_result.json").write_text(
        json.dumps({"baseline": base, "best": final,
                    "variables": dict(zip(NAMES, map(float, final_x))),
                    "evaluations": len(CACHE), "algorithm": "differential_evolution",
                    "constraints": {"swr_max": 2.0, "fb_min_db": 8.0}}, indent=2) + "\n",
        encoding="utf-8",
    )
    (ROOT / "best_candidate.nec").write_bytes(deck(final_x).encode("ascii"))


if __name__ == "__main__":
    main()
