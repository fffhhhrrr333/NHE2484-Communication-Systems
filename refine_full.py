"""Reproducible local refinement of the complete 14-wire NEC2 antenna.

Run after optimize_full.py. Uses a bounded differential-evolution restart around
the previous feasible best and independent direct NEC2 evaluations. All trials
retain the driven element, nine reflector rods, and four directors.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

import optimize_full as full


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "optimization_result.json"
TRACE = ROOT / "refinement.csv"
RESULT = ROOT / "refinement_result.json"


def clamp(x: np.ndarray) -> np.ndarray:
    return np.array([np.clip(a, lo, hi) for a, (lo, hi) in zip(x, full.BOUNDS)])


def score(x: np.ndarray) -> float:
    x = clamp(np.asarray(x, dtype=float))
    if tuple(x) in CACHE:
        return -CACHE[tuple(x)]["score"]
    try:
        metrics = full.evaluate(x)
        penalty = 20 * max(0., metrics["swr"] - 2.)
        penalty += 10 * max(0., 8. - metrics["fb_db"])
        value = metrics["realized_gain_dbi"] - penalty
    except Exception:
        metrics = {key: float("nan") for key in ("swr", "fb_db", "realized_gain_dbi")}
        value = -1000.
    metrics["score"] = value
    CACHE[tuple(x)] = metrics
    with TRACE.open("a", newline="", encoding="utf-8") as fp:
        csv.writer(fp).writerow([len(CACHE), *x, metrics["realized_gain_dbi"], metrics["swr"], metrics["fb_db"], value])
    if metrics["swr"] <= 2. and metrics["fb_db"] >= 8. and value > BEST["score"]:
        BEST.update(x=x.copy(), metrics=metrics, score=value)
        print(f"refine eval={len(CACHE)} gain={value:.4f} SWR={metrics['swr']:.3f} FB={metrics['fb_db']:.2f}", flush=True)
    return -value


CACHE: dict[tuple[float, ...], dict] = {}
BEST: dict = {"score": -np.inf}


def main() -> None:
    prior = json.loads(INPUT.read_text(encoding="utf-8"))
    x0 = np.array([prior["variables"][name] for name in full.NAMES])
    with TRACE.open("w", newline="", encoding="utf-8") as fp:
        csv.writer(fp).writerow(["eval", *full.NAMES, "realized_gain_dbi", "swr", "fb_db", "score"])
    score(x0)
    # Tight global restart makes better use of evaluations after a broad first pass.
    margins = np.array([.035, .025, .015, .015, .015, .015, .04, .05, .05, .05])
    bounds = [(max(lo, xx-dd), min(hi, xx+dd)) for xx, dd, (lo, hi) in zip(x0, margins, full.BOUNDS)]
    result = differential_evolution(score, bounds, x0=x0, seed=20261025,
                                    maxiter=18, popsize=7, mutation=(.4, .9),
                                    recombination=.8, tol=.001, polish=False)
    print("DE restart", result.message, "best", BEST["score"], flush=True)
    # Deterministic coordinate/Nelder-Mead improvement at the best feasible point.
    x1 = np.array(BEST["x"])
    simplex = np.tile(x1, (len(x1)+1, 1))
    steps = [.002, .0015, .001, .001, .001, .001, .003, .003, .003, .003]
    for j, step in enumerate(steps):
        simplex[j+1,j] = np.clip(simplex[j+1,j] + step, *full.BOUNDS[j])
    local = minimize(score, x1, method="Nelder-Mead", bounds=full.BOUNDS,
                     options={"initial_simplex": simplex, "maxfev": 700,
                              "xatol": 2e-6, "fatol": 2e-4})
    print("local", local.message, "best", BEST["score"], flush=True)
    x = np.array(BEST["x"])
    metrics = full.evaluate(x)
    result_json = {"best": metrics, "variables": dict(zip(full.NAMES, x.tolist())),
                   "evaluations": len(CACHE), "method": "bounded DE restart plus Nelder-Mead",
                   "base_file": str(INPUT)}
    RESULT.write_text(json.dumps(result_json, indent=2) + "\n", encoding="utf-8")
    (ROOT / "best_candidate.nec").write_bytes(full.deck(x).encode("ascii"))
    print("FINAL", json.dumps(result_json), flush=True)


if __name__ == "__main__":
    main()
