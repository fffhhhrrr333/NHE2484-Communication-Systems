"""Broader, archive-seeded global search on the complete 14-wire NEC2 model.

The only evaluator is the bundled NEC2 engine. Every population member contains
all 14 physical rods and all 31 segments per rod. No response-surface or
reduced electromagnetic model is evaluated as a candidate.
"""

from __future__ import annotations

import csv
import json
import math
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.stats import qmc

import optimize_full as full


ROOT = Path(__file__).resolve().parent
TRACE = ROOT / "advanced_search.csv"
RESULT = ROOT / "advanced_result.json"
SOURCE = ROOT / "refinement_result.json"
NAMES = full.NAMES
BOUNDS = np.array([
    (.10, .26), (.43, .56),
    (.38, .47), (.36, .47), (.36, .47), (.36, .47),
    (.12, .36), (.17, .50), (.17, .50), (.17, .50),
])
SPAN = BOUNDS[:,1] - BOUNDS[:,0]
LOCK = threading.Lock()
CACHE: dict[tuple[float,...], dict] = {}
BEST: dict = {"gain": -math.inf}


def evaluate_point(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    key = tuple(round(float(t), 7) for t in x)
    with LOCK:
        if key in CACHE:
            return -CACHE[key]["score"]
    try:
        r = full.evaluate(np.asarray(key))
        # Feasibility is maintained as a true selection gate for the archived
        # winner. Penalties only guide population members back to feasibility.
        score = r["realized_gain_dbi"] - 8*max(0.,r["swr"]-2.) - 4*max(0.,8.-r["fb_db"])
    except Exception:
        r = {k:float("nan") for k in ("realized_gain_dbi","gain_dbi","swr","fb_db","r_ohm","x_ohm")}
        score = -1000.
    r["score"] = score
    with LOCK:
        CACHE[key] = r
        with TRACE.open("a", newline="", encoding="utf-8") as fp:
            csv.writer(fp).writerow([len(CACHE), *key, r["realized_gain_dbi"],r["gain_dbi"],r["swr"],r["fb_db"],score])
        if r["swr"] <= 2 and r["fb_db"] >= 8 and r["realized_gain_dbi"] > BEST["gain"]:
            BEST.update(x=np.asarray(key), gain=r["realized_gain_dbi"], metrics=r)
            print(f"best {len(CACHE):5d}: realized={r['realized_gain_dbi']:.5f} "
                  f"SWR={r['swr']:.3f} FB={r['fb_db']:.2f}", flush=True)
    return -score


def initial_population(current: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(20261025)
    sobol = qmc.Sobol(d=10, scramble=True, seed=20261025).random_base2(7)
    pop = qmc.scale(sobol, BOUNDS[:,0], BOUNDS[:,1])
    # Half the population explores the already good basin at several scales;
    # the other half is space filling across all allowed physical dimensions.
    for i in range(64):
        scale = [.025,.06,.12,.20][i % 4]
        pop[i] = np.clip(current + rng.normal(0,scale,size=10)*SPAN, BOUNDS[:,0],BOUNDS[:,1])
    pop[0] = current
    # Seed diverse high-performing archive members where available.
    archive = ROOT/"refinement.csv"
    if archive.exists():
        rows = list(csv.DictReader(archive.open(newline="",encoding="utf-8")))
        rows.sort(key=lambda r: float(r["realized_gain_dbi"]), reverse=True)
        chosen = []
        for row in rows:
            if float(row["swr"])>2 or float(row["fb_db"])<8:
                continue
            x = np.array([float(row[n]) for n in NAMES])
            if np.any(x < BOUNDS[:,0]) or np.any(x > BOUNDS[:,1]):
                continue
            if all(np.linalg.norm((x-y)/SPAN)>.035 for y in chosen):
                chosen.append(x)
            if len(chosen)>=12:
                break
        for i,x in enumerate(chosen):
            pop[i+1] = x
    return pop


def main() -> None:
    current = json.loads(SOURCE.read_text(encoding="utf-8"))
    x0 = np.array([current["variables"][n] for n in NAMES])
    with TRACE.open("w",newline="",encoding="utf-8") as fp:
        csv.writer(fp).writerow(["eval",*NAMES,"realized_gain_dbi","gain_dbi","swr","fb_db","score"])
    evaluate_point(x0)
    with ThreadPoolExecutor(max_workers=4) as pool:
        de = differential_evolution(
            evaluate_point, BOUNDS, init=initial_population(x0), x0=x0,
            strategy="randtobest1bin", mutation=(.35,.95), recombination=.85,
            seed=20261025, maxiter=24, polish=False, tol=.001,
            workers=pool.map, updating="deferred",
            callback=lambda x,convergence: print(
                f"generation: best={BEST['gain']:.5f}, evaluated={len(CACHE)}, convergence={convergence:.4f}",
                flush=True),
        )
    print("global", de.message, "evaluated",len(CACHE),flush=True)
    # A constrained local pass resolves small variations in the strongest basin.
    start = np.asarray(BEST["x"])
    simplex = np.tile(start,(11,1))
    increments = np.array([.003,.002,.0015,.0015,.0015,.0015,.003,.004,.004,.004])
    for j in range(10):
        simplex[j+1,j] = min(BOUNDS[j,1],start[j]+increments[j])
    local = minimize(evaluate_point,start,method="Nelder-Mead",bounds=BOUNDS,
                     options={"initial_simplex":simplex,"maxfev":650,
                              "xatol":2e-6,"fatol":1e-5})
    print("local",local.message,"evaluated",len(CACHE),flush=True)
    x=np.asarray(BEST["x"])
    best=full.evaluate(x)
    result={"best":best,"variables":dict(zip(NAMES,map(float,x))),
            "evaluations":len(CACHE),"method":"archive-seeded Sobol DE (expanded bounds) plus Nelder-Mead",
            "bounds":BOUNDS.tolist(),"reference_realized_gain_dbi":current["best"]["realized_gain_dbi"]}
    RESULT.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    (ROOT/"advanced_candidate.nec").write_bytes(full.deck(x).encode("ascii"))
    print("FINAL",json.dumps(result),flush=True)


if __name__=="__main__":
    main()
