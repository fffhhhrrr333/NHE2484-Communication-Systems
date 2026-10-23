# Complete-model 4NEC2 antenna optimization

This folder contains reproducible optimization code for Objective 3 of the NHE2484 antenna assignment. Every candidate sent to the electromagnetic solver is a **complete 14-wire NEC2 model**. The search does not use an electromagnetic surrogate, a reduced wire model, or a simplified parameter model.

The `models/` directory contains the delivered `E1.nec`, `E2.nec`, and final `E3.nec` files. `results/Experiment3/` contains the final 300 MHz summary, frequency sweep, director parameters, and two complete-model parameter sweeps. These files are snapshots of the delivered result; rerunning `publish_full.py` generates fresh copies in its selected output directory.

## Physical model and performance measure

| Item | Model setting |
|---|---|
| Driven element | Centre-fed cylindrical dipole, length 0.4626 m, on the Z axis at X = Y = 0 |
| Reflector | Nine parallel cylindrical rods forming a 120° corner behind the driven element; adjacent rods are 0.15 m apart along each wing |
| Directors | Four parallel cylindrical rods ahead of the driven element, each with its own length and adjacent gap |
| All rods | Diameter 8 mm (NEC radius 0.004 m); 31 segments per rod |
| Excitation and environment | 1 V source on segment 16 of the driven element; free space; 300 MHz design point |
| Direction and reference | Forward +X, backward −X; 50 Ω reference impedance |

Ten physical dimensions are varied: corner apex distance, common reflector rod length, four director lengths, and four consecutive director gaps. The driven length, rod diameter, nine-rod spacing, corner angle, and segment count stay fixed. The optimization begins from the four-director antenna developed from the assignment's initial 0.42 m director length and 0.25 m spacing.

The target is maximum **realized forward gain at 300 MHz** in dBi, subject to SWR ≤ 2 and forward-to-back ratio ≥ 8 dB. NEC2 reports gain referenced to accepted power. The scripts compute realized gain as `gain_dbi + 10 log10(1 − |Γ|²)`, with `Γ = (Zin − 50)/(Zin + 50)`. The two extra limits are practical quality gates used in this search; they are not stated as numerical limits in the assignment brief. An optimization run is a bounded numerical search, not a proof of a global maximum.

## Requirements

- Windows with 4NEC2 installed; the current engine path is `C:\4nec2\exe\nec2dxs1K5.exe`. Change `ENGINE` in `optimize_full.py` if 4NEC2 is elsewhere.
- Python 3.12 or later with NumPy and SciPy. Verified here with Python 3.12.5, NumPy 2.1.2, and SciPy 1.17.0. Install dependencies with `python -m pip install -r requirements.txt`.
- Each NEC2 engine process receives a unique temporary input/output pair. Four processes may run concurrently in the global searches.

## Reproduce the search

From this directory in PowerShell:

```powershell
python -m pip install -r requirements.txt
python optimize_full.py
python refine_full.py
python advanced_full.py
```

`optimize_full.py` performs a seeded broad differential-evolution search. `refine_full.py` restarts differential evolution around its best feasible point and then performs bounded Nelder–Mead refinement. `advanced_full.py` broadens the physical bounds, seeds the population with Sobol samples and diverse archived candidates, and finishes with another local refinement. All three stages call the full NEC2 solver for each candidate; the population search does not use an interpolated objective.

Each stage writes a CSV evaluation trace and a JSON result file:

| Stage | Trace | Result |
|---|---|---|
| Broad search | `optimization.csv` | `optimization_result.json` |
| Local restart | `refinement.csv` | `refinement_result.json` |
| Expanded search | `advanced_search.csv` | `advanced_result.json` |

The random seeds are fixed in the scripts. Different CPU scheduling can change the order of trace rows, and minor numerical or software-version differences can change the last digits.

To publish a candidate, select its result JSON explicitly:

```powershell
python publish_full.py --source refinement_result.json
# Or, if advanced_result.json improves the verified incumbent:
python publish_full.py --source advanced_result.json
```

The publisher first verifies the **actual delivered deck** independently with NEC2 and 361-point elevation and azimuth requests. It checks that the result improves the original 12.5331 dBi four-director baseline while meeting SWR ≤ 2 and F/B ≥ 8 dB. It then writes `E3.nec`, `Experiment3/results.json`, `frequency_sweep.csv`, `Director_Parameters.csv`, and 41-point `director_length_sweep.csv` and `director_spacing_sweep.csv`. Each sweep point again evaluates all 14 rods. The scan CSV fields are `offset_m, realized_gain_dbi, swr, fb_db`; the common offset is applied to all four director lengths or all four gaps while preserving their individual differences.

When this directory is `E:\nec\code`, the default publication destination is the adjacent `E:\nec` directory. In an independent checkout, the default is `output/` inside the checkout. Use `--output-root PATH` to select another destination. The publisher does not create or upload a GitHub repository.

## Verified result and numerical limits

The current delivered `E3.nec` has been evaluated separately from the two-direction optimization deck. At 300 MHz it gives **13.0500 dBi realized gain**, **1.000013 SWR**, and **10.22 dB F/B**. The 2:1 SWR band containing 300 MHz is **297.15352–301.68335 MHz**, or **4.52982 MHz**. The original E3 result before this search was 12.5331 dBi, 1.17065 SWR, and 8.31 dB F/B. The intermediate 12.9700 dBi design had a 5.91149 MHz 2:1 bandwidth; the extra 0.08 dB of realized gain in the final design comes with a narrower band.

The final corner apex distance is 0.1028049 m and the common reflector length is 0.5103174 m. Director lengths are `[0.4379216, 0.4255981, 0.4206466, 0.4263909]` m; successive gaps are `[0.1826665, 0.2798626, 0.3537916, 0.3556554]` m. `advanced_result.json` records the search result, while `delivery_audit.json` records the independently parsed delivered deck.

The optimizer found an input impedance near `50 + j0 Ω`, so its calculated S11 at exactly 300 MHz is about −104 dB. This is an ideal-model numerical point, sensitive to coordinate rounding, manufacturing tolerances, material properties, and the frequency grid. It should not be interpreted as an achievable measured S11. The finite 2:1 SWR bandwidth and the complete frequency sweep are more useful measures of matching robustness.

The model assumes perfect conductors and free space. Feed hardware, support structure, nearby objects, and fabrication tolerances are outside this simulation. The forward-to-back number compares the exact +X and −X directions, not an angular average over a rear hemisphere.

## License

The code and accompanying model/data files in this repository are released under the MIT License; see `LICENSE`.
