# 2026-Wigner-I-Hardly-Know-Her

Python port and workflow tooling for processing Flash-DSC style exports used in Wigner/indium analysis.

## What is in this repo

1. Legacy MATLAB analysis scripts.
2. New Python pipeline with CLI for parsing, baseline workflow plotting, and indium mass estimation.
3. Example outputs produced from the available baseline dataset.

## File map

### MATLAB source files

1. [find_lefts_rights_mc.m](find_lefts_rights_mc.m): Monte Carlo boundary detection for peak integration.
2. [spline_integral_orig.m](spline_integral_orig.m): baseline-corrected integration for enthalpy-like area.
3. [hook_MC_QAD.m](hook_MC_QAD.m): hook-method tangent/intersection workflow.
4. [indium_interface_experiment.m](indium_interface_experiment.m): multi-file interface plotting with optional baseline subtraction.
5. [indium_example.m](indium_example.m): full indium mass workflow orchestration.

### Python pipeline files

1. [run_pipeline.py](run_pipeline.py): top-level executable wrapper.
2. [wigner_pipeline/cli.py](wigner_pipeline/cli.py): CLI subcommands.
3. [wigner_pipeline/io_parser.py](wigner_pipeline/io_parser.py): parser for multi-segment text exports.
4. [wigner_pipeline/methods.py](wigner_pipeline/methods.py): numerical kernels ported from MATLAB.
5. [wigner_pipeline/models.py](wigner_pipeline/models.py): segment dataclass model.
6. [wigner_pipeline/workflows/interface_experiment.py](wigner_pipeline/workflows/interface_experiment.py): plotting workflow.
7. [wigner_pipeline/workflows/indium_analysis.py](wigner_pipeline/workflows/indium_analysis.py): mass analysis workflow.
8. [requirements.txt](requirements.txt): runtime dependencies.

## Pipeline commands

Run from repository root.

1. Parser smoke test:

```powershell
C:/Users/rianc/AppData/Local/Python/pythoncore-3.14-64/python.exe run_pipeline.py parser-smoke --file cle_ind_baseline_25to650C_200kps_07292026.txt
```

2. Interface plotting workflow:

```powershell
C:/Users/rianc/AppData/Local/Python/pythoncore-3.14-64/python.exe run_pipeline.py run-interface-experiment --workspace .
```

3. Indium analysis workflow:

```powershell
C:/Users/rianc/AppData/Local/Python/pythoncore-3.14-64/python.exe run_pipeline.py run-indium-analysis --file cle_ind_baseline_25to650C_200kps_07292026.txt
```

## Current generated outputs

1. [outputs/interface_experiment/01_cle_ind_baseline_25to650C_200kps_07292026.png](outputs/interface_experiment/01_cle_ind_baseline_25to650C_200kps_07292026.png)
2. [outputs/indium_analysis/indium_analysis_summary.csv](outputs/indium_analysis/indium_analysis_summary.csv)
3. [outputs/indium_analysis/mass_vs_heating_rate.png](outputs/indium_analysis/mass_vs_heating_rate.png)
4. [outputs/indium_analysis/segment_00.png](outputs/indium_analysis/segment_00.png)
5. [outputs/indium_analysis/segment_01.png](outputs/indium_analysis/segment_01.png)
6. [outputs/indium_analysis/segment_02.png](outputs/indium_analysis/segment_02.png)
7. [outputs/indium_analysis/segment_03.png](outputs/indium_analysis/segment_03.png)

## Run in VS Code debug

Use the default interpreter already configured in this workspace.

### Option A: Debug by command line in integrated terminal

Run the same commands listed in "Pipeline commands" from the integrated terminal.

### Option B: Debug with launch profiles

Create or update .vscode/launch.json and add profiles like these:

```json
{
	"version": "0.2.0",
	"configurations": [
		{
			"name": "Parser smoke",
			"type": "python",
			"request": "launch",
			"program": "${workspaceFolder}/run_pipeline.py",
			"console": "integratedTerminal",
			"args": [
				"parser-smoke",
				"--file",
				"cle_ind_baseline_25to650C_200kps_07292026.txt"
			]
		},
		{
			"name": "Interface experiment",
			"type": "python",
			"request": "launch",
			"program": "${workspaceFolder}/run_pipeline.py",
			"console": "integratedTerminal",
			"args": [
				"run-interface-experiment",
				"--workspace",
				"."
			]
		},
		{
			"name": "Indium analysis",
			"type": "python",
			"request": "launch",
			"program": "${workspaceFolder}/run_pipeline.py",
			"console": "integratedTerminal",
			"args": [
				"run-indium-analysis",
				"--file",
				"cle_ind_baseline_25to650C_200kps_07292026.txt"
			]
		}
	]
}
```

Set breakpoints in:

1. [wigner_pipeline/workflows/indium_analysis.py](wigner_pipeline/workflows/indium_analysis.py)
2. [wigner_pipeline/methods.py](wigner_pipeline/methods.py)

Then run the matching launch profile.

## Why enthalpy vs hook mass can diverge (current run)

## Units and Formula Legend (indium_analysis_summary.csv)

The fields below are written by [wigner_pipeline/workflows/indium_analysis.py](wigner_pipeline/workflows/indium_analysis.py) into [outputs/indium_analysis/indium_analysis_summary.csv](outputs/indium_analysis/indium_analysis_summary.csv).

1. segment: zero-based loop index over analyzed segments.
2. seg_id: segment identifier parsed from source export (or fallback index).
3. direction: segment direction label, heating or cooling.
4. n_points: number of valid data points used for that segment.
5. heating_rate_K_per_s: heating rate estimate in K/s.
6. hcrys_mean: Monte Carlo mean integrated area from spline baseline correction, in mW*K.
7. hcrys_std: Monte Carlo standard deviation of hcrys, in mW*K.
8. hook_height: |h_end - h_start| from hook method, in mW.
9. mass_hcrys_ng: enthalpy-based mass estimate in ng.
10. mass_hook_ng: hook-based mass estimate in ng.
11. lr_center_C: detected left-inner transition boundary center in deg C.
12. rl_center_C: detected right-inner transition boundary center in deg C.

Mass conversion equations used in code:

1. Enthalpy path
2. $$m_{hcrys,ng} = \left|\frac{hcrys\_mean\cdot 10^{-3}}{HR\cdot h_{fusion}\cdot 10^{-9}}\right|$$
3. with $HR$ in K/s, $h_{fusion}=28.6$ J/g, and $hcrys\_mean$ in mW*K.

4. Hook path
5. $$m_{hook,ng} = \left|\frac{hook\_height\cdot 10^{-3}}{HR\cdot C_p\cdot 10^{-9}}\right|$$
6. with $C_p=0.23$ J/(g*K) and $hook\_height$ in mW.

Constant reference carried from the MATLAB workflow notes in [indium_example.m](indium_example.m):

1. Enthalpy of fusion: 3283 J/mol (28.6 J/g).
2. Cp at 300 K: 26.98 J/K/mol (0.23 J/K/g).
3. Source: https://akjournals.com/view/journals/10973/13/3/article-p419.xml

The current CSV indicates a large gap for this dataset:

1. mass_hcrys_ng about 16640.
2. mass_hook_ng about 1250.

See [outputs/indium_analysis/indium_analysis_summary.csv](outputs/indium_analysis/indium_analysis_summary.csv).

Without changing code, the most likely contributing sections are:

1. Integration region width and placement:
2. Boundary centers are around lr ≈120 C and rl ≈188 C in the current run, which is a broad window around the expected In peak region.
3. Source: [wigner_pipeline/workflows/indium_analysis.py](wigner_pipeline/workflows/indium_analysis.py) stores lr/rl centers and [wigner_pipeline/methods.py](wigner_pipeline/methods.py) computes bounds.

4. Absolute-area integration inflates area:
5. The integrated signal uses absolute value before integration, so all baseline residual contributes positively.
6. Source: [wigner_pipeline/methods.py](wigner_pipeline/methods.py) in spline_integral_orig where q_corr = abs(qc - baseline).

7. Hook region uses only pre-lr data:
8. Hook estimation is computed on seg.tr_c < mean(lr), so selected pre-peak window can strongly affect hook height.
9. Source: [wigner_pipeline/workflows/indium_analysis.py](wigner_pipeline/workflows/indium_analysis.py).

10. Different physical constants and signal semantics:
11. Enthalpy path uses h_fusion = 28.6, hook path uses heat_capacity = 0.23.
12. These constants are inherited from the MATLAB notes in [indium_example.m](indium_example.m):
13. Enthalpy of fusion 3283 J/mol (28.6 J/g) and Cp at 300 K 26.98 J/K/mol (0.23 J/K/g).
14. Source: https://akjournals.com/view/journals/10973/13/3/article-p419.xml
15. Even with correct constants, mismatched region selection/baseline fit can dominate and create large ratio differences.

## Practical checks to run next (no code change)

1. Verify boundary reasonableness against traces:
2. Open [outputs/indium_analysis/segment_00.png](outputs/indium_analysis/segment_00.png) to [outputs/indium_analysis/segment_03.png](outputs/indium_analysis/segment_03.png) and check whether lr/rl bracket only the intended melting event.

3. Compare LR/RL centers to expected transition window:
4. Review lr_center_C and rl_center_C in [outputs/indium_analysis/indium_analysis_summary.csv](outputs/indium_analysis/indium_analysis_summary.csv).

5. Check sensitivity to Monte Carlo region realization:
6. Re-run with different seeds and observe stability of hcrys_mean vs hook_height.

7. Check baseline fit quality qualitatively:
8. In debugger, inspect left/right subsets feeding polynomial fits in spline_integral_orig and whether they represent true local baselines rather than transition shoulders.

9. Check heating-rate estimate consistency:
10. Confirm each heating segment is truly 200 K/s and no timing irregularity skews conversion from integrated area to mass.

## Notes

1. Current pipeline is parity-first with MATLAB logic where practical.
2. One compatibility update was applied for Python 3.14 runtime support: np.trapezoid is used for numerical integration.
