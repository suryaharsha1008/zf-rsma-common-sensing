# Common-stream-only sensing in ZF-RSMA with imperfect CSIT

MATLAB code that reproduces the numerical results of

> Surya Harsha and Shankar Prakriya, "Feasibility of Common-Stream-Only Sensing in ZF-RSMA with Imperfect CSIT,"
> accepted in *IEEE Global Communications Conference (GLOBECOM)*, 2026.

## Problem
A multi-antenna base station serves K single-antenna users with downlink rate-splitting multiple access (RSMA) and
reuses the RSMA common stream as the sensing waveform for a co-located radar receiver, with no dedicated probing
signal. Private streams use zero-forcing beams built from imperfect CSIT (Gauss–Markov model with error parameter τ),
and successive interference cancellation (SIC) is imperfect (residual δ). With the common beam fixed over a sensing
dwell, the common-power fraction α is the only design variable. It must support common decoding, private decoding and
target detection at once.

The paper shows that the communication-feasible set of α is either a closed interval or empty, localises its endpoints
with analytical outage expressions, derives a high-SNR floor–ceiling law (common decoding sets a floor, imperfect CSIT
and SIC set a ceiling), and shows that a geometry-matched common beam reduces the joint communication–sensing check to
the upper endpoint of the communication-feasible interval. Monte Carlo simulations validate the analysis.

## Scripts
| script | paper figure | what it computes |
|---|---|---|
| `fig2_outage_vs_alpha.m` | Fig. 2 | per-user common- and private-stream outage versus α: analytical expressions against Monte Carlo, with the feasible interval |
| `fig3_tau_Rp_frontier.m` | Fig. 3 | communication-feasibility frontier in the (τ, R_p) plane, with Monte Carlo certification on a 7 × 9 grid |
| `fig4_high_snr_frontier.m` | Fig. 4 | (a) feasible-interval width versus SNR and its asymptotic limit; (b) high-SNR outage under ideal and imperfect SIC |
| `fig5_detection_matched_vs_fixed.m` | Fig. 5 | conditional detection probability over the feasible interval for a geometry-matched common beam versus a fixed benchmark beam |

Fig. 1 is the system-model diagram and has no script.

Each script is self-contained (local functions, fixed random seed) and writes its figure as a vector PDF to the
working directory. Baseline operating point: N_t = 8 antennas, K = 4 users, τ = 0.20, δ = 0.05, R_c = 0.20 and
R_p = 1.00 bit/s/Hz.

## Requirements
- MATLAB R2020a or newer (`exportgraphics`).
- Statistics and Machine Learning Toolbox for `fig5_detection_matched_vs_fixed.m` (`chi2inv`, `ncx2cdf`).

## Run
```matlab
fig2_outage_vs_alpha
fig3_tau_Rp_frontier
fig4_high_snr_frontier
fig5_detection_matched_vs_fixed
```

## Citation
```bibtex
@inproceedings{harsha2026common,
  author    = {Surya Harsha and Shankar Prakriya},
  title     = {Feasibility of Common-Stream-Only Sensing in {ZF-RSMA} with Imperfect {CSIT}},
  booktitle = {Proc. IEEE Global Communications Conference (GLOBECOM)},
  year      = {2026}
}
```
