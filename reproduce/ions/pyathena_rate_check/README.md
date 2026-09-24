# Tigris vs pyathena ion rate coefficients

Ions: He0-He1, N0-N2, O0-O2, S0-S3 (stage q = lower stage). T = 1e3, 3e3, 5e3, 8e3, 1e4, 2e4, 5e4, 1e5, 1e6 K.
Processes: rec = RR+DR of X^(q+1) -> X^q; ci = X^q -> X^(q+1); ctrec = X^(q+1) + H0; ction = X^q + H+.

## Run
Compile: `/opt/homebrew/bin/g++-16 -std=c++17 -O2 -I$HOME/Projects/tigris-gow17/src tigris_rates.cpp -o tigris_rates`
Tigris rates: `./tigris_rates > tigris_rates.csv`
Compare: `/opt/homebrew/Caskroom/miniforge/base/envs/pyathena/bin/python compare_rates.py` (writes rate_comparison.txt)

## Files
tigris_rates.cpp: calls RecombRate::Rate, CollIonRate::Rate, ChargeTransferRate::CtRec/CtIon on inputs/tables/rates/*.dat.
compare_rates.py: pyathena RecRate.get_rec_rate (default caseB=True, only affects H), CollIonRate.get_ci_rate, CT via PhotChem._ct_rate_safe.
rate_comparison.txt: summary (max |pyathena/tigris - 1| over T, rates > 1e-15 cm^3 s^-1) and full table; py_raw = ChargeTransferRate.get_ct_*_rate without the PhotChem wrapper.

## Findings
Data files: Badnell RR/DR 2023, Voronov coll_ion, Kingdon-Ferland CT ion/rec are identical in content between the two repos (comments/whitespace aside).
rec, ci: agree to < 4e-11 relative for all ions (limit set by 11-digit CSV output).
CT for N, S, O1, O2: agree to < 4e-11 relative.
O0 ctrec (O+ + H0): pyathena uses Draine (2011) eqs 14.24-14.26 (Stancil 1999), sum of 3 J channels; Tigris uses KF96 table with T clamped to [Tmin, Tmax] = fit range (flat 1.04e-9 above 1e4 K). pyathena/tigris - 1 = +0.17 at 1e3 K, +0.93 at 1e4 K, +19.8 at 1e6 K.
O0 ction (O + H+): pyathena PhotChem uses Draine k0i (J=2 only); Tigris uses KF96. pyathena/tigris - 1 = +0.35 at 1e3 K, +0.95 at 1e4 K, +17.0 at 1e6 K.
Draine fits are power laws in T4 with no upper clamp; above ~1e4 K the pyathena O0 CT rates are extrapolations.
He ctrec: PhotChem skips CT for He (0); Tigris CtRec gives He+ + H0 1.3e-15 to 2.1e-13 and He2+ + H0 1.0e-14 cm^3 s^-1.
He+ + H0 raw pyathena get_ct_rec_rate is 3.99x Tigris: pyathena overrides a = 7.47e-6 (KF96 Table 1) in place of the file's 1.87e-6.
He ction: 0 in both.
pyathena ct_rate.py special-case loop zips 10 elements with 9 a1 values, so the S I + H+ override is never applied; S0 ction therefore equals the file fit in both codes.
