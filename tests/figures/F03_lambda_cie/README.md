# F03 CIE cooling Lambda(T)
Harness: lambda_curve.cpp (build line in its header; built against a source snapshot with
GOW17_ENABLED 1, PHOTCHEM_IONS 1), tigris-gow17 at d340425e2 plus the uncommitted
higher-ion steps 3-4 and returned C, O, Si. Run: ./lambda_curve 1.0 > lambda_curve_n1.txt
(and 0.01 > lambda_curve_n001.txt; within 0.9% of n_H = 1). Plot: python plot.py.
