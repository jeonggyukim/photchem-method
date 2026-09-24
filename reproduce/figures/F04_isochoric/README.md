# F04 isochoric cooling (E1)
Harness: isochoric.cpp (build line in its header), tigris-gow17 at d340425e2 plus the
uncommitted higher-ion steps 3-4 and returned C, O, Si. Runs: ./isochoric 1.0 >
isochoric_n1.txt; ./isochoric 0.01 > isochoric_n001.txt (about 1 s each). Next: full-ladder
python reference along the same T(t), and plot.py.
Reference: fullladder.cpp (every stage, same rate files, along the run's T(t), n_e, x_H):
./fullladder isochoric_n1.txt 1.0 > fullladder_n1.txt (same for 0.01). Columns: t, T, then
C(0..6) Si(0..14) O(0..8) S(0..16) N(0..7).
