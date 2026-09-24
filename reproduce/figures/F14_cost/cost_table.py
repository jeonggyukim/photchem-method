"""F14 data: per-cycle cost from the loop_time files of the runs behind F10, F12 and F13.
Each loop_time row holds the seconds rank 0 spent over the last ncycle_out cycles, per
task list. Columns, per cycle: All, Photchem, hydro (TimeIntegratorTaskList), OpSplit,
RaytPoint, and other = All minus the listed tasks, which is mostly output (the GOW17
rad_snr runs write HDF5 every 0.005 code, the classic and NCR ones none). Writes
cost_table.txt."""
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.expanduser('~/Documents/tigris-photchem-gow17-multi-ion')
RUNS = [('rad_snr', 'classic', 'M6_rad_snr/F13_series/run_classic'),
        ('rad_snr', 'NCR', 'M6_rad_snr/F13_series/run_ncr'),
        ('rad_snr', 'GOW17 core', 'M6_rad_snr/F13_series/run_core'),
        ('rad_snr', 'GOW17 + O3,S3,N3', 'M6_rad_snr/F13_series/run_ions'),
        ('hii_dtype', 'GOW17 + O3,S3,N3', 'M5_hii_dtype_ions/F12_run/run')]
TASKS = ('All', 'Photchem', 'TimeIntegratorTaskList', 'OpSplit', 'RaytPoint', 'RaytDiffuse',
         'Userwork', 'NewDt', 'Feedback', 'SelfGravity', 'AMR', 'ProcessNewParticles',
         'MassReturn')


def parse(fn):
    rows = []
    for ln in open(fn):
        kv = dict(p.split('=') for p in ln.strip().split(',') if '=' in p)
        if 'ncycle' in kv:
            rows.append(kv)
    ncycle = int(rows[-1]['ncycle'])
    tot = {k: sum(float(r.get(k, 0.0)) for r in rows)/ncycle for k in TASKS}
    tot['other'] = tot['All'] - sum(v for k, v in tot.items() if k != 'All')
    return ncycle, tot


out = ['# problem mode cycles; per cycle [s]: all photchem hydro opsplit rayt other']
for prob, mode, d in RUNS:
    ncycle, t = parse(glob.glob(os.path.join(BASE, d, '*.loop_time.txt'))[0])
    line = '%-9s %-18s %5d %.4e %.4e %.4e %.4e %.4e %.4e' % (
        prob, mode.replace(' ', '_'), ncycle, t['All'], t['Photchem'],
        t['TimeIntegratorTaskList'], t['OpSplit'], t['RaytPoint'], t['other'])
    out.append(line)
    print(line)
open(os.path.join(HERE, 'cost_table.txt'), 'w').write('\n'.join(out) + '\n')
