"""Regenerate rounded manuscript-cell metadata from the full CF campaign."""
from pathlib import Path
import json
import mpmath as mp
p=Path(__file__).resolve().parent
mp.mp.dps=80
data=json.loads((p/'schwarzschild_cf_results.json').read_text())
rows=data['rows'];assert len(rows)==19
num=lambda s:mp.mpf(s)
string=lambda v:mp.nstr(v,36)
report={'scope':data['scope'],'mpmath_version':data['mpmath_version'],
 'normalization':data['normalization'],'truncations':[256,512,1024],
 'decimal_precisions':[40,70],'asymptotic_tail_orders':[5,3],
 'max_depth512_to1024_change':string(max(num(r['steps'][2]['absolute_change']) for r in rows)),
 'max_40_to70_digit_change':string(max(num(r['steps'][3]['absolute_change']) for r in rows)),
 'max_tail5_to3_change_at1024':string(max(num(r['tail_order_3_change']) for r in rows)),
 'max_inversion_n_to_nplus1_change_at1024':string(max(num(r['inversion_plus_1_change']) for r in rows)),
 'max_final_fraction_residual':string(max(num(r['steps'][3]['fraction_residual']) for r in rows)),
 'important':'These numerical differences are convergence diagnostics, not rigorous error bounds.',
 'cells':[{'s':r['s'],'ell':r['ell'],'n':r['n'],'rounded_4dp':r['rounded_4dp'],
 'real':r['final_real'],'imag':r['final_imag']} for r in rows]}
(p/'schwarzschild_cf_table_cells.json').write_text(json.dumps(report,indent=2)+'\n')
print('Summarized 19 continued-fraction roots and convergence diagnostics.')
