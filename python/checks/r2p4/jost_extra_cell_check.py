"""Bounded independent check of the scalar ell2 n1 table cell at Lambda=.09.
Existing 24-anchor output files are neither read as seeds nor overwritten.
"""
from pathlib import Path
import json,math,hashlib
from jost_check import Matcher
B=Path(__file__).resolve().parent
paths=[B/'jost_check_results.json',B/'jost_check_results.csv']
before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
seed=complex(math.sqrt(5.75),-1.5)
continuation=[]
for lam in [.1111,.11,.105,.10,.095,.09]:
 m=Matcher(lam,0,2)
 w,res,ok,msg=m.solve(seed)
 continuation.append(dict(Lambda_tilde=lam,varpi_re=w.real,varpi_im=w.imag,residual=res,success=ok))
 seed=w
fine=Matcher(.09,0,2,order=32,eps=.025,mid=.57,rtol=8e-14,atol=8e-15)
w,res,ok,msg=fine.solve(seed)
third=Matcher(.09,0,2,order=36,eps=.06,mid=.37,rtol=8e-14,atol=8e-15)
thirdres=abs(third.determinant(w))
Omega=fine.G['S']*w
result={'purpose':'Independent physical-Frobenius-sector verification of scalar ell2 n1 at Lambda_tilde=.09, using Nariai-seeded continuation and original-master-equation shooting.','continuation':continuation,'target':{'Lambda_tilde':.09,'spin':0,'ell':2,'n':1,'varpi_re':w.real,'varpi_im':w.imag,'Omega_re':Omega.real,'Omega_im':Omega.imag,'coarse_refined_difference':abs(w-seed),'refined_normalized_match_residual':res,'third_configuration_match_residual':thirdres,'four_decimal_Omega':f'{Omega.real:.4f}{Omega.imag:+.4f}i'},'existing_anchor_sha256':before}
assert all(hashlib.sha256(p.read_bytes()).hexdigest()==before[p.name] for p in paths)
(B/'jost_extra_cell_check.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result['target'],indent=2))
