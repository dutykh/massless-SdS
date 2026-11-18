from pathlib import Path
import json,csv,math
B=Path(__file__).resolve().parent
j=json.loads((B/'jost_check_results.json').read_text())['rows']
with (B/'near_nariai_anchor_results.csv').open() as f:s=list(csv.DictReader(f))
sectors={0:'scalar',1:'em',2:'odd'}
comp=[]
for a in j:
 matches=[b for b in s if b['sector']==sectors[a['spin']] and int(b['ell'])==a['ell'] and int(b['n'])==a['n'] and abs(float(b['Lambda_tilde'])-a['lam'])<1e-13]
 assert len(matches)==1,(a,matches)
 b=matches[0];jw=complex(a['varpi_re'],a['varpi_im']);sw=complex(float(b['varpi_re']),float(b['varpi_im']))
 comp.append({'Lambda_tilde':a['lam'],'sector':b['sector'],'ell':a['ell'],'n':a['n'],'jost_varpi_re':jw.real,'jost_varpi_im':jw.imag,'spectral_varpi_re':sw.real,'spectral_varpi_im':sw.imag,'absolute_difference':abs(jw-sw),'jost_refinement_difference':a['refinement_difference'],'spectral_max_resolution_spread':float(b['max_resolution_spread']),'jost_match_residual':a['refined_normalized_residual'],'third_setting_jost_match_residual':a['independent_match_residual']})
with (B/'jost_spectral_comparison.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(comp[0]));w.writeheader();w.writerows(comp)
summary={'count':len(comp),'max_absolute_difference':max(a['absolute_difference'] for a in comp),'max_spectral_resolution_spread':max(a['spectral_max_resolution_spread'] for a in comp),'max_jost_refinement_difference':max(a['jost_refinement_difference'] for a in comp),'rows':comp}
(B/'jost_spectral_comparison.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2))
