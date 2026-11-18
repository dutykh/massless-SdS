"""Independent precision/representation check of the *same* finite pencil.

The differential operators preserve polynomials of degree <=N-1. Writing
these polynomials in the monomial basis gives exact integer coefficient
matrices, below. Their evaluation at the N distinct Chebyshev roots is an
invertible Vandermonde left factor; changing from monomials to Chebyshev
polynomials is an invertible right factor. Thus this pencil and the fresh
roots-grid collocation pencil have the same exact finite eigenvalues.

The entire companion matrix is diagonalized with mpmath.eig at 60 decimal
digits. No analytic frequency list, factorized characteristic polynomial,
or selected diagonal-block eigenvalues are used. All 2N roots are retained.
"""
from pathlib import Path
import csv,json,hashlib,time
import mpmath as mp
from nariai_unfiltered_check import classify,in_window,WINDOW,max_available_metric

BASE=Path(__file__).resolve().parent

def exact_companion(A,N):
    # Interleaved state [b0,w*b0,b1,w*b1,...] is a permutation of the
    # usual companion state [b;w*b]. M2=I in this exact coefficient basis.
    C=mp.matrix(2*N)
    for k in range(N):
        C[2*k,2*k+1]=1
        C[2*k+1,2*k]=k*(k+1)+A
        C[2*k+1,2*k+1]=mp.mpc(0,-(1+2*k))
        if k>=2:C[2*(k-2)+1,2*k]=-k*(k-1)
    return C

def run():
    rows=[];summary=[]
    with mp.workdps(60):
        for A in (0,2):
            for N in (16,24,32):
                start=time.time()
                vals=mp.eig(exact_companion(A,N),left=False,right=False)
                vals=sorted(vals,key=lambda w:(-w.imag,w.real))
                case=[]
                for k,w in enumerate(vals,1):
                    r=dict(A=A,N=N,raw_index=k,varpi_re=mp.nstr(w.real,60),varpi_im=mp.nstr(w.imag,60),in_audited_window=in_window(w),**classify(A,w))
                    case.append(r);rows.append(r)
                audited=[r for r in case if r['in_audited_window']]
                counts={status:sum(r['status']==status for r in audited) for status in sorted(set(r['status'] for r in audited))}
                summary.append(dict(A=A,N=N,total_eigenvalues=len(vals),finite_eigenvalues=len(vals),window_eigenvalues=len(audited),status_counts=counts,max_Jost_residual=max_available_metric(audited,'normalized_Jost_residual'),max_precision_difference=max_available_metric(audited,'Jost_precision_difference'),max_gamma_identity_error=max_available_metric(audited,'gamma_identity_relative_error'),elapsed_seconds=time.time()-start))
    fields=list(dict.fromkeys(k for r in rows for k in r))
    for name,data in [('nariai_polynomial_raw_all.csv',rows),('nariai_polynomial_raw_window.csv',[r for r in rows if r['in_audited_window']])]:
        with (BASE/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(data)
    result=dict(method='Full exact-matrix monomial-coefficient companion eigensolve, mpmath.eig, 60 decimal digits; independently classified by normalized physical Jost determinant, without frequency-list selection.',equivalence='Each differential operator preserves P_(N-1). Evaluation at distinct collocation nodes left-multiplies by an invertible Vandermonde matrix, and basis changes right-multiply by an invertible matrix. This is the same exact finite pencil as roots-grid Chebyshev collocation.',coefficient_action={'L0_on_y^k':'k(k-1)y^(k-2) - [k(k+1)+A] y^k','L1_on_y^k':'(2k+1)y^k','L2_on_y^k':'y^k'},pencil='L0+i*varpi*L1+varpi^2*L2',window=WINDOW,all_raw_multiplicities_retained=True,Jost_decimal_precisions=[80,120],Jost_tolerance='1e-8',scope='Representative finite pencils only; the repeated exact free-wave polynomial roots are retained, not inferred to be scattering QNMs. All roots outside the displayed rectangle remain in the full CSV, unclassified.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),summaries=summary)
    (BASE/'nariai_polynomial_precision_results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':run()
