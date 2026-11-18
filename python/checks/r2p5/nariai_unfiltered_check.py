"""Unfiltered finite Nariai pencils and independently normalized Jost tests.

Run: python nariai_unfiltered_check.py
No analytic QNM frequency list is used in assembly, selection, or classification.
The N roots grid and degree N-1 Chebyshev expansion reproduce the uploaded
Maple assembler's discretization, here with a fresh binary64 implementation.
Every finite eigenvalue is saved, with algebraic multiplicity. Only the stated
rectangle is audited physically; the remaining roots are retained unclassified.
"""
from pathlib import Path
import csv, json, platform, hashlib
import numpy as np
import scipy
from scipy.linalg import eigvals
import mpmath as mp

BASE = Path(__file__).resolve().parent
WINDOW = dict(real_min=-3., real_max=3., imag_min=-3.75, imag_max=.25)
JOST_TOL = mp.mpf('1e-8')

def pencil(A, N):
    x = np.cos((2*np.arange(1,N+1)-1)*np.pi/(2*N))
    T = np.polynomial.chebyshev.chebvander(x,N-1)
    D = np.column_stack([np.polynomial.Chebyshev.basis(k).deriv()(x) for k in range(N)])
    D2 = np.column_stack([np.polynomial.Chebyshev.basis(k).deriv(2)(x) for k in range(N)])
    M0 = (1-x*x)[:,None]*D2 - 2*x[:,None]*D - A*T
    M1 = T + 2*x[:,None]*D
    M2 = T
    Z = np.zeros((N,N)); I = np.eye(N)
    AA = np.block([[Z,I],[-M0,-1j*M1]])
    BB = np.block([[I,Z],[Z,M2]])
    return eigvals(AA,BB)

def in_window(w):
    return WINDOW['real_min'] <= w.real <= WINDOW['real_max'] and WINDOW['imag_min'] <= w.imag <= WINDOW['imag_max']

def midpoint_series(A, w, dps=80):
    """Physical Jost solutions at chi*=0, direct 2F1 series, no gamma factors.

    J_-=exp(-i*w*chi*) F(a,b;c;q), q=(1+tanh chi*)/2.
    At q=1/2, J_-=F; J_-'=-i*w*F+F_q/2. Reflection gives J_+.
    A=0 is simplified exactly before any singular gamma/2F1 parameter.
    """
    with mp.workdps(dps):
        w = mp.mpc(w.real,w.imag)
        if A == 0:
            F=mp.mpc(1); P=-1j*w
            W=2j*w
            return dict(W=W,D=W/(abs(F)**2+abs(P)**2),F=F,P=P,terms=0)
        q=mp.sqrt(mp.mpf(1)/4-A)
        a,b,c=mp.mpf('.5')+q,mp.mpf('.5')-q,1-1j*w
        # Normalized Jost solutions have parameter poles at c=0,-1,... .
        if abs(c.imag) < mp.mpf('1e-20') and c.real <= 0 and abs(c.real-mp.nint(c.real)) < mp.mpf('1e-20'):
            raise ValueError('Jost normalization pole: requires separate Laurent analysis')
        term=mp.mpc(1); F=term; Fq=mp.mpc(0)
        small=0
        for n in range(1,10001):
            term *= (a+n-1)*(b+n-1)/(c+n-1)/n/2
            F += term
            derivative_term=2*n*term
            Fq += derivative_term
            bound=mp.power(10,-dps+15)*max(1,abs(F),abs(Fq))
            small=small+1 if max(abs(term),abs(derivative_term)) < bound else 0
            if small >= 8: break
        else: raise RuntimeError('Hypergeometric series failed to converge')
        P=-1j*w*F+Fq/2
        W=-2*F*P
        den=abs(F)**2+abs(P)**2
        return dict(W=+W,D=+(W/den),F=+F,P=+P,terms=n)

def gamma_wronskian(A,w,dps=100):
    """Uncancelled physical W, simplified only with gamma recurrence.

    W=2iw Gamma(1-iw)Gamma(-iw)/(Gamma(a-iw)Gamma(b-iw))
     =-2 Gamma(1-iw)^2 rgamma(a-iw)rgamma(b-iw).
    For A=0 all frequency-dependent gamma factors cancel: W=2iw.
    The formula at Jost normalization poles is not evaluated as a zero.
    """
    with mp.workdps(dps):
        w=mp.mpc(w.real,w.imag)
        if A == 0: return 2j*w
        q=mp.sqrt(mp.mpf(1)/4-A)
        return -2*mp.gamma(1-1j*w)**2*mp.rgamma(mp.mpf('.5')-1j*w+q)*mp.rgamma(mp.mpf('.5')-1j*w-q)

def classify(A,w):
    if not in_window(w): return dict(status='outside audited window')
    try:
        low=midpoint_series(A,w,80)
        high=midpoint_series(A,w,120)
        gamma=gamma_wronskian(A,w,120)
    except ValueError as exc:
        if not str(exc).startswith('Jost normalization pole:'):
            raise
        # The normalized Jost basis is singular here. Preserve the candidate
        # as unresolved; neither clear denominators nor manufacture a zero.
        return dict(status='unresolved: Jost normalization pole',
                    normalized_Jost_residual=None,W_re=None,W_im=None,
                    Jost_precision_difference=None,gamma_identity_relative_error=None,
                    series_terms=None)
    with mp.workdps(120):
        precision_delta=abs(low['D']-high['D'])
        gamma_delta=abs(high['W']-gamma)/max(1,abs(gamma))
        residual=abs(high['D'])
        if A == 0:
            status='static threshold candidate' if abs(w)<1e-8 else 'excluded: nonzero free-wave W'
        elif precision_delta > mp.mpf('1e-40') or gamma_delta > mp.mpf('1e-60'):
            status='unresolved: evaluation disagreement'
        elif residual <= JOST_TOL:
            status='passes physical Jost check'
        else:
            status='fails physical Jost tolerance'
        return dict(status=status,normalized_Jost_residual=float(residual),W_re=float(high['W'].real),W_im=float(high['W'].imag),Jost_precision_difference=float(precision_delta),gamma_identity_relative_error=float(gamma_delta),series_terms=high['terms'])

def max_available_metric(rows,key):
    """An unresolved pole has no numerical residual, not a zero residual."""
    return max((r[key] for r in rows if r.get(key) is not None),default=None)

def run():
    rows=[]; summaries=[]
    for A in (0,2):
        for N in (16,24,32):
            roots=pencil(A,N)
            finite=sorted([complex(w) for w in roots if np.isfinite(w.real) and np.isfinite(w.imag)],key=lambda w:(-w.imag,w.real))
            case=[]
            for k,w in enumerate(finite,1):
                row=dict(A=A,N=N,raw_index=k,varpi_re=w.real,varpi_im=w.imag,in_audited_window=in_window(w),**classify(A,w))
                rows.append(row);case.append(row)
            audited=[r for r in case if r['in_audited_window']]
            counts={status:sum(r['status']==status for r in audited) for status in sorted(set(r['status'] for r in audited))}
            summaries.append(dict(A=A,N=N,total_eigenvalues=len(roots),finite_eigenvalues=len(finite),window_eigenvalues=len(audited),status_counts=counts,max_Jost_residual=max_available_metric(audited,'normalized_Jost_residual'),max_precision_difference=max_available_metric(audited,'Jost_precision_difference'),max_gamma_identity_error=max_available_metric(audited,'gamma_identity_relative_error')))
    fields=list(dict.fromkeys(k for row in rows for k in row))
    for name,data in [('nariai_raw_all.csv',rows),('nariai_raw_window.csv',[r for r in rows if r['in_audited_window']])]:
        with (BASE/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(data)
    metadata=dict(method='Fresh roots-grid Chebyshev coefficient pencil and independent physical midpoint Jost determinant; no analytic-frequency proximity criterion.',window=WINDOW,raw_multiplicities_retained=True,N_is_number_of_nodes=True,polynomial_degree='N-1',Jost_tolerance=str(JOST_TOL),Jost_evaluation_decimal_precisions=[80,120],normalization='D=W[J_-,J_+]/(||(J_-,dJ_-/dchi*)||_2 ||(J_+,dJ_+/dchi*)||_2), evaluated at chi*=0. Unit-amplitude physical horizon Jost functions. Reflection makes denominator |F|^2+|P|^2.',gamma_identity='W=-2 Gamma(1-i varpi)^2/[Gamma(1/2-i varpi+sqrt(1/4-A)) Gamma(1/2-i varpi-sqrt(1/4-A))]. A=0 is simplified exactly to W=2i varpi before evaluation. Normalization poles are not classified by a regularized determinant zero.',scope='Finite binary64 pencils in the specified rectangle only; not a rerun of historical Maple/Julia matrices, not interval certification or all-mode completeness.',runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),summaries=summaries)
    (BASE/'nariai_unfiltered_results.json').write_text(json.dumps(metadata,indent=2))
    print(json.dumps(summaries,indent=2))

if __name__=='__main__': run()
