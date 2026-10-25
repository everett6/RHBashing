"""Actual-prime obstruction to separate, cutoff-uniform Type II bounds."""
import numpy as np
from check_fourier_reduction import mangoldt
from check_global_vaughan import convolution, mobius

def windows(a,n):
    c=np.cumsum(a)
    x=np.arange(n,2*n)
    return c[x]-c[x-n//2]

def check(n):
    lam=mangoldt(2*n)
    logs=np.zeros(2*n)
    logs[1:]=np.log(np.arange(1,2*n))
    mu=mobius(2*n); mu[1]=0
    # U=V=1: b_U=log and B=mu_{>1}*log.
    b=convolution(mu,logs)
    assert np.allclose(b,lam-logs,atol=1e-10,rtol=1e-10)
    t=windows(logs,n)
    q=windows(b,n)
    s=windows(lam,n)-n//2
    p=t-n//2
    assert np.allclose(p+q,s,atol=1e-8,rtol=1e-9)
    lhs=float(s@s)
    rhs=float(p@p+q@q+2*(p@q))
    assert np.isclose(lhs,rhs,atol=1e-5,rtol=1e-7)
    return float(q@q)/(n**3*np.log(n)**2),lhs,rhs

if __name__=='__main__':
    print('N,TypeII_energy_over_N3log2N,full_energy,covariance_reconstruction')
    for n in (64,256,1024,4096):
        print(n,*check(n),sep=',')
