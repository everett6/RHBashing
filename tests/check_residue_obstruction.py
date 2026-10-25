"""Check residue-class energy identity and exact modulus-six centering.

Finite numerical checks support implementation, not the asymptotic theorem.
"""
import math
import numpy as np
from check_fourier_reduction import mangoldt

def windows(a,n):
    p=np.concatenate(([0.0],np.cumsum(a)))
    x=np.arange(n,2*n)
    return p[x+1]-p[x-n//2+1]

def check(n):
    h=n//2
    idx=np.arange(2*n)
    lam=mangoldt(2*n)
    odd=idx%2==1
    codd=np.where(odd,lam-2,0.0)
    weight=np.maximum(0,np.minimum(2*n-1,idx+h-1)-np.maximum(n,idx)+1)
    diagonal=float(np.sum(codd*codd*weight))
    classes=np.asarray([windows(np.where(idx%6==a,codd,0.0),n) for a in range(6)])
    grouped=float(np.sum(classes*classes)-diagonal)
    # Class 3 has only powers of 3 as nonzero Mangoldt contributions.
    bound=2+math.log(3)
    assert np.max(np.abs(classes[3]+h/3))<=bound+1e-9
    lower=n*max(0,h/3-bound)**2-diagonal
    assert grouped>=lower-1e-7*max(1,abs(lower))
    if n<=64:
        direct=0.0
        for shift in range(6,h,6):
            j=np.arange(1,2*n-shift)
            k=np.maximum(0,np.minimum(2*n-1,j+h-1)-np.maximum(n,j+shift)+1)
            direct+=2*float(np.sum(codd[j]*codd[j+shift]*k))
        assert np.isclose(grouped,direct,rtol=1e-10,atol=1e-8)

    reduced=np.gcd(idx,6)==1
    c6=np.where(reduced,lam-3,0.0)
    delta=3*windows(reduced.astype(float),n)-h
    powers=windows(np.where(reduced,0.0,lam),n)
    s=windows(lam,n)-h
    w6=windows(c6,n)
    assert np.max(abs(delta))<=2
    assert np.max(powers)<=math.log(6)+1e-9
    assert np.allclose(s-w6,powers+delta,atol=1e-9,rtol=1e-10)
    assert np.max(abs(s-w6))<=2+math.log(6)+1e-9
    return grouped,grouped/n**3,lower,float(np.max(abs(s-w6)))

if __name__=='__main__':
    print('N,multiples6_signed_offdiag,offdiag_over_N3,proved_lower_formula,max_mod6_window_error')
    for n in (16,64,256,1024,4096):
        print(n,*check(n),sep=',')
