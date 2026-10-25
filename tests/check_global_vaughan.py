"""Check exact algebra and three independent energy formulas. Not an RH test."""
import math
import numpy as np
from check_fourier_reduction import mangoldt

def convolution(a,b):
    out = np.zeros(len(a))
    for d in range(1,len(a)):
        for m in range(1,(len(a)-1)//d+1):
            out[d*m] += a[d]*b[m]
    return out

def mobius(limit):
    # Independent divisor recurrence: sum_{d|n} mu(d)=1_{n=1}.
    mu = np.zeros(limit)
    mu[1] = 1
    for n in range(1,limit):
        for m in range(2*n,limit,n):
            mu[m] -= mu[n]
    return mu

def check_vaughan(n,u,v):
    limit=2*n
    lam=mangoldt(limit)
    mu=mobius(limit)
    one=np.ones(limit); one[0]=0
    log=np.zeros(limit); log[1:]=np.log(np.arange(1,limit))
    small_mu=mu.copy(); small_mu[v+1:]=0
    small_lam=lam.copy(); small_lam[u+1:]=0
    b=convolution(one,lam-small_lam)
    term1=convolution(small_mu,log)
    term3=convolution(convolution(small_mu,one),small_lam)
    term4=convolution(mu-small_mu,b)
    reconstructed=term1+small_lam-term3+term4
    error=float(np.max(np.abs(lam-reconstructed)))
    assert np.allclose(lam,reconstructed,atol=1e-11,rtol=1e-11)
    return error

def check_global(n):
    h=n//2
    a=mangoldt(2*n)[1:]-1
    padded=np.zeros(4*n); padded[1:2*n]=a
    w=np.zeros(4*n); w[:h]=1
    # Linear convolution has no wrap: its support fits inside the period.
    g=np.convolve(padded,np.ones(h))[:4*n]
    direct=float(g@g)
    spectral=float(np.sum(abs(np.fft.fft(padded)*np.fft.fft(w))**2)/(4*n))
    diagonal=h*float(a@a)
    correlations=diagonal+2*sum((h-d)*float(a[:-d]@a[d:]) for d in range(1,h))
    assert np.isclose(direct,spectral,atol=1e-8,rtol=1e-11)
    assert np.isclose(direct,correlations,atol=1e-8,rtol=1e-11)
    assert float(g[n:2*n]@g[n:2*n]) <= direct+1e-8
    return direct,diagonal,abs(direct-spectral),abs(direct-correlations)

if __name__=='__main__':
    errors=[check_vaughan(n,u,v) for n in (16,64,128)
            for u,v in ((1,1),(3,5),(8,4),(2*n,2*n))]
    print(f'12 Vaughan cases passed; max coefficient error={max(errors):.3g}')
    print('N,global_energy,diagonal_energy,parseval_error,correlation_error')
    for n in (16,64,128,512,2048):
        print(n,*check_global(n),sep=',')
