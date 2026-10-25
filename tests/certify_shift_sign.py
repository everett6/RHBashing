"""Exact rational enclosures for finite parity-adjusted correlations.

No floating-point arithmetic enters the sign certificates. For 0<=t<1,
log((1+t)/(1-t))=2 sum_{j>=0}t^(2j+1)/(2j+1). After m terms the
positive tail is at most 2*t^(2m+1)/((2m+1)*(1-t*t)).
This certifies finite counterexamples to termwise nonpositivity, not RH.
"""
from fractions import Fraction as F
from functools import lru_cache

def atanh_log_interval(t, m=20):
    assert F(0) <= t <= F(1, 3)
    lower = sum((2*t**(2*j+1)/F(2*j+1) for j in range(m)), F(0))
    tail = 2*t**(2*m+1)/(F(2*m+1)*(1-t*t))
    return lower, lower+tail

@lru_cache(None)
def log_interval(p):
    k = p.bit_length()-1
    ratio = F(p, 2**k)
    t = (ratio-1)/(ratio+1)
    l2,u2 = atanh_log_interval(F(1,3))
    lr,ur = atanh_log_interval(t)
    return k*l2+lr, k*u2+ur

def isprime(p):
    return p>=2 and all(p%d for d in range(2,p) if d*d<=p)

def mangoldt_interval(n):
    for p in range(2,n+1):
        if isprime(p):
            power=p
            while power<n:
                power*=p
            if power==n:
                return log_interval(p)
    return F(0),F(0)

def coefficient(n):
    if n%2==0:
        return F(0),F(0)
    lo,hi=mangoldt_interval(n)
    return lo-2,hi-2

def multiply(a,b):
    products=[x*y for x in a for y in b]
    return min(products),max(products)

def correlation(n,shift):
    h=n//2
    lower=upper=F(0)
    for j in range(1,2*n-shift):
        weight=max(0,min(2*n-1,j+h-1)-max(n,j+shift)+1)
        lo,hi=multiply(coefficient(j),coefficient(j+shift))
        lower+=2*weight*lo
        upper+=2*weight*hi
    return lower,upper

if __name__=='__main__':
    neg=correlation(16,2)
    pos=correlation(16,6)
    assert F(-77)<neg[0]<=neg[1]<F(-76)
    assert F(37)<pos[0]<=pos[1]<F(38)
    print('CERTIFIED using exact rational arithmetic:')
    print('N=16, shift=2: -77 < weighted correlation < -76')
    print('N=16, shift=6:  37 < weighted correlation <  38')
