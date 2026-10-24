"""Finite numerical consistency checks; these do not prove RH. Requires numpy."""
import math
import numpy as np

def mangoldt(limit):
    a = np.zeros(limit)
    prime = np.ones(limit, dtype=bool)
    prime[:2] = False
    for p in range(2, limit):
        if prime[p]:
            prime[p*p::p] = False
            q = p
            while q < limit:
                a[q] = math.log(p)
                q *= p
    return a

def check(n, synthetic=False):
    length, h = 4*n, n//2
    a = np.zeros(length)
    if synthetic:
        a[1:n+1] = 1
        a[n+1:2*n] = -1
    else:
        a[1:2*n] = mangoldt(2*n)[1:] - 1
    window = np.zeros(length)
    window[:h] = 1
    spectrum = np.fft.fft(a)*np.fft.fft(window)
    k = np.arange(length)
    low = np.minimum(k, length-k) < math.ceil(math.sqrt(n*math.log(2*n)))
    high = np.fft.ifft(spectrum*(~low)).real
    even = np.fft.ifft(spectrum*low*(k%2 == 0)).real
    direct = np.array([sum(a[x-r] for r in range(h)) for x in range(n,2*n)])
    reconstructed = 2*even[n:2*n]+high[n:2*n]+high[3*n:4*n]
    error = np.max(np.abs(direct-reconstructed))
    assert np.allclose(direct,reconstructed,atol=1e-9,rtol=1e-10)
    f, q = float(direct@direct), float(even[n:2*n]@even[n:2*n])
    bound = 4*n*n / math.ceil(math.sqrt(n*math.log(2*n)))**2 * float(a@a)
    assert float(high@high) <= bound*(1+1e-10)
    if synthetic:
        assert abs(f-(n**3/6+n/3)) < 1e-7
    print(f'{n},{synthetic},{error:.3g},{f:.12g},{q:.12g},{q/n**3:.12g}')

if __name__ == '__main__':
    print('N,synthetic,max_identity_error,F,Q,Q_over_N_cubed')
    for n in (16,64,128,512,2048):
        for synthetic in (False,True):
            check(n,synthetic)
