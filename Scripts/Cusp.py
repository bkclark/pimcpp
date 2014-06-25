
tau = 0.125
lam = 1.0
Z1Z2 = 1.0
gam = tau*(Z1Z2*Z1Z2)/lam

P = [1.772453851,-0.074137740,0.005834805,-0.000382686,0.000008738,0.000002138,-0.000000356,0.000000021]

u00 = 0.
for j in range(1,len(P)+1):
    u00 += ((1)**j) * P[j-1] * (gam**(j/2.))
    print u00
#print u00/tau

#du00 = -Z1Z2/lam
#print du00
#
#from scipy.misc import factorial, comb
#from scipy.special import zeta
#
#def kappa(n):
#    tot = mu(n)
#    for m in range(1,n):
#        tot -= comb(n-1,m-1)*kappa(m)*mu(n-m)
#    return tot
#
#def mu(n):
#    if(n % 2 == 0):
#        return (((-1)**n)/factorial(n,exact=1)) * factorial(n/2.,exact=1) * zeta(n+2,1)
#    else:
#        
#        return (((-1)**n)/factorial(n,exact=1)) * factorial(n/2.,exact=1) * zeta(n+2,1)
#
#print kappa(1)
