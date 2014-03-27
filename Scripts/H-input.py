import genDM

###### Inputs ######

## Hydrogen Atom
unit1 = 'H' # energy unit
unit2 = 'A' # length unit
tau = 0.125 # time step
pairs = [['e','e',0.5,0.5,1.0],['e','p',0.5,0.0002723089072243553,-1.0],['p','p',0.0002723089072243553,0.0002723089072243553,1.0]]
suffix = '' # suffix to input files

## Potgen
OurEwald = False # use our Ken's ewald (True) or David's (False)
nGrid = 100 # number of grid points
rMin = 0.00001 # first grid point
L = 20.0 # length of box
breakup = 1 # 2 - Short-ranged only, 1 - Optimized breakup, 0 - Classical Ewald breakup
rCut = 5. # r cutoff for ewald
nCut = 4 # n cutoff for ewald (kCut = nCut*2.0*math.pi/L)
gridType = "LOG" # LOG/LINEAR

## Squarer
nTemp = 1 # number of temperatures at which to compute the density matrix (must be <= 8)
D = 3 # dimension

## Run through pairs
for [type1,type2,lam1,lam2,Z1Z2] in pairs:
    genDM.GenDM(suffix,unit1,unit2,type1,type2,lam1,lam2,Z1Z2,L,D,tau,gridType,nGrid,rMin,rCut,nCut,nTemp,breakup,OurEwald)
