import HEGInput as hi

prefix = 'heg'
Ns = [33]
D = 3
rss = [4.0]
ts = [1.0]
pols = [1]
nMax = 4
nTemp = 1
RefToTF = 1.0
RefMs = { 4.0:64 }
MinM = 16
ParticleType = 1 # 2 - boltzmannon, 1 - fermion, 0 - boson
TrackSign = 0
restart = 0
PPC = 1
PPNODE = 16

for N in Ns:
  for rs in rss:
    for t in ts:
      for pol in pols:
        RefM = RefMs[rs]
        hi.main(['nada',prefix,N,D,rs,t,pol,nMax,nTemp,RefToTF,RefM,MinM,ParticleType,TrackSign,restart,PPC,PPNODE])
