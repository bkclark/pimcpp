---
layout: base
title: Function
---

`[python,N] def VMDOut(coordList):   numPtcls=len(coordList[0]) outFile=open("myTrajectory.xyz","w")   outFile.write(str(numPtcls)+"\n") for coord in coordList: for i   in range(0,len(coord)): outFile.write(str(i)+"   "+str(coord[i][0])+" "+str(coord[i][1])+"   "+str(coord[i][2])+"\n") outFile.close()`
