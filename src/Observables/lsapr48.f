      subroutine lsapr48(n,c,iperm)
c real version of linear sum assignment problem from Derigs book
c input c(mnp,n) is n by n cost matrix: on output  iperm is the optimal 
c permutation: z the cost function; ys and yt are the optimal row and column
c dual values
c minimize z=sum_i c(iperm(i),i)
      implicit none
      integer mnp,n,i,j
      parameter (mnp=48) !you must change this if it changes in calling routine
      integer j0,iu,ius,index,iw,iws,ind,ize(mnp),iperm(n),ivor(mnp)
      real*8 z,eps,sup,cc, ui,vj,d,vgl,c(mnp,n),ys(mnp),yt(mnp),dm(mnp)
     & ,dp(mnp)
      logical la(mnp)
      data eps,sup/1.d-15,1.d15/
      do  i=1,n
      ize(i)=0
      iperm(i)=0
      ivor(i)=0
      ys(i)=0.d0
      yt(i)=0.d0
      enddo

      do 2 i=1,n
      do 3 j=1,n
      cc=c(j,i)
      if(j.eq.1) goto 4
      if(cc-ui.ge.eps) goto 3
4     ui=cc
      j0=j 
3     continue
      ys(i)=ui
      if(ize(j0).ne.0) go to 2
      ize(j0)=i
      iperm(i)=j0
2     continue

      do j=1,n
       if(ize(j).eq.0) yt(j)=sup
      enddo

      do 6 i=1,n
      ui=ys(i)
      do 7 j=1,n
      vj=yt(j)
      if(vj.le.eps) go to 7
      cc=c(j,i)-ui
      if(cc+eps.ge.vj) go to 7
      yt(j)=cc
      ivor(j)=i
7     continue
6     continue

      do 8 j=1,n
      i=ivor(j)
      if(i.eq.0) go to 8
      if(iperm(i).ne.0) go to 8
      iperm(i)=j
      ize(j)=i
8     continue
    
      do 9 i=1,n
      if(iperm(i).ne.0) go to 9
      ui=ys(i)
      do 10 j=1,n
      if(ize(j).ne.0) go to 10
      cc=c(j,i)
      if(cc-ui-yt(j)+eps.gt.0.) go to 10
      iperm(i)=j
      ize(j)=i
      go to 9
10    continue
9     continue

      do 1000 iu=1,n
      if(iperm(iu).gt.0) go to 1000

      ius=(iu-1)*n
      do 100 i=1,n
      ivor(i)=iu
      la(i)=.false.
      dp(i)=sup
100   dm(i)=c(i,iu)-ys(iu)-yt(i)
      dp(iu)=0.d0
105   d=sup
      do 110 i=1,n
      if(la(i)) go to 110
      if(dm(i)+eps.ge.d) go to 110
      d=dm(i)
      index=i
110   continue
      if(ize(index).le.0) go to 400
      la(index)=.true.
      iw=ize(index)
      iws=(iw-1)*n
      dp(iw)=d
      do 130 i=1,n
      if(la(i)) go to 130
      vgl=d+c(i,iw)-ys(iw)-yt(i)
      if(dm(i).le.vgl+eps) go to 130
      dm(i)=vgl
      ivor(i)=iw
130   continue
      go to 105

400   iw=ivor(index)
      ize(index)=iw
      ind=iperm(iw)
      iperm(iw)=index
      if(iw.eq.iu) go to 500

      index=ind
      go to 400

500   do 510 i=1,n
      if(dp(i).eq.sup) go to 505
      ys(i)=ys(i)+d-dp(i)
505   if(dm(i)+eps.ge.d) go to 510
      yt(i)=yt(i)+dm(i)-d
510   continue
1000  continue
   
      z=0.d0
      do i=1,n
      z=z+c(iperm(i),i)
      enddo
      end
