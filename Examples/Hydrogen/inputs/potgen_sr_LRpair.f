      program potgen_sr
!version of potgen for short ranged potenials only. No k-space part 5.24.07
      implicit none
      integer mn,mx
      parameter (mn=20,mx=3000)
      real*8 pott(mx),rv(mx),diff(mx),cutr,pot0,a1,a2,a3,a4,x,dp
     &,rlread,potx,expon(20),cexpon(20)
      character fname*14,p(mn)*28,eunit*3,lunit*3
     +,pot(mn)*28,grid(mn)*28
       integer ntypes,nx,ifpair,j,ipickoff,n,ng,ix,i,k,npts
     &,intread,ln,n1,ntail

      write (*,*)' input prefix of files '
      read (*,99) fname
99    format(a14)
      ln=index(fname,' ')-1
      open(1,file=fname(1:ln)//'.in')
c output file
      open (3,file=fname(1:ln)//'.dm',status='unknown')

      ntypes=0
      nx=0
      ifpair=0
c loop thru reading the input
1     continue
      j=ipickoff(1,p,n,mn)
      if(j.eq.1) go to 2
c repeat input on output
      call echo(3,p,n)

c input energy units and length units
      if(p(1).eq.'UNITS') then
      eunit=p(2)
      lunit=p(3)
 
      elseif(p(1).eq.'GRID')then
c input grid parameters
        do i=1,n
        grid(i)=p(i)
        enddo
        ng=n
        nx=intread(p(2))
        call setgrid(mx,nx,rv,grid(3))

      elseif(p(1).eq.'POT') then
c  input type of potential
        ifpair=ifpair+1
c store away the potential specification
        do i=2,n
         pot(i-1)=p(i)
        enddo
        n1=n-1

      endif
      go to 1

2     continue
c check if potential and grid have been specified.
      if(ifpair.ne.1)stop
      if(nx.le.0)stop
c setup potential
      do i=1,nx
      call ipot(pot,n1,pott(i),rv(i),eunit,lunit,ntail,expon,cexpon)
      enddo

c now check grid
      npts=10
      do i=1,nx-1
       diff(i)=0.
       do k=1,npts
         x=rv(i)+(rv(i+1)-rv(i))*k/(npts+1.d0)
         call ipot(pot,n1,potx,x,eunit,lunit,ntail,expon,cexpon)   
         call interp(x,ix,a1,a2,a3,a4)
         dp=a1*pott(ix-1)+a2*pott(ix)+a3*pott(ix+1)+a4*pott(ix+2)-potx
         diff(i)=max(diff(i),abs(dp))
       enddo
         diff(nx)=diff(nx-1)
      enddo

c zero potential at cutr
!! 7/27/12 WARNING: cutoff turned off 
      cutr=rlread(pot(2))
!       call ipot(pot,n1,pot0,cutr,eunit,lunit,ntail,expon,cexpon)
      pot0=0.d0  ! just for testing wall
      write (3,33)pot0,ntail,(expon(i),cexpon(i),i=1,ntail)
33    format(' POTTAIL ',e14.5,i3,5(f7.2,e13.5))
       do i=1,nx
!       if(rv(i).le.cutr) then
         pott(i)=pott(i)-pot0
!       else
!         pott(i)=0
!       endif
      enddo
 
      grid(2)='1 '
      write (3,*)'RANK 2 ',nx,1
      call echo(3,grid,ng)
      write (3,*)'LABEL 1 r'
      write (3,*)'BEGIN potential 0'
      write (3,969) (pott(i),i=1,nx)
969   format(5e18.10)
      close(3)

      open (9,file=fname(1:ln)//'.dg',status='unknown')
      do i=1,nx
      write (9,969) rv(i),pott(i),diff(i)
      enddo
      close(9)
      end

      subroutine ipot(p,n,pott,rv,eunit,lunit,ntail,expon,cexpon)
c interprets potential parameters and sets up table.
      implicit real*8 (a-h,o-z)
      real*8 bb(6), cent,g
      character p(n)*(*),eunit*(*),lunit*(*),cite*80
      dimension expon(10),cexpon(10)
      save icall,nderv,cutr,b,bb,eps,sigma,rc,potmin
     & ,c6,c8,c9,c10,alpha,beta,gamma,abohr,hartok,rm,d
     &,beta1,g,evtok
      data icall/0/
      data cite/' '/

      icall=icall+1
      if(icall.eq.1)ntail=0

      if(p(1).eq.'HEDF2') then
c helium-helium (aziz-3) potential.
c       data eps,rm/10.948,2.963/
c       data a,alpha,beta,d,c6,c8,c10/184431.01,10.43329537,-2.27965105
c    +  ,1.4826,1.36745214,.42123807,.17473318/
c helium-helium (aziz-2) potential
!     data d,eps,a,alpha,c6,c8,c10/1.241314d0,10.8d0,544850.4d0,13.353384d0
!    +,1.3732412d0,.4253785d0,.178100d0/
!     data rm,beta/2.9673d0,0.d0/
!     DATA beta1,bb/3.32316d0,1.263030d6,1.399649d6,-8.389601d5
!    +,7.426020d6,-2.006420d6,8.570426d5/
c helium-helium (aziz-1992) potential.
c      data eps,rm/10.94,2.970/
c      data a,alpha,beta,d,c6,c8,c10/192215.29,10.73520708,-1.89296514
c    +  ,1.4135,1.34920045,.41365922,.17078164/
c He-He potential Aziz, PRL 74, 1586 (1995) ab initio HFD-B3-FCI1
       data eps,rm/10.956d0,2.9683d0/
       data a,alpha,beta,d,c6,c8,c10/1.86924404d5,10.5717543d0
     +,-2.07758779d0
     +,1.438d0,1.35186623d0,.41495143d0,.17151143d0/
       cite=
     +'R. A. Aziz, A. R. Janzen, M. R. Moldover, Phys. Rev. Letts. 74,
     +1586 (1995)'
       if(icall.eq.1) then
c      xe=convert('energy','meV',eunit)
c      xl=convert('length','au',lunit)
        ntail=3
        expon(1)=-6.
        cexpon(1)=-c6*eps*rm**6
        expon(2)=-8.
        cexpon(2)=-c8*eps*rm**8
        expon(3)=-10.
        cexpon(3)=-c10*eps*rm**10
        write (3,33)cite
33      format('CITATION:',a80)
        
       endif
 
!      if(rv.gt.1.828d0) then
        x=rv/rm
        a2i=1.d0/(x*x)
        a6i=a2i*a2i*a2i
        if(x.lt.d) then
          f=exp(-(d/x-1.d0)**2)
        else
          f=1.d0
        endif
         pott=eps*(a*exp(-alpha*x+beta*x*x)-f*a6i*
     +  (c6+a2i*(c8+a2i*c10)))
!       else
C THIS IS THE CEPERLEY-PARTRIDGE FORM FOR SMALL R
!       x=rv*1.8897266d0
!       sum=0.d0
!       DO  i=1,6
!         sum=bb(7-i)+x*sum
!       enddo
!       pott=sum*exp(-beta1*x)/x
!     endif


      elseif(p(1).eq.'HeNe') then
      echg=1.6021917e-19
      dkbolt=1.380622e-23
      evtokb=echg/dkbolt
      rm=3.029
      sig=2.699
      eps=0.001827*evtokb
      alfa=4.031
      beta=0.0987
      a=797.*evtokb
      d=1.28
      c6=1.810*evtokb
      c8=5.503*evtokb
      c10=20.5*evtokb
      c12=96.*evtokb
      c14=554.*evtokb
      drm=d*rm
       if(icall.eq.1) then
        ntail=5
        expon(1)=-6.
        cexpon(1)=-c6
        expon(2)=-8.
        cexpon(2)=-c8
        expon(3)=-10.
        cexpon(3)=-c10
        expon(4)=-12.
        cexpon(4)=-c12
        expon(5)=-14.
        cexpon(5)=-c14
       endif
 
      r=rv
      vhf=a*exp(-r*(alfa+r*beta))
      ri2=1./(r*r)
      ri6=ri2*ri2*ri2
      vd=-ri6*(c6+ri2*(c8+ri2*(c10+ri2*(c12+ri2*c14))))
      if(r.lt.drm)then
        f=drm/r-1.
        vd=vd*exp(-f*f)
      endif
      pott=vhf+vd

      elseif(p(1).eq.'WALL') then
        rc=rlread(p(2))
        if(n.ge.3) then
          cent=rlread(p(3))
        else
          cent=0.
        endif
        potmin=1000.
        ntail=0
        x=(rc-rv)/rm
         if(x.le.0)then
          pott=potmin
         else
        a2i=1./(x*x)
        a6i=a2i*a2i*a2i
        if(x.lt.d) then
          f=exp(-(d/x-1.)**2)
        else
          f=1
        endif
         pott=eps*(a*exp(-alpha*x+beta*x*x)-f*a6i*
     +  (c6+a2i*(c8+a2i*c10)))
         if(pott.lt.0)pott=0.
         if(pott.gt.potmin)pott=potmin
         endif
         pott=pott+cent/rv**2
         
 
      elseif(p(1).eq.'COUL') then
      nexp=1
      cutr=rlread(p(2))
      accuracy=rlread(p(4))
      eps=rlread(p(3))


      if(accuracy.gt.0) then
        alpha=-.5*log(cutr*accuracy/abs(eps))/cutr**2
        xpo=.5d0*nexp
        call gammi(gr,xpo,alpha*rv**2,gr0)
        pott=eps*gr/(gr0*rv**nexp)
       else
        pott=eps/rv**nexp
       endif

      elseif(p(1).eq.'COULLR') then
      nexp=1
      cutr=rlread(p(2))
!      accuracy=rlread(p(4))
      eps=rlread(p(3))
      pott=eps/rv**nexp


      elseif(p(1).eq.'SCCOUL') then
      if(icall.eq.1)then
       eps=rlread(p(3))
       rs=rlread(p(4))
      endif
       call epp(eps,rs,rv,pott)

      elseif(p(1).eq.'LJ') then
      if(icall.eq.1)then
       eps=rlread(p(3))
       sigma=rlread(p(4))
       ntail=2
       expon(1)=-6.
       cexpon(1)=-4*eps*sigma**6
       expon(2)=-12.
       cexpon(2)=4*eps*sigma**12
      endif
        x=rv/sigma
        a2i=1./(x*x)
        a6i=a2i*a2i*a2i
        pott=4*eps*a6i*(a6i-1.)

      elseif(p(1).eq.'SECH2') then
      if(icall.eq.1)then
       eps=rlread(p(3))
       sigma=rlread(p(4))
       ntail=0
      endif
        pott=eps/cosh(rv/sigma)**2
      elseif(p(1).eq.'GAUSS') then
      if(icall.eq.1)then
       eps=rlread(p(3))
       sigma=rlread(p(4))
       ntail=0
      endif
        pott=eps*exp(-(rv/sigma)**2)
      elseif(p(1).eq.'IWALL') then
      if(icall.eq.1)then
       eps=rlread(p(3))
       sigma=4.d0*rlread(p(4))
       npower=intread(p(5))
       ntail=0
       potmin=1200.d0
      endif
        x=rv/sigma
        pott=eps*x**npower
         pott=min(pott,potmin)

      elseif(p(1).eq.'LJ39WALL') then
      if(icall.eq.1)then
       eps=rlread(p(3))
       sigma=rlread(p(4))
       rc=rlread(p(5))
       ntail=0
       potmin=200.d0
      endif
        x=(rc-rv)/sigma
        if(x.gt.0.d0) then
        ai=1.d0/x
        a3i=ai*ai*ai
        pott=eps*a3i*(a3i*a3i-1.d0)
        pott=min(pott,potmin)
        else
          pott=potmin
        endif

      elseif(p(1).eq.'POWER') then
      if(icall.eq.1)then
       eps=rlread(p(3))
       sigma=rlread(p(4))
       expon(1)=-rlread(p(5))
       ntail=1
       cexpon(1)=4*eps*sigma**(-expon(1))
      endif
        x=rv/sigma
        pott=4*eps*x**expon(1)

      elseif(p(1).eq.'H2H2') then
      if(icall.eq.1) then
c potential parameters in atomic units
      abohr=0.529177249d00
c conversion factor for energy from hartree to K
      hartok=3.15777321d05
      alpha= 1.713d00
      beta= 1.5671d00
      gamma=0.00993d00
      c6=12.14d00
      c8=215.2d00
      c10=4813.9d00
      c9=143.1d00
      cite='I. F. Silvera and V. V. Goldman, J. Chem. Phys. 69, 42
     +09(1978)'
!     evtok=1.d0/8.617342d-5
!     abohr=1.d0
!     alpha=log(101.4)
!     beta=2.779
!     gamma=0.08
!     g=5.102
!     c6=7.254
!     c8=36.008
!     c9=0.d0
!     c10=225.56
!     cite='Buck et al, J. Chem. Phys 78,4439 (1983)'
        ntail=4
        expon(1)=-6.
        cexpon(1)=-c6*abohr**6
        expon(2)=-8.
        cexpon(2)=-c8*abohr**8
        expon(3)=-9.
        cexpon(3)=+c9*abohr**9
        expon(4)=-10.
        cexpon(4)=-c10*abohr**10
        write (3,33)cite
      endif

c distance in abohr
      rmm=3.41d0/abohr
      x=rv / abohr
      g=1.28d0*rmm

      ri=g / x
      if(x .lt. g) then
        f=exp(-( ri - 1.d00)**2) 
      else
        f=1.d00
      endif
      ri=1.d0 / x
      r2i=ri * ri
      r6i=r2i * r2i * r2i

      pott=hartok*( exp(alpha - x * (beta + gamma * x))
     &           - f*r6i*(c6 + r2i*(c8 + r2i*c10 - c9*ri)) )

      else
         write (*,*)' this potential not defined ',p(1),n
         stop

      endif

      end

