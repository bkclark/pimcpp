      subroutine evtricub(xget,yget,zget,x,nx,y,ny,z,nz,
     >                   ilinx,iliny,ilinz,
     >                   f,inf2,inf3,ict,fval,ier)
      implicit real*8 (A-H,O-Z)
C
C  use mktricub to set up spline coefficients...
C
C  evaluate a 3d cubic Spline interpolant on a rectilinear
C  grid -- this is C2 in all directions.
C
C  this subroutine calls two subroutines:
C     herm3xyz  -- find cell containing (xget,yget,zget)
C     fvtricub  -- evaluate the spline function (w/derivatives if req.)
C
C  input arguments:
C  ================
C
      real*8 xget,yget,zget               ! target of this interpolation
      real*8 x(nx)                        ! ordered x grid
      real*8 y(ny)                        ! ordered y grid
      real*8 z(nz)                        ! ordered z grid
      integer ilinx                     ! ilinx=1 => assume x evenly spaced
      integer iliny                     ! iliny=1 => assume y evenly spaced
      integer ilinz                     ! ilinz=1 => assume z evenly spaced
C
      real*8 f(0:7,inf2,inf3,nz)          ! function data
C
C       f 2nd dimension inf2 must be .ge. nx; 3rd dim inf3 .ge. ny
C       contents of f:
C
C  f(0,i,j,k) = f @ x(i),y(j),z(k)
C  f(1,i,j,k) = d2f/dx2 @ x(i),y(j),z(k)
C  f(2,i,j,k) = d2f/dy2 @ x(i),y(j),z(k)
C  f(3,i,j,k) = d2f/dz2 @ x(i),y(j),z(k)
C  f(4,i,j,k) = d4f/dx2dy2 @ x(i),y(j),z(k)
C  f(5,i,j,k) = d4f/dx2dz2 @ x(i),y(j),z(k)
C  f(6,i,j,k) = d4f/dy2dz2 @ x(i),y(j),z(k)
C  f(7,i,j,k) = d6f/dx2dy2dz2 @ x(i),y(j),z(k)
C
      integer ict(10)                   ! code specifying output desired
C
C  ict(1)=1 -- return f  (0, don't)
C  ict(2)=1 -- return df/dx  (0, don't)
C  ict(3)=1 -- return df/dy  (0, don't)
C  ict(4)=1 -- return df/dz  (0, don't)
C  ict(5)=1 -- return d2f/dx2  (0, don't)
C  ict(6)=1 -- return d2f/dy2  (0, don't)
C  ict(7)=1 -- return d2f/dz2  (0, don't)
C  ict(8)=1 -- return d2f/dxdy (0, don't)
C  ict(9)=1 -- return d2f/dxdz (0, don't)
C  ict(10)=1-- return d2f/dydz (0, don't)
C
C output arguments:
C =================
C
      real*8 fval(10)                     ! output data
      integer ier                       ! error code =0 ==> no error
C
C  fval(1) receives the first output (depends on ict(...) spec)
C  fval(2) receives the second output (depends on ict(...) spec)
C  fval(3) receives the third output (depends on ict(...) spec)
C  fval(4) receives the 4th output (depends on ict(...) spec)
C  fval(5-10) receive 5th thru 10th outputs (if required by ict(...) spec)
C
C  examples:
C    on input ict = [1,1,1,1,0,0,0,0,0,0,0]
C   on output fval= [f,df/dx,df/dy,df/dz]
C
C    on input ict = [1,0,0,0,0,0,0,0,0,0,0]
C   on output fval= [f] ... elements 2-10 never referenced
C
C    on input ict = [0,1,1,0,0,0,0,0,0,0,0]
C   on output fval= [df/dx,df/dy] ... elements 3-10 never referenced
C
C    on input ict = [0,0,0,0,1,0,0,0,0,0,0]
C   on output fval= [d2f/dx2] ... elements 2-10 never referenced.
C
C  ier -- completion code:  0 means OK
C-------------------
C  local:
C
      integer i,j,k                     ! cell indices
C
C  normalized displacement from (x(i),y(j),z(k)) corner of cell.
C    xparam=0 @x(i)  xparam=1 @x(i+1)
C    yparam=0 @y(j)  yparam=1 @y(j+1)
C    zparam=0 @z(k)  zparam=1 @z(k+1)
C
      real*8 xparam,yparam,zparam
C
C  cell dimensions and
C  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
C
      real*8 hx,hy,hz
      real*8 hxi,hyi,hzi
C
C  0 .le. xparam .le. 1
C  0 .le. yparam .le. 1
C  0 .le. zparam .le. 1
C
C---------------------------------------------------------------------
C  use lookup routine as in Hermite interpolation
C
      call herm3xyz(xget,yget,zget,x,nx,y,ny,z,nz,ilinx,iliny,ilinz,
     >   i,j,k,xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi,ier)
      if(ier.ne.0) return
c
      call fvtricub(ict,1,1,
     >   fval,i,j,k,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   f,inf2,inf3,nz)
C
      return
      end
C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 3d fcn
C   --vectorized-- dmc 10 Feb 1999
C
      subroutine fvtricub(ict,ivec,ivecd,
     >   fval,ii,jj,kk,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   fin,inf2,inf3,nz)
      implicit real*8 (A-H,O-Z)
C
C  use mktricub to set up spline coefficients...
C
      integer ict(10)                   ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)
C
      integer ii(ivec),jj(ivec),kk(ivec) ! target cells (i,j,k)
      real*8 xparam(ivec),yparam(ivec),zparam(ivec)
                          ! normalized displacements from (i,j,k) corners
C
      real*8 hx(ivec),hy(ivec),hz(ivec)   ! grid spacing, and
      real*8 hxi(ivec),hyi(ivec),hzi(ivec) ! inverse grid spacing
           ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))
C
      real*8 fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "evtricub")
C
      real*8 fval(ivecd,10)               ! output returned
C
C  for detailed description of fin, ict and fval see subroutine evtricub
C  comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.
C
C  note that the index inputs ii,jj,kk and parameter inputs
C     xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument
C
C  to use this routine in scalar mode, pass in ivec=ivecd=1
C
C---------------
C
      integer v
C
      real*8 sum
      real*8 sixth
C
      data sixth/0.166666666666666667/
C
C---------------
C
      z36th=sixth*sixth
      z216th=sixth*sixth*sixth
C
C  prepare useful parameters...
C
      do v=1,ivec
         i=ii(v)
         j=jj(v)
         k=kk(v)
C
C   ...in x direction
C
         xp=xparam(v)
         xpi=1.0-xp
         xp2=xp*xp
         xpi2=xpi*xpi
C
         if((ict(1).eq.1).or.(ict(3).eq.1).or.(ict(4).eq.1).or.
     >      (ict(6).eq.1).or.(ict(7).eq.1).or.(ict(10).eq.1)) then
            cx=xp*(xp2-1.0)
            cxi=xpi*(xpi2-1.0)
            hx2=hx(v)*hx(v)
         endif
         if((ict(2).eq.1).or.(ict(8).eq.1).or.(ict(9).eq.1)) then
            cxd=3.0*xp2-1.0
            cxdi=-3.0*xpi2+1.0
         endif
C
C   ...and in y direction
C
         yp=yparam(v)
         ypi=1.0-yp
         yp2=yp*yp
         ypi2=ypi*ypi
C
         if((ict(1).eq.1).or.(ict(2).eq.1).or.(ict(4).eq.1).or.
     >      (ict(5).eq.1).or.(ict(7).eq.1).or.(ict(9).eq.1)) then
            cy=yp*(yp2-1.0)
            cyi=ypi*(ypi2-1.0)
            hy2=hy(v)*hy(v)
         endif
         if((ict(3).eq.1).or.(ict(8).eq.1).or.(ict(10).eq.1)) then
            cyd=3.0*yp2-1.0
            cydi=-3.0*ypi2+1.0
         endif
C
C   ...and in z direction
C
         zp=zparam(v)
         zpi=1.0-zp
         zp2=zp*zp
         zpi2=zpi*zpi
C
         if((ict(1).eq.1).or.(ict(2).eq.1).or.(ict(3).eq.1).or.
     >      (ict(5).eq.1).or.(ict(6).eq.1).or.(ict(8).eq.1)) then
            cz=zp*(zp2-1.0)
            czi=zpi*(zpi2-1.0)
            hz2=hz(v)*hz(v)
         endif
         if((ict(4).eq.1).or.(ict(9).eq.1).or.(ict(10).eq.1)) then
            czd=3.0*zp2-1.0
            czdi=-3.0*zpi2+1.0
         endif
C
         iadr=0
C
C  get desired values:
C
         if(ict(1).eq.1) then
C
C  function value:
C
            iadr=iadr+1
            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+
     >            xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+
     >            xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
            sum=sum+sixth*hx2*(
     >         zpi*(
     >           cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
            sum=sum+sixth*hy2*(
     >         zpi*(
     >           xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+
     >            xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+
     >            xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
            sum=sum+sixth*hz2*(
     >         czi*(
     >           xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >            xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >            xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx2*hy2*(
     >         zpi*(
     >           cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx2*hz2*(
     >         czi*(
     >           cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
            sum=sum+z36th*hy2*hz2*(
     >         czi*(
     >           xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >            xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >            xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
            sum=sum+z216th*hx2*hy2*hz2*(
     >         czi*(
     >           cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(2).eq.1) then
C
C  df/dx:
C
            iadr=iadr+1
C
            sum=hxi(v)*(
     >         zpi*(
     >              -(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))
     >              +(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >        +zp*(
     >              -(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))
     >              +(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
            sum=sum+sixth*hx(v)*(
     >         zpi*(
     >           cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >        +zp*(
     >           cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
            sum=sum+sixth*hxi(v)*hy2*(
     >         zpi*(
     >              -(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))
     >              +(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >        +zp*(
     >              -(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))
     >              +(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
            sum=sum+sixth*hxi(v)*hz2*(
     >         czi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +cz*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx(v)*hy2*(
     >         zpi*(
     >           cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx(v)*hz2*(
     >         czi*(
     >           cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +cz*(
     >           cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
            sum=sum+z36th*hxi(v)*hy2*hz2*(
     >         czi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +cz*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
            sum=sum+z216th*hx(v)*hy2*hz2*(
     >         czi*(
     >           cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(3).eq.1) then
C
C  df/dy:
C
            iadr=iadr+1
C
            sum=hyi(v)*(
     >         zpi*(
     >           xpi*(-fin(0,i,j,k)  +fin(0,i,j+1,k))+
     >            xp*(-fin(0,i+1,j,k)+fin(0,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(-fin(0,i,j,k+1)  +fin(0,i,j+1,k+1))+
     >            xp*(-fin(0,i+1,j,k+1)+fin(0,i+1,j+1,k+1))))
C
            sum=sum+sixth*hyi(v)*hx2*(
     >         zpi*(
     >           cxi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >            cx*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >            cx*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
            sum=sum+sixth*hy(v)*(
     >         zpi*(
     >           xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+
     >            xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+
     >            xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
            sum=sum+sixth*hyi(v)*hz2*(
     >         czi*(
     >           xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >            xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >            xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx2*hy(v)*(
     >         zpi*(
     >           cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >            cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >            cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
            sum=sum+z36th*hyi(v)*hx2*hz2*(
     >         czi*(
     >           cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >            cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >            cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
            sum=sum+z36th*hy(v)*hz2*(
     >         czi*(
     >           xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+
     >            xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+
     >            xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
            sum=sum+z216th*hx2*hy(v)*hz2*(
     >         czi*(
     >           cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >            cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >            cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(4).eq.1) then
C
C  df/dz:
C
            iadr=iadr+1
C
            sum=hzi(v)*(
     >           -(
     >           xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+
     >            xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >           +(
     >           xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+
     >            xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
            sum=sum+sixth*hx2*hzi(v)*(
     >           -(
     >           cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >           +(
     >           cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
            sum=sum+sixth*hy2*hzi(v)*(
     >           -(
     >           xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+
     >            xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >           +(
     >           xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+
     >            xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
            sum=sum+sixth*hz(v)*(
     >         czdi*(
     >           xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >            xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +czd*(
     >           xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >            xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx2*hy2*hzi(v)*(
     >           -(
     >           cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >           +(
     >           cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx2*hz(v)*(
     >         czdi*(
     >           cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +czd*(
     >           cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
            sum=sum+z36th*hy2*hz(v)*(
     >         czdi*(
     >           xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >            xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +czd*(
     >           xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >            xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
            sum=sum+z216th*hx2*hy2*hz(v)*(
     >         czdi*(
     >           cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +czd*(
     >           cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(5).eq.1) then
C
C  d2f/dx2:
C
            iadr=iadr+1
C
            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
            sum=sum+sixth*hy2*(
     >         zpi*(
     >           xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
            sum=sum+sixth*hz2*(
     >         czi*(
     >           xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
            sum=sum+z36th*hy2*hz2*(
     >         czi*(
     >           xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(6).eq.1) then
C
C  d2f/dy2:
C
            iadr=iadr+1
C
            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+
     >            xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+
     >            xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))
C
            sum=sum+sixth*hx2*(
     >         zpi*(
     >           cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >            cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >            cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))
C
            sum=sum+sixth*hz2*(
     >         czi*(
     >           xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >            xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+
     >            xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx2*hz2*(
     >         czi*(
     >           cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >            cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >            cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(7).eq.1) then
C
C  d2f/dz2:
C
            iadr=iadr+1
C
            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >            xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >            xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
            sum=sum+sixth*hx2*(
     >         zpi*(
     >           cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
            sum=sum+sixth*hy2*(
     >         zpi*(
     >           xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >            xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >            xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx2*hy2*(
     >         zpi*(
     >           cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(8).eq.1) then
C
C  d2f/dxdy:
C
            iadr=iadr+1
C
            sum=hxi(v)*hyi(v)*(
     >         zpi*(
     >               (fin(0,i,j,k)  -fin(0,i,j+1,k))-
     >               (fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >        +zp*(
     >               (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))-
     >               (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
C
            sum=sum+sixth*hyi(v)*hx(v)*(
     >         zpi*(
     >           cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >            cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >        +zp*(
     >           cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >            cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))
C
            sum=sum+sixth*hxi(v)*hy(v)*(
     >         zpi*(
     >              -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))
     >              +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >        +zp*(
     >              -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))
     >              +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
            sum=sum+sixth*hxi(v)*hyi(v)*hz2*(
     >         czi*(
     >               (fin(3,i,j,k)  -fin(3,i,j+1,k))-
     >               (fin(3,i+1,j,k)-fin(3,i+1,j+1,k)))
     >        +cz*(
     >               (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))-
     >               (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx(v)*hy(v)*(
     >         zpi*(
     >           cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >            cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >            cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
            sum=sum+z36th*hyi(v)*hx(v)*hz2*(
     >         czi*(
     >           cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >            cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >        +cz*(
     >           cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >            cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
            sum=sum+z36th*hxi(v)*hy(v)*hz2*(
     >         czi*(
     >               -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))
     >               +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >        +cz*(
     >               -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))
     >               +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
            sum=sum+z216th*hx(v)*hy(v)*hz2*(
     >         czi*(
     >           cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >            cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >            cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(9).eq.1) then
C
C  d2f/dxdz:
C
            iadr=iadr+1
C
            sum=hxi(v)*hzi(v)*(
     >            (
     >              (ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k)) -
     >              (ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >           -(
     >              (ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1)) -
     >              (ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))
C
            sum=sum+sixth*hx(v)*hzi(v)*(
     >           -(
     >           cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >           +(
     >           cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))
C
            sum=sum+sixth*hxi(v)*hy2*hzi(v)*(
     >            (
     >              (cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k)) -
     >              (cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >           -(
     >              (cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1)) -
     >              (cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))
C
            sum=sum+sixth*hxi(v)*hz(v)*(
     >         czdi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +czd*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx(v)*hy2*hzi(v)*(
     >           -(
     >           cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >           +(
     >           cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx(v)*hz(v)*(
     >         czdi*(
     >           cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +czd*(
     >           cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))
C
            sum=sum+z36th*hxi(v)*hy2*hz(v)*(
     >         czdi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +czd*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))
C
            sum=sum+z216th*hx(v)*hy2*hz(v)*(
     >         czdi*(
     >           cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +czd*(
     >           cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
         if(ict(10).eq.1) then
C
C  d2f/dydz:
C
            iadr=iadr+1
C
            sum=hyi(v)*hzi(v)*(
     >            (
     >           xpi*(fin(0,i,j,k)  -fin(0,i,j+1,k))+
     >            xp*(fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >           -(
     >           xpi*(fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))+
     >            xp*(fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))
C
            sum=sum+sixth*hyi(v)*hx2*hzi(v)*(
     >            (
     >           cxi*(fin(1,i,j,k)  -fin(1,i,j+1,k))+
     >            cx*(fin(1,i+1,j,k)-fin(1,i+1,j+1,k)))
     >           -(
     >           cxi*(fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+
     >            cx*(fin(1,i+1,j,k+1)-fin(1,i+1,j+1,k+1))))
C
            sum=sum+sixth*hy(v)*hzi(v)*(
     >           -(
     >           xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+
     >            xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >           +(
     >           xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+
     >            xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))
C
            sum=sum+sixth*hyi(v)*hz(v)*(
     >         czdi*(
     >           xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >            xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >        +czd*(
     >           xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >            xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))
C
            sum=sum+z36th*hx2*hy(v)*hzi(v)*(
     >           -(
     >           cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >            cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >           +(
     >           cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >            cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))
C
            sum=sum+z36th*hyi(v)*hx2*hz(v)*(
     >         czdi*(
     >           cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >            cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >        +czd*(
     >           cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >            cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))
C
            sum=sum+z36th*hy(v)*hz(v)*(
     >         czdi*(
     >           xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+
     >            xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >        +czd*(
     >           xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+
     >            xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))
C
            sum=sum+z216th*hx2*hy(v)*hz(v)*(
     >         czdi*(
     >           cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >            cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >        +czd*(
     >           cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >            cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))
C
            fval(v,iadr)=sum
         endif
C
      enddo                             ! vector loop
C
      return
      end
