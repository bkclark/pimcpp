      function p1 (t)
      implicit none
      real*8 p1,t
      p1 = ((t-1.0)*(t-1.0)*(1.0+2.0*t))
      end function p1

      function p2 (t)
      implicit none
      real*8 p2,t
      p2 = (t*t*(3.0-2.0*t))
      end function p2

      function q1 (t)
      implicit none
      real*8 q1,t
      q1 = (t*(t-1.0)*(t-1.0))
      end function q1

      function q2 (t)
      implicit none
      real*8 q2,t
      q2 = (t*t*(t-1.0))
      end function q2

      function dp1 (t)
      implicit none
      real *8 dp1, t
      dp1 = 6.0*t*(t-1.0)
      end function dp1

      function dp2 (t)
      implicit none
      real *8 dp1,dp2, t
      dp2 = -dp1(t)
      end function dp2

      function dq1 (t)
      implicit none
      real *8 dq1, t
      dq1 = (t-1.0)*(3.0*t-1.0)
      end function dq1

      function dq2 (t)
      implicit none
      real *8 dq2, t
      dq2 = (3.0*t*t - 2.0*t)
      end function dq2



      subroutine r3spline (x,y,z,x0,dx,nx,y0,dy,ny,z0,dz,nz,F,num,vals)
      implicit none
      real*8 x,x0,dx,y,y0,dy,z,z0,dz,xlo,ylo,zlo,u,v,w,
     +       a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3
      real*8 p1,p2,q1,q2
      integer nx,ny,nz, ixl,ixh,iyl,iyh,izl,izh,num,i
      real*8 F(8,num,nz,ny,nx)
      real*8 vals(num)
      
      ixl = int((x-x0)/dx)+1
      ixh = ixl+1
      xlo = dx*(ixl-1)
      iyl = int((y-y0)/dy)+1
      iyh = iyl+1
      ylo = dy*(iyl-1)
      izl = int((z-z0)/dx)+1
      izh = izl + 1
      zlo = dz*(izl-1)

      u = (x-xlo)/dx
      v = (y-ylo)/dy
      w = (z-zlo)/dz

      a0 = p1(u)
      a1 = p2(u)
      a2 = dx*q1(u)
      a3 = dx*q2(u)

      b0 = p1(v)
      b1 = p2(v)
      b2 = dy*q1(v)
      b3 = dy*q2(v)

      c0 = p1(w)
      c1 = p2(w)
      c2 = dz*q1(w)
      c3 = dz*q2(w)

      do i = 1,num
         vals(i) =a0*(b0*(c0*F(1,i,izl,iyl,ixl)+c1*F(1,i,izh,iyl,ixl)+
     +                    c2*F(4,i,izl,iyl,ixl)+c3*F(4,i,izh,iyl,ixl))+
     +                b1*(c0*F(1,i,izl,iyh,ixl)+c1*F(1,i,izh,iyl,ixl)+
     +                    c1*F(4,i,izl,iyh,ixl)+c3*F(4,i,izh,iyh,ixl))+
     +                b2*(c0*F(3,i,izl,iyl,ixl)+c1*F(3,i,izh,iyl,ixl)+
     +                    c2*F(7,i,izl,iyl,ixl)+c3*F(7,i,izh,iyl,ixl))+
     +                b3*(c0*F(3,i,izl,iyh,ixl)+c1*F(3,i,izh,iyl,ixl)+
     +                    c1*F(7,i,izl,iyh,ixl)+c3*F(7,i,izh,iyh,ixl)))+
     +            a1*(b0*(c0*F(1,i,izl,iyl,ixh)+c1*F(1,i,izh,iyl,ixh)+
     +                    c2*F(4,i,izl,iyl,ixh)+c3*F(4,i,izh,iyl,ixh))+
     +                b1*(c0*F(1,i,izl,iyh,ixh)+c1*F(1,i,izh,iyl,ixh)+
     +                    c1*F(4,i,izl,iyh,ixh)+c3*F(4,i,izh,iyh,ixh))+
     +                b2*(c0*F(3,i,izl,iyl,ixh)+c1*F(3,i,izh,iyl,ixh)+
     +                    c2*F(7,i,izl,iyl,ixh)+c3*F(7,i,izh,iyl,ixh))+
     +                b3*(c0*F(3,i,izl,iyh,ixh)+c1*F(3,i,izh,iyl,ixh)+
     +                    c1*F(7,i,izl,iyh,ixh)+c3*F(7,i,izh,iyh,ixh)))+
     +            a2*(b0*(c0*F(2,i,izl,iyl,ixl)+c1*F(2,i,izh,iyl,ixl)+
     +                    c2*F(6,i,izl,iyl,ixl)+c3*F(6,i,izh,iyl,ixl))+
     +                b1*(c0*F(2,i,izl,iyh,ixl)+c1*F(2,i,izh,iyl,ixl)+
     +                    c1*F(6,i,izl,iyh,ixl)+c3*F(6,i,izh,iyh,ixl))+
     +                b2*(c0*F(5,i,izl,iyl,ixl)+c1*F(5,i,izh,iyl,ixl)+
     +                    c2*F(8,i,izl,iyl,ixl)+c3*F(8,i,izh,iyl,ixl))+
     +                b3*(c0*F(5,i,izl,iyh,ixl)+c1*F(5,i,izh,iyl,ixl)+
     +                    c1*F(8,i,izl,iyh,ixl)+c3*F(8,i,izh,iyh,ixl)))+
     +            a3*(b0*(c0*F(2,i,izl,iyl,ixh)+c1*F(2,i,izh,iyl,ixh)+
     +                    c2*F(6,i,izl,iyl,ixh)+c3*F(6,i,izh,iyl,ixh))+
     +                b1*(c0*F(2,i,izl,iyh,ixh)+c1*F(2,i,izh,iyl,ixh)+
     +                    c1*F(6,i,izl,iyh,ixh)+c3*F(6,i,izh,iyh,ixh))+
     +                b2*(c0*F(5,i,izl,iyl,ixh)+c1*F(5,i,izh,iyl,ixh)+
     +                    c2*F(8,i,izl,iyl,ixh)+c3*F(8,i,izh,iyl,ixh))+
     +                b3*(c0*F(5,i,izl,iyh,ixh)+c1*F(5,i,izh,iyl,ixh)+
     +                    c1*F(8,i,izl,iyh,ixh)+c3*F(8,i,izh,iyh,ixh)))

      end do
     

      end subroutine r3spline



      subroutine r3valgrad (x,y,z,x0,dx,nx,y0,dy,ny,z0,dz,nz,F,num,vals,
     +                      grads)
      implicit none
      real*8 x,x0,dx,y,y0,dy,z,z0,dz,xlo,ylo,zlo,u,v,w,
     +       a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,
     +       da0,da1,da2,da3,db0,db1,db2,db3,dc0,dc1,dc2,dc3,
     +       dxInv,dyInv,dzInv
      real*8 p1,p2,q1,q2,dp1,dp2,dq1,dq2
      integer nx,ny,nz, ixl,ixh,iyl,iyh,izl,izh,num,i
      real*8 F(8,num,nz,ny,nx), vals(num), grads(3,num)
      
      dxInv = 1.0/dx
      dyInv = 1.0/dy
      dzInv = 1.0/dz

      ixl = int((x-x0)*dxInv)+1
      ixh = ixl+1
      xlo = dx*(ixl-1)
      iyl = int((y-y0)*dyInv)+1
      iyh = iyl+1
      ylo = dy*(iyl-1)
      izl = int((z-z0)*dzInv)+1
      izh = izl + 1
      zlo = dz*(izl-1)

      u = (x-xlo)/dx
      v = (y-ylo)/dy
      w = (z-zlo)/dz

c      print *,"u = ",u," v = ",v, " w = ", w
c      print *,"F77: ixl=",ixl, " iyl=",iyl, " izl=",izl
      
      a0 = p1(u)
      a1 = p2(u)
      a2 = dx*q1(u)
      a3 = dx*q2(u)

      b0 = p1(v)
      b1 = p2(v)
      b2 = dy*q1(v)
      b3 = dy*q2(v)

      c0 = p1(w)
      c1 = p2(w)
      c2 = dz*q1(w)
      c3 = dz*q2(w)

c      print *,"a0=",a0,"a1=",a1,"a2=",a2,"a3=",a3
c      print *,"b0=",b0,"b1=",b1,"b2=",b2,"b3=",b3
c      print *,"c0=",c0,"c1=",c1,"c2=",c2,"c3=",c3
c      print *,"nx=",nx, "ny=", ny, "nz=", nz, "num=",num


      da0 = dxInv*dp1(u)
      da1 = dxInv*dp2(u)
      da2 = dq1(u)
      da3 = dq2(u)

      db0 = dyInv*dp1(v)
      db1 = dyInv*dp2(v)
      db2 = dq1(v)
      db3 = dq2(v)

      dc0 = dzInv*dp1(w)
      dc1 = dzInv*dp2(w)
      dc2 = dq1(w)
      dc3 = dq2(w)

c      print *,"F(10,10,10,10,10)=",F(8,10,10,10,10)

      do i = 1,num
         vals(i) =a0*(b0*(c0*F(1,i,izl,iyl,ixl)+c1*F(1,i,izh,iyl,ixl)  +
     +                    c2*F(4,i,izl,iyl,ixl)+c3*F(4,i,izh,iyl,ixl)) +
     +                b1*(c0*F(1,i,izl,iyh,ixl)+c1*F(1,i,izh,iyh,ixl)  +
     +                    c2*F(4,i,izl,iyh,ixl)+c3*F(4,i,izh,iyh,ixl)) +
     +                b2*(c0*F(3,i,izl,iyl,ixl)+c1*F(3,i,izh,iyl,ixl)  +
     +                    c2*F(7,i,izl,iyl,ixl)+c3*F(7,i,izh,iyl,ixl)) +
     +                b3*(c0*F(3,i,izl,iyh,ixl)+c1*F(3,i,izh,iyh,ixl)  +
     +                    c2*F(7,i,izl,iyh,ixl)+c3*F(7,i,izh,iyh,ixl)))+
     +            a1*(b0*(c0*F(1,i,izl,iyl,ixh)+c1*F(1,i,izh,iyl,ixh)  +
     +                    c2*F(4,i,izl,iyl,ixh)+c3*F(4,i,izh,iyl,ixh)) +
     +                b1*(c0*F(1,i,izl,iyh,ixh)+c1*F(1,i,izh,iyh,ixh)  +
     +                    c2*F(4,i,izl,iyh,ixh)+c3*F(4,i,izh,iyh,ixh)) +
     +                b2*(c0*F(3,i,izl,iyl,ixh)+c1*F(3,i,izh,iyl,ixh)  +
     +                    c2*F(7,i,izl,iyl,ixh)+c3*F(7,i,izh,iyl,ixh)) +
     +                b3*(c0*F(3,i,izl,iyh,ixh)+c1*F(3,i,izh,iyh,ixh)  +
     +                    c2*F(7,i,izl,iyh,ixh)+c3*F(7,i,izh,iyh,ixh)))+
     +            a2*(b0*(c0*F(2,i,izl,iyl,ixl)+c1*F(2,i,izh,iyl,ixl)  +
     +                    c2*F(6,i,izl,iyl,ixl)+c3*F(6,i,izh,iyl,ixl)) +
     +                b1*(c0*F(2,i,izl,iyh,ixl)+c1*F(2,i,izh,iyh,ixl)  +
     +                    c2*F(6,i,izl,iyh,ixl)+c3*F(6,i,izh,iyh,ixl)) +
     +                b2*(c0*F(5,i,izl,iyl,ixl)+c1*F(5,i,izh,iyl,ixl)  +
     +                    c2*F(8,i,izl,iyl,ixl)+c3*F(8,i,izh,iyl,ixl)) +
     +                b3*(c0*F(5,i,izl,iyh,ixl)+c1*F(5,i,izh,iyh,ixl)  +
     +                    c2*F(8,i,izl,iyh,ixl)+c3*F(8,i,izh,iyh,ixl)))+
     +            a3*(b0*(c0*F(2,i,izl,iyl,ixh)+c1*F(2,i,izh,iyl,ixh)  +
     +                    c2*F(6,i,izl,iyl,ixh)+c3*F(6,i,izh,iyl,ixh)) +
     +                b1*(c0*F(2,i,izl,iyh,ixh)+c1*F(2,i,izh,iyh,ixh)  +
     +                    c2*F(6,i,izl,iyh,ixh)+c3*F(6,i,izh,iyh,ixh)) +
     +                b2*(c0*F(5,i,izl,iyl,ixh)+c1*F(5,i,izh,iyl,ixh)  +
     +                    c2*F(8,i,izl,iyl,ixh)+c3*F(8,i,izh,iyl,ixh)) +
     +                b3*(c0*F(5,i,izl,iyh,ixh)+c1*F(5,i,izh,iyh,ixh)  +
     +                    c2*F(8,i,izl,iyh,ixh)+c3*F(8,i,izh,iyh,ixh)))
         grads(1,i) =
     +           da0*(b0*(c0*F(1,i,izl,iyl,ixl)+c1*F(1,i,izh,iyl,ixl)  +
     +                    c2*F(4,i,izl,iyl,ixl)+c3*F(4,i,izh,iyl,ixl)) +
     +                b1*(c0*F(1,i,izl,iyh,ixl)+c1*F(1,i,izh,iyh,ixl)  +
     +                    c2*F(4,i,izl,iyh,ixl)+c3*F(4,i,izh,iyh,ixl)) +
     +                b2*(c0*F(3,i,izl,iyl,ixl)+c1*F(3,i,izh,iyl,ixl)  +
     +                    c2*F(7,i,izl,iyl,ixl)+c3*F(7,i,izh,iyl,ixl)) +
     +                b3*(c0*F(3,i,izl,iyh,ixl)+c1*F(3,i,izh,iyh,ixl)  +
     +                    c2*F(7,i,izl,iyh,ixl)+c3*F(7,i,izh,iyh,ixl)))+
     +           da1*(b0*(c0*F(1,i,izl,iyl,ixh)+c1*F(1,i,izh,iyl,ixh)  +
     +                    c2*F(4,i,izl,iyl,ixh)+c3*F(4,i,izh,iyl,ixh)) +
     +                b1*(c0*F(1,i,izl,iyh,ixh)+c1*F(1,i,izh,iyh,ixh)  +
     +                    c2*F(4,i,izl,iyh,ixh)+c3*F(4,i,izh,iyh,ixh)) +
     +                b2*(c0*F(3,i,izl,iyl,ixh)+c1*F(3,i,izh,iyl,ixh)  +
     +                    c2*F(7,i,izl,iyl,ixh)+c3*F(7,i,izh,iyl,ixh)) +
     +                b3*(c0*F(3,i,izl,iyh,ixh)+c1*F(3,i,izh,iyh,ixh)  +
     +                    c2*F(7,i,izl,iyh,ixh)+c3*F(7,i,izh,iyh,ixh)))+
     +           da2*(b0*(c0*F(2,i,izl,iyl,ixl)+c1*F(2,i,izh,iyl,ixl)  +
     +                    c2*F(6,i,izl,iyl,ixl)+c3*F(6,i,izh,iyl,ixl)) +
     +                b1*(c0*F(2,i,izl,iyh,ixl)+c1*F(2,i,izh,iyh,ixl)  +
     +                    c2*F(6,i,izl,iyh,ixl)+c3*F(6,i,izh,iyh,ixl)) +
     +                b2*(c0*F(5,i,izl,iyl,ixl)+c1*F(5,i,izh,iyl,ixl)  +
     +                    c2*F(8,i,izl,iyl,ixl)+c3*F(8,i,izh,iyl,ixl)) +
     +                b3*(c0*F(5,i,izl,iyh,ixl)+c1*F(5,i,izh,iyh,ixl)  +
     +                    c2*F(8,i,izl,iyh,ixl)+c3*F(8,i,izh,iyh,ixl)))+
     +           da3*(b0*(c0*F(2,i,izl,iyl,ixh)+c1*F(2,i,izh,iyl,ixh)  +
     +                    c2*F(6,i,izl,iyl,ixh)+c3*F(6,i,izh,iyl,ixh)) +
     +                b1*(c0*F(2,i,izl,iyh,ixh)+c1*F(2,i,izh,iyh,ixh)  +
     +                    c2*F(6,i,izl,iyh,ixh)+c3*F(6,i,izh,iyh,ixh)) +
     +                b2*(c0*F(5,i,izl,iyl,ixh)+c1*F(5,i,izh,iyl,ixh)  +
     +                    c2*F(8,i,izl,iyl,ixh)+c3*F(8,i,izh,iyl,ixh)) +
     +                b3*(c0*F(5,i,izl,iyh,ixh)+c1*F(5,i,izh,iyh,ixh)  +
     +                    c2*F(8,i,izl,iyh,ixh)+c3*F(8,i,izh,iyh,ixh)))

         grads(2,i) = 
     +           a0*(db0*(c0*F(1,i,izl,iyl,ixl)+c1*F(1,i,izh,iyl,ixl)  +
     +                    c2*F(4,i,izl,iyl,ixl)+c3*F(4,i,izh,iyl,ixl)) +
     +               db1*(c0*F(1,i,izl,iyh,ixl)+c1*F(1,i,izh,iyh,ixl) +
     +                    c2*F(4,i,izl,iyh,ixl)+c3*F(4,i,izh,iyh,ixl)) +
     +               db2*(c0*F(3,i,izl,iyl,ixl)+c1*F(3,i,izh,iyl,ixl) +
     +                    c2*F(7,i,izl,iyl,ixl)+c3*F(7,i,izh,iyl,ixl)) +
     +               db3*(c0*F(3,i,izl,iyh,ixl)+c1*F(3,i,izh,iyh,ixl) +
     +                    c2*F(7,i,izl,iyh,ixl)+c3*F(7,i,izh,iyh,ixl)))+
     +           a1*(db0*(c0*F(1,i,izl,iyl,ixh)+c1*F(1,i,izh,iyl,ixh)  +
     +                    c2*F(4,i,izl,iyl,ixh)+c3*F(4,i,izh,iyl,ixh)) +
     +               db1*(c0*F(1,i,izl,iyh,ixh)+c1*F(1,i,izh,iyh,ixh) +
     +                    c2*F(4,i,izl,iyh,ixh)+c3*F(4,i,izh,iyh,ixh)) +
     +               db2*(c0*F(3,i,izl,iyl,ixh)+c1*F(3,i,izh,iyl,ixh) +
     +                    c2*F(7,i,izl,iyl,ixh)+c3*F(7,i,izh,iyl,ixh)) +
     +               db3*(c0*F(3,i,izl,iyh,ixh)+c1*F(3,i,izh,iyh,ixh) +
     +                    c2*F(7,i,izl,iyh,ixh)+c3*F(7,i,izh,iyh,ixh)))+
     +           a2*(db0*(c0*F(2,i,izl,iyl,ixl)+c1*F(2,i,izh,iyl,ixl)  +
     +                    c2*F(6,i,izl,iyl,ixl)+c3*F(6,i,izh,iyl,ixl)) +
     +               db1*(c0*F(2,i,izl,iyh,ixl)+c1*F(2,i,izh,iyh,ixl) +
     +                    c2*F(6,i,izl,iyh,ixl)+c3*F(6,i,izh,iyh,ixl)) +
     +               db2*(c0*F(5,i,izl,iyl,ixl)+c1*F(5,i,izh,iyl,ixl) +
     +                    c2*F(8,i,izl,iyl,ixl)+c3*F(8,i,izh,iyl,ixl)) +
     +               db3*(c0*F(5,i,izl,iyh,ixl)+c1*F(5,i,izh,iyh,ixl) +
     +                    c2*F(8,i,izl,iyh,ixl)+c3*F(8,i,izh,iyh,ixl)))+
     +           a3*(db0*(c0*F(2,i,izl,iyl,ixh)+c1*F(2,i,izh,iyl,ixh)  +
     +                    c2*F(6,i,izl,iyl,ixh)+c3*F(6,i,izh,iyl,ixh)) +
     +               db1*(c0*F(2,i,izl,iyh,ixh)+c1*F(2,i,izh,iyh,ixh) +
     +                    c2*F(6,i,izl,iyh,ixh)+c3*F(6,i,izh,iyh,ixh)) +
     +               db2*(c0*F(5,i,izl,iyl,ixh)+c1*F(5,i,izh,iyl,ixh) +
     +                    c2*F(8,i,izl,iyl,ixh)+c3*F(8,i,izh,iyl,ixh)) +
     +               db3*(c0*F(5,i,izl,iyh,ixh)+c1*F(5,i,izh,iyh,ixh) +
     +                    c2*F(8,i,izl,iyh,ixh)+c3*F(8,i,izh,iyh,ixh)))

         grads(3,i) = 
     +          a0*(b0*(dc0*F(1,i,izl,iyl,ixl)+dc1*F(1,i,izh,iyl,ixl)  +
     +                  dc2*F(4,i,izl,iyl,ixl)+dc3*F(4,i,izh,iyl,ixl)) +
     +              b1*(dc0*F(1,i,izl,iyh,ixl)+dc1*F(1,i,izh,iyh,ixl)  +
     +                  dc2*F(4,i,izl,iyh,ixl)+dc3*F(4,i,izh,iyh,ixl)) +
     +              b2*(dc0*F(3,i,izl,iyl,ixl)+dc1*F(3,i,izh,iyl,ixl)  +
     +                  dc2*F(7,i,izl,iyl,ixl)+dc3*F(7,i,izh,iyl,ixl)) +
     +              b3*(dc0*F(3,i,izl,iyh,ixl)+dc1*F(3,i,izh,iyh,ixl)  +
     +                  dc2*F(7,i,izl,iyh,ixl)+dc3*F(7,i,izh,iyh,ixl)))+
     +          a1*(b0*(dc0*F(1,i,izl,iyl,ixh)+dc1*F(1,i,izh,iyl,ixh)  +
     +                  dc2*F(4,i,izl,iyl,ixh)+dc3*F(4,i,izh,iyl,ixh)) +
     +              b1*(dc0*F(1,i,izl,iyh,ixh)+dc1*F(1,i,izh,iyh,ixh)  +
     +                  dc2*F(4,i,izl,iyh,ixh)+dc3*F(4,i,izh,iyh,ixh)) +
     +              b2*(dc0*F(3,i,izl,iyl,ixh)+dc1*F(3,i,izh,iyl,ixh)  +
     +                  dc2*F(7,i,izl,iyl,ixh)+dc3*F(7,i,izh,iyl,ixh)) +
     +              b3*(dc0*F(3,i,izl,iyh,ixh)+dc1*F(3,i,izh,iyh,ixh)  +
     +                  dc2*F(7,i,izl,iyh,ixh)+dc3*F(7,i,izh,iyh,ixh)))+
     +          a2*(b0*(dc0*F(2,i,izl,iyl,ixl)+dc1*F(2,i,izh,iyl,ixl)  +
     +                  dc2*F(6,i,izl,iyl,ixl)+dc3*F(6,i,izh,iyl,ixl)) +
     +              b1*(dc0*F(2,i,izl,iyh,ixl)+dc1*F(2,i,izh,iyh,ixl)  +
     +                  dc2*F(6,i,izl,iyh,ixl)+dc3*F(6,i,izh,iyh,ixl)) +
     +              b2*(dc0*F(5,i,izl,iyl,ixl)+dc1*F(5,i,izh,iyl,ixl)  +
     +                  dc2*F(8,i,izl,iyl,ixl)+dc3*F(8,i,izh,iyl,ixl)) +
     +              b3*(dc0*F(5,i,izl,iyh,ixl)+dc1*F(5,i,izh,iyh,ixl)  +
     +                  dc2*F(8,i,izl,iyh,ixl)+dc3*F(8,i,izh,iyh,ixl)))+
     +          a3*(b0*(dc0*F(2,i,izl,iyl,ixh)+dc1*F(2,i,izh,iyl,ixh)  +
     +                  dc2*F(6,i,izl,iyl,ixh)+dc3*F(6,i,izh,iyl,ixh)) +
     +              b1*(dc0*F(2,i,izl,iyh,ixh)+dc1*F(2,i,izh,iyh,ixh)  +
     +                  dc2*F(6,i,izl,iyh,ixh)+dc3*F(6,i,izh,iyh,ixh)) +
     +              b2*(dc0*F(5,i,izl,iyl,ixh)+dc1*F(5,i,izh,iyl,ixh)  +
     +                  dc2*F(8,i,izl,iyl,ixh)+dc3*F(8,i,izh,iyl,ixh)) +
     +              b3*(dc0*F(5,i,izl,iyh,ixh)+dc1*F(5,i,izh,iyh,ixh)  +
     +                  dc2*F(8,i,izl,iyh,ixh)+dc3*F(8,i,izh,iyh,ixh)))

      end do
     

      end subroutine r3valgrad



