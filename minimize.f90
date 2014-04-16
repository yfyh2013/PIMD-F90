subroutine minimize(RR, tolg, tolx, stpmx, fret)
use system_mod
use consts
implicit none
double precision, dimension(3, Natoms) :: RR
double precision :: tolg, tolx, stpmx, fret
integer :: Neq, iter

Neq = 3*Natoms
write(*,'(80("-"),/,80("-"),/,20x,a,/,80("-"),/,80("-"))')"E N E R G Y    M I N I M I Z A T I O N"
print*,'Total number of eqs to minimize = ', Neq
call dfpmin(Neq, RR, tolg, tolx, stpmx, iter, fret)

print*,'Miinimized Energy = ', fret

end subroutine minimize
!!!!!!!!!!!!!!!!!
subroutine func(n, p, fp)
use system_mod
use consts
implicit none
double precision :: fp
integer          :: n
double precision, dimension(n) :: p
double precision, dimension(3, Natoms) :: RR, dRR
double precision, dimension(3,3) :: virt
double precision, dimension(3, Nwaters)   :: dip_momI
double precision, dimension(Natoms)   :: chg

RR(1:3, 1:Natoms) = reshape(p(1:3*Natoms), (/3, Natoms/))
call potential(RR, fp, dRR, virt, dip_momI, chg)

end subroutine func
!!!!!!!!!!!!!!!!!
subroutine dfunc(n, p, g, fp)
use system_mod
use consts
implicit none
double precision :: fp
integer          :: n
double precision, dimension(n) :: p, g
double precision, dimension(3, Natoms) :: RR, dRR
double precision, dimension(3,3) :: virt
double precision, dimension(3, Nwaters)   :: dip_momI
double precision, dimension(Natoms)   :: chg

RR(1:3, 1:Natoms) = reshape(p(1:3*Natoms), (/3, Natoms/))
call potential(RR, fp, dRR, virt, dip_momI, chg)
g(1:3*Natoms) = reshape(dRR(1:3, 1:Natoms), (/3*Natoms/))

end subroutine dfunc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE dfpmin(n,p,gtol,tolx, stpmx,iter,fret)
!.......................................................................
      implicit none
      INTEGER iter,n,nwat
      double precision :: fret,gtol,tolx,stpmx,p(n),g(n)
      integer, parameter :: ITMAX=9999
      INTEGER i,j,its,jn
      LOGICAL check
      double precision, dimension(:,:), allocatable :: hessin
      double precision :: den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp,test, &
      dg(N),tmp, hdg(N),pnew(N),xi(N)
      double precision :: testx, testg,EPS
      external func, dfunc
      integer :: iu

      allocate(hessin(N,N))
      EPS=TOLX/4.d0
      call dfunc(n, p, g, fp)
      sum=0.d0 
      do 12 i=1,n
        do 11 j=1,n
          hessin(i,j)=0.D0
11      continue
        hessin(i,i)=1.D0
        xi(i)=-g(i)
        sum=sum+p(i)**2
12    continue
      stpmax=STPMX*dmax1(dsqrt(sum),dble(n))
      do 27 its=1,ITMAX
        iter=its
        call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check)
        fp=fret
        do 13 i=1,n
          xi(i)=pnew(i)-p(i)
          p(i)=pnew(i)
13      continue
        test=0.d0
        do 14 i=1,n
          temp=dabs(xi(i))/dmax1(dabs(p(i)),1.d0)
          if(temp.gt.test)test=temp
14      continue
        if(test.lt.TOLX) goto 1111
        testx=test
        do 15 i=1,n
          dg(i)=g(i)
15      continue
        call dfunc(n, p, g, tmp)
        test=0.d0
        den=dmax1(fret,1.d0)
        do 16 i=1,n
          temp=dabs(g(i))*dmax1(dabs(p(i)),1.d0)/den
          if(temp.gt.test)test=temp
16      continue
        if(test.lt.gtol) goto 1111
        testg=test
         write(*,776)iter,fp, testx, tolx,testg, gtol
        do 17 i=1,n
          dg(i)=g(i)-dg(i)
17      continue
        do 19 i=1,n
          hdg(i)=0.d0
          do 18 j=1,n
            hdg(i)=hdg(i)+hessin(i,j)*dg(j)
18        continue
19      continue
        fac=0.d0
        fae=0.d0
        sumdg=0.d0
        sumxi=0.d0
        do 21 i=1,n
          fac=fac+dg(i)*xi(i)
          fae=fae+dg(i)*hdg(i)
          sumdg=sumdg+dg(i)**2
          sumxi=sumxi+xi(i)**2
21      continue
        if(fac.gt.dsqrt(EPS*sumdg*sumxi))then
          fac=1.d0/fac
          fad=1.d0/fae
          do 22 i=1,n
            dg(i)=fac*xi(i)-fad*hdg(i)
22        continue
          do 24 i=1,n
            do 23 j=i,n
              hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+  &
                                              fae*dg(i)*dg(j)
              hessin(j,i)=hessin(i,j)
23          continue
24        continue
        endif
        do 26 i=1,n
          xi(i)=0.d0
          do 25 j=1,n
            xi(i)=xi(i)-hessin(i,j)*g(j)
25        continue
26      continue
27    continue
      stop 'too many iterations in dfpmin'
1111  deallocate(hessin)
      return
776              format("%Miter=",i6,1x,"fp=",f16.8,1x,&
               "tolx=",d10.5,"(",1d6.1,")",1x, &
               "tolg=",d9.4,"(",1d6.1,")")
777              format("%Miter=",i6,1x,"fp=",f16.8,1x,&
               "tolx=",d10.5,"(",1d6.1,")",1x, &
               "tolg=",d9.4,"(",1d6.1,")",  &
            /,1x,"box:",3(f11.6,1x), /,1x,"ang:",3(f11.6,1x))
      END
!.......................................................................
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check)
!.......................................................................
      implicit none
      LOGICAL check
      integer :: n
      double precision ::  f,fold,stpmax,g(n),p(n), x(n),xold(n),POT
!      double precision, parameter :: ALF=1.d-4,TOLX=1.d-7
      double precision, parameter :: ALF=1.d-4,TOLX=1.d-7
      EXTERNAL func,pot
!U    USES func,pot
      INTEGER i, nwat
      double precision :: a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum, &
                          temp, test,tmplam

      check=.false.
      sum=0.d0
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=dsqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.d0
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      if(slope.ge.0.d0) stop 'roundoff problem in lnsrch'
      test=0.d0
      do 14 i=1,n
        temp=dabs(p(i))/dmax1(dabs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.d0
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        call func(n, x, f)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          goto 1111
        else if(f.le.fold+ALF*alam*slope)then
          goto 1111
        else
          if(alam.eq.1.d0)then
            tmplam=-slope/(2.d0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.d0)then
              tmplam=-slope/(2.d0*b)
            else
              disc=b*b-3.d0*a*slope
              if(disc.lt.0.d0)then
                tmplam=.5d0*alam
              else if(b.le.0.d0)then
                tmplam=(-b+dsqrt(disc))/(3.d0*a)
              else
                tmplam=-slope/(b+dsqrt(disc))
              endif
            endif
            if(tmplam.gt.0.5d0*alam)tmplam=0.5d0*alam
          endif
        endif
        alam2=alam
        f2=f
        alam=dmax1(tmplam,0.1d0*alam)
      goto 1
1111  return
      END
