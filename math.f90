!*****************************************************************************
!*** The following subroutines used for the calculation of the Gamma function
!*** have been taken from the "Numerical Recipes"   
!*****************************************************************************
!--------------------------------------------------
function gammln(XX)

!***   taken from the "Numerical Recipes"

!--------------------------------------------------
implicit none
double precision :: XX
double precision :: GAMMLN
integer :: J
double precision :: SER,STP,TMP,X,Y
double precision, DIMENSION(6) :: COF
SAVE cof,stp
data COF,STP/76.18009172947146D0,-86.50532032941677d0,24.01409824083091D0,  &
            -1.231739572450155D0,.1208650973866179D-2,-.5395239384953D-5, &
            2.5066282746310005D0/
!
!-------- Executable code
!
X=XX
Y=X
TMP=X+5.5D0
TMP=(X+0.5D0)*DLOG(TMP)-TMP
SER=1.000000000190015D0
do j=1,6
   Y=Y+1.D0
   SER=SER+COF(J)/Y
enddo
GAMMLN=TMP+DLOG(STP*SER/X)
RETURN
END FUNCTION gammln

!--------------------------------------------------
function gammq2_3(x)

!***   taken from the "Numerical Recipes"

!--------------------------------------------------
implicit none
double precision :: gammq2_3
double precision :: a, x
double precision :: gamser, gammcf, gln

a=2.d0/3.d0
if (x<0.d0 .or. a <= 0.d0) stop 'bad arguments in gammq'
if (x<a+1.d0) then
   call gser2_3(gamser, x)
   gammq2_3 = 1.d0 - gamser
else
   call gcf2_3(gammcf, x)
   gammq2_3 = gammcf
endif
return
end function gammq2_3

!--------------------------------------------------
subroutine gser2_3(gamser, x)

!***   taken from the "Numerical Recipes"

!--------------------------------------------------
implicit none
double precision :: gamser
double precision :: a, x, gln
integer          :: n
integer, parameter :: ITMAX = 800
double precision, parameter :: EPS = 3.d-7
double precision :: ap, ssum, del, gammln

a=2.d0/3.d0
!gln = gammln(a)
gln=0.303150275147547d0
!print*,'gammln=',gammln(a)
if (x <= 0.d0) then
   if (x < 0.d0) stop 'x<0 in gser'
   gamser = 0.d0
   return
endif

ap = a
ssum = 1.d0 / a
del =ssum

do n=1, ITMAX
   ap = ap+1.d0
   del = del*x/ap
   ssum = ssum+del
   if (dabs(del) .lt. dabs(ssum)*EPS) goto 1
enddo
stop 'a too large, ITMAX too small in gser'
1 gamser = ssum*dexp(-x+a*dlog(x)-gln)
return
end subroutine gser2_3

!--------------------------------------------------
subroutine gcf2_3(gammcf, x)

!***   taken from the "Numerical Recipes"

!--------------------------------------------------
implicit none
double precision :: gammcf
double precision :: a, x, gln
integer          :: i
double precision :: an, b, c, d, del, h, gammln
integer,parameter :: ITMAX = 100
double precision, parameter :: EPS = 3.d-7, FPMIN=1.d-30

a=2.d0/3.d0
!gln = gammln(a)
gln=0.303150275147547d0
b = x+1.d0-a
c = 1.d0/FPMIN
d = 1.d0/b
h=d
do i=1, ITMAX
   an = -i*(i-a)
   b = b+2.d0
   d = an*d+b
   if (dabs(d) < FPMIN) d=FPMIN
   c = b+an/c
   if (dabs(c) < FPMIN) c=FPMIN
   d = 1.d0/d
   del = d*c
   h = h*del
   if (dabs(del-1.d0) < EPS) goto 1
enddo
stop ' a too large ITMAX too small in gcf'
1 gammcf=dexp(-x+a*dlog(x)-gln)*h
return
end subroutine gcf2_3
!............................................................
FUNCTION ERFC ( X )
implicit none
!    *******************************************************************
!    ** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
!    **                                                               **
!    ** REFERENCE:                                                    **
!    **                                                               **
!    ** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
!    **    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
!    *******************************************************************

double precision :: erfc
double precision, parameter ::  A1 = 0.254829592d0, A2 = -0.284496736d0, &
                                A3 = 1.421413741d0, A4 = -1.453152027d0, &
                                A5 = 1.061405429d0, P  =  0.3275911d0

double precision :: T, X, XSQ, TP


T  = 1.d0 / ( 1.d0 + P * X )
XSQ = X * X

TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) )

ERFC = TP * dEXP ( -XSQ )

RETURN
END function erfc
!.........................................................
DOUBLE PRECISION FUNCTION ran3(idum)
!.........................................................
INTEGER                     :: idum
INTEGER,          PARAMETER :: MBIG  = 1000000000
INTEGER,          PARAMETER :: DBIG  = 1000000000.d0
INTEGER,          PARAMETER :: MSEED = 161803398
INTEGER,          PARAMETER :: MZ    = 0
DOUBLE PRECISION, PARAMETER :: FAC   = 1.d0/DBIG
!DOUBLE PRECISION            :: ran3
INTEGER                     :: i,iff,ii,inext,inextp,k
INTEGER                     :: mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
if (idum.lt.0.or.iff.eq.0)then
   iff=1
   mj=abs(MSEED-abs(idum))
   mj=mod(mj,MBIG)
   ma(55)=mj
   mk=1
   do i=1,54
      ii=mod(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      if(mk.lt.MZ)mk=mk+MBIG
      mj=ma(ii)
   enddo
   do k=1,4
      do i=1,55
         ma(i)=ma(i)-ma(1+mod(i+30,55))
         if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
      enddo
   enddo
   inext=0
   inextp=31
   idum=1
endif
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=dble(mj)*FAC
END FUNCTION ran3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE CALERF(ARG,RESULT,JINT)
DOUBLE PRECISION FUNCTION ERFC2(ARG)
    DOUBLE PRECISION                            :: ARG, RESULT
    INTEGER                                  :: JINT

    INTEGER                                  :: I
    DOUBLE PRECISION :: A, B, C, D, DEL, FOUR, HALF, ONE, P, Q, SIXTEN, SQRPI, &
      THRESH, TWO, X, XBIG, XDEN, XHUGE, XINF, XMAX, XNEG, XNUM, XSMALL, Y, &
      YSQ, ZERO

  DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
!------------------------------------------------------------------
!  Mathematical constants
!------------------------------------------------------------------
  DATA FOUR,ONE,HALF,TWO,ZERO/4.d0,1.d0,0.5d0,2.d0,0.d0/, &
       SQRPI/5.6418958354775628695D-1/,THRESH/0.46875d0/, &
       SIXTEN/16.d0/
!------------------------------------------------------------------
!  Machine-dependent constants
!------------------------------------------------------------------
  DATA XINF,XNEG,XSMALL/1.79D+308,-26.628,1.11D-16/, &
       XBIG,XHUGD,XMAX/26.543,6.71D+7,2.53D+307/
!------------------------------------------------------------------
!  Coefficients for approximation to  erf  in first interval
!------------------------------------------------------------------
  DATA A/3.16112374387056560D+00,1.13864154151050156D+02, &
       3.77485237685302021D+02,3.20937758913846947D+03, &
       1.85777706184603153D-1/
  DATA B/2.36012909523441209D+01,2.44024637934444173D+02, &
       1.28261652607737228D+03,2.84423683343917062D+03/
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in second interval
!------------------------------------------------------------------
  DATA C/5.64188496988670089D-1,8.88314979438837594D+0, &
       6.61191906371416295D+01,2.98635138197400131D+02, &
       8.81952221241769090D+02,1.71204761263407058D+03, &
       2.05107837782607147D+03,1.23033935479799725D+03, &
       2.15311535474403846D-8/
  DATA D/1.57449261107098347D+01,1.17693950891312499D+02, &
       5.37181101862009858D+02,1.62138957456669019D+03, &
       3.29079923573345963D+03,4.36261909014324716D+03, &
       3.43936767414372164D+03,1.23033935480374942D+03/
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in third interval
!------------------------------------------------------------------
  DATA P /3.05326634961232344D-1,3.60344899949804439D-1, &
       1.25781726111229246D-1,1.60837851487422766D-2, &
       6.58749161529837803D-4,1.63153871373020978D-2/
  DATA Q /2.56852019228982242D+00,1.87295284992346047D+00, &
       5.27905102951428412D-1,6.05183413124413191D-2, &
       2.33520497626869185D-3/

!------------------------------------------------------------------
  JINT= 1

  X = ARG
  Y = ABS(X)
  IF (Y <= THRESH) THEN
!------------------------------------------------------------------
!  Evaluate  erf  for  |X| <= 0.46875
!------------------------------------------------------------------
     YSQ = ZERO
     IF (Y > XSMALL) YSQ = Y * Y
     XNUM = A(5)*YSQ
     XDEN = YSQ
     DO I = 1, 3
        XNUM = (XNUM + A(I)) * YSQ
        XDEN = (XDEN + B(I)) * YSQ
     END DO
     RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
     IF (JINT /= 0) RESULT = ONE - RESULT
     IF (JINT == 2) RESULT = EXP(YSQ) * RESULT
     GO TO 800
!------------------------------------------------------------------
!  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
!------------------------------------------------------------------
  ELSE IF (Y <= FOUR) THEN
     XNUM = C(9)*Y
     XDEN = Y
     DO I = 1, 7
        XNUM = (XNUM + C(I)) * Y
        XDEN = (XDEN + D(I)) * Y
     END DO
     RESULT = (XNUM + C(8)) / (XDEN + D(8))
     IF (JINT /= 2) THEN
        YSQ = AINT(Y*SIXTEN)/SIXTEN
        DEL = (Y-YSQ)*(Y+YSQ)
        RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
     END IF
!------------------------------------------------------------------
!  Evaluate  erfc  for |X| > 4.0
!------------------------------------------------------------------
  ELSE
     RESULT = ZERO
     IF (Y >= XBIG) THEN
        IF ((JINT /= 2) .OR. (Y >= XMAX)) GO TO 300
        IF (Y >= XHUGE) THEN
           RESULT = SQRPI / Y
           GO TO 300
        END IF
     END IF
     YSQ = ONE / (Y * Y)
     XNUM = P(6)*YSQ
     XDEN = YSQ
     DO I = 1, 4
        XNUM = (XNUM + P(I)) * YSQ
        XDEN = (XDEN + Q(I)) * YSQ
     END DO
     RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
     RESULT = (SQRPI -  RESULT) / Y
     IF (JINT /= 2) THEN
        YSQ = AINT(Y*SIXTEN)/SIXTEN
        DEL = (Y-YSQ)*(Y+YSQ)
        RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
     END IF
  END IF
!------------------------------------------------------------------
!  Fix up for negative argument, erf, etc.
!------------------------------------------------------------------
300 IF (JINT == 0) THEN
     RESULT = (HALF - RESULT) + HALF
     IF (X < ZERO) RESULT = -RESULT
  ELSE IF (JINT == 1) THEN
     IF (X < ZERO) RESULT = TWO - RESULT
  ELSE
     IF (X < ZERO) THEN
        IF (X < XNEG) THEN
           RESULT = XINF
        ELSE
           YSQ = AINT(X*SIXTEN)/SIXTEN
           DEL = (X-YSQ)*(X+YSQ)
           Y = EXP(YSQ*YSQ) * EXP(DEL)
           RESULT = (Y+Y) - RESULT
        END IF
     END IF
  END IF
800 ERFC2 = RESULT
RETURN
!---------- Last card of CALERF ----------
END FUNCTION ERFC2

