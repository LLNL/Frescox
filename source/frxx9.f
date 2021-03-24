!***********************************************************************
! 
!    Copyright 2018, I.J. Thompson
!
!    This file is part of FRESCOX.
!
!    FRESCOX is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License
!    as published by the Free Software Foundation; either version 2
!    of the License, or (at your option) any later version.
!    
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!    
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
!    MA  02110-1301, USA.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************
*****FRXX9*************************************************************
      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP,
     *                  MODE1,KFN,IFAIL,M1)
C  REVISED IJT WITH L-T ALGORITHM FOR CONTINUED FRACTIONS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  A. R. BARNETT           MANCHESTER  MARCH   1981                    C
C                                                                      C
C  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           C
C                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           C
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           C
C  THIS VERSION WRITTEN UP       IN    CPC 27 (1982) 147-166           C
C                                                                      C
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    C
C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   C
C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     C
C   THE DIRAC EQUATION ,ALSO SPHERICAL & CYLINDRICAL BESSEL EQUATIONS  C
C                                                                      C
C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    C
C  STARTING ARRAY ELEMENT IS M1 = MAX0(int(XLMIN+ACCUR),0) + 1       C
C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES                 C
C                                                                      C
C  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES     C
C            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    C
C            = 3      F               CALL TO AT LEAST LENGTH (1)      C
C  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED            C
C            = 1 SPHERICAL   BESSEL      "      "     "                C
C            = 2 CYLINDRICAL BESSEL      "      "     "                C
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
C                                                                      C
C  PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    C
C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     C
C   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   C
C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   C
C   FOR MANTISSAS OF 56 & 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) C
C   REASSIGN sqrt=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	use io
	use drier
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION    FC(*),GC(*),FCP(*),GCP(*)
      LOGICAL      ETANE0,XLTURN
      COMMON       /STEE / PACCQ,NFP,NPQ,IEXP,M11
C***  COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE
C***  COULFG HAS CALLS TO: sqrt,abs,mod,int,sign,real,min
      DATA ZERO,ONE,TWO,TEN2,ABORT /0.0D0, 1.0D0, 2.0D0, 1.0D2, 1.0D7/
      DATA HALF,TM30 / 0.5D0, 1.0D-30 /
      DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 D0/
C *** THIS CONSTANT IS  sqrt(TWO/PI):  USE Q0 FOR IBM REAL*16: D0 FOR
C ***  REAL*8 & CDC DOUBLE P:  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.
C
C                        ACCUR = 1.0D-16
C ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
      ACCUR = ACC8
C ***
      MODE  = 1
      IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      ETA   = ETA1
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO
      ACC   = ACCUR * 10D0
      ACC2  = ACC*TEN2
      ACC4  = ACC*TEN2*TEN2
!     ACC4  = sqrt(ACC)
      ACCH  = sqrt(ACC)
C ***    TEST RANGE OF XX, EXIT IF.LE.sqrt(ACCUR) OR IF NEGATIVE
C
      IF(XX .LE. ACCH)                          GO TO 100
      X     = XX
      XLM   = XLMIN
      IF(KFN .EQ. 2)  XLM = XLM - HALF
      IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 105
      E2MM1 = ETA*ETA + XLM*XLM + XLM
      XLTURN= X*(X - TWO*ETA) .LT. XLM*XLM + XLM
      DELL  = XLMAX - XLMIN + ACC
C     IF(abs(mod(DELL,ONE)) .GT. ACC) WRITE(KO,2040)XLMAX,XLMIN,DELL
      LXTRA = int(DELL)
      XLL   = XLM + real(LXTRA)
C ***       LXTRA IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
C ***       XLL  IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
C ***         DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN
      M1  = MAX0(int(XLMIN + ACC),0) + 1
      L1  = M1 + LXTRA
C
C ***    EVALUATE CF1  =  F   =  FPRIME(XL,ETA,X)/F(XL,ETA,X)
C
      XI  = ONE/X
      FCL = ONE
      PK  = XLL + ONE
      PX  = PK  + ABORT
      F   =  ETA/PK + PK*XI
         IF(abs(F).LT.TM30) F = TM30
         D = ZERO
         C = F
C
C ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
C
    4 PK1   = PK + ONE
        EK  = ETA / PK
        RK2 = ONE + EK*EK
        TK  = (PK + PK1)*(XI + EK/PK1)
        D   =  TK - RK2 * D
        C   =  TK - RK2 / C
         IF(abs(C).LT.TM30) C = TM30
         IF(abs(D).LT.TM30) D = TM30
         D = ONE/D
         DF = D * C
         F  = F * DF
            IF(D .LT. ZERO) FCL = - FCL
         PK = PK1
                          IF( PK .GT. PX ) GO TO 110
      IF(abs(DF-ONE) .GE. ACC)             GO TO 4
                  NFP = PK - XLL - 1
      IF(LXTRA .EQ. 0)                          GO TO 7
C
C *** DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
C
      FCL = FCL*TM30
      FPL = FCL*F
      IF(MODE .EQ. 1) FCP(L1) = FPL
                      FC (L1) = FCL
      XL  = XLL
      RL  = ONE
      EL  = ZERO
      DO 6  LP = 1,LXTRA
         IF(ETANE0) EL = ETA/XL
         IF(ETANE0) RL = sqrt(ONE + EL*EL)
         SL    =  EL  + XL*XI
         L     =  L1  - LP
         FCL1  = (FCL *SL + FPL)/RL
         FPL   =  FCL1*SL - FCL *RL
         FCL   =  FCL1
         FC(L) =  FCL
         IF(MODE .EQ. 1) FCP(L)  = FPL
         IF(MODE .NE. 3 .AND. ETANE0) GC(L+1) = RL
    6 XL = XL - ONE
      IF(FCL .EQ. ZERO) FCL = ACC
      F  = FPL/FCL
C ***    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
C ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
C ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM
C
    7 IF( XLTURN ) CALL JWKB(X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP)
      IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE/(ACCH*TEN2))  GO TO 9
          XLTURN = .FALSE.
      TA =  TWO*ABORT
      PK =  ZERO
      WI =  ETA + ETA
      P  =  ZERO
      Q  =  ONE - ETA*XI
      AR = -E2MM1
      AI =  ETA
      BR =  TWO*(X - ETA)
      BI =  TWO
      DR =  BR/(BR*BR + BI*BI)
      DI = -BI/(BR*BR + BI*BI)
      DP = -XI*(AR*DI + AI*DR)
      DQ =  XI*(AR*DR - AI*DI)
    8 P     = P  + DP
         Q  = Q  + DQ
         PK = PK + TWO
         AR = AR + PK
         AI = AI + WI
         BI = BI + TWO
         D  = AR*DR - AI*DI + BR
         DI = AI*DR + AR*DI + BI
         C  = ONE/(D*D + DI*DI)
         DR =  C*D
         DI = -C*DI
         A  = BR*DR - BI*DI - ONE
         B  = BI*DR + BR*DI
         C  = DP*A  - DQ*B
         DQ = DP*B  + DQ*A
         DP = C
         IF(PK .GT. TA)                         GO TO 120
      IF(abs(DP)+abs(DQ).GE.(abs(P)+abs(Q))*ACC)   GO TO 8
                      NPQ   = PK/TWO
                      PACCQ = HALF*ACC/min(abs(Q),ONE)
                      IF(abs(P) .GT. abs(Q)) PACCQ = PACCQ*abs(P)
C
C *** SOLVE FOR FCM = F AT LAMBDA = XLM,THEN FIND NORM FACTOR W=W/FCM
C
      GAM = (F - P)/Q
!            IF(Q .LE. ACC2*abs(P)) GO TO 130 
       IF(Q .LE. ACC2*abs(P)) WRITE(KO,2030) P,Q,ACC2,DELL,LXTRA,M1


      W   = ONE/sqrt((F - P)*GAM + Q)
            GO TO 10
C *** ARRIVE HERE IF G(XLM) .GT. 10**6 OR IEXP .GT. 70 & XLTURN = .TRUE.
    9 W   = FJWKB
      GAM = GJWKB*W
      P   = F
      Q   = ONE
C
C *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C
   10                     ALPHA = ZERO
          IF(KFN  .EQ. 1) ALPHA = XI
          IF(KFN  .EQ. 2) ALPHA = XI*HALF
                          BETA  = ONE
          IF(KFN  .EQ. 1) BETA  = XI
          IF(KFN  .EQ. 2) BETA  = sqrt(XI)*RT2DPI
      FCM  = sign(W,FCL)*BETA
           FC(M1)  = FCM
                      IF(MODE .EQ. 3)           GO TO 11
           IF(.NOT. XLTURN)   GCL =  FCM*GAM
           IF(      XLTURN)   GCL =  GJWKB*BETA
           IF( KFN .NE. 0 )   GCL = -GCL
           GC(M1)  = GCL
           GPL =  GCL*(P - Q/GAM) - ALPHA*GCL
                      IF(MODE .EQ. 2)           GO TO 11
           GCP(M1) = GPL
           FCP(M1) = FCM*(F - ALPHA)
   11 IF(LXTRA .EQ. 0 ) RETURN
C *** UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
C *** RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
         W    = BETA*W/abs(FCL)
         MAXL = L1 - 1
      DO 13 L = M1,MAXL
                      IF(MODE .EQ. 3)           GO TO 12
                      XL = XL + ONE
         IF(ETANE0)   EL = ETA/XL
         IF(ETANE0)   RL = GC(L+1)
                      SL = EL + XL*XI
         GCL1     = ((SL - ALPHA)*GCL - GPL)/RL
         GPL      =   RL*GCL -  (SL + ALPHA)*GCL1
         GCL      = GCL1
         GC(L+1)  = GCL1
                      IF(MODE .EQ. 2)           GO TO 12
         GCP(L+1) = GPL
         FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
   12 FC(L+1)     = W* FC(L+1)
!      write(201,*) 'COULFG =',XX,ETA1,XL,FC(L+1),GC(L+1)
   13 continue
      RETURN
C
C ***    ERROR MESSAGES
C
  100 IFAIL = -1
      WRITE(KO,2000) XX,ACCH
 2000 FORMAT(' FOR XX = ',1P,D12.3,' TRY SMALL-X  SOLUTIONS',
     *' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',D12.3/)
      RETURN
  105 IFAIL = -2
      WRITE(KO,2005) XLMAX,XLMIN,XLM
 2005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES:XLMAX,XLMIN,XLM = ',
     *1P,3D15.6/)
      RETURN
  110 IFAIL =  1
      WRITE(KO,2010) ABORT,F ,DF,PK,XX,ACC
 2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
     *' F,DF,PK,XX,ACCUR =  ',1P,5D12.3//)
      RETURN
  120 IFAIL =  2
      WRITE(KO,2020) ABORT,P,Q,DP,DQ,ACC
 2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
     *' P,Q,DP,DQ,ACCUR =  ',1P,4D17.7,D12.3//)
      RETURN
  130 IFAIL =  3
      WRITE(KO,2030) P,Q,ACC,DELL,LXTRA,M1
 2030 FORMAT(' FINAL Q.LE.abs(P)*ACC2 , P,Q,ACC2 = ',1P,3D12.3,4X,
     *' DELL,LXTRA,M1 = ',D12.3,2I5 /)
      RETURN
 2040 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P,3D20.10/)
      END
C
      SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS    FOR XL.GE. 0
C *** AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
C *** CALLS DMAX1,SQRT,DLOG,EXP,ATAN2,FLOAT,INT        BARNETT FEB 1981
      DATA   ZERO,HALF,ONE,SIX,TEN/ 0.0D0, 0.5D0, 1.0D0, 6.0D0, 10.0D0 /
      DATA  DZERO, RL35, DLOGE  /0.0D0, 35.0D0, 0.43429 45 D0 /
      X     = XX
      ETA   = ETA1
      GH2   = X*(ETA + ETA - X)
      XLL1  = DMAX1(XL*XL + XL,DZERO)
      IF(GH2 + XLL1 .LE. ZERO) RETURN
       HLL  = XLL1 + SIX/RL35
       HL   = SQRT(HLL)
       SL   = ETA/HL + HL/X
       RL2  = ONE + ETA*ETA/HLL
       GH   = SQRT(GH2 + HLL)/X
       PHI  = X*GH - HALF*( HL*DLOG((GH + SL)**2/RL2) - DLOG(GH) )
          IF(ETA .NE. ZERO) PHI = PHI - ETA*ATAN2(X*GH,X - ETA)
      PHI10 = -PHI*DLOGE
      IEXP  =  INT(PHI10)
      IF(IEXP .GT. 70) GJWKB = TEN**(PHI10 - FLOAT(IEXP))
      IF(IEXP .LE. 70) GJWKB = EXP(-PHI)
      IF(IEXP .LE. 70) IEXP  = 0
      FJWKB = HALF/(GH*GJWKB)
      RETURN
      END

      SUBROUTINE WHIT(HETA,R,XK,E,LL,F,FD,IE)
C
C     CALCULATES  WHITTAKER  FUNCTION  WL(K,R)  WITH
C     ASYMPTOTIC  FORM  EXP(-(KR + ETA(LOG(2KR)))
C     E  IS  NEGATIVE
C     If IE = 0, allowed to return result e**IE larger than Whittaker,
C                for the IE value returned.
C     If IE > 0, must scale results by that amount.
C
	use drier
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(LL+1),FD(LL+1) ,T(12),S(7)
      L = LL+1
C              NOW L = NO. OF VALUES TO FIND
      EE=-1.0
      AK=XK
      ETA=HETA
      LP1=L+1
      RHO=AK*R
	S(:) = 0
      IF(L-50)1,1,2
    1 LM=60
      GO TO 3
    2 LM=L+10
    3 LMP1=LM+1
      IS=7
      PJE=30.0*RHO+1.0
      H=max(INT(PJE),4)
      H=RHO/H
      RHOA=10.0*(ETA+1.0)
      IF(RHOA-RHO)13,13,14
   13 IFEQL=1
      RHOA=RHO
      GO TO 15
   14 IFEQL=0
   15 PJE=RHOA/H+0.5
      RHOA=H*INT(PJE)
      IF(IFEQL)16,16,18
   16 IF(RHOA-RHO-1.5*H)17,18,18
   17 RHOA=RHO+2.0*H
   18 IF(EE)55,55,19
   19 STOP 'WHIT'
   27 A=2.0-10.0/12.0*H*H*EE
      B=1.0/6.0*H*ETA
      C=1.0+1.0/12.0*H*H*EE
      M1=INT(RHOA/H-0.5)
      M2=INT(RHO/H-1.5)
      T(2)=B/FLOAT(M1+1)
      T(3)=B/FLOAT(M1)
      JS=M1
      DO 29 IS=M2,M1
      DO 28 I=1,6
      S(I)=S(I+1)
   28 CONTINUE
      T(1)=T(2)
      T(2)=T(3)
      T(3)=B/FLOAT(JS-1)
      S(7)=((A+10.0*T(2))*S(6)-(C-T(1))*S(5))/(C-T(3))
      JS=JS-1
      IF(ABS(S(7)).LE.FPMAX) GO TO 29
       DO 285 I=2,7
  285   S(I) = S(I) / FPMAX
   29 CONTINUE
      T(1)=S(4)
      T(2)=(1.0/60.0*(S(1)-S(7))+0.15*(S(6)-S(2))+0.75*(S(3)-S(5)))/H
      GO TO 60
   55 C=1.0/RHOA
      A=1.0
      B=1.0-C*ETA
      F(1)=A
      FD(1)=B
      DO 56 M=1,26
      D=0.5*(ETA+FLOAT(M-1))*(ETA+FLOAT(M))*C/FLOAT(M)
      A=-A*D
      B=-B*D-A*C
      F(1)=F(1)+A
      FD(1)=FD(1)+B
   56 CONTINUE
      A=-ETA*LOG(2.0*RHOA)-RHOA
      FPMINL = -LOG(FPMAX)
      if(IE.eq.0.and.A.LT.FPMINL) IE = INT(FPMINL-A)
      A=EXP(A+IE)
      F(1)=A*F(1)
c      FD(1)=A*FD(1)
      FD(1)=A*FD(1) * (-1d0 - 2*ETA/(RHOA))
      IF(IFEQL)57,57,61
   57 S(IS)=F(1)
      IF(IS-7)27,58,27
   58 IS=6
      RHOA=RHOA+H
      GO TO 55
   60 F(1)=T(1)
      FD(1)=T(2)
   61 C=1.0/RHO
      DO 63 M=1,L-1
      A=ETA/FLOAT(M)
      B=A+C*FLOAT(M)
      F(M+1)=(B*F(M)-FD(M))/(A+1.0)
      FD(M+1)=(A-1.0)*F(M)-B*F(M+1)
   63 CONTINUE
      DO 65 M=1,L
      FD(M)=AK*FD(M)
   65 CONTINUE
      RETURN
      END
      SUBROUTINE PHASES(ETA,LLMAX,SIGMA)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 ETA,SIGMA(LLMAX)
      INTEGER MAX
      COMPLEX*16 S,CLGAMM,Z
      Z=CMPLX( 1D0,ETA)
      S=CLGAMM(Z)
      SIGMA(1)=S*(0D0,-1D0)
      IF (LLMAX.EQ.1) RETURN
      DO 1 MAX=2,LLMAX
    1 SIGMA(MAX)= ATAN2(ETA,(MAX-1D0))+SIGMA(MAX-1)
      RETURN
      END
      FUNCTION CLGAMM (Z)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 S,U,V,ZP,Z,CLGAMM
      DIMENSION A(6)
      DATA C,A,D/4.1893853321D-1,-1.3346722722D-2,8.4175084173D-4,     -
     15.9523809524D-4,7.9365079365D-4,-2.777777777D-3,8.33333333D-2,0D0/
      ZRE = Z
      I=9- int(ZRE)
      IF(I.LT.0)I=0
      Z = Z +CMPLX( real(I),D)
      ZP=Z**2
      V=(1.,0.)
      S=(0.,0.)
      DO 1 J=1,6
      S=V*A(J)+S
    1 V=V*ZP
      ZP=S*Z/V
      U=Z-(0.5,0.)
      CLGAMM = LOG(Z)*U-U+ZP+C
      IF(I.EQ.0)RETURN
      DO 2 J=1,I
      Z=Z-(1.,0.)
    2 CLGAMM=CLGAMM - LOG(Z)
      RETURN
      END
!	include 'cf2.f'
      FUNCTION CF2(X,ETA,ZL,PM,EPS,ERR,LIMIT,ACC8)
      IMPLICIT COMPLEX*16(A-H,O-Z)
      REAL*8 EPS,ERR,ACC8,TA,RK,ABSC,ZERO,HALF,ONE,TWO
      DATA ZERO,HALF,ONE,TWO / 0D+0, .5D+0, 1D+0, 2D+0 /
      ABSC(W) = ABS(DREAL(W)) + ABS(AIMAG(W))
C
C                                    (omega)        (omega)
C *** Evaluate  CF2  = p + PM.q  =  H   (ETA,X)' / H   (ETA,X)
C                                    ZL             ZL
C     where PM = omega.i
C
	ERR = 1.0
	CF2 = 0.0
      TA = TWO*LIMIT
      E2MM1 = ETA*ETA + ZL*ZL + ZL
      ETAP = ETA * PM
      XI = ONE/X
      WI = TWO*ETAP
      RK = ZERO
      PQ = (ONE - ETA*XI) * PM
      AA = -E2MM1 + ETAP
      BB = TWO*(X - ETA + PM)
         RL = XI * PM
      IF(ABSC(BB).LT.EPS) THEN
         IF(ABSC(AA + RK + WI).LT.EPS) return
         RL = RL * AA / (AA + RK + WI)
         PQ = PQ + RL * (BB + TWO*PM)
            AA = AA + TWO*(RK+ONE+WI)
            BB = BB + (TWO+TWO)*PM
            RK = RK + (TWO+TWO)
         ENDIF
      IF(ABSC(BB).LT.EPS) return
      DD = ONE/BB
      DL = AA*DD* RL
   10 PQ    = PQ + DL
         RK = RK + TWO
         AA = AA + RK + WI
         BB = BB + TWO*PM
         DD = ONE/(AA*DD + BB)
         DL = DL*(BB*DD - ONE)
            ERR = ABSC(DL)/ABSC(PQ)
         IF(ERR.GE.MAX(EPS,ACC8*RK*HALF) .AND. RK.LE.TA) GO TO 10
C
         NPQ   = RK/TWO
         PQ    = PQ + DL
           IF(NPQ.GE.LIMIT-1)
     X             WRITE(6,1000) INT(AIMAG(PM)),NPQ,ERR,ZL
 1000 FORMAT(' CF2(',I2,') NOT CONVERGED FULLY IN ',I7,
     X' ITERATIONS, So error in irregular solution =',1p,d11.2,' at zl
     X=', 0P,2F8.3)
      CF2 = PQ
      RETURN
      END
