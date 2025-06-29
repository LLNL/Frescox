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
*****FRXX3***************************************************************
      SUBROUTINE KERNEL(KCOEF,NLL,ALOSS,ALOSSM,RER,
     &                  C1,C2,IC1,IC2,REPEAT, LVAL,ICTO,ICFROM,REV,
     &                  PART,DNL,A,B,NK,FPT,FI,CHNO,CUTOFF,CP,KIND,NIB,
     &                  NUMLT,LTMIN,P,Q,LTRANS,QNF,IC7,HF,HT,FFR,EXCIT,
     &                  JTOTAL,PARITY,JVAL)
	use parameters
	use io
	use factorials
	use drier
	use kcom
	use kcom2
	use trace
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 KCOEF(MAXQRN,NUMLT),ALOSS(MAXQRN),DNL(NLO),JTOTAL,
     x       JVAL(MAXCH)
      INTEGER LVAL(MAXCH),D,C,PART(MAXCH),FPT(2,NK),FI,PARITY,
     &    CP,QNF(19,MSP),CUTOFF,CHNO(MFNL,6),C1,C2,EXCIT(MAXCH,3)
      LOGICAL REV,PRES,ODD,THERE,FFR,REPEAT,C1FR,LTRANS(MAXQRN)
      CHARACTER*1 PSIGN(3)
      DATA PSIGN / '-','?','+' /

      DATA ONE,TWO /1D0, 2D0/
      ODD(I) = I.NE.(I/2)*2
C
      IF(REPEAT) GO TO 7
      if(FFR) then
	if(allocated(QRLN)) deallocate(QRLN,FNL)
	allocate (QRLN(NLO,NLL,NIB),FNL(NLL,NLO))
	else
	if(allocated(QERN)) deallocate(QERN,FNC)
      	allocate (QERN(NLO,NLL,NIB),FNC(NLL,NLO))
	endif
	if(allocated(WHERE)) deallocate(WHERE,WHOL,WHOI)
	allocate (WHERE(MAXMUL,MAXQRN),WHOL(NIB),WHOI(NIB))
      Z = 0.0
      EPS = 1E-14
      THMN = 0.1 *PI/180.  ! at least 0.1 degres, and moreover non-zero!
      DO 1    J=1,NLO
1     DNL(J)  = (J - NLC - ONE) * MLT * HF
      HFHT = HF/HT
      RINTO = HT * MR
      MAXL1 = MAXL + 1
      MINL1 = MINL + 1
!	write(48,*) 'Kernel: MINL,MAXL,FFR =',MINL,MAXL,FFR
      MINLRQ = MAXL1
      MAXLRQ = -1
      DO 5 IB=1,NIB
      WHOL(IB) = 0
5     WHOI(IB) = 0
      DO 6 IN=1,NK
      DO 6 L1=MINL1,MAXL1
6     WHERE(L1,IN) = 0
      REPEAT = .TRUE.
      IB = 0
7     C1FR = ICFROM.EQ.IC1
      IF(C1FR) THEN
         C = C1
         D = C2
      ELSE
         C = C2
         D = C1
      ENDIF
C     IF(ICTO.NE.PART(D).OR.ICFROM.NE.PART(C)) STOP 7
      IF(ICTO.NE.PART(D).OR.ICFROM.NE.PART(C)) CALL ABEND(32)
C       SO TO 'D' FROM 'C' CHANNEL
C  AND WITH NAME(1,IC1) LIKE D & NAME(1,IC2) LIKE P IN (D,P) REACTION
      LD = LVAL(D)
      LC = LVAL(C)
      R7 = SQRT((TWO*LD+ONE)*(TWO*LC+ONE))
      THERE = .FALSE.
      DO 8 IN=1,NK
      DO 8 ILT=1,NUMLT
 8     IF(ABS(KCOEF(IN,ILT)).GT.EPS.AND.LTRANS(IN)) THERE = .TRUE.
      IF(.NOT.THERE) GO TO 100
C
      IF(NL.EQ.0) REWIND 12
      NL = NL + 1
      DRY = DRY .OR. NL.GT.MFNL
9     FORMAT(//' ****** NOT ENOUGH ROOM FOR',I4,' NON-LOCAL FORM',
     &' FACTORS  IN MFNL ARRAY OF',I4,' *****')
	if(FFR.and..not.allocated(FNL)) write(6,*) 'FNL fail!'
	if(.not.FFR.and..not.allocated(FNC)) write(6,*) 'FNC fail!'
      DO 1111 I=1,NLL
         IF(FFR) THEN
            DO 10 J=1,NLO
10            FNL(I,J) = 0
      ELSE
             DO 11 J=1,NLO
11                 FNC(I,J) = 0
      ENDIF
1111   CONTINUE

C CHNO(NL,i) meanings:
C
C  CHNO(NL,1) =  final channel. If negative then no reverse coupling
C  CHNO(NL,2) =  initial channel
C  CHNO(NL,3) =  0 : post
C                1 : prior
C                2 : nonorthogonality coupling (KIND=8)
C                10,11,12:  complex couplings
C	         13 : complex numerical nonlocal potential
C  CHNO(NL,4) =  NLL: the number of radial interpolation points for first radial index.
C  CHNO(NL,5) =  0 : operate on wf
C                1 : operate on derivative (wf' - (L+1)/r * wf)
C  CHNO(NL,6) =  0 : Moshinsky transformation method
C                1 : m-state summation method
C


      CHNO(NL,1) = D
      IF(.NOT.REV) CHNO(NL,1) = -D
      CHNO(NL,2) = C
      CHNO(NL,3) = 2
      IF(KIND.EQ.8) CHNO(NL,3) = IC7
      IF(.NOT.FFR) CHNO(NL,3) = CHNO(NL,3) + 10
      CHNO(NL,4) = NLL
      CHNO(NL,5) = 0
      CHNO(NL,6) = 0
      PRES = .FALSE.
C
      DO 80 IN=1,NK
       IF(.NOT.LTRANS(IN)) GO TO 80
C
C   HERE, TO START :    FPT(1  IS NO. OF B.S. IN CHANNEL 'C1' ("DEUT")
C                   AND FPT(2  IS NO. OF B.S. IN CHANNEL 'C2' ("PROT")
         IF(C1FR) THEN
           KNP = FPT(1,IN)
           KN = FPT(2,IN)
         ELSE
           KNP = FPT(2,IN)
           KN  = FPT(1,IN)
         ENDIF
C
C   NOW, MORE PROPERLY, KNP IS NO. OF B.S. IN CHANNEL 'C' (FROM)
C                   AND KN  IS NO. OF B.S. IN CHANNEL 'D' (TO)
         LN = QNF(9,KN)
         LNP= QNF(9,KNP)
         VINT = 0.0
         THERE = .FALSE.
         DO 16 ILT=1,NUMLT
 16        IF(ABS(KCOEF(IN,ILT)).GT.EPS) THERE = .TRUE.
         IF(.NOT.THERE) GO TO 80
      IF(LISTCC.GT.1) WRITE(KO,15) NL,D,C,IN,KNP,KN,NUMLT,LTMIN
15    FORMAT('0NL. interaction #',I3,' to',I3,' from',I3,' :',5I3)
      LN1 = LN + 1
      LNP1 = LNP + 1
      DO 60 NU1=1,LN1
         NU = NU1 - 1
         LNMNU = LN-NU
      DO 60 NUP1=1,LNP1
         NUP = NUP1 - 1
         LNPMNU = LNP-NUP
         LAMIN1 = ABS(LNPMNU-NU)+ 1
         LAMAX1 =      LNPMNU+NU + 1
         LAMIN2 = ABS(LNMNU-NUP)+ 1
         LAMAX2 =      LNMNU+NUP + 1
         LMAX = MIN(LC+LAMAX1-1, LD+LAMAX2-1)
         LMAXP = LMAX + 1
         LMINP= MAX(LC-(LAMAX1-1), LD-(LAMAX2-1)) + 1
      R2 =       FACT(2*LN+1 +1)-FACT(2*NU +1)-FACT(2*LNMNU+1 +1)
     &         + FACT(2*LNP+1+1)-FACT(2*NUP+1)-FACT(2*LNPMNU+1+1)
      R2 = EXP( R2 * 0.5D0 )
      R3 = SQRT(ONE*(2*LN+1)*(2*LNMNU+1)*(2*LNP+1)*(2*LNPMNU+1))
      DO 55 L1=LMINP,LMAXP
         L = L1-1
C
      RL = 0.0
      DO 25 LAM11=LAMIN1,LAMAX1
         LAM1 = LAM11-1
      DO 25 LAM21=LAMIN2,LAMAX2
         LAM2 = LAM21-1
      IF(ODD(LAM1+LC+L).OR.ODD(LAM2+LD+L).OR.ODD(LNPMNU+NU+LAM1)
     &   .OR. ODD(LNMNU+NUP+LAM2)) GO TO 25
      RL1 = (TWO*LAM1+ONE) * (TWO*LAM2+ONE)
      RL2 = WIG3J(LAM1+Z,LC+Z,L+Z,Z,Z,Z)
      RL3 = WIG3J(LAM2+Z,LD+Z,L+Z,Z,Z,Z)
      RL4 = WIG3J(LNPMNU+Z,NU+Z,LAM1+Z,Z,Z,Z)
      RL5 = WIG3J(LNMNU+Z,NUP+Z,LAM2+Z,Z,Z,Z)
      RLP =       RL1*RL2*RL3*RL4*RL5
      IF(ABS (RLP).LT.EPS) GO TO 25
      IQMIN = MAX(ABS(LN-LNP),ABS(LAM1-LAM2))  + 1
      IQMAX = MIN(    (LN+LNP),    (LAM1+LAM2))  + 1
      RQ = 0.0
      DO 22 IQ1=IQMIN,IQMAX
         IQ = IQ1 - 1
      RQ1 = (2*IQ+1)*(2*L+1d0)
      RQ2 = RACAH(LAM1+Z,LC+Z,LAM2+Z,LD+Z,L+Z,IQ+Z)
      RQ3 = WIG9J(LN+Z    ,IQ+Z,    LNP+Z,
     &            NU+Z,    LAM1+Z,  LNPMNU+Z,
     &            LNMNU+Z, LAM2+Z,  NUP+Z)
      RQP =       RQ1 * RQ2 * RQ3
      IF(ABS (RQP).LT.EPS) GO TO 22
      RT = 0.0
         DO 20 ILT=1,NUMLT
            LTOTAL =LTMIN + ILT - 1
         RT1 = KCOEF(IN,ILT)
         IF(ABS(RT1).LT.EPS) GO TO 20
         RT2 = RACAH(LNP+Z,LC+Z,LN+Z,LD+Z,LTOTAL+Z,IQ+Z)
         RT3 = (-1) ** (LTOTAL + L + LC - LD)
      IF(LISTCC.GT.2) THEN
      REST = RL1*RL4*RL5 * (TWO*IQ+ONE)*RQ3 * RT2 * SQRT(TWO*LC+ONE)
                   WRITE(KO,18) D,C,KN,NU,NUP,L,LD,LC,LN,LNP,LAM1,LAM2,
     &   IQ,LTOTAL,RT1,R2,RL3,RL2,R3         ,RQ2              ,REST
C    #   IQ,LTOTAL,RT1,R2,RL3,RL2,R3*(2*L+1.),RQ2*SQRT(2*LD+1.),REST
18    FORMAT('0INNER :',14I3,7F9.4/)
      IF(LISTCC.GT.3) WRITE(KO,*)  RL1,RL4,RL5,RQ3 , RT2
         ENDIF
         RT = RT +        RT1 * RT2 * RT3
20       CONTINUE
      RQ = RQ + RQP * RT
22    CONTINUE
      RL = RL + RLP * RQ
25    CONTINUE
      T = RL * R2 * R3 * R7
      IF(ABS(T).LT.EPS) GO TO 55
      PRES = .TRUE.
      IF(C1FR) THEN
      TC = T * A**NU * B**LNMNU * Q**NUP * P**LNPMNU
C    TO GET  (A*RF)**NU      * (B*RT)**(LN-NU)     IF C1FR
C          * (P*RF)**(LNP-NUP)*(Q*RT)**NUP
         NRF = LNPMNU + NU
         NRT = LNMNU + NUP
      ELSE
      TC = T * P**NU * Q**LNMNU * B**NUP * A**LNPMNU
C    TO GET  (B*RT)**NUP     * (A*RF)**(LNP-NUP)     IF NOT C1FR
C          * (Q*RT)**(LN-NU)  *(P*RF)**NU
         NRF = LNPMNU + NU
         NRT = LNMNU + NUP
      ENDIF
C
      IF(LISTCC.GT.1)WRITE(KO,30) D,C,KN,NU,NUP,
     & L,LD,LC,LN,LNP,          RL,R2,R3,R7,    T,TC
30    FORMAT('0' ,10I3,E13.6,3F10.5,   1P,2E20.9)
         IF(MAXLRQ.LT.L) MAXLRQ = L
         IF(MINLRQ.GT.L) MINLRQ = L
      IF(L .LE. MAXL.AND.L .GE. MINL) GO TO 40
      NL = NL - 1
      GO TO 100
35    FORMAT(' OVERLAP MULTIPOLE ORDER',I4,' IS REQUIRED, AND ',A3,'IMUM
     & PRE-CALCULATED IS',I4,' SO SOME CHANNELS ARE UNCOUPLED.')
40    IF(WHERE(L1,IN).NE.0) GO TO 47
         IB = IB + 1
         IF(IB.GT.NIB) IB = 1
         WHERE(L1,IN) = IB
         IF(WHOL(IB).NE.0) then
 		I = WHOL(IB); J=WHOI(IB)
 		WHERE(I,J) = 0
 		endif
         WHOL(IB) = L1
         WHOI(IB) = IN
      IF(LISTCC.GT.3) WRITE(KO,46) L1,IN,IB,NIB
46    FORMAT(' Read multipole #',I3,' of pair',I2,' to',I3,' in block of
     &',I4)
         IF(DRY) GO TO 55
C
      NREC =  (IN-1)*(MAXL1-MINL1+1) + L1 - MINL1    + 1    +  FI

!!!	write(6,*) 'NREC,IB,NLO,NLL = ',NREC,IB,NLO,NLL
      IF(     FFR) READ(11,REC=NREC) ((QRLN(J,I,IB),J=1,NLO),I=2,NLL)
!     IF(     FFR) CALL FIOR(11,QRLN(1,1,IB),NLO*(NLL  ),NREC,ISTAT)
      IF(.NOT.FFR) READ(9 ,REC=NREC) ((QERN(J,I,IB),J=1,NLO),I=2,NLL)
!     IF(.NOT.FFR) CALL FIORC( 9,QERN(1,1,IB),NLO*(NLL  )*2,NREC,ISTAT)
      IF(LISTCC.GT.2.and..not.FFR) then
!!!	write(6,*) 'QERN from NREC = ',NREC,' @ ',IN,MINL1,MAXL1,L1,FI
         	CALL DISPLR(QERN(1,1,IB),NLO,NLL,NLO,SCALR)
		call flush(6)
           	SCALI = SCALR * 1e-12
         	CALL DISPLI(QERN(1,1,IB),NLO,NLL,NLO,SCALI)
	endif
C
 47      RQ  = ONE/(MLT * HF)
         CR  = (CUTOFF-1) * HF
      DO 54 I=2,NLL
      RT = (I-1)*RINTO
         RF = RT * HFHT
      TT = TC * RT**NRT
         JMIN = MAX( NLC+2 + INT((CR-RF)*RQ),1)
      IF(FFR) THEN
         DO 50 J=JMIN,NLO
   50    FNL(I,J) = FNL(I,J) + TT * QRLN(J,I,WHERE(L1,IN))
     &                     * (RF + DNL(J))**NRF
      ELSE
         DO 51 J=JMIN,NLO
   51    FNC(I,J) = FNC(I,J) +      QERN(J,I,WHERE(L1,IN))
     &                     *(TT* (RF + DNL(J))**NRF )
      ENDIF
      IF(NLPL.LE.0.AND.JTEST.GT.4) GO TO 54
        IF(FFR) THEN
         DO 52 J=JMIN,NLO
   52    VINT = MAX(VINT,ABS(FNL(I,J)))
        ELSE
         DO 53 J=JMIN,NLO
   53    VINT = MAX(VINT,ABS(FNC(I,J)))
        ENDIF
54    CONTINUE
55    CONTINUE
60    CONTINUE
      VINT = VINT * 3.0
      ALOSS(IN) = MAX(ALOSS(IN),VINT)
      ALOSSM = MAX(ALOSSM,VINT)
80    CONTINUE
!      WRITE(48,90) D,C,CP
!      WRITE(48,99) ALOSSM,ALOSSM*RER*100.
!      written(48) = .true.
      IF(DRY) GO TO 100
	if(TRENEG>2) then
 	 IA = EXCIT(C2,1); IB = EXCIT(C1,1)
!! 	if(C1FR) then
!! 	 IA = EXCIT(C1,1); IB = EXCIT(C2,1)
!!	 endif
!        WRITE(87,355) NLL,RINTP,0,0,0,IB,IA,' NL transfer potential'
355      FORMAT(I4,F8.5,'      0.      1.',5I4,1x,A)
         WRITE(87,356) NLL,RINTP,NLO,MLT * HF,NLC
356      FORMAT(I4,F8.5,I4,f8.5,I4)
!	 write(87,358) C1,C2,LVAL(C1),JVAL(C1),LVAL(C2),JVAL(C2),IB,IA,
!    x               REV,JTOTAL,PSIGN(PARITY+2)
358	 format('# Nonlocal transfer coupling for chs',2i4,' (L,Jp=',
     x   I4,f5.1,'<',I4,F5.1,') to ex# ',i4,' from',i4,1x,L1,
     x           ' in J/pi =', F5.1,A1)
	 write(87,359) JTOTAL,PARITY,IB,LVAL(C1),JVAL(C1),
     x               IA,LVAL(C2),JVAL(C2),
     x               REV,C1,C2,JTOTAL,PSIGN(PARITY+2)
359	 format(f8.1,i4,2(2i4,f6.1),L2,
     x       ': Nonlocal couplings for chs',2i4,' in J/pi =', F5.1,A1)
232	   format(1p,6e12.4)
	endif

      IF(.NOT.PRES) NL = NL - 1
      IF(.NOT.PRES) GO TO 100
      IF(FFR) THEN
             WRITE(12) FNL
	  do II=1,NLL
	   if(TRENEG>2) write(87,232) (FNL(II,J)+(0.,0.),J=1,NLO)
          enddo
      ELSE
              WRITE(12) FNC
	  do II=1,NLL
	   if(TRENEG>2) write(87,232) (FNC(II,J),J=1,NLO)
          enddo
      ENDIF
C
C     IF(.NOT.(NLPL.GT.0.AND.NL.LE.1)) GO TO 100
      IF(NLPL.LE.0) GO TO 100
      WRITE(KO,90) D,C,CP
  90  FORMAT(' NL interaction V-c(',I4,') c(',I4,') of Cplg',I3,' is'/)
      IF(FFR) CALL DISPLY(FNL,NLL,NLO,NLL,SCALE)
      IF(.NOT.FFR) THEN
         CALL DISPLR(FNC,NLL,NLO,NLL,SCALR)
           SCALI = SCALR * RER * 1E3
         CALL DISPLI(FNC,NLL,NLO,NLL,SCALI)
         SCALE = MAX(SCALR,SCALI)
         ENDIF
      IF(SCALE.ne.0.) then
      NLPL = NLPL - 1
       IF(FFR) THEN
         DO 95 I=2,NLL
         FNL(I,1) = 0
         DO 95 J=1,NLO
95       FNL(I,1) = FNL(I,1) + FNL(I,J) * HF *MLT
         WRITE(KO,98) (FNL(I,1),I=2,NLL)
98       FORMAT(1X,18F7.2)
       ENDIF
       endif
      WRITE(KO,99) ALOSSM,ALOSSM*RER*100.
99    FORMAT(1X,'Largest intermediate sum =',1P,E12.4,
     X ', so expect errors of',F10.4,' %')
100   CONTINUE
      IF(NL.GT.MFNL)  WRITE(KO,9) NL,MFNL
      IF(MAXLRQ.GT.MAXL) WRITE(KO,35) MAXLRQ,'MAX',MAXL
      IF(MINLRQ.LT.MINL) WRITE(KO,35) MINLRQ,'MIN',MINL
!     if(FFR) then
!	deallocate (QRLN,FNL)
!	else
!     	deallocate (QERN,FNC)
!	endif
      RETURN
      END
      SUBROUTINE QERNEL(IFL,NLN,NLM,NLO,MAXL1,NK,A,B,P,Q,FPT,NG,
     &             IC7,ICV,IP1,IREM,NNT,MINL1,RINTO,EPC,HNL,NLT,
     &             VFOLD,VCORE,RINS,RINC,NLL,NIBL,cxwf,FORML,VSP,
     &             FORMC,CSP,MINT,RIN,NLC,WID,CENTR,CENTRE,NONO,HFHT,
     &             FFREAL,JACOBN,LTRANS)
	use parameters
	use io
	use drier
	use trace
	use fresco1, only:rnl
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 JACOBN,COFFIN(4),XG(6),WG(6),C2,CI,WD,TH
      REAL*8,allocatable:: QRLN(:,:,:)
      COMPLEX*16,allocatable:: QERN(:,:,:)
      INTEGER FPT(7,NK),KLIM(MAXNLN)
      COMPLEX*16 VCOR(MAXNNU),VSUBI,DV,FFC,WC,SUM(MAXL1),WRT,WRP
      LOGICAL FFREAL,LTRANS(NG),THREV,VCC,cxwf
      REAL*8 RP4,NLOC(NLM),VR,FFR,FFR4,RCOR2,VR1,VR2,
     X       RUM(MAXL1),VSP(MAXNLR,MSP),FORML(MAXNLR,MSP)
      COMPLEX*16 CSP(MAXNLC,MSP),FORMC(MAXNLC,MSP),FFC4
      COMPLEX*16 VFOLD(MAXNLN),VCORE(MAXNLN),DV0,DV1,VSUBF
      REAL*8 PLEG(MAXNNU,MAXMUL),UK(MAXNNU),GK(MAXNNU),THM(MAXNLN),
     &       WT(MAXNNU),RN(MAXNNU),RDINT(MAXNNU),RCOR(MAXNNU)
      DATA NWW,XG(4),XG(5),XG(6),WG(4),WG(5),WG(6)/3,
     1   .2386191861D0,.6612093865D0,.9324695142D0,
     2   .4679139346D0,.3607615730D0,.1713244924D0/, PI /3.14159D0/
      CALL CCTIME(IPTIME)
      if(FFREAL) then
	allocate (QRLN(NLO,NLL,NIBL))
	else
      	allocate (QERN(NLO,NLL,NIBL))
	endif
      THMN = 0.1 *PI/180.  ! at least 0.1 degres, and moreover non-zero!
      EP = EPC * 0.01
      IF(EP.EQ.0.) EP = 1.D-5
      HIRS = 1.0/RINS
      HIRC = 1.0/RINC
      RDRPM = 1.0/((NLL-1)*RINTO)**2
      VR0 = 0.; VR1 = 0.

      NW=2*NWW
      DO 1 N=1,NWW
      NN=NW-N+1
      XG(N)=-XG(NN)
  1   WG(N)=WG(NN)
C
           MAXL = MAXL1 - 1
           IC3 = NNT/NW
           NNU = IC3 * NW
            CI = 1D0/IC3
            C1 = 0.5D0 * CI
      DO 10 I=1,NLM
10    NLOC(I) = 0.0
      THREV = IP1.LT.-1
      VCC   = IP1.LE.-3
      PA = P - A
      QB = Q - B
      PQ = P*HFHT + Q
      AB = A*HFHT + B
      PAQB2 = PA*QB*2.
      PAQB  = PA*HFHT+QB
      AB2 = A*B*2
      PQ2 = P*Q*2
      C = 0.5 * JACOBN
      NREC = IFL
      NBLOCK = (MAXL1-MINL1+1 + NIBL-1) / NIBL
      DO 150 IGG=1,NG
      IF(.NOT.LTRANS(IGG)) GO TO 150
        IMINL1 = MINL1
         THMAX = PI
         XTI = 0.0
      DO 120 IBL=1,NBLOCK
        IMAXL1 = MIN(IMINL1 + NIBL - 1, MAXL1)
        IMAXL  = IMAXL1 - 1
        IMINL  = IMINL1 - 1
      DO 100 I=2,NLL
         RP = (I-1) * RINTO
       IF(IBL.EQ.1) THEN
         IF(I.LE.3) GO TO 11
         IF(KLIM(I-1).EQ.NNU) THMAX = MIN(THMAX*1.333,PI)
         IF(KLIM(I-1).LT.NNU.AND..NOT.THREV)
     X                             THMAX = ACOS(UK(KLIM(I-1))+1.)
         IF(KLIM(I-1).LT.NNU.AND.THREV)
     X                        THMAX = PI - ACOS(UK(KLIM(I-1))+1.)
	 THMAX = max(THMAX,THMN)  ! at least non-zero!
 11      KLIM(I) = NNU/2
         THM(I) = THMAX
         XTL = XTI
         XTI = 0.0
       ELSE
        THMAX = THM(I)
       ENDIF
            C2 = C1
               K = 0
            DO 14  J=1,IC3
            DO 13 NN=1,NW
               K = K + 1
               X = C1 * XG(NN) + C2
               TH = (3D0*X*X + 1D0)*X*THMAX*0.25D0
               IF(THREV) TH = PI - TH
               UK(K) =  COS(TH) - 1D0
               GK(K) = C1 * WG(NN)
               WT(K) =  SIN(TH) * (9D0 *X*X+1D0 )*THMAX*0.25D0
               PLEG(K,1) = 1D0
               if(MAXMUL>1) PLEG(K,2) = UK(K) + 1D0
13             CONTINUE
14           C2 = C2 + CI
               DO 15 L=2,IMAXL
               DO 15 K=1,NNU
               WD = UK(K)  + 1D0
15         PLEG(K,L+1) = ((2*L-1)*WD*PLEG(K,L-1+1) -(L-1)*PLEG(K,L-2+1))
     &                       / DBLE(L)
      VSUBF = 0.0
      DO 25 K=1,NNU
25    VCOR(K) = 0.0
         IF(IREM.NE.0.AND.IC7.ne.1) then
	   if(.NOT.FFREAL) VSUBF = FFC(RP*HIRS,VFOLD,NLN)   ! Vopt subtracted in post
	   if(     FFREAL) VSUBF = FFR(RP*HIRS,VFOLD,NLN)
	 endif
         DO 29 JJ=1,NLO
      IF(FFREAL) THEN
         DO 27 L1=IMINL1,IMAXL1
27       QRLN(JJ,I,L1-IMINL) = 0.0
      ELSE
         DO 28 L1=IMINL1,IMAXL1
28       QERN(JJ,I,L1-IMINL) = 0.0
      ENDIF
29    CONTINUE
C
         G = 4./3.
         G = 1.
         DO 85 J=1,NLM
            DNL =     (J - NLM/2 - 1) * HNL  +  CENTRE
            RT = RP*HFHT
            RDCM = RT    + DNL
            PQR = (P * DNL + PQ * RP)**2
            ABR = (A * DNL + AB * RP)**2
            PAQBR = (PA * DNL + (PAQB)*RP)**2
            RDRP = RDCM*RP
            DO 30 L1=IMINL1,IMAXL1
               IF(FFREAL) RUM(L1)      = 0.
30             IF(.NOT.FFREAL) SUM(L1) = 0
            IF(RDCM.LE.0.01) GO TO 60
            DO 31 K=1,NNU
               RDINT2 =     PQR + PQ2*RDRP*UK(K)
               RN2  =       ABR + AB2*RDRP*UK(K)
               RDINT(K) = SQRT(ABS(RDINT2)) * RIN
31             RN(K) =    SQRT(ABS(RN2))    * RIN
               IF(IREM.EQ.0) GO TO 32
            IF(IC7.ne.0.AND..NOT.FFREAL) VSUBI= FFC(RDCM*HIRS,VFOLD,NLN)   ! prior
            IF(IC7.ne.0.AND.     FFREAL) VSUBI= FFR(RDCM*HIRS,VFOLD,NLN)
            IF(IC7.eq.0)                 VSUBI= VSUBF			   ! post
	    if(IC7.eq.4) VSUBI = VSUBI - VSUBF			   	   ! prior - post
               IF(VCC) VSUBI = 0.0
                  DO 315 K=1,NNU
                  RCOR2 = PAQBR + PAQB2*RDRP * UK(K)
 315              RCOR(K) = SQRT(ABS(RCOR2)) * HIRC
32             DO 50 K=1,NNU
		   VCOR(K) = 0.0
                 IF(IREM.ne.0) then
 		   if(IC7.ne.4) VCOR(K)=FFC(RCOR(K),VCORE,NLN)
                   VCOR(K)=VCOR(K) - VSUBI
	         endif
               DO 38 IN=1,NK
                   IG  = FPT(3,IN)
                   IF(IG.NE.IGG) GO TO 38
                   IFT = FPT(2,IN)
                   IFP = FPT(1,IN)
		   if(.not.cxwf) then
                   WRT= FFR4(RN(K),FORML(1,IFT),MINT)
                   WRP= FFR4(RDINT(K),FORML(1,IFP),MINT)
		   else
                   WRT= FFC4(RN(K),FORMC(1,IFT),MINT)
                   WRP= FFC4(RDINT(K),FORMC(1,IFP),MINT)
		   endif
         IF(.NOT.FFREAL) THEN
C                            Complex form factors
                   IF(NONO.GT.0) THEN
                      DV = WRT * WRP
                      ELSE
		  DV0=0.; DV1=0.;
	 if(.not.cxwf) then
            IF(ICV.ne.1) DV0=FFR4(RDINT(K),VSP(1,IFP),MINT) * WRT  ! with post vertex function
            IF(ICV.ne.0) DV1=FFR4(RN(K)   ,VSP(1,IFT),MINT) * WRP  ! with prior
	 else
            IF(ICV.ne.1) DV0=FFC4(RDINT(K),CSP(1,IFP),MINT) * WRT  ! with post vertex function
            IF(ICV.ne.0) DV1=FFC4(RN(K)   ,CSP(1,IFT),MINT) * WRP  ! with prior
	 endif
		if(IC7<4)  DV = DV1 + DV0    ! only 1 in non-zero
		if(IC7==4) DV = DV1 - DV0    ! prior-post
                IF(VCC) DV = 0.0
                IF(IREM.ne.0) DV = DV + VCOR(K) *WRT*WRP
                      ENDIF
               WC= DV * WT(K)
                 IF(IBL.EQ.1) THEN
                   W = WC
                   IF(ABS(W).GT.EP*XTL) KLIM(I) = MAX(KLIM(I),K)
                   XTI = MAX(XTI,ABS(W))
                 ENDIF
                WC = WC * GK(K)
C
               DO 35 L1=IMINL1,IMAXL1
               SUM(L1) = SUM(L1) + WC * PLEG(K,L1)
35             CONTINUE
            ELSE
C                        Real form factors
                   IF(NONO.GT.0) THEN
                      VR = WRT * WRP
                      ELSE
         IF(ICV.ne.1) VR0=FFR4(RDINT(K),VSP(1,IFP),MINT) * WRT
         IF(ICV.ne.0) VR1=FFR4(RN(K)   ,VSP(1,IFT),MINT) * WRP
                if(IC7<4)  VR = VR1 + VR0    ! only 1 in non-zero
                if(IC7==4) VR = VR1 - VR0    ! prior-post
                IF(VCC) VR = 0.0
                IF(IREM.NE.0) VR = VR + DBLE(VCOR(K)) *WRT*WRP
                      ENDIF
               WD= VR * WT(K)
                IF(IBL.EQ.1) THEN
                   W = WD
                   IF(ABS(W).GT.EP*XTL) KLIM(I) = MAX(KLIM(I),K)
                   XTI = MAX(XTI,ABS(W))
                 ENDIF
                WD = WD * GK(K)
               DO 36 L1=IMINL1,IMAXL1
36             RUM(L1) = RUM(L1) + WD * PLEG(K,L1)
            ENDIF
38      CONTINUE
C
50          CONTINUE
              C2 = RDRP * RDRPM
               IF(FFREAL) THEN
                 DO 57 L1=IMINL1,IMAXL1
57                NLOC(J)=NLOC(J) +    (RUM(L1-IMINL))**2 * C2
               ELSE
                 DO 58 L1=IMINL1,IMAXL1
58                NLOC(J)=NLOC(J) + ABS(SUM(L1-IMINL))**2 * C2
               ENDIF
60          CALL SPLINT( DNL/(HNL*NLT) + NLC    ,NLO,IJ,NP,COFFIN)
            G = 2. - G
            DO 82 M=1,NP
            WD  = (COFFIN(M) * G / NLT) * C * RDRP
            JJ = IJ + M - 1
         IF(FFREAL) THEN
            DO 79 L1=IMINL1,IMAXL1
79          QRLN(JJ,I,L1-IMINL) = QRLN(JJ,I,L1-IMINL) + RUM(L1) * WD
         ELSE
            DO 80 L1=IMINL1,IMAXL1
80          QERN(JJ,I,L1-IMINL) = QERN(JJ,I,L1-IMINL) + SUM(L1) * WD
         ENDIF
82       CONTINUE
85            CONTINUE
100      CONTINUE
              DO 90 L1=IMINL1,IMAXL1
              NREC = (IGG-1)*(MAXL1-MINL1+1) + L1 - MINL1  + 1 + IFL
                 IF(DRY) GO TO 90
            IF(FFREAL) THEN
             WRITE(11,REC=NREC) ((QRLN(JJ,I,L1-IMINL),JJ=1,NLO),I=2,NLL)
!           CALL FIOW(11,QRLN(1,1,L1-IMINL),NLO*(NLL  ),NREC,ISTAT)
            ELSE
             WRITE(9 ,REC=NREC) ((QERN(JJ,I,L1-IMINL),JJ=1,NLO),I=2,NLL)
!           CALL FIOWC( 9,QERN(1,1,L1-IMINL),NLO*(NLL  )*2,NREC,ISTAT)
!!!	   write(6,*) 'WRITE 9, NREC,L1,NLO,NLL = ',NREC,L1,NLO,NLL
!!!         	CALL DISPLR(QERN(1,1,L1-IMINL),NLO,NLL,NLO,SCALR)
!!!           	SCALI = SCALR * 1e-12
!!!         	CALL DISPLI(QERN(1,1,L1-IMINL),NLO,NLL,NLO,SCALI)
            ENDIF
90         CONTINUE
120    IMINL1 = IMAXL1 + 1
150    CONTINUE
      IFL = IFL + NG * (MAXL1-MINL1+1)
      CALL CCTIME(IPTIME)
       PT = IPTIME * 0.01
      IF(LISTCC.GE.3) WRITE(KO,992) PT,(THM(I)*180/PI,I=2,NLL)
992   FORMAT(/' Theta - maxima (after',F6.2,' secs) are'/(1X,20F6.1))
!992   FORMAT(/' Theta - maxima (after',F6.2,' secs) are'/(1X,1p,10e10.2))
      RP4 = 0.0
      DO 210 J=1,NLM
210   RP4 = MAX(NLOC(J),RP4)
      VR = 100.*EP * 0.1
      IF(ABS(RP4).LT.1.E-30) WRITE(KO,992) PT,(THM(I)*180/PI,I=2,NLL)
      IF(ABS(RP4).LT.1.E-30) GO TO 219
      IFT = 10000
      IFP = -10000
      DO 215 J=1,NLM
            DNL = (J - NLC - 1) * HNL
      NLOC(J) = SQRT( NLOC(J) / RP4) * 100.
      IF(NLOC(J).LE.VR) GO TO 215
         IFP = MAX(IFP,J)
         IFT = MIN(IFT,J)
215   CONTINUE
      IF(NLOC(1).GT.VR)
     X IFT = 1 - MAX(INT(LOG(NLOC(1)/VR)/LOG(NLOC(2)/NLOC(1))),0)
      IF(NLOC(NLM).GT.VR)
     X IFP =NLM+MAX(INT(LOG(NLOC(NLM)/VR)/LOG(NLOC(NLM-1)/NLOC(NLM))),0)
      WID  = (IFP - IFT) * HNL
      CENTR = ((IFP + IFT) * 0.5 - NLM/2-1) * HNL + CENTRE
219   WRITE(KO,220) WID,CENTR,PT,(NLOC(J),J=1,NLM)
220   FORMAT('0RECOMMENDED NON-LOCAL WIDTH IS GREATER THAN',F8.2,' FM.',
     &     ',  RECOMMENDED CENTRATION ',F6.2,',',
     &' after',F8.2,' secs.'/'0Relative non-local usages (rms) are'
     &                       /(1X,15F8.3))
      if(WID>rnl*1.5) then
	write(KO,221)  rnl,WID
  221 	format(' Non local width RNL',f8.3,' is too small: INCREASE TO',
     x              ' >',f8.3,', so STOP!')
!	  stop
	endif
      if(FFREAL) then
	deallocate (QRLN)
	else
      	deallocate (QERN)
	endif
      RETURN
      END
      SUBROUTINE SOURCE(INHOMG,PSI,N,H,NCH,IEX,FORMF,NF,FORMFR,CUTOFF,
     &  ICUTC,SIMPLE,ITNL,NLN,NLO,MR,NL,EMPTY,SAME,SH,WAVES,LVAL,
     &  MLT,SKIP,FED,NLC, CHNO,MLM,NICH,NAXICH,PTYPE,LOCFIL,
     &  CLIST,NCLIST,NFLIST)
	use parameters
	use io
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 PSI(MAXN,NICH),INHOMG(MAXN,MAXCH),FORMF(MAXM,MLOC),
     &      FORMFR(NF),PSID(MAXN),
     &      CLIST(MAXCH,MAXCH,MCLIST),S,ECF(MAXNLN),FNC(NLN,NLO),CI,PH
      INTEGER C,D,WAVES,CUTOFF,CHNO(MFNL,6),LVAL(MAXCH),PTYPE(12,NF),
     & 	    NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH)
      REAL*8 EPS,EXF(NLN),FNL(NLN,NLO)
      REAL*8 H(NCH)
      LOGICAL ITNL,EMPTY(MAXCH),SAME(MAXCH),SKIP(MAXCH),NREV,NFOR,SH,
     &        SIMPLE(MAXCH),FFR,FED(MAXCH),CP,LOCFIL,DERIV
      DATA EPS / 1E-12 /, CI / (0.0,1.0) /
C
      DO 3 C=1,NCH
        SKIP(C) = SIMPLE(C)
        FED(C) = .FALSE.
        DO 3 D=1,NCH
        DO 3 NC=1,NCLIST(C,D)
	 IF = NFLIST(C,D,NC)
         IF((C.EQ.D .OR. C.LE.IEX .AND. D.LE.IEX)
     X      .AND.mod(PTYPE(6,IF),2).EQ.0) GO TO 3
            CP = ABS(CLIST(C,D,NC)).GT.EPS
        SKIP(C) = SKIP(C) .AND. (SAME(D).OR..NOT.CP)
        FED(C) = FED(C) .OR. CP.AND..NOT.EMPTY(D)
3     CONTINUE
      ICH = 0
      DO 4 INL=1,NL
            IF(MOD(CHNO(INL,3),10).LE.1) GO TO 4
         D = ABS(CHNO(INL,1))
         C = CHNO(INL,2)
         SKIP(D) = SKIP(D) .AND. SAME(C)
         SKIP(C) = SKIP(C) .AND.(SAME(D) .OR. CHNO(INL,1).LT.0)
         FED(D) = FED(D) .OR. .NOT.EMPTY(C)
         FED(C) = FED(C) .OR. .NOT.EMPTY(D)  .AND.CHNO(INL,1).GT.0
        IF(.NOT.EMPTY(C)) ICH = MAX(C,ICH)
        IF(.NOT.EMPTY(D).AND.CHNO(INL,1).GT.0) ICH = MAX(D,ICH)
4     CONTINUE
      NAXICH = MAX(ICH,NAXICH)
!      CALL CHECK(ICH,MAXICH,-25)
	if(ICH>MAXICH) then
	write(KO,7) ICH,MAXICH+1,ICH
7	format(//'  ***** Channels up to',i4,' have unexpected flux,'/
     X		 '        probably from long-range Coulomb couplings.'/
     X		 '        FLUX IN CHANNELS ',i4,' TO',i4,' IS IGNORED'/)
	SH=.true.
	endif
      IF(SH)WRITE(KO,8) (C,EMPTY(C),SAME(C),SIMPLE(C),
     X                 SKIP(C),FED(C),C=1,NCH)
8     FORMAT(' #,EMPTY,SAME,SIMPLE,SKIP,FED =',12(1X,I2,5L1,';') )
C
      DO 10 C=1,NCH
      IF(SKIP(C)) GO TO 10
      DO 9 I=1,N
9     INHOMG(I,C) = 0.0
10    CONTINUE
C
      if(.not.LOCFIL) then
      DO 15 D=1,NCH
      DO 15 C=1,NCH
      DO 15 NC=1,NCLIST(C,D)
	 JF = NFLIST(C,D,NC)
C        IF(C.EQ.D .OR. C.LE.IEX .AND. D.LE.IEX) GO TO 15
         IF((C.EQ.D .OR. C.LE.IEX .AND. D.LE.IEX)
     X      .AND.mod(PTYPE(6,JF),2).EQ.0) GO TO 15
         IF(D.gt.NICH) go to 15
         S = CLIST(C,D,NC)
         IF(ABS(S).LT.EPS.OR.SKIP(C).OR.EMPTY(D)) GO TO 15
C     IF(SH) WRITE(KO,11) C,D,JF,S,IEX
      IF(ABS(WAVES).GE.2) WRITE(41,11) C,D,JF,S,IEX
11    FORMAT('#ZR to',I3,' fr',I3,' by',I3,' of',2E14.3,'(IEX=',I3,')')
         DO 12 I=max(1,ICUTC),N
      INHOMG(I,C) = INHOMG(I,C) + S * FORMF(I,JF) * PSI(I,D)
      IF(ABS(WAVES).GE.2.AND.ABS(WAVES).le.3)
     X   WRITE(41,13) (I-1)*H(D) , INHOMG(I,C),PSI(I,D),FORMF(I,JF)
12    CONTINUE
13     FORMAT(F8.4,1P,6E12.3)
      IF(ABS(WAVES).GE.2.AND.ABS(WAVES).le.3) write(41,*) '&'
15    CONTINUE
	else
      rewind 19
      IS=max(1,ICUTC)
	 do I=1,IS-1
	 read(19)
	 enddo
      DO 26 I=IS,N
	 read(19) FORMFR
      DO 25 D=1,NCH
      DO 25 C=1,NCH
      DO 25 NC=1,NCLIST(C,D)
	 JF = NFLIST(C,D,NC)
         IF((C.EQ.D .OR. C.LE.IEX .AND. D.LE.IEX)
     X      .AND.mod(PTYPE(6,JF),2).EQ.0) GO TO 25
         IF(D.gt.NICH) go to 25
         S = CLIST(C,D,NC)
         IF(ABS(S).LT.EPS.OR.SKIP(C).OR.EMPTY(D)) GO TO 25
         IF(ABS(WAVES).GE.2.and.I==IS) WRITE(41,11) C,D,JF,S,IEX
      INHOMG(I,C) = INHOMG(I,C) + S * FORMFR(JF) * PSI(I,D)
25    CONTINUE
26    CONTINUE
	endif
C
      IF(.NOT.ITNL.OR.NL.EQ.0) GO TO 51
      HI = 1.0 / DBLE(MR)
       REWIND 12
      DO 50 INL=1,NL
         FFR = CHNO(INL,3).LT.10
       NLL = CHNO(INL,4)
          IF(FFR) READ(12) ((FNL(I,J),I=1,NLL),J=1,NLO)
          IF(.NOT.FFR) READ(12) ((FNC(I,J),I=1,NLL),J=1,NLO)
         D = ABS(CHNO(INL,1))
         C = CHNO(INL,2)
         DERIV = CHNO(INL,5)>0
         PH = CI ** (LVAL(C) - LVAL(D))
            IF(MOD(CHNO(INL,3),10).LE.1) GO TO 50
         NFOR = SKIP(D) .OR. EMPTY(C)
         NREV = SKIP(C) .OR. EMPTY(D) .OR. CHNO(INL,1).LT.0 .OR. C.EQ.D
         IF(SH) WRITE(KO,28) INL,D,C,NICH,.NOT.NFOR,.NOT.NREV,FFR,NLL,
     &                       DERIV
28       FORMAT(' NL coupling #',I4,' to',I4,' from Ch.',I3,'<=',i3,
     &         ', Forw =',L2,', Reverse=',L2,', Real =',L2,' NLL=',i4,
     &         ', Deriv=',L1)
         IF(NFOR.AND.NREV) GO TO 50
            SQH = SQRT(H(C)/H(D))
C        IF(SH.AND.FFR) CALL DISPLY(FNL,NLL,NLO,NLN,SCALE)
!          CALL DISPLR(FNC,NLL,NLO,NLL,SCALR)
         IF(NFOR .OR. C.GT.NICH) GO TO 32
         
         if(DERIV) then
!          PSID = deriv(PSI(*,C)) - (L+1)/r*PSI(*,C)
          CALL DERIVC(PSI(:,C),PSID,H(C),N)   
           do I=2,N
            PSID(I) = PSID(I) - (LVAL(C)+1)/((I-1)*H(C)) * PSI(I,C)
            enddo        
         endif
         
         DO 31 J=1,MLM  ! do FORWARD
            JJ = (J - NLC*MLT  - 1)
            IMIN = 1 + MAX(-JJ ,0)  +  MAX(CUTOFF,ICUTC)
            IMAX = N - MAX(JJ ,0)   -  5
            KMIN = MAX((IMIN-1)/MR+2,2)
            KMAX = MIN(IMAX/MR,NLL-2)
            DO 31 II=1,MR
               P = (II-1)*HI
              P1 = P - 1.
              P2 = P - 2.
              Q  = P + 1.
              X = P * P1 / 6.0  * H(C) * SQH
              Y = Q * P2 * 0.5  * H(C) * SQH
         IF(.not.DERIV) then
         IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 29 K=KMIN,KMAX
             V = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
            INHOMG(I,D) = INHOMG(I,D) + (V*PH) * PSI(I+JJ,C)
29           I = I + MR
          ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 30 K=KMIN,KMAX
             S = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
            INHOMG(I,D) = INHOMG(I,D) + (S*PH) * PSI(I+JJ,C)
30           I = I + MR
          ENDIF
          else  ! act on (wf' - (L+1)/r*wf)
          
         IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 291 K=KMIN,KMAX
             V = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
            INHOMG(I,D) = INHOMG(I,D) + (V*PH) * PSID(I+JJ)
291           I = I + MR
          ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 301 K=KMIN,KMAX
             S = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
            INHOMG(I,D) = INHOMG(I,D) + (S*PH) * PSID(I+JJ)
301           I = I + MR
          ENDIF        
          
          endif
31      CONTINUE
32       IF(NREV .OR. D.GT.NICH) GO TO 50
	   if(DERIV) then
	     write(0,*) ' Reverse derivative on ch',D,' not implemented!!'
	      go to 50
	    endif
         DO 42 J=1,MLM ! do REVERSE
            JJ =-(J - NLC*MLT - 1)
            IMIN = 1 + MAX(-JJ ,0)  +  MAX(CUTOFF,ICUTC)
            IMAX = N - MAX(JJ ,0)   -  5
            KMIN = MAX((IMIN-1+JJ)/MR+2,2)
            KMAX = MIN((IMAX+JJ)/MR,NLL-2)
            DO 42 II=1,MR
               P = (II-1)*HI
              P1 = P - 1.
              P2 = P - 2.
              Q  = P + 1.
              X = P * P1 / 6.0  * H(D) / SQH
              Y = Q * P2 * 0.5  * H(D) / SQH
          IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II - JJ
            DO 38 K=KMIN,KMAX
             V = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
            INHOMG(I,C) = INHOMG(I,C) + (V*CONJG(PH)) * PSI(I+JJ,D)
38           I = I + MR
          ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II - JJ
            DO 40 K=KMIN,KMAX
             S = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
            INHOMG(I,C) = INHOMG(I,C) + S * CONJG(PH) * PSI(I+JJ,D)
40           I = I + MR
         ENDIF
42    CONTINUE
50    CONTINUE
51    IF(ABS(WAVES).GE.2.AND.ABS(WAVES).le.3) then
      DO 80 C=1,NCH
80     if(FED(C)) WRITE(KO,90) C, (INHOMG(I,C),I=1,N-1)
90    FORMAT(//' For channel no.',I3,', The Inhomogeneous driving',
     & ' Terms are'// 1P,(5(E12.3,' +i*',E9.2,',') ) )
	C = 1
      IF(.false..and.FED(C)) WRITE(331,91) 
     x    C, ((I-1)*H(C),PSI(I,C),INHOMG(I,C)/PSI(I,C),I=2,N-1)
91    FORMAT(' For channel no.',I3,', TELP  ',
     & ' is'/ (f8.3,2(2f12.6,',')) )

	DO 100 C=1,NCH
        if(FED(C)) WRITE(290,94) C, ((I-1)*H(C),INHOMG(I,C),I=1,N-3)
94    format('# Source in channel ',i4,':'/(f8.3,2e14.5))
100	WRITE(190,*) '&'
      endif
	
      RETURN
      END
      SUBROUTINE RENO(W,PHI,N,H,NCH,CUTOFF,NLN,NLO,MR,MLM,IC7, NF,
     &      LB,NL,EMPTY,MLT,NLC,SH,CLIST,NCLIST,NFLIST,CHNO,NPOST,
     &   NPRIOR,INHOMG,COEF,SIMPLE,IEX,FORMF,FORMFR,LOCFIL,ECMR,LVAL,
     &   PTYPE,ICUTC,NICH)
	use parameters
	use io
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 W(MAXN,MAXCH),PHI(N,2),INHOMG(N),S,ENS,FNC(NLN,NLO),
     &          ECF(NLN),FORMF(MAXM,NF),CI,PH,CLIST(MAXCH,MAXCH,MCLIST),
     & 		FORMFR(NF)
      REAL*8 EPS,FNL(NLN,NLO),EXF(NLN)
      REAL*8 COEF(MAXCH),H(MAXCH),ECMR(MAXCH)
      LOGICAL EMPTY(MAXCH),NFOR,NREV,SIMPLE(MAXCH),NPRIOR,FFR,NPOST,
     & 		LOCFIL
      INTEGER CHNO(MFNL,6),CUTOFF,D,C,SH,LVAL(MAXCH),C2,PTYPE(12,NF),
     & 	    NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH)
C

      EPS = 1.E-20
      CI = (0.0,1.0)
      HI = 1.0 / DBLE(MR)
      LAST1 = 0
      LAST2 = 0
      DO 10 C=1,NCH
      SIMPLE(C) = .TRUE.
10    EMPTY(C) = .FALSE.
      IF(NL.EQ.0) GO TO 51
      REWIND 12
      DO 50 INL=1,NL
         FFR = CHNO(INL,3).LT.10
       NLL = CHNO(INL,4)
          IF(FFR) READ(12) ((FNL(I,J),I=1,NLL),J=1,NLO)
          IF(.NOT.FFR) READ(12) ((FNC(I,J),I=1,NLL),J=1,NLO)
         D = ABS(CHNO(INL,1))
         C = CHNO(INL,2)
         PH = CI ** (LVAL(C) - LVAL(D))
            IF(MOD(CHNO(INL,3),10).EQ.2) GO TO 50
         NFOR = EMPTY(C) .OR. CHNO(INL,3).EQ.1-IC7
         NREV = EMPTY(D) .OR. CHNO(INL,3).EQ.IC7 .OR. CHNO(INL,1) .LT. 0
         IF(SH.GE.7) WRITE(KO,22) INL,CHNO(INL,1),C,.NOT.NFOR,
     X                                          .NOT.NREV,CHNO(INL,3)
22       FORMAT(' NONO NL coupling #',I4,' to',I4,' from Ch.',I3,', Forw
     & =',      L2,', Reverse=',L2,' of CHNO =',I3)
         IF(NFOR.AND.NREV) GO TO 50
            SQH = SQRT(H(C)/H(D))
         IF(SH.GE.8.AND.FFR) CALL DISPLY(FNL(1,1),NLL,NLO,NLN,SCALE)
         IF(NFOR) GO TO 32
            SIMPLE(D) = .FALSE.
                  IF(LAST1.NE.LB+C) READ(8,REC=LB+C) (PHI(I,1),I=1,N)
                  LAST1 = LB+C

c	add in close-coupled contributions to source term, for prior nono corrections:
          DO 23 C2=1,NCH
          DO 23 NC=1,NCLIST(C,C2)
             JF = NFLIST(C,C2,NC)
             IF(C.EQ.C2) GO TO 23
             IF(C2.gt.NICH) go to 23
             S = CLIST(C,C2,NC)
             IF(ABS(S).LT.EPS) GO TO 23
           if(C.LE.IEX .AND. C2.LE.IEX.AND.mod(PTYPE(6,JF),2).EQ.0) then
              READ(8,REC=C2) (PHI(I,2),I=1,N)
!	    T = sum(abs(PHI(:,1)))
             DO  I=max(1,ICUTC),N
              PHI(I,1) = PHI(I,1) + S * FORMF(I,JF) * PHI(I,2)
             ENDDO
!	   if(SH>7)write(ko,2299) C,C2,S,
!	           write(ko,2299) C,C2,S,
!     x             T,sum(abs(PHI(:,1))),sum(abs(PHI(:,2)))
!2299	  format(' PR-NO corr',2i4,' @',2f10.5,:,' o,n,c=',3f10.5)
           endif
23        CONTINUE


	    T = sum(abs(W(:,D)))
             JJ = N/10
           DO 24 I=1,N,JJ
24         IF(ABS(DBLE(PHI(I,1)))+ABS(AIMAG(PHI(I,1))).GT.EPS) GO TO 25
           EMPTY(C) = .TRUE.
            GO TO 32
25       DO 31 J=1,MLM
            JJ = (J - NLC*MLT  - 1)
            IMIN = 1 + MAX(-JJ ,0)  +  CUTOFF
            IMAX = N - MAX(JJ ,0)   -  5
            KMIN = MAX((IMIN-1)/MR+2,2)
            KMAX = MIN(IMAX/MR,NLL-2)
            DO 31 II=1,MR
               P = (II-1)*HI
              P1 = P - 1.
              P2 = P - 2.
              Q  = P + 1.
              X = P * P1 / 6.0  * H(C) * SQH
              Y = Q * P2 * 0.5  * H(C) * SQH
         IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 29 K=KMIN,KMAX
             V = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
            W(I,D) = W(I,D) - V * PH * PHI(I+JJ,1)
29           I = I + MR
          ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 30 K=KMIN,KMAX
             S = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
            W(I,D) = W(I,D) - S * PH * PHI(I+JJ,1)
30           I = I + MR
         ENDIF
31    CONTINUE
           write(ko,3199) D,C,T,sum(abs(W(:,D)))
3199	  format(' Src-osrc corr',2i4,' o,n=',3f10.5)
32       IF(NREV) GO TO 50
            SIMPLE(C) = .FALSE.
            IF(LAST2.NE.LB+D) READ(8,REC=LB+D) (PHI(I,2),I=1,N)
            LAST2 = LB+D
             JJ = N/10
           DO 33 I=1,N,JJ
33         IF(ABS(DBLE(PHI(I,2)))+ABS(AIMAG(PHI(I,2))).GT.EPS) GO TO 35
           EMPTY(D) = .TRUE.
            GO TO 50
35       DO 42 J=1,MLM
            JJ =-(J - NLC*MLT - 1)
            IMIN = 1 + MAX(-JJ ,0)  +  CUTOFF
            IMAX = N - MAX(JJ ,0)   -  5
            KMIN = MAX((IMIN-1+JJ)/MR+2,2)
            KMAX = MIN((IMAX+JJ)/MR,NLL-2)
            DO 42 II=1,MR
               P = (II-1)*HI
              P1 = P - 1.
              P2 = P - 2.
              Q  = P + 1.
              X = P * P1 / 6.0  * H(D) / SQH
              Y = Q * P2 * 0.5  * H(D) / SQH
         IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II - JJ
            DO 39 K=KMIN,KMAX
             V = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
            W(I,C) = W(I,C) - V * CONJG(PH) * PHI(I+JJ,2)
39           I = I + MR
         ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II - JJ
            DO 40 K=KMIN,KMAX
             S = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
            W(I,C) = W(I,C) - CONJG(S * PH)* PHI(I+JJ,2)
40           I = I + MR
         ENDIF
42    CONTINUE
50    CONTINUE
51    DO 53 C=1,NCH
53    IF(IC7.EQ.1.AND.NPOST.AND..NOT.SIMPLE(C))
     #    WRITE(8,REC=NCH+C) (W(I,C),I=1,N)
      IF(NPRIOR.AND.IC7.EQ.0) then
      IF(SH.GE.6) WRITE(KO,55) IC7,(SIMPLE(C),C=1,NCH)
55    FORMAT('0RENO',I3,' : sources still simple are ',   40L2)
C56    FORMAT('0RENO',I3,' : empty channels were      ',   40L2)
      DO 70 C=1,NCH
         IF(SIMPLE(C)) GO TO 70
C		Mistake in original theory:
C	   	Should be previous PHI, not current W
         READ(8,REC=C) (PHI(I,1),I=1,N)
      INHOMG(:) = 0.0
      IMIN = 1+ CUTOFF + 3 + MR
      IMAX = N - 5 - 3        - MR
      T = COEF(C) /(H(C)**2  * 12.0)
       if(LOCFIL) then
	 rewind 19
	 do i=1,IMIN-1
	 read(19) 
	 enddo
	 endif
      DO 60 I=IMIN,IMAX
         ENS = - LVAL(C)*(LVAL(C)+1)*COEF(C)/(H(C)*(I-1))**2
     X         - ECMR(C)
	if(LOCFIL) then
	 read(19) FORMFR
         DO 58 NC=1,NCLIST(C,C)
	   JF = NFLIST(C,C,NC)
58       ENS = ENS + CLIST(C,C,NC) * FORMFR(JF)
	else
         DO 59 NC=1,NCLIST(C,C)
	   JF = NFLIST(C,C,NC)
59       ENS = ENS + CLIST(C,C,NC) * FORMF(I,JF)
	endif
60    INHOMG(I) = - (ENS * PHI(I,1) +
     & T * (-PHI(I+2,1)+16.*PHI(I+1,1)-30.*PHI(I,1)+
     & 	                16.*PHI(I-1,1)-PHI(I-2,1)))
      WRITE(8,REC=NCH+C) (INHOMG(I),I=1,N)
70    CONTINUE
	endif
      RETURN
      END
      SUBROUTINE CHECK(I,LIM,KIND)
	use io
	use drier
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*19  LI(30)
      CHARACTER*6 WHO(30)
      DATA LI /
     &'EXCITATION PAIRS   ','MASS PARTITIONS    ','BOUND STATES       ',
     &'COUPLING TYPES     ','CHANNELS           ','NON-LOCAL FORMS.   ',
     &'C.M. RADIAL POINTS ','CM INTERPOLATN PTS.','NL INTERPOLATN PTS ',
     &'PW COUPLING SLOTS  ','L-VALS & MULTIPOLES','ANGULAR QUADRTR PTS',
     &'KERNELS/COUPLING   ','MATRIX ROWS:ERWIN  ','BLOCKS OF EQUATIONS',
     &'COUPLED PAIRS OF WF','                   ','P(L,M) VALUES.     ',
     &'                   ','PLACES IN FAM()    ','TOTAL NO. STATES   ',
     &'PROJECTILE M VALUES','IN YSIG ARRAY,     ','LOCAL FORMS        ',
     &'INTERMEDIATE CHANS.','MRXY in EXTERN2    ','MM in EXTERN2      ',
     &'arrays in EXTERN1  ','                   ','COUPLINGS/PW-PAIR  '/
      DATA WHO /'MXX   ','MXP   ','MSP   ','MAXCPL','MAXCH ','MFNL  ',
     &          'MAXN  ','MAXNLN','MAXNLO','MPWCOU','LMAX1 ','MAXNNU',
     &          'MAXQRN','MAXNR ','MAXB  ','MPAIR ','NFDEC ','PL-DEC',
     &          '      ','MAXF  ','MXPEX ','MSPIN ','MXYSIG','MLOC  ',
     &          'MAXICH','MRXY  ','MM    ','there ','      ','MCLIST'/
      IF(I.LE.LIM) RETURN
      WRITE(KO,833) LIM,LI(ABS(KIND)),I,WHO(ABS(KIND))
833   FORMAT(' ****** THERE IS ONLY ROOM FOR',I15,1X,A19,' BUT',I15,
     & ' ARE REQUIRED,  SO INCREASE PARAMETER ',A6,' !!')
      DRY = .TRUE.
	write(KO,20)
20	format(/'  *** INTERNAL ERROR IN ALGORITHM TO CALCULATE',
     X ' PARAMETERS!! ***',/,'  *** Please report this (along with',
     X ' a copy of the input file)',/,
     X '      to Ian@kernz.org ***')
      IF(KIND.LT.0) RETURN
C     STOP
      CALL ABEND(8)
      END
      SUBROUTINE CCTIME(I)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /TIMER/ START
      T = SECOND()*100
C
C
      I = T - START
      START = T
      RETURN
      END
      SUBROUTINE DISPLR(A,M,N,MA,DMAX)
      use io
      REAL*8 EMAX
      COMPLEX*16 A(MA,N)
      REAL*8 DMAX
      CHARACTER CHARS(11),SIGNS(3),LINE(132),PLUS,CAPI
      LOGICAL P
      DATA  CHARS     / ' ','1',' ','3',' ','5',' ','7',' ','9','A' /
      DATA  SIGNS    / '.','0',' ' /
      DATA PLUS,CAPI / '+','I' /
      MP= MIN(M,130)
      EMAX = 0
      DO 10 I=1,N
      DO 10 J=1,MP
10    EMAX = MAX(EMAX, ABS(DBLE(A(J,I))) )
      DMAX = EMAX
      IF(EMAX.EQ.0) GOTO 100
      SCALE = 10/EMAX
C20   FORMAT(// 20('*') //)
      IF(ABS(LOG10(EMAX)).LE.4) THEN
      WRITE(KO,12) EMAX
       ELSE
      WRITE(KO,13) EMAX
       ENDIF
12    FORMAT(' Real part: full scale =',F12.5/ )
13    FORMAT(' Real part: full scale =',1P,E12.4/ )
      DO 50 I=1,N
      DO 30 J=1,MP
         ANO= SCALE * ABS(DBLE(A(J,I)))
         NO = ANO
         NO = max(NO,0)
         LINE(J) = CHARS(NO+1)
         IF(NO.LT.1 .AND. ANO .GE. 0.10) LINE(J) = PLUS
30    CONTINUE
      WRITE(KO,35) (LINE(K),K=1,MP),CAPI
35    FORMAT(' I', 131A1)
      P = .FALSE.
      DO 40 J=1,MP
         IS = SIGN(1D0,DBLE(A(J,I)))
         P = P .OR. IS.LT.0
         P = P .OR. IS.EQ.0
40       LINE(J) = SIGNS(IS+2)
!     IF(P) WRITE(KO,45)   (LINE(J),J=1,MP)
45    FORMAT('+ ',130A1)
50    CONTINUE
C     WRITE(KO,20)
      RETURN
100   WRITE(KO,105)
105   FORMAT(/' Real part everywhere zero???'/)
      RETURN
      END
      SUBROUTINE DISPLI(A,M,N,MA,DMAX)
      use io
      REAL*8 EMAX
      COMPLEX*16 A(MA,N)
      REAL*8 DMAX
      CHARACTER CHARS(11),SIGNS(3),LINE(132),PLUS,CAPI
      LOGICAL P
      DATA  CHARS     / ' ','1',' ','3',' ','5',' ','7',' ','9','A' /
      DATA  SIGNS    / '.','0',' ' /
      DATA PLUS,CAPI / '+','I' /
      MP= MIN(M,130)
      EMAX = 0
      DO 10 I=1,N
      DO 10 J=1,MP
10    EMAX = MAX(EMAX, ABS(AIMAG(A(J,I))) )
        RIMAX= MAXVAL(abs(IMAG(A(:,:))))

C      DMAX (ON INPUT) IS A SMALL FACTOR OF SCALE OF REAL PART
      IF(EMAX.LE.DMAX) GOTO 100
      DMAX = EMAX
      SCALE = 10/EMAX
C20   FORMAT(// 20('*') //)
      IF(ABS(LOG10(EMAX)).LE.4) THEN
      WRITE(KO,12) EMAX !,RIMAX
       ELSE
      WRITE(KO,13) EMAX !,RIMAX
       ENDIF
12    FORMAT(' Imaginary part: full scale =',2F12.5/ )
13    FORMAT(' Imaginary part: full scale =',1P,2E12.4/ )
      DO 50 I=1,N
      DO 30 J=1,MP
         ANO= SCALE * ABS(AIMAG(A(J,I)))
         NO = ANO
         NO = max(NO,0)
         LINE(J) = CHARS(NO+1)
         IF(NO.LT.1 .AND. ANO .GE. 0.10) LINE(J) = PLUS
30    CONTINUE
      WRITE(KO,35) (LINE(K),K=1,MP),CAPI
35    FORMAT(' I', 131A1)
      P = .FALSE.
      DO 40 J=1,MP
         IS = SIGN(1D0,DBLE(AIMAG(A(J,I))))
         P = P .OR. IS.LT.0
         P = P .OR. IS.EQ.0
40       LINE(J) = SIGNS(IS+2)
!     IF(P) WRITE(KO,45)   (LINE(J),J=1,MP)
45    FORMAT('+ ',130A1)
50    CONTINUE
C     WRITE(KO,20)
      RETURN
100   WRITE(KO,13) EMAX
      RETURN
      END
****DISPX**************************************************************
      SUBROUTINE DISPX(SMAT,NCH,JIN,EL,JCOEF,XS,SIGJ,PART,
     X           EXCIT,JEX,ITC,K,RMASS,PEL,EXL,LVAL,FUSL,OUTJ,
     X           JVAL,JPROJ,JTARG,JTOTAL,PARITY,
     X           SMATS,CHANS,JEND,SCALE,ITCM,IF1,IF2)
	use io
	use parameters
	use searchpar, only: final
	use trace, only: smatl
	use parallel, only: mpisets
      IMPLICIT REAL*8(A-H,O-Z)
C
      INTEGER PEL,EXL,PARITY,SMATS,CHANS,C,ITC(MXP,MXX),EL,
     X        PART(NCH),EXCIT(NCH),LVAL(NCH)
      COMPLEX*16 SMAT(NCH)
      REAL*8 XS(MAXCH),SIGJ(0:MXPEX),K(MXP,MXX),SCSIG(0:MXPEX)
      REAL*8 JCOEF,FUSL(1+NFUS)
      REAL*8 RMASS(MXP),JEX(6,MXP,MXX),JTOTAL
      LOGICAL JEND
      REAL*8 JVAL(MAXCH),JPROJ(MAXCH),JTARG(MAXCH)
      CHARACTER*1 BLANK,SLASH,HATCH,PSIGN(3),SCALE(0:MXPEX),
     x		CHF(1:1+NFUS)
C
      DATA  SLASH,HATCH /'/','#'/
      DATA PSIGN / '-','?','+' /, BLANK / ' ' /, Z / 0D0 /
C
      CHF(1) = 'f'
	i0 = min(NFUS,24)
      do i=1,i0
      CHF(1+i) = char(ichar('c')+i-1)
      enddo
      do i=i0+1,NFUS
      CHF(1+i) = char(ichar('A')+i-i0-1)
      enddo
      FUSL(1) = 0.0
      SIGJ(0) = 0.0
      IFT = 1
      DO 690 C=1,NCH
         IC = PART(C)
         IA = EXCIT(C)
         IT = ITC(IC,IA)
        IF(DBLE(K(IC,IA)**2) .LT. 0.0) GOTO 690
      AMDSQS = ABS(SMAT(C))**2
      IF(C.EQ.EL) AMDSQS = 1 - AMDSQS
!@@
         IF(RMASS(IC).lt.1e-5) then
             S =          RMASS(PEL)*amu/ (HBC*K(PEL,EXL))
!@@@         S = K(IC,IA)*RMASS(PEL)*amu/ (HBC*K(PEL,EXL))
         ELSE IF(RMASS(PEL).lt.1e-5) then
             S = HBC*K(IC,IA)/(RMASS(IC)*amu)
         ELSE
             S = K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
         ENDIF
!	 S is v_out/v_in but for photons this is just
!	      c/v_in for emission so since
!	      v_in/c = h*K_in/(RMASS_in*c) 
!	             = (hc)*K_in/(RMASS_in*c^2)
!		     = (hc)*K_in/(RMASS_in*amu)
!	 S = RMASS(PEL)*amu/ (HBC*K(PEL,EXL))	or as before
!	 
!        S = K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
!@@
      XS(C) = JCOEF* AMDSQS * S
      IF(C.EQ.EL) THEN
         SIGJ(0) = SIGJ(0) + XS(C)
         FUSL(1) = FUSL(1) + XS(C)
        ELSE
         SIGJ(IT) = SIGJ(IT) + XS(C)
         FUSL(1) = FUSL(1) - XS(C)
        ENDIF
      IFT = IF2
  690 CONTINUE
C
	if(SMATL>=2) 
     X WRITE(38,1445) JTOTAL,PSIGN(PARITY+2),LVAL(EL),JIN
     X       , (SIGJ(IT),IT=0,ITCM),(FUSL(I),I=IF1,IFT)
 1445 FORMAT(F7.1,A1,I6,I2,10G12.4,:,/(16X,10G12.4))
	OUTJ = SIGJ(0)
	if(SMATL>=2) then
       	written(38) = .true.
       	call flush(38)
	endif
	SCSIG(:) = SIGJ(:)
      IF(final.and.(max(0,CHANS)+SMATS.GE.1.OR.JEND)) THEN
      DO 710 IT=0,ITCM
      SCALE(IT) = BLANK
      IF(ABS(SCSIG(IT)).GT.0.010) GO TO 710
         SCSIG(IT) = SCSIG(IT) * 1000.0
         SCALE(IT)= SLASH
      IF(ABS(SCSIG(IT)).GT.0.010) GO TO 710
         SCSIG(IT) = SCSIG(IT) * 1000.0
         SCALE(IT)= HATCH
  710 CONTINUE
       WRITE(KO,1450) JTOTAL,PSIGN(PARITY+2),MOD(LVAL(EL),100),JIN
     X    , (SCSIG(IT),SCALE(IT),IT=0,ITCM),(FUSL(I),CHF(I),I=IF1,IFT)
	if(mpisets.gt.1) 
     X WRITE(48,1450) JTOTAL,PSIGN(PARITY+2),MOD(LVAL(EL),100),JIN
     X    , (SCSIG(IT),SCALE(IT),IT=0,ITCM),(FUSL(I),CHF(I),I=IF1,IFT)
       ENDIF
 1450 FORMAT(' Reaction Xsec',F8.1,A1,'/',I2,' @',I2,' =',F10.3,A1,',',
     X   ' Out:', 9(F8.3,A1),:,/(11X,'Xsec',23X,10(F8.3,A1)))
      if(SMATL>1) write(KO,*)
      CALL FLUSH(KO)
!     DO 720 C=1,NCH
!       IF(ABS(SMAT(C)).GT.0E-10)
!    X WRITE(77,750) SMAT(C),LVAL(C),JVAL(C),JTOTAL,PART(C),EXCIT(C),
!    X    LVAL(EL),JVAL(EL)
!    *  ,jproj(c),jtarg(c),jproj(el),jtarg(el)
!call flush(77)
!720    CONTINUE
!750  FORMAT(2F15.10,I6,2F6.1,2I3,I6,F6.1,4F4.1)
C
      RETURN
      END
      function lnbl(s)
      character*70 s
      l=len(s)
      do 1 i=l,1,-1
      if(s(i:i).ne.' ') go to 5
1     continue
      i=0
      i=1
5     lnbl=i
      end
      SUBROUTINE FILEMV(KIN,KOUT)
      CHARACTER*200 ALINE
	write(51,*) ' Append all file ',KIN,' to ',KOUT
      WRITE(KIN,713)
  713    FORMAT('EOF')
      REWIND KIN
  777 READ(KIN,'(a)',END=778) ALINE
      IF(ALINE(1:3).eq.'EOF') GO TO 778
      do 1 l=200,1,-1
      if(ALINE(l:l).ne.' ') go to 5
1     continue
      l=0
5     WRITE(KOUT,'(a)') ALINE(1:l)
      GO TO 777
  778 CALL FLUSH(KOUT)
      REWIND KIN
      end
	subroutine rewop(ifile)
	integer ifile
	logical op
	character*8 name
	inquire(ifile,opened=op)
	if(.not.op) then
	   if(ifile<10) then
		write(name,'(''fort.'',i1)') ifile
	   else if(ifile<100) then
		write(name,'(''fort.'',i2)') ifile
	   else if(ifile<1000) then
		write(name,'(''fort.'',i3)') ifile
	   endif
	   open(ifile,form='formatted',file=name)
	endif
	rewind ifile
	return
	end


	subroutine openif(ifile)
	integer ifile
	logical op
	character*8 name
	inquire(ifile,opened=op)
	if(.not.op) then
	   if(ifile<10) then
		write(name,'(''fort.'',i1)') ifile
	   else if(ifile<100) then
		write(name,'(''fort.'',i2)') ifile
	   else if(ifile<1000) then
		write(name,'(''fort.'',i3)') ifile
	   endif
	   open(ifile,form='formatted',file=name)
	endif
!!!!!!	rewind ifile
	return
	end

	subroutine openuf(ifile)
	integer ifile
	logical op
	character*8 name
	inquire(ifile,opened=op)
	if(.not.op) then
	   if(ifile<10) then
		write(name,'(''fort.'',i1)') ifile
	   else if(ifile<100) then
		write(name,'(''fort.'',i2)') ifile
	   else if(ifile<1000) then
		write(name,'(''fort.'',i3)') ifile
	   endif
	   open(ifile,form='unformatted',file=name)
	endif
!!!!!!	rewind ifile
	return
	end

