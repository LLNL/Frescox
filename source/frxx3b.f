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
*****FRXX3B***************************************************************
      SUBROUTINE ERWIN(W,ECM,COEF,IEX,FORMF,NF,FORMFR,INHOMG,CH,
     $  EL,SMAT,L,JVAL,CORESP,LL1,NEQS,NICH,N,H,M,MD,SCL,RENORM,
     $  REPEAT,WREQ,BLOCKD,AL,RE,IC,SHOW,SHSMAT,SMATEL, 
     $  CUTOFF,SKIP,FED,PTYPE,F,FIM,FIMD,LOCFIL,
     $  CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)
	use io
	use factorials
	use parameters
	use drier
	use searchpar, only: final
      IMPLICIT REAL*8(A-H,O-Z)
C
C        SOLVE 'NEQ' COUPLED SCHROEDINGERS EQUATIONS
C        SOLVE   (FOR NOT BLOCKD(K))
C             (COEF(K). D2/DR2 + EN(R,K)).W(R,K)
C                  + SUM(J): COUPL(R,K,J).W(R,J) + INHOMG(R,K) = 0
C           WHERE COUPL(R,K,J) = SUM(NC): FORMF(R,JF)*CLIST(K,J,NC)
C           WHERE    JF = NFLIST(K,J,NC)
C         FOR K & J <= IEX (ALL OTHER COUPLINGS MUST BE ITERATED)
C       &   WHERE EN(R,K) IS THE DIAGONAL PART OF COUPL
C
C     ASSUMED FOR ARRAYS HERE THAT N<=MAXN AND NEQS <= MAXCH
C
C  Using 'Enhanced Numerov' of Thorlacius & Cooper (JCP 72(1987) 70)
C   with 5 terms in cosh(sqrt(12T)) expansion, but only diagonal potl.
C
      PARAMETER(MINVEC=5)
C               MINIMUM VECTOR LENGTH FOR EFFICIENT OPERATION.
      SAVE IS
      COMPLEX*16 F(NFDEC),WVD(MAXN),FIMD(MAXB,MAXCH),FIM(MAXB,MAXCH)
      COMPLEX*16 MAT(2*NEQS+1,2*NEQS+2),S,ZIV,ZPV,FI(MAXB,MAXCH),
     &           V(MAXB,MAXCH),ZI(MAXB,MAXCH),ZM(MAXB,MAXCH),CZ,
     &           SE,CH(2,MAXCH),SMAT(MAXCH),ONEC,ZI2,CI
      COMPLEX*16 FORMF(MAXM,NF),INHOMG(MAXN,MAXCH),F8,FORMFR(NF),
     &       C,W(MAXN,NICH),CLIST(MAXCH,MAXCH,MCLIST),COUPL(MAXB,MAXB)
      COMPLEX*16 ZDL,ZHPLUS,ZHMINUS
      REAL*8 COEF(MAXCH),JVAL(MAXCH),CORESP(MAXCH),H2C(MAXCH),MAGN(MAXN)
     &          ,ECM(MAXCH),H(MAXCH),LL1(MAXCH),
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),CRCRAT(MAXCH)
      INTEGER EL,L(MAXCH),FIT,CUTOFF,SHOW,SKIPL,SKIPU,FEDL,FEDU,
     X        CUTVAL(MAXCH),NFLIST(MAXCH,MAXCH,MCLIST),
     X        NCLIST(MAXCH,MAXCH),HOLES,SKFED,PTYPE(12,NF)
      LOGICAL WREQ,SING,REPEAT,SHSMAT,SMATEL,BLOCKD(MAXCH),CRAY,ITV,
     X           SKIP(MAXCH),FED(MAXCH),DO1,DO2,DO3,BFEED,BSKIP,LOCFIL,
     X           FCWFN,FJSWTCH
      AMD1(C) = ABS(DBLE(C)) + ABS(AIMAG(C))
C
      CRAY = MACH .EQ. 9 
      NM1 = N-1
      NR = 1 + 2*NEQS
      NP = NR + 1
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
      ZI2 = (0.0D0,0.5D0)
      CI = (0.0D0,1.0d0)
      CZ = (0D0,0D0)
      ONEC = (1D0,0D0)
C     AL = 0.0
      ALUN = 0.0
      NQ1 = 1+NEQS
      IEX1= IEX + 1
      IIX = IEX
      NITS = 2+IEX
      FIT = 2
      NLIF = NITS*NICH

       IF(REPEAT) THEN
         NITS = 1
         FIT = 1
       ENDIF
      IF(.NOT.REPEAT .AND. IEX.GE.NEQS) FIT = 3

      DO1 = REPEAT
      DO2 = .NOT.REPEAT .AND. IEX.LT.NEQS
      DO3 = .NOT.REPEAT .AND. IEX.GT.0
      IF(.NOT.REPEAT) IS =  CUTOFF
      IF(.NOT.REPEAT) IS = MIN(1*N/2, MAX(2, IS ))
		do I=1,NEQS
		CUTVAL(I) = max(CUTVAL(I),IS)
		enddo
       if(SHOW.ge.3) 
     *       write(6,*) 'ERWIN: CUTOFF,IS =',CUTOFF,IS
       if(SHOW.ge.3) 
     *       write(6,*) 'ERWIN: CUTVAL =',(CUTVAL(I),I=1,NEQS)
      NN = NM1 - IS + 1
      TMAX = 20.
      TMIN = -125.
        SMALL = 1.0/FPMAX
        EPS = SQRT(SMALL)
        BIG = 1./EPS
      NVREQ =  NN*NLIF
!      write(48,*) 'NICH,IEX,NITS,NLIF,NN =',NICH,IEX,NITS,NLIF,NN
!#      write(48,*)  'NN,NLIF,NVREQ,NFDEC',NN,NLIF,NVREQ,NFDEC
!      write(48,*) 'ERWIN solutions, with IEX=',IEX,NEQS
      call flush(48)
      CALL CHECK(NVREQ,NFDEC,17)
      if(SHOW>1) then
	do K=1,min(3,NEQS)
	write(KO,*) 'Ch ',K,' ECM/COEF=',real(ECM(K)/COEF(K))
	enddo
      endif

C
      IF (FJSWTCH)  THEN
c                              skip nmrov , match CRCWFN to zero
761      MAT(:,:) = 0.0
         MAT(1,1) = 1
         DO 780 IT=1,NEQS
               MAT(IT+NQ1,IT+1) = ONEC
               MAT(IT+1,IT+1)   = ONEC*CRCRAT(it)
            MAT(1 ,IT+1) = 1
            MAT(1 ,IT+NQ1) = 0
780       CONTINUE
C
      ELSE
C
      IF(SHOW.GE.2) WRITE(KO, 5) IEX,NEQS,NF,M,REPEAT,NICH,NQ1,NR,
     X                          NVREQ
      IF(SHOW.GE.2) WRITE(KO, 4) (H(K),K=1,NEQS)
4     FORMAT(' ERWIN step sizes are',10F8.5)
5     FORMAT(' ERWIN given',4I4,L4,' & requires',3I6,' i.e.',I8,'/',I8)
      CALL CHECK(NITS,MAXB,15)
C
C   If SKIP(K), don't integrate channel K, as same as last iteration!
C   If not FED(K) too, set channel to zero to start with.
C   NB. Some channels are non-zero if FED(K) & SKIP(K).
C
C    KFIRST - KLAST (incl) = range of channels integrated
C    FEDl   - FEDL  (incl) = range of channels non-zero
C     (may be larger than KFIRST->KLAST, if some SKIPed)
C
       BFEED = .FALSE.
       BSKIP = .TRUE.
       DO 5001 K=1,IEX
        BSKIP = BSKIP .AND. SKIP(K)
5001    IF(FED(K)) BFEED = .TRUE.
C     in a closely-coupled set,
C     all FED if any one is,
C and only allow any SKIP if all are.
C
       DO 5002 K=1,IEX
       IF(BFEED) FED(K) = .TRUE.
5002   IF(.NOT.BSKIP) SKIP(K) = .FALSE.
      FEDU = 0
      FEDL = NEQS+1
      SKIPU = 0
      SKIPL = NEQS+1
          DO 7 K=1,NEQS
            SMAT(K) = CZ
             IF(.NOT.REPEAT) THEN
               SKIP(K) = .FALSE.
               FED(K)  = .FALSE.
             ELSE
              FED(K) = FED(K) .AND..NOT. BLOCKD(K)
              SKIP(K) = SKIP(K) .OR. .NOT.FED(K)
               IF(FED(K)) FEDL = MIN(FEDL,K)
               IF(FED(K)) FEDU = MAX(FEDU,K)
               IF(.NOT.SKIP(K)) SKIPL = MIN(SKIPL,K)
               IF(.NOT.SKIP(K)) SKIPU = MAX(SKIPU,K)
             ENDIF
          H2C(K) = H(K)**2 / COEF(K)
C7        LL1(K) = -L(K)*(L(K)+1)/H2C(K)
C7        LL1(K) = -LL1(K)/H2C(K)
7        CONTINUE
         HOLES = 0
         SKFED = 0
      IF(REPEAT) THEN
          KFIRST = MAX(FEDL,SKIPL)
          KLAST  = MIN(FEDU,SKIPU)
            DO 701 K=KFIRST+1,KLAST-1
             IF(SKIP(K).AND.FED(K))  SKFED = SKFED + 1
701          IF(SKIP(K))  HOLES = HOLES + 1
            DO 702 K=FEDL,FEDU
702          IF(SKIP(K).AND.FED(K))  SKFED = SKFED + 1
       ELSE
           FEDL   = 1
           FEDU   = NEQS
           KFIRST = 1
           KLAST  = NEQS
         ENDIF
        ITV = IIX-KFIRST+1 .GE. MINVEC
C
         IF(SHOW.GE.2.OR.KFIRST.GT.1.AND.KFIRST.LE.IEX) THEN
          WRITE(KO,*) 'FEDL,FEDU,SKIPL,SKIPU,KFIRST,KLAST,HOLES,SKFED ='
     X            ,FEDL,FEDU,SKIPL,SKIPU,KFIRST,KLAST,HOLES,SKFED
            WRITE(KO,*) 'FED(*) =',(FED(K),K=1,NEQS)
            WRITE(KO,*) 'SKIP(*) =',(SKIP(K),K=1,NEQS)
            ENDIF
         DO 705 I=1,N
705           MAGN(I) = 1.0
C
      DO 11 IT=1,NITS
           KMAX = NEQS
           IF(IT.GE.3) KMAX = IEX
          KIMAX = MIN(KMAX,NICH)
           L1 = (IT-1) * NICH
!           IF(IT.GE.3) L1 = (IT-3)*NICH + 2*NICH
          DO 10 K=1,NEQS
            IF(WREQ.AND.K.LE.KIMAX.AND..NOT.(SKIP(K).AND.FED(K))) THEN
C                            THIS VECTORISES OK.
            II = L1+K - IS*NLIF
             DO 8 I=IS,NM1
8           F(II+I*NLIF) = CZ
            ENDIF
         V(IT,K)   = 0
         FI(IT,K) = 0
         ZM(IT,K)  = 0
         ZI(IT,K)  = 0
          IF(REPEAT) GO TO 10
            FIM(IT,K) = 0
            FIMD(IT,K)= 0
10       CONTINUE
11     CONTINUE

C
        DO 141 K=1,FEDU
141      IF(SHOW.GE.6) WRITE(KO,142) K,H(K),COEF(K),(CH(II,K),II=1,2)
142      FORMAT(' FOR CH.',I4,' H,COEF,CH(1,2)=',10F10.5)
C-----------------------------------------------------------------------
          IF(REPEAT .AND. KLAST.LT.KFIRST) GO TO 61
	  if(LOCFIL) then
	  rewind 19
	  do I=1,IS-1
	   read(19) 
	  enddo
	  endif
         DO 60 I=IS,NM1
          RI2 = 1D0 / DBLE(I-1)**2
	  if(LOCFIL) read(19) FORMFR
          DO 148 K=KFIRST,KLAST
           IF(SKIP(K)) GO TO 148
            C = CZ
	   if(LOCFIL) then
CDIR$                NOVECTOR
            DO NC=1,NCLIST(K,K)
		JF = NFLIST(K,K,NC)
            IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,K,NC) * FORMFR(JF)
	    enddo
	   else
CDIR$                NOVECTOR
            DO NC=1,NCLIST(K,K)
		JF = NFLIST(K,K,NC)
            IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,K,NC) * FORMF(I,JF)
	    enddo
	   endif
CDIR$ VECTOR
C          SMAT(K) = (LL1(K)*RI2 - ECM(K) + C) * H2C(K)
           SMAT(K) = -LL1(K)*RI2 + (-ECM(K) + C) * H2C(K)
            F8 = SMAT(K)
            IF(ABS(AIMAG(F8)).LT.1E-20*ABS(DBLE(F8)))
     X         SMAT(K) = DBLE(SMAT(K))
148        CONTINUE
C
       if(DO2) then
         DO 12 K=IEX1,NEQS
           IF(SKIP(K)) GO TO 12
           IF(I.ne.CUTVAL(K)) GO TO 12
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
               ZI(2,K) = SCL * EXP(T)
       	IF(SHOW.GE.3) write(KO,1302) I,2,K,ZI(2,K)
12       CONTINUE
       endif
       if(DO3) then
          DO 13 IT=3,NITS
            K = IT-2
              IF(SKIP(K)) GO TO 13
              IF(I.ne.CUTVAL(K)) GO TO 13
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
              ZI(IT,K) = SCL * EXP(T)
       	IF(SHOW.GE.3) write(KO,1302) I,IT,K,ZI(IT,K)
13         CONTINUE
       endif
1302         FORMAT(' AT I=',I3,' ZI(',I5,',',I5,')=',1P,2E12.2)

         DO 16 IT=FIT,NITS
            IF(SKFED.GT.0) THEN
C       Approx. ZI(IT,J), if J is FED and SKIPed,
C         by F(L1+J,I) already stored, for use in DO 24 loop.
              KMAX = FEDU
              IF(IT.GE.3) KMAX = MIN(IEX,FEDU)
              KIMAX = MIN(KMAX,NICH)
               L1 = (IT-1) * NICH
!                IF(IT.GE.3) L1 = (IT-3)*NICH + 2*NICH
                DO 157 K=FEDL,KIMAX
 157             IF(SKIP(K).AND.FED(K)) ZI(IT,K) = F(L1+K+(I-IS)*NLIF)
             ENDIF
           KMAX = KLAST
           IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
        IF(KMAX-KFIRST.LT.MINVEC-1) THEN
CDIR$ NOVECTOR
         DO 158 K=KFIRST,KMAX
158       IF(.NOT.SKIP(K)) FI(IT,K) = ZI(IT,K) *
C    X      (ONEC - SMAT(K) * R12 + ENA2*SMAT(K)**2 + ENA3*SMAT(K)**3)
     X      (ONEC - SMAT(K) * (R12 - SMAT(K)*(ENA2 + SMAT(K)*ENA3)))
CDIR$ VECTOR
        ELSE
         DO 159 K=KFIRST,KMAX
159       IF(.NOT.SKIP(K)) FI(IT,K) = ZI(IT,K) *
     X      (ONEC - SMAT(K) * (R12 - SMAT(K)*(ENA2 + SMAT(K)*ENA3)))
        ENDIF
16      CONTINUE
         IF(SHOW.GE.7) THEN
          DO 161 K=KFIRST,KLAST
161          WRITE(KO,165) I,K,H(K)*(I-1),L(K),
     &    SMAT(K)/H2C(K)+ECM(K),FI(FIT,K)
165      FORMAT(' AT I,K,R =',I6,I4,F7.3,
     &          ' L,PE,PSI =',I4,2F9.3,1P,2E10.1)
          ENDIF
C
         IF(DRY) GO TO 40
         IF(KFIRST.GT.IEX) GO TO 28
         DO 25 K=KFIRST,IEX
            DO 17 J=FEDL,IIX
17          COUPL(K,J) = 0.0
         IF(SKIP(K)) GO TO 25
            H2 = H2C(K) * R12
           IF(I.LT.IC) GO TO 21
           DO 19 J=FEDL,IIX
            C = CZ
	   if(LOCFIL) then
CDIR$                NOVECTOR
            DO 18 NC=1,NCLIST(K,J)
		JF = NFLIST(K,J,NC)
18          IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,J,NC) * FORMFR(JF)
	   else
CDIR$                NOVECTOR
            DO 181 NC=1,NCLIST(K,J)
		JF = NFLIST(K,J,NC)
181         IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,J,NC) * FORMF(I,JF)
	   endif
CDIR$ VECTOR
19         COUPL(K,J) = C
C (THE DIAGONAL PART COUPL(K,K) IS NOT USED, AS CCU(K,K)=0 ALWAYS
C
21       DO 24 J=FEDL,IIX
            IF(K==J.or.NCLIST(K,J)==0) GO TO 24
            C = COUPL(K,J) * H2
c    if(i.gt.200.and.i.lt.220) write(6,*) i,k,j,c
          IF(DO2) FI(2,K) = FI(2,K) - C * ZI(2,J)
          IF(DO1) FI(1,K) = FI(1,K) - C * ZI(1,J)
            IF(DO3) THEN
            IF(.NOT.ITV) THEN
CDIR$                      NOVECTOR
             DO 22 IT=3,NITS
 22          FI(IT,K) = FI(IT,K) - C * ZI(IT,J)
CDIR$ VECTOR
            ELSE
             DO 23 IT=3,NITS
 23          FI(IT,K) = FI(IT,K) - C * ZI(IT,J)
            ENDIF
            ENDIF
24       CONTINUE
25    CONTINUE
C
28     IF(.NOT.REPEAT) THEN
         DO 29 IT=FIT,NITS
         DO 29 K=KFIRST,KLAST
29       V(IT,K) = 0.0
       ELSE
C           REPEAT: IT = 1 ONLY
         IF(KLAST-KFIRST.LT.MINVEC-1) THEN
CDIR$             NOVECTOR
         DO 314 K=KFIRST,KLAST
          IF(.NOT.SKIP(K)) V(1,K) = INHOMG(I,K) * H2C(K)
314      FI(1,K) = FI(1,K) - V(1,K) * R12
CDIR$ VECTOR
          ELSE
         DO 315 K=KFIRST,KLAST
C                                 THIS STILL VECTORISES!
          IF(.NOT.SKIP(K)) V(1,K) = INHOMG(I,K) * H2C(K)
315      FI(1,K) = FI(1,K) - V(1,K) * R12
         ENDIF
        ENDIF
         IF(KFIRST.GT.IEX) GO TO 40
         DO 35 K=KFIRST,IEX
         IF(SKIP(K)) GO TO 35
         DO 34 J=FEDL,IIX
           IF(K==J.or.NCLIST(K,J)==0) GO TO 34
            C = COUPL(K,J) * H2C(K)
          IF(DO2) V(2,K) = V(2,K) + C * FI(2,J)
          IF(DO1) V(1,K) = V(1,K) + C * FI(1,J)
          IF(DO3) THEN
            IF(.NOT.ITV) THEN
CDIR$                      NOVECTOR
            DO 33 IT=3,NITS
33          V(IT,K) = V(IT,K) + C * FI(IT,J)
CDIR$ VECTOR
            ELSE
            DO 331 IT=3,NITS
331          V(IT,K) = V(IT,K) + C * FI(IT,J)
            ENDIF
           ENDIF
34       CONTINUE
35       CONTINUE
40       DO 46 IT=FIT,NITS
           KMAX = KLAST
           IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
         K = KMAX-KFIRST+1
         IF(K.LT.MINVEC.OR..NOT.CRAY.OR.HOLES*2.GT.K) THEN
CDIR$              NOVECTOR
         DO 43 K=KFIRST,KMAX
            IF(SKIP(K)) GO TO 43
            ZIV = ZI(IT,K)
            ZPV = ZIV + ZIV - ZM(IT,K) - V(IT,K) - SMAT(K) * FI(IT,K)
            ZM(IT,K) = ZIV
            ZI(IT,K) = ZPV
43            CONTINUE
CDIR$ VECTOR
           ELSE
         DO 44 K=KFIRST,KMAX
C           IF(SKIP(K)) GO TO 44
C                   OMIT SKIP TO GET VECTORISING; RESET IN DO 45 LOOP.
            ZIV = ZI(IT,K)
            ZPV = ZIV + ZIV - ZM(IT,K) - V(IT,K) - SMAT(K) * FI(IT,K)
            ZM(IT,K) = ZIV
            ZI(IT,K) = ZPV
44           CONTINUE
                  IF(HOLES.GT.0) THEN
                    DO 45 K=KFIRST,KMAX
                     IF(SKIP(K)) ZI(IT,K) = CZ
45                   IF(SKIP(K).AND.FED(K)) ZM(IT,K) = CZ
                   ENDIF
       ENDIF
46    CONTINUE
      KMAX = KLAST
      DO 468 IT=FIT,NITS
       IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
       IF(I.EQ.M) THEN
        DO 463 K=KFIRST,KMAX
463     IF(.NOT.SKIP(K)) FIM(IT,K) = FI(IT,K)
       ELSE IF(I.EQ.MD) THEN
        DO 465 K=KFIRST,KMAX
465     IF(.NOT.SKIP(K)) FIMD(IT,K) = FI(IT,K)
       ENDIF
      IF(.NOT.WREQ) GO TO 468
       KIMAX = MIN(KMAX,NICH)
        L1 = (IT-1) * NICH
!        IF(IT.GE.3) L1 = (IT-3)*NICH + 2*NICH
        DO 466 K=KFIRST,KIMAX
466     IF(.NOT.SKIP(K)) F(L1+K + (I-IS)*NLIF) = FI(IT,K)
468   CONTINUE
            T = 0.0
         IF(REPEAT)  GO TO 60
        IF(CRAY) THEN
         K = ICAMAX(MAXB*MAXCH,FI,1)
         F8 = FI(K,1)
         T = AMD1(F8)
         IF(T.GT.BIG) GO TO 48
        ELSE
         DO 47 IT=FIT,NITS
           KMAX = KLAST
           IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
            DO 47 K=KFIRST,KMAX
            F8 = FI(IT,K)
            T = MAX(T, AMD1(F8))
   47       IF(T .GT. BIG ) GO TO 48
         ENDIF
            MAGN(I) = T
         IF(I.NE.NM1) GO TO 60
 48      IF(T.EQ.0.0) GO TO 60
            MAGN(I) = T
         T = RENORM/T
         IF(SHOW.GE.1) WRITE(KO,49) T,I,RENORM
49       FORMAT(' Renormalising by',E12.4,' at',I4,' to max. WF.',E12.4)
         if(t.eq.0d0) then
          write(KO,*) 'RENORMALISING TO 0! at I=',I,' as lWF=',MAGN(I)
          write(KO,*) ' FI:',FI
          T = SMALL
          endif
         FMIN = small/T
         DO 55 IT=FIT,NITS
           KMAX = KLAST
           IF(IT.GE.3) KMAX = MIN(IEX,KLAST)
          L1 = (IT-1) * NICH
!           IF(IT.GE.3) L1 = (IT-3)*NICH + 2*NICH
         DO 51 K=KFIRST,KMAX
            IF(BLOCKD(K)) GO TO 51
            ZI(IT,K) = ZI(IT,K) * T
            ZM(IT,K) = ZM(IT,K) * T
            FIM(IT,K) = FIM(IT,K) * T
            FIMD(IT,K) = FIMD(IT,K) * T
           IF(.NOT.WREQ.OR.K.GT.NICH) GO TO 51
           IF(CRAY) THEN
            II = L1+K - IS*NLIF
            DO 499 J=IS,I
499           F(II+J*NLIF) = F(II+J*NLIF) * T
           ELSE
CDIR$             NOVECTOR
            DO 50 J=IS,I
             JJ = L1+K+(J-IS)*NLIF
             F8 = F(JJ)
50          IF(AMD1(F8).GT.FMIN) F(JJ) = F(JJ) * T
CDIR$ VECTOR
           ENDIF
51       CONTINUE
55       CONTINUE
C
60       CONTINUE
C-----------------------------------------------------------------------
c
61    DO 65 I=1,NR
      DO 65 J=1,NP
65     MAT(I,J) = 0
         DO 80 IT=1,NQ1
         IF(IT.EQ.1 .AND..NOT.REPEAT) GO TO 75
         J = IT-1
         IF(J.GT.IEX) I = 2
         IF(J.LE.IEX) I = 2 + J
         IF(IT.EQ.1) I = 1
           KMAX = NEQS
           IF(I.GE.3) KMAX = IEX
         DO 70 K=1,KMAX
            IF(BLOCKD(K)) GO TO 70
            IF(IT.EQ.1.AND..NOT.FED(K)) GO TO 70
         IF(I.EQ.2 .AND. K.NE.IT-1) GO TO 70
            MAT(K+NQ1,IT) = FIMD(I,K)
            MAT(K+1,IT)   = FIM(I,K)
70          CONTINUE
75       MAT(1 ,IT) = 1
         IF(IT .GT. 1 .AND. REPEAT) MAT(1,IT) = 0
C           IF(IT.GT.K) GO TO 80
         MAT(1 ,IT+NQ1) = 0
80       CONTINUE
C
c
      ENDIF
c
C                    MATCH TO REQUIRED ASYMPTOTIC CONDITIONS
C
      IF (FCWFN) THEN
C                           match Numerov to CRCWFN
CDIR$   NOVECTOR
       DO 82 K=1,NEQS
       DO 82 K2=1,NEQS
           IF(SHOW.ge.4.and.ABS(CGMAT(K2,K,2)).gt.0.0) 
     *       WRITE(6,716) K,K2,CGMAT(K2,K,2),CFMAT(K2,K,2),
     *           CGMAT(K2,K,1),CFMAT(K2,K,1)
716            FORMAT(1X,2I3,': CRCWFN= ',4(D15.7,3X))
c
c james
         ZDL = CI**(L(K)-L(K2))
         ZHPLUS = CGMAT(K2,K,2)+CI*CFMAT(K2,K,2)
         MAT(K2+1,K+NQ1)=   ZI2*ZDL*ZHPLUS
         ZHPLUS = CGMAT(K2,K,1)+CI*CFMAT(K2,K,1)
         MAT(K2+NQ1,K+NQ1)= ZI2*ZDL*ZHPLUS
c
c        MAT(K2+1,K+NQ1)=   ZI2*(CGMAT(K2,K,2)+CI*CFMAT(K2,K,2))
c        MAT(K2+NQ1,K+NQ1)= ZI2*(CGMAT(K2,K,1)+CI*CFMAT(K2,K,1))
82       CONTINUE
CDIR$ VECTOR
         MAT(1,NP)=1
         DO 85 K=1,NEQS
c james
         ZDL = CI**(L(EL)-L(K))
         ZHMINUS = CGMAT(K,EL,2)-CI*CFMAT(K,EL,2)
         MAT(K+1,NP)=   ZI2*ZDL*ZHMINUS
         ZHMINUS = CGMAT(K,EL,1)-CI*CFMAT(K,EL,1)
         MAT(K+NQ1,NP)= ZI2*ZDL*ZHMINUS
c
c          MAT(K+1,NP)=   ZI2*(CGMAT(K,EL,2)-CI*CFMAT(K,EL,2))
c          MAT(K+NQ1,NP)= ZI2*(CGMAT(K,EL,1)-CI*CFMAT(K,EL,1))
85       CONTINUE
      ELSE
c                         match Numerov to Uncoupled Coulomb wfns
       DO 90  K=1,NEQS
         MAT(K+1,K+NQ1) = CH(1,K)
         MAT(K+NQ1 ,K+NQ1) = CH(2,K)
C                                      C(K) == MAT(K,NP) IS THE RHS.
         MAT(K+1,NP) = 0
         MAT(K+NQ1 ,NP) = 0
90       CONTINUE
         MAT(1 ,NP) = 1
c                                         note minus as ch = i/2*H+
         MAT(EL+1,NP) = - CONJG(CH(1,EL))
         MAT(EL+NQ1 ,NP) = - CONJG(CH(2,EL))
      ENDIF
         DO 100 I=1,NR
100        IF (SHOW.GE.5) WRITE(KO,105) I,(MAT(I,J),J=1,NP)
105      FORMAT(' MAT(',I2,',*) =', /(6(1P,E10.2,E9.1,1X)))
C
      CALL GAUSS(NR,NR,MAT,SING,S,SMALL,.FALSE.)
         IF(SING) GO TO 600
C
      T = EXP(DBLE(S) * 0.5/NEQS)
      IF(.NOT.REPEAT) RENORM = RENORM / T
       IF(RENORM.GT.BIG) RENORM = BIG
       IF(RENORM.LT.EPS) RENORM = EPS
      IF(SHOW.GE.3) WRITE(KO,130) (MAT(I,NP),I=1,NR)
130   FORMAT(/' The combining coefficients are',1P,/(2(E20.8,E15.8)))
C
c     ENDIF
c                  end of FJSWITCH loop matching CRCWFN to zero
      DO 150 K=1,NEQS
         DO 138 I=1,N
138         WVD(I) = CZ
      SMAT(K) = MAT(K+NQ1,NP)
C
      IF(ABS(SMAT(K)).GT.1E-10) THEN
        T = 0.0
       DO 145 IT=1,NQ1
         J = IT-1
         IF(J.GT.IEX) I = 2
         IF(J.LE.IEX) I = 2 + J
         IF(IT.EQ.1) I = 1
           KMAX = NEQS
           IF(I.GE.3) KMAX = IEX
         KIMAX = MIN(KMAX,NICH)
          L1 = (I-1) * NICH
!           IF(I.GE.3) L1 = (I-3)*NICH + 2*NICH
         IF(I.EQ.2 .AND. K.NE.IT-1 .OR. MAT(IT,NP).EQ.CZ
     &                 .OR. K.GT.KMAX) GO TO 145
           IF(WREQ.AND.K.LE.KIMAX) THEN
                II = L1+K - IS*NLIF
         DO 140 J=IS,NM1
140         WVD(J) = WVD(J) + MAT(IT,NP)*F(II+J*NLIF)
          ELSE
            WVD(MD) = WVD(MD) + MAT(IT,NP)*FIMD(I,K)
          ENDIF
            T = MAX(T,ABS(WVD(MD)))
145      CONTINUE
             H2C(K) = T / (ABS(WVD(MD))+1D-20)
	   if(ECM(K)>0.0) then
             AL = MAX(AL, H2C(K))
             ALUN = MAX(ALUN, H2C(K)*ABS(SMAT(K)))
	   endif
!	     write(106,*) ALUN,K,WVD(MD),H2C(k),smat(K),AL
      ENDIF
         IF(.NOT.WREQ.OR.K.GT.NICH)  GO TO 150
         DO 149 I=1,NM1
149         W(I,K) = WVD(I)
150      CONTINUE
      SING = SHOW.GE.1.AND..NOT.REPEAT
      IF(SING) WRITE(KO,1157) S,AL
1157  FORMAT(' ERWIN: LOG10(DETERMINANT) =',2F11.3,',  ACCURACY LOSS ='
     X  , 1P,E10.2)
C                       NOW CHECK THE LOGARITHMIC DERIVS. MATCH REQD.
      DO 200 K=1,NEQS
        IF(SHOW.LT.3.OR..NOT.WREQ.OR.K.GT.NICH) GO TO 185
         IF(ABS(W(M,K)) .EQ. 0) W(M,K) = (1.0E-20,0.0)
        SE = W(MD,K) / W(M,K)
        INIEX = 0
        IF(K.EQ.EL) INIEX = 1
        S  =            ( -INIEX*CONJG(CH(2,K)) - SMAT(K)*CH(2,K))
     &                 /( -INIEX*CONJG(CH(1,K)) - SMAT(K)*CH(1,K)+1D-20)
      C  = S - SE
      WRITE(KO,180) SE,C
180   FORMAT(/'  ',2F10.6,' matched within ',2E12.4)
c185   IF((SHSMAT.OR.K.EQ.EL.AND.SMATEL).AND.ABS(SMAT(K)).GT.1.D-6)
185   IF((SHSMAT.OR.K.EQ.EL.AND.SMATEL).and.final  )
     & WRITE(KO,190) K,SMAT(K),L(K),JVAL(K),CORESP(K),L(EL),
     &           LOG10(MAX(H2C(K),1D0))
190   FORMAT( '   S-matrix',I5,' = '
     & , 2F10.5,' for L=',I5,', J=',F7.1,' channel on core I =',F6.1,
     & ' from L=',I5,',  Acc. loss =',F5.1,' D.')
      IF(K.EQ.EL.AND.SMATEL.AND.ABS(SMAT(K)).GT.1.D-6.and.final) THEN
         SE = LOG(SMAT(EL))
         WRITE(KO,195) K,SE*(0.,-.5)*180./3.14159,L(K),JVAL(K)
195      FORMAT( ' Elastic phase shift ',I3,' = '
     &     , 2F8.3,' deg. for the L =',I5,', J =',F7.1,' channel.')
         WRITE(45,196) ECM(EL),SE*(0.,-.5)*180./3.14159,L(K),JVAL(K)
196        format(f10.3,2f9.3,' for LJin =',i6,f6.1)
	 written(45) = .true.
         ENDIF
200   CONTINUE
      if(.not.FJSWTCH) then
C     IF(AL.GT.1D10.OR.AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
      IF(AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
     X          3.*ALUN*RE * 100.
205   FORMAT('0****** WARNING : ACCURACY LOSS IN ERWIN =',F5.1,
     X ' DIGITS, SO EXPECT ERRORS OF',F9.4,' % OF UNITARITY'/)
      IF(AL*RE.GT..03) WRITE(KO,*) AL,ALUN,RE
	endif
      RETURN
600   IF(REPEAT) CALL ABEND(32)
C600   IF(REPEAT) STOP
CDIR$   NOVECTOR
      DO 605 I=IS,NM1
605   IF(MAGN(I).NE.0.0) MAGN(I) = LOG10(MAGN(I))
CDIR$ VECTOR
      WRITE(KO,610) (MAGN(I),I=IS,NM1)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
C     STOP
      CALL ABEND(4)
      END
      SUBROUTINE ERWINCC(ECM,COEF,FORMF,NF,FORMFR,CH,REPEAT,
     $  EL,SMAT,L,JVAL,CORESP,LL1,NEQS,N,H,M,MD,SCL,RENORM,
     $  AL,RE,IC,SHOW,SHSMAT,SMATEL,CUTOFF,PTYPE,FIM,FIMD,LOCFIL,
     $  CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)
	use io
	use factorials
	use parameters
	use drier
	use searchpar, only: final
      IMPLICIT REAL*8(A-H,O-Z)
C
C        SOLVE 'NEQ' COUPLED SCHROEDINGERS EQUATIONS BY EXACT CC
C        SOLVE  
C             (COEF(K). D2/DR2 + EN(R,K)).W(R,K)
C                  + SUM(J): COUPL(R,K,J).W(R,J) + INHOMG(R,K) = 0
C           WHERE COUPL(R,K,J) = SUM(NC): FORMF(R,JF)*CLIST(K,J,NC)
C           WHERE    JF = NFLIST(K,J,NC)
C       &   WHERE EN(R,K) IS THE DIAGONAL PART OF COUPL
C
C     ASSUMED FOR ARRAYS HERE THAT N<=MAXN AND NEQS <= MAXCH
C
C  Using 'Enhanced Numerov' of Thorlacius & Cooper (JCP 72(1987) 70)
C   with 5 terms in cosh(sqrt(12T)) expansion, but only diagonal potl.
C
      COMPLEX*16 FIMD(MAXB,MAXCH),FIM(MAXB,MAXCH)
      COMPLEX*16 MAT(2*NEQS,2*NEQS+1),S,ZIV,ZPV,FI(MAXB,MAXCH),
     &           V(MAXB,MAXCH),ZI(MAXB,MAXCH),ZM(MAXB,MAXCH),CZ,SE,
     &           CH(2,MAXCH),SMAT(MAXCH),ONEC,ZI2,CI,ZDL,ZHPLUS,ZHMINUS
      COMPLEX*16 FORMF(MAXM,NF),FORMFR(NF),F8,
     &       C,CLIST(MAXCH,MAXCH,MCLIST),COUPL(MAXB,MAXB)
      REAL*8 COEF(MAXCH),JVAL(MAXCH),CORESP(MAXCH),H2C(MAXCH),MAGN(MAXN)
     &          ,ECM(MAXCH),H(MAXCH),LL1(MAXCH),
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),CRCRAT(MAXCH)
      INTEGER EL,L(MAXCH),CUTOFF,SHOW,
     X        CUTVAL(MAXCH),NFLIST(MAXCH,MAXCH,MCLIST),
     X        NCLIST(MAXCH,MAXCH),PTYPE(12,NF)
      LOGICAL SING,SHSMAT,SMATEL,FCWFN,FJSWTCH,BLAS,REPEAT,LOCFIL
      AMD1(C) = ABS(DBLE(C)) + ABS(AIMAG(C))
C
!	BLAS = NEQS.ge.20 .and. (MACH==2 .or. MACH==6)
	BLAS = NEQS.ge.20 .and. (MACH==2 .or. MACH==3 .or. MACH==4)
!	BLAS = NEQS.ge.20 .and. MACH==2
	if(NEQS>0.and.(SHSMAT.and.SMATEL))
     X   write(6,*) 'ERWINCC solutions, with BLAS=',BLAS,NEQS
     X   , NEQS,SHSMAT,SMATEL
      write(48,*) 'ERWINCC solutions, with BLAS=',BLAS,NEQS
      call flush(48)
      NM1 = N-1
      NR = 2*NEQS
      NP = NR + 1
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
      ZI2 = (0.0D0,0.5D0)
      CI = (0.0D0,1.0d0)
      CZ = (0D0,0D0)
      ONEC = (1D0,0D0)
      RAD = 45d0/atan(1d0)
C     AL = 0.0
      ALUN = 0.0
      IS =  CUTOFF
      IS = MIN(1*N/2, MAX(2, IS ))
       if(SHOW.ge.3) write(6,*) 'ERWINCC: CUTOFF,IS =',CUTOFF,IS
       if(SHOW.ge.3) write(6,*) 'ERWINCC: CUTVAL =',(CUTVAL(I),I=1,NEQS)
      NN = NM1 - IS + 1
      TMAX = 20.
      TMIN = -125.
        SMALL = 1.0/FPMAX
        EPS = SQRT(SMALL)
        BIG = 1./EPS

C
      IF (FJSWTCH)  THEN
c                              skip nmrov , match CRCWFN to zero
c
761    DO 765 I=1,NR
       DO 765 J=1,NP
765       MAT(I,J) = 0
c
         DO 780 IT=1,NEQS
               MAT(IT+NEQS,IT) = ONEC
               MAT(IT,IT)   = ONEC*CRCRAT(it)
780       CONTINUE
C
      ELSE
C
      IF(SHOW.GE.2) WRITE(KO, 5) NF,M,NEQS,NR,REPEAT
      IF(SHOW.GE.2) WRITE(KO, 4) (H(K),K=1,NEQS)
4     FORMAT(' ERWINCC step sizes are',10F8.5)
5     FORMAT(' ERWINCC given',4I4,L4,', requires',3I6,';',I8,'/',I8,L3)
      CALL CHECK(NEQS,MAXB,15)
	if(REPEAT) go to 61
C
C
          DO 7 K=1,NEQS
            SMAT(K) = CZ
          H2C(K) = H(K)**2 / COEF(K)
7        CONTINUE
C
         DO 705 I=1,N
705           MAGN(I) = 1.0
C
      DO 11 IT=1,NEQS
          DO 10 K=1,NEQS
         V(IT,K)   = 0
         FI(IT,K) = 0
         ZM(IT,K)  = 0
         ZI(IT,K)  = 0
            FIM(IT,K) = 0
            FIMD(IT,K)= 0
10       CONTINUE
11     CONTINUE

C
        DO 141 K=1,NEQS
141      IF(SHOW.GE.9) WRITE(KO,142) K,H(K),COEF(K),(CH(II,K),II=1,2)
142      FORMAT(' FOR CH.',I4,' H,COEF,CH(1,2)=',10F10.5)
C-----------------------------------------------------------------------
	  if(LOCFIL) then
	  rewind 19
	  do I=1,IS-1
	   read(19) 
	  enddo
	  endif
         DO 60 I=IS,NM1
          RI2 = 1D0 / DBLE(I-1)**2
	  if(LOCFIL) read(19) FORMFR
          DO 144 K=1,NEQS
             C = CZ
	     if(LOCFIL) then
             DO NC=1,NCLIST(K,K)
	 	  JF = NFLIST(K,K,NC)
             IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,K,NC) * FORMFR(JF)
	     ENDDO
	     else
             DO NC=1,NCLIST(K,K)
	 	  JF = NFLIST(K,K,NC)
             IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,K,NC) * FORMF(I,JF)
	     ENDDO
	     endif
           SMAT(K) = -LL1(K)*RI2 + (-ECM(K) + C) * H2C(K)
144        CONTINUE
C
          DO 13 IT=1,NEQS
            K = IT
              IF(I.ne.CUTVAL(K)) GO TO 13
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
              ZI(IT,K) = SCL * EXP(T)
       	IF(SHOW.GE.3) write(KO,1302) I,IT,K,ZI(IT,K)
13         CONTINUE
1302         FORMAT(' AT I=',I3,' ZI(',I5,',',I5,')=',1P,2E12.2)

         DO 158 IT=1,NEQS
         DO 158 K=1,NEQS
158       FI(IT,K) = ZI(IT,K) *
     X      (ONEC - SMAT(K) * (R12 - SMAT(K)*(ENA2 + SMAT(K)*ENA3)))

         IF(SHOW.GE.7) THEN
          DO 161 K=1,NEQS
161          WRITE(KO,165) I,K,H(K)*(I-1),L(K),
     &    SMAT(K)/H2C(K)+ECM(K),FI(1,K)
165      FORMAT(' AT I,K,R =',I6,I4,F7.3,
     &          ' L,PE,PSI =',I4,2F9.3,1P,2E10.1)
          ENDIF
C
         DO 20 K=1,NEQS
            DO 17 J=1,NEQS
17          COUPL(K,J) = 0.0
           IF(I.LT.IC) GO TO 20
           DO 19 J=1,NEQS
            C = CZ
	    if(K==J) go to 19
	   if(LOCFIL) then
            DO 18 NC=1,NCLIST(K,J)
		JF = NFLIST(K,J,NC)
18          IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,J,NC) * FORMFR(JF)
	   else
            DO 181 NC=1,NCLIST(K,J)
		JF = NFLIST(K,J,NC)
181          IF(mod(PTYPE(6,JF),2).EQ.0) C=C+CLIST(K,J,NC) * FORMF(I,JF)
	   endif
19         COUPL(K,J) = C * H2C(K)
20	 CONTINUE
C THE DIAGONAL PART COUPL(K,K) IS ZERO
C
	if(BLAS) then
	 CALL ZGEMM('N','T',NEQS,NEQS,NEQS,-ONEC*R12,ZI,MAXB,COUPL,MAXB,
     X		    ONEC,FI,MAXB)
C
	 CALL ZGEMM('N','T',NEQS,NEQS,NEQS,ONEC,FI,MAXB,COUPL,MAXB,
     X		    CZ,V,MAXB)
	else
!	 CALL POTWF(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)
	 CALL POTWF2(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)   ! works better even on SUNs!
	endif

         DO 46 IT=1,NEQS
         DO 43 K=1,NEQS
            ZIV = ZI(IT,K)
            ZPV = ZIV + ZIV - ZM(IT,K) - V(IT,K) - SMAT(K) * FI(IT,K)
            ZM(IT,K) = ZIV
            ZI(IT,K) = ZPV
43            CONTINUE
46    CONTINUE
      DO 468 IT=1,NEQS
       IF(I.EQ.M) THEN
        DO 463 K=1,NEQS
463     FIM(IT,K) = FI(IT,K)
       ELSE IF(I.EQ.MD) THEN
        DO 465 K=1,NEQS
465     FIMD(IT,K) = FI(IT,K)
       ENDIF
468   CONTINUE
            T = 0.0
         DO 47 IT=1,NEQS
            DO 47 K=1,NEQS
            F8 = FI(IT,K)
            T = MAX(T, AMD1(F8))
   47       IF(T .GT. BIG ) GO TO 48
            MAGN(I) = T
         IF(I.NE.NM1) GO TO 60
 48      IF(T.EQ.0.0) GO TO 60
            MAGN(I) = T
         T = RENORM/T
         IF(SHOW.GE.1) WRITE(KO,49) T,I,RENORM
49       FORMAT(' Renormalising by',E12.4,' at',I4,' to max. WF.',E12.4)
         if(t.eq.0d0) then
          write(KO,*) 'RENORMALISING TO 0! at I=',I,' as lWF=',MAGN(I)
          write(KO,*) ' FI:',FI
          T = SMALL
          endif
         FMIN = small/T
         DO 55 IT=1,NEQS
         DO 51 K=1,NEQS
            ZI(IT,K) = ZI(IT,K) * T
            ZM(IT,K) = ZM(IT,K) * T
            FIM(IT,K) = FIM(IT,K) * T
            FIMD(IT,K) = FIMD(IT,K) * T
51       CONTINUE
55       CONTINUE
C
60       CONTINUE
C-----------------------------------------------------------------------
c
61    DO 65 I=1,NR
      DO 65 J=1,NP
65     MAT(I,J) = 0d0
         DO 70 IT=1,NEQS
         DO 70 K=1,NEQS
            MAT(K+NEQS,IT) = FIMD(IT,K)
            MAT(K,IT)   = FIM(IT,K)
70          CONTINUE
C
	ENDIF
C                    MATCH TO REQUIRED ASYMPTOTIC CONDITIONS
C
      IF (FCWFN) THEN
C                           match Numerov to CRCWFN
           IF(SHOW.ge.4) write(6,*) ' CG2,CF2(M) CG1,CF1(MD) ='
       DO 82 K=1,NEQS
       DO 82 K2=1,NEQS
           IF(SHOW.ge.4.and.ABS(CGMAT(K2,K,2)).gt.0.0) 
     *       WRITE(6,716) K,K2,CGMAT(K2,K,2),CFMAT(K2,K,2),
     *           CGMAT(K2,K,1),CFMAT(K2,K,1)
716            FORMAT(1X,2I3,': CRCWFN= ',4(D15.7,3X))
c
c james
         ZDL = CI**(L(K)-L(K2))
         ZHPLUS = CGMAT(K2,K,2)+CI*CFMAT(K2,K,2)	! M
         MAT(K2,K+NEQS)=   ZI2*ZDL*ZHPLUS
         ZHPLUS = CGMAT(K2,K,1)+CI*CFMAT(K2,K,1)   	! MD
         MAT(K2+NEQS,K+NEQS)= ZI2*ZDL*ZHPLUS
c
c        MAT(K2,K+NEQS)=   ZI2*(CGMAT(K2,K,2)+CI*CFMAT(K2,K,2))
c        MAT(K2+NEQS,K+NEQS)= ZI2*(CGMAT(K2,K,1)+CI*CFMAT(K2,K,1))
82       CONTINUE
         DO 85 K=1,NEQS
c james
         ZDL = CI**(L(EL)-L(K))
         ZHMINUS = CGMAT(K,EL,2)-CI*CFMAT(K,EL,2)
         MAT(K,NP)=   ZI2*ZDL*ZHMINUS
         ZHMINUS = CGMAT(K,EL,1)-CI*CFMAT(K,EL,1)
         MAT(K+NEQS,NP)= ZI2*ZDL*ZHMINUS
c
c          MAT(K,NP)=   ZI2*(CGMAT(K,EL,2)-CI*CFMAT(K,EL,2))
c          MAT(K+NEQS,NP)= ZI2*(CGMAT(K,EL,1)-CI*CFMAT(K,EL,1))
85       CONTINUE
      ELSE
c                         match Numerov to Uncoupled Coulomb wfns
       DO 90  K=1,NEQS
         MAT(K,K+NEQS) = CH(1,K)
         MAT(K+NEQS ,K+NEQS) = CH(2,K)
C                                      C(K) == MAT(K,NP) IS THE RHS.
         MAT(K,NP) = 0d0
         MAT(K+NEQS ,NP) = 0d0
90       CONTINUE
c                                         note minus as ch = i/2*H+
         MAT(EL,NP) = - CONJG(CH(1,EL))
         MAT(EL+NEQS ,NP) = - CONJG(CH(2,EL))
      ENDIF
         DO 100 I=1,NR
100        IF (SHOW.GE.5) WRITE(KO,105) I,(MAT(I,J),J=1,NP)
105      FORMAT(' MAT(',I2,',*) =', /(6(1P,E10.2,E9.1,1X)))
C
      CALL GAUSS(NR,NR,MAT,SING,S,SMALL,.FALSE.)
         IF(SING) GO TO 600
C
      T = EXP(DBLE(S) * 0.5/NEQS)
      IF(.NOT.REPEAT) RENORM = RENORM / T
       IF(RENORM.GT.BIG) RENORM = BIG
       IF(RENORM.LT.EPS) RENORM = EPS
      IF(SHOW.GE.3) WRITE(KO,130) (MAT(I,NP),I=1,NEQS)
130   FORMAT(/' The combining coefficients are',1P,/(2(E20.8,E15.8)))
C
c     ENDIF
c                  end of FJSWITCH loop matching CRCWFN to zero
      DO 150 K=1,NEQS
      SMAT(K) = MAT(K+NEQS,NP)
C
      IF(ABS(SMAT(K)).GT.1E-10) THEN
        T = 0.0
	WVD = 0.0
       DO 149 IT=1,NEQS
         IF(MAT(IT,NP).EQ.CZ) GO TO 149
            WVD = WVD + MAT(IT,NP)*FIMD(IT,K)
            T = MAX(T,ABS(WVD))
149      CONTINUE
             H2C(K) = T / (ABS(WVD)+1D-20)
           if(ECM(K)>0.0) then
             AL = MAX(AL, H2C(K))
             ALUN = MAX(ALUN, H2C(K)*ABS(SMAT(K)))
	   endif
!	     write(106,*) ALUN,K,WVD,H2C(k),smat(K),AL
      ENDIF
150      CONTINUE
      IF(SHOW.ge.1) WRITE(KO,1157) S,AL
1157  FORMAT(' ERWINCC: LOG10(DETERMINANT) =',2F11.3,', ACCURACY LOSS ='
     X  , 1P,E10.2)
C                       NOW CHECK THE LOGARITHMIC DERIVS. MATCH REQD.
      DO 200 K=1,NEQS
      IF((SHSMAT.OR.K.EQ.EL.AND.SMATEL).and.final)
     & WRITE(KO,190) K,SMAT(K),L(K),JVAL(K),CORESP(K),L(EL),
     &           LOG10(MAX(H2C(K),1D0))
190   FORMAT( '   S-matrix',I5,' = '
     & , 2F10.5,' for L=',I5,', J=',F7.1,' channel on core I =',F6.1,
     & ' from L=',I5,',  Acc. loss =',F5.1,' D.')
      IF(K.EQ.EL.AND.SMATEL.AND.ABS(SMAT(K)).GT.1.D-6.and.final) THEN
         SE = LOG(SMAT(EL))
         WRITE(KO,195) K,SE*(0.,-.5)*RAD,L(K),JVAL(K)
195      FORMAT( ' Elastic phase shift ',I3,' = '
     &     , 2F8.3,' deg. for the L =',I5,', J =',F7.1,' channel.')
         WRITE(45,196) ECM(EL),SE*(0.,-.5)*RAD,L(K),JVAL(K)
196        format(f10.3,2f9.3,' for LJin =',i6,f6.1)
	 written(45) = .true.
         ENDIF
200   CONTINUE
      if(.not.FJSWTCH) then
C     IF(AL.GT.1D10.OR.AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
      IF(AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
     X          3.*ALUN*RE * 100.
205   FORMAT('0****** WARNING : ACCURACY LOSS IN ERWINCC =',F5.1,
     X ' DIGITS, SO EXPECT ERRORS OF',F9.4,' % OF UNITARITY'/)
      IF(AL*RE.GT..03) WRITE(KO,*) AL,ALUN,RE
	endif
      RETURN
600   continue
      DO 605 I=IS,NM1
605   IF(MAGN(I).NE.0.0) MAGN(I) = LOG10(MAGN(I))
      WRITE(KO,610) (MAGN(I),I=IS,NM1)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      CALL ABEND(4)
      END
	SUBROUTINE POTWF(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)
	IMPLICIT NONE
	REAL*8 R12
	INTEGER NEQS,MAXB,K,J,IT,MAXCH
	COMPLEX*16 COUPL(MAXB,MAXB),ZI(MAXB,MAXCH),FI(MAXB,MAXCH),
     X		   V(MAXB,MAXCH),SE,ZIV
C					SCALAR MATRIX MULTIPLIES
C					NO SKIPS IF ZERO COUPLING

         DO 24 K=1,NEQS
             DO 22 IT=1,NEQS
	     SE = 0.0
	     DO J=1,NEQS
             SE = SE + COUPL(K,J) * ZI(IT,J)
	     ENDDO
 22          FI(IT,K) = FI(IT,K) - SE * R12
 24       CONTINUE
         DO 34 K=1,NEQS
            DO 33 IT=1,NEQS
	    ZIV = 0.
            DO J=1,NEQS
		ZIV = ZIV + COUPL(K,J) * FI(IT,J)
	    ENDDO
33          V(IT,K) = ZIV
34       CONTINUE
	return
	end
	SUBROUTINE POTWF2(COUPL,ZI,FI,V,NEQS,MAXB,MAXCH,R12)
	IMPLICIT NONE
	REAL*8 R12
	INTEGER NEQS,MAXB,K,J,IT,MAXCH
	COMPLEX*16 COUPL(MAXB,MAXB),ZI(MAXB,MAXCH),FI(MAXB,MAXCH),
     X		   V(MAXB,MAXCH),C,ZERO
	PARAMETER(ZERO = (0d0,0d0))
C					VECTORISED MATRIX MULTIPLIES
C					ALLOWS SKIPS IF ZERO COUPLING

         DO 24 K=1,NEQS
	     DO 24 J=1,NEQS
             C = COUPL(K,J) * R12
	      if(C/=ZERO) FI(1:NEQS,K) = FI(1:NEQS,K) - C * ZI(1:NEQS,J)
 24       CONTINUE
         DO 34 K=1,NEQS
	    V(:,K)  = ZERO
            DO 34 J=1,NEQS
	      C = COUPL(K,J)
	      if(C/=ZERO) V(1:NEQS,K) = V(1:NEQS,K) + C * FI(1:NEQS,J)
34       CONTINUE
	return
	end
****ERWINRC**************************************************************
      SUBROUTINE ERWINRC(ECM,COEF,FORMF,NF,FORMFR,CH,REPEAT,
     $  EL,SMAT,L,JVAL,CORESP,LL1,NEQS,N,H,M,MD,SCL,RENORM,
     $  AL,RE,IC,SHOW,SHSMAT,SMATEL,CUTOFF,PTYPE,FIM,FIMD,LOCFIL,
     $  CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)
	use io
	use factorials
	use parameters
	use drier
	use searchpar, only: final
!$      use omp_lib
#ifdef MPI
        use mpi
#endif /* MPI */
      IMPLICIT REAL*8(A-H,O-Z)
C
C        SOLVE 'NEQ' COUPLED SCHROEDINGERS EQUATIONS BY EXACT CC, REAL COUPLINGS
C        SOLVE  
C             (COEF(K). D2/DR2 + EN(R,K)).W(R,K)
C                  + SUM(J): COUPL(R,K,J).W(R,J) + INHOMG(R,K) = 0
C           WHERE COUPL(R,K,J) = SUM(NC): FORMF(R,JF)*CLIST(K,J,NC)
C           WHERE    JF = NFLIST(K,J,NC)
C       &   WHERE EN(R,K) IS THE DIAGONAL PART OF COUPL
C
C     ASSUMED FOR ARRAYS HERE THAT N<=MAXN AND NEQS <= MAXCH
C
C  Using 'Enhanced Numerov' of Thorlacius & Cooper (JCP 72(1987) 70)
C   with 5 terms in cosh(sqrt(12T)) expansion, but only diagonal potl.
C
      COMPLEX*16 MAT(2*NEQS,2*NEQS+1),S,CZ,SE,
     &           CH(2,MAXCH),SMAT(MAXCH),ONEC,ZI2,CI,ZDL,ZHPLUS,ZHMINUS,
     &	     FORMF(MAXM,NF),FORMFR(NF),CLIST(MAXCH,MAXCH,MCLIST)
      REAL*8     FIMD(NEQS,NEQS),FIM(NEQS,NEQS),C,ZP(NEQS),
     &           SR,ZIV,ZPV,FI(NEQS,NEQS),RLIST(MAXCH,MAXCH,MCLIST),
     &           ZI(NEQS,NEQS),ZM(NEQS,NEQS),F8
      REAL*8, allocatable:: COUPL(:,:,:),DIAG(:,:)
      REAL*8 COEF(MAXCH),JVAL(MAXCH),CORESP(MAXCH),H2C(MAXCH),MAGN(MAXN)
     &          ,ECM(MAXCH),H(MAXCH),LL1(MAXCH),
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),CRCRAT(MAXCH)
      INTEGER EL,L(MAXCH),CUTOFF,SHOW,
     X        CUTVAL(MAXCH),NFLIST(MAXCH,MAXCH,MCLIST),
     X        NCLIST(MAXCH,MAXCH),PTYPE(12,NF)
      LOGICAL SING,SHSMAT,SMATEL,FCWFN,FJSWTCH,BLAS,REPEAT,LOCFIL
      AMD1(C) = ABS(C) 
C
        numthread = 1
!$      numthread = OMP_GET_MAX_THREADS()
	if(NEQS>0.and.(SHSMAT.and.SMATEL))
     X   write(6,*) 'ERWINRC2 solutions =',NEQS,numthread
     X   , NEQS,SHSMAT,SMATEL
      write(48,*) 'ERWINRC2 solutions =',NEQS,numthread
	allocate(COUPL(NEQS,NEQS,numthread),DIAG(NEQS,numthread))
      call flush(48)
      NM1 = N-1
      NR = 2*NEQS
      NP = NR + 1
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
      ZI2 = (0.0D0,0.5D0)
      CI = (0.0D0,1.0d0)
      CZ = (0D0,0D0)
      ONEC = (1D0,0D0)
      RAD = 45d0/atan(1d0)
C     AL = 0.0
      ALUN = 0.0
      IS =  CUTOFF
      IS = MIN(1*N/2, MAX(2, IS ))
       if(SHOW.ge.3) write(6,*) 'ERWINRC: CUTOFF,IS =',CUTOFF,IS
       if(SHOW.ge.3) write(6,*) 'ERWINRC: CUTVAL =',(CUTVAL(I),I=1,NEQS)
      NN = NM1 - IS + 1
      TMAX = 20.
      TMIN = -125.
        SMALL = 1.0/FPMAX
        EPS = SQRT(SMALL)
        BIG = 1./EPS

C
      IF (FJSWTCH)  THEN
c                              skip nmrov , match CRCWFN to zero
  	 MAT(:,:) = 0.0
         DO 780 IT=1,NEQS
               MAT(IT+NEQS,IT) = ONEC
               MAT(IT,IT)   = ONEC*CRCRAT(IT)
 780       CONTINUE
C
      ELSE
C
      IF(SHOW.GE.2) WRITE(KO, 5) NF,M,NEQS,NR,REPEAT
      IF(SHOW.GE.2) WRITE(KO, 4) (H(K),K=1,NEQS)
4     FORMAT(' ERWINRC step sizes are',10F8.5)
5     FORMAT(' ERWINRC given',4I6,L4)
      CALL CHECK(NEQS,MAXB,15)
	if(REPEAT) go to 61
C
	MAGN(:) = 0.0
C
	do K=1,NEQS
	    DIAG(K,:) = 0d0
          H2C(K) = H(K)**2 / COEF(K)

         FI(:,K) = 0; ZM(:,K)  = 0; ZI(:,K)  = 0
	   FIM(:,K) = 0; FIMD(:,K)= 0
	 do J=1,NEQS
	 do nc=1,NCLIST(J,K)
	   RLIST(J,K,nc) = CLIST(J,K,nc) *(0.,1.)**(L(J)-L(K)) * H2C(K)
	   IF(mod(PTYPE(6,NFLIST(J,K,NC)),2)/=0) RLIST(J,K,nc) = 0.
	 enddo
	 enddo
	enddo
C
        DO 141 K=1,NEQS
 141     IF(SHOW.GE.9) WRITE(KO,142) K,H(K),COEF(K),(CH(II,K),II=1,2)
 142     FORMAT(' FOR CH.',I4,' H,COEF,CH(1,2)=',10F10.5)
C-----------------------------------------------------------------------

	 NBL = (NEQS-1)/numthread + 1

!$OMP  PARALLEL DO  
!$OMP&  PRIVATE (ITHR,I,K,J,IT,I1,I2,JF,R,T,F8,RI2,C,ZIV,ZPV,NC)

	DO 60 ITHR=1,numthread
	I1 = (ITHR-1)*NBL + 1
	I2 = min(ITHR*NBL,NEQS)
		write(48,*) 'Thread ',ITHR,': ',I1,I2
	  if(LOCFIL) then
	  rewind 19
	  do I=1,IS-1
	   read(19) 
	  enddo
	  endif
	  
         DO 55 I=IS,NM1
          RI2 = 1D0 / DBLE(I-1)**2
	    if(LOCFIL) read(19) FORMFR
	  
          DO 13 K=1,NEQS
             C = 0d0
	     if(LOCFIL) then
             DO NC=1,NCLIST(K,K)
		   C=C+RLIST(K,K,NC) * real(FORMFR(NFLIST(K,K,NC)))
	       ENDDO
	     else
             DO NC=1,NCLIST(K,K)
	 	   C=C+RLIST(K,K,NC) * real(FORMF(I,NFLIST(K,K,NC)))
	       ENDDO
	     endif
           DIAG(K,ITHR) = -LL1(K)*RI2 - ECM(K) * H2C(K) + C
	     
           IT = K
              IF(I==CUTVAL(K).and.I1<=IT.and.IT<=I2) then
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
              ZI(IT,K) = SCL * EXP(T)
!		  IF(SHOW.GE.3) write(KO,1302) I,IT,K,ZI(IT,K),0.0
		  endif
 13        CONTINUE

		IF(SHOW.GE.3) then
		do 14 K=I1,I2
		IT = K
 14		 if(I==CUTVAL(K)) write(KO,1302) I,IT,K,ZI(IT,K)
 1302         FORMAT(' AT I=',I3,' ZI(',I5,',',I5,')=',1P,2E12.2)
		endif

         DO 20 K=1,NEQS
           COUPL(K,:,ITHR) = 0.0
           IF(I.LT.IC) GO TO 20
           DO 19 J=1,NEQS
            C = 0d0
	    if(K==J) go to 19
	   if(LOCFIL) then
            DO 18 NC=1,NCLIST(K,J)
 18          C=C+RLIST(K,J,NC) * real(FORMFR(NFLIST(K,J,NC)))
	   else
            DO 181 NC=1,NCLIST(K,J)
 181          C=C+RLIST(K,J,NC) * real(FORMF(I,NFLIST(K,J,NC)))
	   endif
 19         COUPL(K,J,ITHR) = C 
 20	 CONTINUE
C THE DIAGONAL PART COUPL(K,K,ITHR) IS ZERO

         DO 158 K=1,NEQS
	     T = (1d0 - DIAG(K,ITHR) * 
     x                  (R12 - DIAG(K,ITHR)*(ENA2 + DIAG(K,ITHR)*ENA3)))
 158       FI(I1:I2,K) = ZI(I1:I2,K) * T

	   IF(SHOW.GE.7.and.I1==1) THEN
          DO K=1,NEQS
		WRITE(KO,165) I,K,H(K)*(I-1),L(K),
     &    DIAG(K,ITHR)/H2C(K)+ECM(K),0.0,FI(1,K),ZI(1,K)!*(0.,1.)**(L(EL)-L(K))
          ENDDO
 165      FORMAT(' AT I,K,R =',I6,I4,F7.3,
     &          ' L,PE,PSI =',I4,2F9.3,1P,2E10.1)
          ENDIF
	    
	    DO 24 K=1,NEQS
	     DO 24 J=1,NEQS
             C = COUPL(K,J,ITHR) * R12
	      if(C/=0d0) FI(I1:I2,K) = FI(I1:I2,K) - C * ZI(I1:I2,J)
 24       CONTINUE
	    
         DO 43 K=1,NEQS
	   
	    ZP(I1:I2) = 2d0*ZI(I1:I2,K)-ZM(I1:I2,K) - 
     x                      DIAG(K,ITHR)*FI(I1:I2,K)
	    
           DO 34 J=1,NEQS
	      C = COUPL(K,J,ITHR)
 34	      if(C/=0d0) ZP(I1:I2) = ZP(I1:I2) - C * FI(I1:I2,J)

            ZM(I1:I2,K) = ZI(I1:I2,K)
            ZI(I1:I2,K) = ZP(I1:I2)
 43       CONTINUE
	    
       IF(I.EQ.M) THEN
        FIM(I1:I2,:) = FI(I1:I2,:)
       ELSE IF(I.EQ.MD) THEN
        FIMD(I1:I2,:) = FI(I1:I2,:)
       ENDIF

            T = 0.0
         DO 47 IT=I1,I2
            DO 47 K=1,NEQS
            F8 = FI(IT,K)
            T = MAX(T, AMD1(F8))
   47       IF(T .GT. BIG ) GO TO 48
            MAGN(I) = T
         IF(I.NE.NM1) GO TO 51
 48      IF(T.EQ.0.0) GO TO 51
            MAGN(I) = T
         T = RENORM/T
         IF(SHOW.GE.1) WRITE(KO,49) T,I,RENORM
 49      FORMAT(' Renormalising by',E12.4,' at',I4,' to max. WF.',E12.4)
         if(t.eq.0d0) then
          write(KO,*) 'RENORMALISING TO 0! at I=',I,' as lWF=',MAGN(I)
          write(KO,*) ' FI:',FI
          T = SMALL
          endif
         FMIN = small/T
         DO 50 IT=I1,I2
            ZI(IT,:) = ZI(IT,:) * T
            ZM(IT,:) = ZM(IT,:) * T
            FIM(IT,:) = FIM(IT,:) * T
            FIMD(IT,:) = FIMD(IT,:) * T
 50       CONTINUE
C
C
 51      CONTINUE

 55       CONTINUE
 60       CONTINUE
!$OMP END PARALLEL DO
	deallocate (COUPL,DIAG)
C-----------------------------------------------------------------------
c
 61    DO 65 I=1,NR
 65    MAT(I,1:NP) = 0d0
         DO 70 IT=1,NEQS
         DO 70 K=1,NEQS
            MAT(K+NEQS,IT) = FIMD(IT,K)
            MAT(K,IT)   = FIM(IT,K)
 70          CONTINUE
C
	ENDIF
C                    MATCH TO REQUIRED ASYMPTOTIC CONDITIONS
C
      IF (FCWFN) THEN
C                           match Numerov to CRCWFN
           IF(SHOW.ge.4) write(6,*) ' CG2,CF2(M) CG1,CF1(MD) ='
       DO 82 K=1,NEQS
       DO 82 K2=1,NEQS
           IF(SHOW.ge.4.and.ABS(CGMAT(K2,K,2)).gt.0.0) 
     *       WRITE(6,716) K,K2,CGMAT(K2,K,2),CFMAT(K2,K,2),
     *           CGMAT(K2,K,1),CFMAT(K2,K,1)
716            FORMAT(1X,2I3,': CRCWFN= ',4(D15.7,3X))
c
c james
         ZDL = CI**(L(K)-L(K2))
         ZHPLUS = CGMAT(K2,K,2)+CI*CFMAT(K2,K,2)	! M
         MAT(K2,K+NEQS)=   ZI2*ZDL*ZHPLUS
         ZHPLUS = CGMAT(K2,K,1)+CI*CFMAT(K2,K,1)   	! MD
         MAT(K2+NEQS,K+NEQS)= ZI2*ZDL*ZHPLUS
82       CONTINUE
         DO 85 K=1,NEQS
c james
         ZDL = CI**(L(EL)-L(K))
         ZHMINUS = CGMAT(K,EL,2)-CI*CFMAT(K,EL,2)
         MAT(K,NP)=   ZI2*ZDL*ZHMINUS
         ZHMINUS = CGMAT(K,EL,1)-CI*CFMAT(K,EL,1)
         MAT(K+NEQS,NP)= ZI2*ZDL*ZHMINUS
c
85       CONTINUE
      ELSE
c                         match Numerov to Uncoupled Coulomb wfns
       DO 90  K=1,NEQS
         MAT(K,K+NEQS) = CH(1,K)
         MAT(K+NEQS ,K+NEQS) = CH(2,K)
C                                      C(K) == MAT(K,NP) IS THE RHS.
         MAT(K,NP) = 0d0
         MAT(K+NEQS ,NP) = 0d0
90       CONTINUE
c                                         note minus as ch = i/2*H+
         MAT(EL,NP) = - CONJG(CH(1,EL))
         MAT(EL+NEQS ,NP) = - CONJG(CH(2,EL))
      ENDIF
         DO 100 I=1,NR
100        IF (SHOW.GE.5) WRITE(KO,105) I,(MAT(I,J),J=1,NP)
105      FORMAT(' MAT(',I2,',*) =', /(6(1P,E10.2,E9.1,1X)))
C
      CALL GAUSS(NR,NR,MAT,SING,S,SMALL,.FALSE.)
         IF(SING) GO TO 600
C
      T = EXP(DBLE(S) * 0.5/NEQS)
      IF(.NOT.REPEAT) RENORM = RENORM / T
       IF(RENORM.GT.BIG) RENORM = BIG
       IF(RENORM.LT.EPS) RENORM = EPS
      IF(SHOW.GE.3) WRITE(KO,130) (MAT(I,NP),I=1,NEQS)
130   FORMAT(/' The combining coefficients are',1P,/(2(E20.8,E15.8)))
C
c     ENDIF
c                  end of FJSWITCH loop matching CRCWFN to zero
      DO 150 K=1,NEQS
      SMAT(K) = MAT(K+NEQS,NP) *(0.,1.)**(L(EL)-L(K))
C
      IF(ABS(SMAT(K)).GT.1E-10) THEN
        T = 0.0
	WVD = 0.0
       DO 149 IT=1,NEQS
         IF(MAT(IT,NP).EQ.CZ) GO TO 149
            WVD = WVD + MAT(IT,NP)*FIMD(IT,K)
            T = MAX(T,ABS(WVD))
149      CONTINUE
             H2C(K) = T / (ABS(WVD)+1D-20)
           if(ECM(K)>0.0) then
             AL = MAX(AL, H2C(K))
             ALUN = MAX(ALUN, H2C(K)*ABS(SMAT(K)))
	   endif
!	     write(106,*) ALUN,K,WVD,H2C(k),smat(K),AL
      ENDIF
150      CONTINUE
      IF(SHOW.ge.1) WRITE(KO,1157) S,AL
1157  FORMAT(' ERWINRC: LOG10(DETERMINANT) =',2F11.3,', ACCURACY LOSS ='
     X  , 1P,E10.2)
C                       NOW CHECK THE LOGARITHMIC DERIVS. MATCH REQD.
      DO 200 K=1,NEQS
      IF((SHSMAT.OR.K.EQ.EL.AND.SMATEL).and.final)
     & WRITE(KO,190) K,SMAT(K),L(K),JVAL(K),CORESP(K),L(EL),
     &           LOG10(MAX(H2C(K),1D0))
190   FORMAT( '   S-matrix',I5,' = '
     & , 2F10.5,' for L=',I5,', J=',F7.1,' channel on core I =',F6.1,
     & ' from L=',I5,',  Acc. loss =',F5.1,' D.')
      IF(K.EQ.EL.AND.SMATEL.AND.ABS(SMAT(K)).GT.1.D-6.and.final) THEN
         SE = LOG(SMAT(EL))
         WRITE(KO,195) K,SE*(0.,-.5)*RAD,L(K),JVAL(K)
195      FORMAT( ' Elastic phase shift ',I3,' = '
     &     , 2F8.3,' deg. for the L =',I5,', J =',F7.1,' channel.')
         WRITE(45,196) ECM(EL),SE*(0.,-.5)*RAD,L(K),JVAL(K)
196        format(f10.3,2f9.3,' for LJin =',i6,f6.1)
	 written(45) = .true.
         ENDIF
200   CONTINUE
      if(.not.FJSWTCH) then
C     IF(AL.GT.1D10.OR.AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
      IF(AL*RE.GT..03) WRITE(KO,205) LOG10(MAX(AL,1D0)),
     X          3.*ALUN*RE * 100.
205   FORMAT('0****** WARNING : ACCURACY LOSS IN ERWINRC =',F5.1,
     X ' DIGITS, SO EXPECT ERRORS OF',F9.4,' % OF UNITARITY'/)
      IF(AL*RE.GT..03) WRITE(KO,*) AL,ALUN,RE
	endif
      RETURN
600   continue
      DO 605 I=IS,NM1
605   IF(MAGN(I).NE.0.0) MAGN(I) = LOG10(MAGN(I))
      WRITE(KO,610) (MAGN(I),I=IS,NM1)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      CALL ABEND(4)
      END
	
