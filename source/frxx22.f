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
*****FRXX10*************************************************************
      SUBROUTINE KERNES(FFR,KCOEF,NLL,NONO,IREM,ICV,NM,LDMIN,LDMAX,
     &               C1,C2,IC1,IC2,REPEAT,  LVAL,ICTO,ICFROM,REV,PART,
     &               NK,NG,FPT,GPT,CHNO,CP, MCG,MAXMV,MKNL,NL0,QNF,
     &               KNL,NUMLT,LTMIN,IC7,LCMIN,LCMAX,MXLN1,MXLNP1,MMX1,
     X               MXMV1,MXMVP1)
	use parameters
	use io
	use kcom
	use drier
	use trace
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 KCOEF(MAXQRN,NUMLT),MCG(0:3,MAXMV,NG,MKNL)
      INTEGER LVAL(MAXCH),D,C,PART(MAXCH),GPT(2,NG),CP,C1,C2,
     &    QNF(19,MSP),CHNO(MFNL,6),FPT(7,NK),NM(MAXQRN,MFNL)
      LOGICAL REV,PRES,IFAIL3,THERE,FFR,REPEAT,C1FR,NONO,
     X        LTRANS(MAXQRN)
      DATA ONE/1D0/
      IFAIL3(I,J,K) = I.GT.J+K .OR. I.LT.ABS(J-K)
C
!	write(0,*) 'KERNES started'
      Z = 0.0
      EPS = 1E-14
      PI = 4.0*ATAN(ONE)
      PISQ8 = 8.0*PI**2
      REPEAT = .TRUE.
      C1FR = ICFROM.EQ.IC1
      IF(C1FR) THEN
         C = C1
         D = C2
      ELSE
         C = C2
         D = C1
      ENDIF
      IF(ICTO.NE.PART(D).OR.ICFROM.NE.PART(C)) STOP 7
C       SO TO 'D' FROM 'C' CHANNEL
C  AND WITH NAME(1,IC1) LIKE D & NAME(1,IC2) LIKE P IN (D,P) REACTION
      LD = LVAL(D)
      LC = LVAL(C)
      LDMIN = MIN(LD,LDMIN);  LDMAX = MAX(LD,LDMAX)
      LCMIN = MIN(LC,LCMIN);  LCMAX = MAX(LC,LCMAX)
      THERE = .FALSE.
      DO 8 IN=1,NG
      DO 8 ILT=1,NUMLT
 8     IF(ABS(KCOEF(IN,ILT)).GT.EPS) THERE = .TRUE.
      IF(.NOT.THERE) GO TO 100

      NL = NL + 1
      NLD = NL+1  ! derivative term
      DRY = DRY .OR. NLD.GT.MFNL
      IF(NLD.GT.MFNL)  WRITE(KO,9) NLD,MFNL
9     FORMAT(//' ****** NOT ENOUGH ROOM FOR',I4,' NON-LOCAL FORM',
     &' FACTORS  IN MFNL ARRAY OF',I4,' *****')
      KNL = NL-NL0
      CHNO(NL,1) = D
      IF(.NOT.REV) CHNO(NL,1) = -D
      CHNO(NL,2) = C
      CHNO(NL,3) = 2
      IF(NONO) CHNO(NL,3) = IC7
      IF(.NOT.FFR) CHNO(NL,3) = CHNO(NL,3) + 10
      CHNO(NL,4) = NLL
      CHNO(NL,5) = 0  ! X term first
      CHNO(NL,6) = 1
      
      CHNO(NLD,1) = D
      IF(.NOT.REV) CHNO(NLD,1) = -D
      CHNO(NLD,2) = C
      CHNO(NLD,3) = 2
      IF(NONO) CHNO(NLD,3) = IC7
      IF(.NOT.FFR) CHNO(NLD,3) = CHNO(NLD,3) + 10
      CHNO(NLD,4) = NLL
      CHNO(NLD,5) = 1  ! Y term second
      CHNO(NLD,6) = 1
     
      
      PRES = .FALSE.
!      MCG(:,:,:,KNL) = 0.0
C
      DO 40 IN=1,NG
       NM(IN,KNL) = 0
       NM(IN,KNL+1) = 0
C
C   HERE, TO START :    GPT(1  IS NO. OF B.S. IN CHANNEL 'C1' ("DEUT")
C                   AND GPT(2  IS NO. OF B.S. IN CHANNEL 'C2' ("PROT")
         IF(C1FR) THEN
           KNP = GPT(1,IN)
           KN = GPT(2,IN)
         ELSE
           KNP = GPT(2,IN)
           KN  = GPT(1,IN)
         ENDIF
C
C   NOW, MORE PROPERLY, KNP IS NO. OF B.S. IN CHANNEL 'C' (FROM)
C                   AND KN  IS NO. OF B.S. IN CHANNEL 'D' (TO)

C  NOTE that the primes (P) are here in the 'FROM' channel: right rather than left!!!
C        (opposite of the 'Surface' notes!)

         LN = QNF(9,KN)
         LNP= QNF(9,KNP)
         MXLN1 = MAX(MXLN1,LN+1,2)
         MXLNP1= MAX(MXLNP1,LNP+1,2)
         THERE = .FALSE.
         DO 16 ILT=1,NUMLT
 16       IF(ABS(KCOEF(IN,ILT)).GT.EPS) THERE = .TRUE.
         IF(.NOT.THERE) GO TO 40
      IF(LISTCC.GT.1) WRITE(KO,18) NL,D,C,IN,KNP,KN,NUMLT,LTMIN
 18   FORMAT('0NL. S-interaction #',I3,' to',I3,' from',I3,' :',5I3)
C
      IM = 0
      DO 30 M=-1,1
      DO 30 MVP=-LNP,LNP
      DO 30 MV =-LN,LN
        MP = M + MVP - MV
        IF(ABS(M).GT.LC) GO TO 30
        IF(ABS(MP).GT.LD) GO TO 30
!        IF(MP.LT.0) GO TO 30
!          R0 = 2.0
!          IF(MP.EQ.0 .AND. MV.EQ.0) R0 = 1.0
!          IF(MP.EQ.0 .AND. MV.LT.0) GO TO 30
	  R0 = 1.0
        IM = IM + 1
        MCG(0,IM,IN,KNL:KNL+1) = M
        MCG(1,IM,IN,KNL:KNL+1) = MV
        MCG(2,IM,IN,KNL:KNL+1) = MVP
        MMX1 = MAX(ABS(MP)+1,MMX1)
        MXMV1 = MAX(ABS(MV)+1,MXMV1,2)
        MXMVP1 = MAX(ABS(MVP)+1,MXMVP1,2)
        MCG(3,IM,IN,KNL:KNL+1) = 0.0
        R1 = PISQ8 * YLMC(LD,MP) * YLMC(LN,MV) 
!					omit initials! :       * YLMC(LC,M) *  YLMC(LNP,MVP)
       if(C1FR) then
           R1 = R1 * (-1)**iabs(MP)
           else
           endif
         DO 25 ILT=1,NUMLT
            LTOTAL = LTMIN + ILT - 1
            IF(LTOTAL.LT.0.OR.ABS(MP+MV).GT.LTOTAL) GO TO 25
            IF(IFAIL3(LD,LN,LTOTAL).OR.IFAIL3(LC,LNP,LTOTAL)) GO TO 25
          R2 = KCOEF(IN,ILT)
          R3 = WIG3J(LD+Z,LN+Z, LTOTAL+Z,MP+Z,MV+Z, Z-MP-MV)
          R4 = WIG3J(LC+Z,LNP+Z,LTOTAL+Z,M+Z ,MVP+Z,Z-M-MVP)
          R = R0 * R1 * R2 * R3 * R4
          IF(LISTCC.GT.2) WRITE(KO,23) D,C,IN,M,MVP,MV,LTOTAL,
     X          R0,R1,R2,R3,R4,R
23        FORMAT('0KERNES1:',2i6,5I3,3x,5F8.4,F9.6)
          PRES = PRES .OR. ABS(R) .GT. EPS
         MCG(3,IM,IN,KNL) = MCG(3,IM,IN,KNL) + R
25       CONTINUE
         MCG(3,IM,IN,KNL+1) = MCG(3,IM,IN,KNL)  !  same KCOEF and final CG coefficients
          IF(LISTCC.GT.1) then
             WRITE(KO,26) D,C,IN,M,MVP,MV,IM,NL,MCG(3,IM,IN,KNL)
             WRITE(KO,26) D,C,IN,M,MVP,MV,IM,NLD,MCG(3,IM,IN,KNL+1),' D'
26        FORMAT('0KERNES2:',2i6,6I3, F8.4,a2)
                call flush(KO)
             endif
30    CONTINUE
      NM(IN,KNL) = IM
      NM(IN,KNL+1) = IM
40    CONTINUE
      NL = NLD ! include both! 
      KNL = KNL+1
      IF(.NOT.PRES) NL = NL - 2
          IF(LISTCC.GT.1) WRITE(KO,*) 'KERNES done to ',NL,KNL
          call flush(KO)      
100   RETURN
      END
      SUBROUTINE QERNES(NLN,NLM,NLO,NK,NG,XA,XB,XP,XQ,FPT,CUTOFF,
     &             MCG,MAXMV,NL0,KNL,NM,MMX1,C1FR,cxwf,RMAS,
     &             RINTO,EPC,HNL,NLT,LCMAX,LDMAX,CP,
     &             RHO,THM,NLL,QNF,NLOC,LRANGF1,
     &             FORML,FORMC,MINT,RIN,NLC,CENTRE,HF,HT,
     &             FFREAL,JACOBN,LISTCC,CHNO,RERR,NCH,LVAL,NLPLOT,
     X             MXMV1,MXMVP1,MXLN1,MXLNP1)
	use parameters
	use io
	use drier
	use fresco1, only:cutl
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MCG(0:3,MAXMV,NG,KNL),ALT(LDMAX+1,MMX1),AL(LCMAX+1,1),
     X       ALN(MXLN1,MXMV1),ALNP(MXLNP1,MXMVP1)
      REAL*8 JACOBN,COFFIN(4),TH,THM(MAXNLN),RUM(KNL),RMAS(MSP)
      INTEGER FPT(7,NK),QNF(19,MSP),CUTOFF,NM(MAXQRN,KNL),CHNO(MFNL,6),
     X        CP,LVAL(NCH)
      COMPLEX*16 FFC,SUMT,FNC(NLL,NLO,KNL), SUM(KNL),S,S1,S2,S3,CG,
     X          FORMC(MAXNLC,MSP),VSC(MAXNLC,MSP),WF1C,WF2C(NK),FFC4,
     X          FORMCFD(MAXNLR,NK),FORMCTD(MAXNLR,NK),WF1CD,WF2CD(NK)
      LOGICAL FFREAL,C1FR,cxwf,PR,PPR
      REAL*8 FNL(NLL,NLO,KNL),NLOC(NLM),VR,WF2,WF1,FFR,FFR4,
     X  FORML(MAXNLR,MSP),FORMLFD(MAXNLR,NK),FORMLTD(MAXNLR,NK),Z,
     x  cleba(MXLNP1-1,-(MXLNP1-1):MXLNP1-1,-1:1),  		! cleb6(lnp-1,mvp-lam,1,lam,lnp,mvp) = cleba(lnp,mvp,lam)
     x  clebb(LCMAX,-1:1,-1:1)			! cleb6(LF-1,M-lam,1,lam,LF,M) 	 = clebb(LF,M,lam)
      REAL*8,allocatable:: YLMCF(:,:)
      DATA  PI /3.14159D0/
C
	   PPR=.false.  ! do not print non-local couplings at all

          IF(LISTCC.GT.1) WRITE(KO,*) 'QERNES starting, rho=',RHO,KO
!          call flush(KO)      
      CALL CCTIME(IPTIME)
      EP = EPC * 0.01
      IF(EP.lt.1e-5) EP = 1E-5
      Z = 0d0

      RFRTM = 1.0/((NLL-1)*RINTO)**2
      MMX = MMX1 - 1
      MXLN = MXLN1 - 1
      MXLNP= MXLNP1 -1
      MXMV = MXMV1 - 1
      MXMVP= MXMVP1 -1
      MX = 0 !  as hat(R)=hat(z)
C
C   In QERNES, we use a pure initial/final set of parameters,
C      as we are not concerned with any projectile/target asymmetries.
C      (These exist in the coupling order, but the transfer kernel
C
C       <LT,LN; LTOTAL | Surf | LF,LNP; LTOTAL>
C
C  ??     being calculated here IS symmetric)
C
C    FROM channel LF, radius RF=ri, M-projection = 0,
C         internal b.s radius R1, l,M = LNP,MVP, state KNF
C
C    TO   channel LT, radius RT=rf, M-projection = MM,
C         internal b.s radius R2, l,M = LN,MV,   state KNT
C
C   NOTE THAT A,B,P,Q ADHERE TO THIS CONVENTION:
C       R1 = P * RF + Q * RT = P * ri + Q * rf     (ri & rf in Buttle's
C       R2 = A * RF + B * RT = A * ri + B * rf       DAISY)
C   AS, IN QERNEL, A,B, P,Q ARE DEFINED FOR:
C       RDINT(proj) = XP * RF + XQ * RT
C       RN(targ)    = XA * RF + XB * RT
C
      IF(C1FR) THEN
C                  from: R1=RDINT, to: R2=RN
        A = XA
        B = XB
        P = XP
        Q = XQ
       ELSE
C                  from: R1=RN, to: R2=RDINT
        A = XP
        B = XQ
        P = XA
        Q = XB
       ENDIF
C
C	Now calculate coefficients for the Surface Operator expressions:
C
	RHO2 = RHO**2

	SA = P; SB = Q; SAP = A; SBP = B
	
	SP = P/A ; SQ = Q - P*B/A ; SPP = 1d0/A ; SQQ = -B/A      
	
	CON1 = RHO/(SAP*SBP) 
	CON2 = - JACOBN
	
	
	 do INF=1,NK
	   IF(C1FR) THEN
	     KNF = FPT(1,INF)
	     KNT = FPT(2,INF)
	   ELSE
	     KNF = FPT(2,INF)
	     KNT = FPT(1,INF)
	   ENDIF
	 LNP= QNF(9,KNF)
	 LN = QNF(9,KNT)
	 if(.not.cxwf) then
        CALL DERIV(FORML(1,KNF),FORMLFD(1,INF),1d0/RIN,MINT)           
        CALL DERIV(FORML(1,KNT),FORMLTD(1,INF),1d0/RIN,MINT)  
		do I=1,MINT
	  	  RF = (I-1)/RIN
		  RT = RF
!		write(300+KNT,*) RF,RT*FORML(I,KNT),
!    x            RT*(FORMLTD(I,INF)+FORML(I,KNT)/RF),
!    x            RT* FORMLTD(I,INF)
!		write(320+KNT,*) RF,FORML(I,KNT)
		enddo
          else
        CALL DERIVC(FORMC(1,KNF),FORMCFD(1,INF),1d0/RIN,MINT)           
        CALL DERIVC(FORMC(1,KNT),FORMCTD(1,INF),1d0/RIN,MINT)  
          endif         
   	enddo ! INF derivatives                

      PA = P - A
      QB = Q - B
      HFHT = HF/HT
      PQ = P*HFHT + Q
      AB = A*HFHT + B
      PAQB2 = PA*QB*2.
      PAQB  = PA*HFHT+QB
      AB2 = A*B*2
      PQ2 = P*Q*2
      IF(LISTCC.GE.3) WRITE(KO,10) ' A,B,P,Q,SAP,SBP,JACOBN',
     x                               A,B,P,Q,SAP,SBP,JACOBN
 10    FORMAT(1X,'QERNES: ', A20 / (1X,10F8.4))

          IF(FFREAL)      FNL(1:NLL,1:NLO,1:KNL) = 0.0
          IF(.NOT.FFREAL) FNC(1:NLL,1:NLO,1:KNL) = 0.0
        THM(:) = 0.0
        CALL PLM(1d0,LCMAX,MX,LCMAX+1,AL)      ! Y_{L_f}(R_f)
        
        LYMAX = max(LCMAX,LDMAX,MXLN,MXLNP)
        MYMAX = max(MX,MMX,MXMV,MXMVP)
        allocate(YLMCF(0:LYMAX,-MYMAX:MYMAX))
         YLMCF(:,:) = 0.0
        do 12 L=0,LYMAX
         do 12 M=-min(L,MYMAX),min(MYMAX,L)
         YLMCF(L,M) = YLMC(L,M)
12         continue
	do 14 LF=1,LCMAX
	do 14 M=-1,1
	do 14 lam=-1,1
!     x  clebb(0:LCMAX,-1:1,-1:1)			! cleb6(LF-1,M-lam,1,lam,LF,M) 	 = clebb(LF,M,lam)
14	clebb(LF,M,lam) = cleb6(LF-1d0,M-lam+Z,1d0,lam+Z,LF+Z,M+Z)
	do 16 lnp=1,MXLNP1-1
	do 16 mvp=-lnp,lnp
	do 16 lam=-1,1
16	cleba(lnp,mvp,lam) = 
     x        cleb6(lnp-1d0,mvp-lam+Z,1d0,lam+Z,lnp+Z,mvp+Z)
        
          CALL PLM(1d0,LCMAX,MX,LCMAX+1,AL)      ! Y_{L_f}(R_f)
C
	    R2 = RHO
               DO 21 INF=1,NK
                     IF(C1FR) THEN
                       KNT = FPT(2,INF)
                     ELSE
                       KNT = FPT(1,INF)
                     ENDIF
                if(.not.cxwf) then  
                 WF2C(INF)= FFR4(R2*RIN,FORML(1,KNT),MINT)
                 WF2CD(INF)= FFR4(R2*RIN,FORMLTD(1,INF),MINT)
                 else  ! cxwf
                 WF2C(INF)= FFC4(R2*RIN,FORMC(1,KNT),MINT) 
                 WF2CD(INF)= FFC4(R2*RIN,FORMCTD(1,INF),MINT)
                 endif ! ~ cxwf        
21		 CONTINUE        
        
C------------------------------------------START R-TO LOOP
      DO 100 I=2,NLL
         RT = (I-1) * RINTO
         IF(RT.LT.(CUTOFF-1)*HT) GO TO 100
          RFMAX = (RHO-B*RT)/A
          RFMIN = (RHO+B*RT)/A
         if(LISTCC>4) write(KO,23) RFMIN,RFMAX, abs(Q*RT+P*RFMAX)
23		format(/f8.3,'< RF < ',f8.3,',  r > ',f8.3)

C------------------------------------------START R-FROM LOOP
         DO 85 J=1,NLM
            DNL = (J - NLM/2 - 1) * HNL  +  CENTRE
            RT0 = RT*HFHT
            RF  = RT0    + DNL
            IF(RF.LT.(CUTOFF-1)*HF .OR. RF.LE.0.5.or.
     x         RF < abs(cutl)*HF ) GO TO 85
            PQR = P * DNL + PQ * RT
            ABR = A * DNL + AB * RT
            PQRS = PQR*PQR
            ABRS = ABR*ABR
            PAQBR = (PA * DNL + PAQB*RT)**2
            RFRT = RF*RT
!		from 					R2SQ  =  ABRS + AB2*RFRT*UK           
            UK = (RHO2 - ABRS)/(AB2*RFRT)
            U = UK + 1.0
            if(ABS(U)>1d0) go to 85
            TH = ACOS(U)
            THM(I) = max(THM(I),TH)
		SINTH = sin(TH)
  
                
               R1SQ  =  PQRS + PQ2*RFRT*UK
               R1 = SQRT(ABS(R1SQ))
               R2 = RHO
               COST = (ABR + B * RT * UK)/R2
               COSF = (PQR + Q * RT * UK)/R1
  	if(LISTCC>6) write(KO,25) RT,RF,U,TH*180/pi,U
25 	 format(' RT,RF,TH =',3f8.3,f8.1,f8.3)

		 	SINF = SBP*SINTH*RT/R1
		 	SINT = SB *SINTH*RT/R2
               
!        	   CALL PLM(1d0,LCMAX,MX,LCMAX+1,AL)      ! Y_{L_f}(R_f)
               CALL PLM(U,LDMAX,MMX,LDMAX+1,ALT)       ! Y_{L_t}(R_t)
               CALL PLM(COST,MXLN,MXMV,MXLN1,ALN)     ! Y_{l_t}(r_t)
               CALL PLM(COSF,MXLNP,MXMVP,MXLNP1,ALNP) ! Y_{l_f}(r_f)
               
                IF(FFREAL)      RUM(1:KNL) = 0.0
                IF(.NOT.FFREAL) SUM(1:KNL) = 0.0
C
               DO 50 INF=1,NK
                   IN = FPT(3,INF)
C  Sum all the INF with the coefficient MCG(3,IM,IN,IL) (i.e. using IN)
                     IF(C1FR) THEN
                       KNF = FPT(1,INF)
                       KNT = FPT(2,INF)
                     ELSE
                       KNF = FPT(2,INF)
                       KNT = FPT(1,INF)
                     ENDIF
                   LNP= QNF(9,KNF)
                   LN = QNF(9,KNT)
                   
			RMASS2 =   RMAS(KNT) ! reduced mass of final single particle
			CON3 = 1d0/(FMSCAL*RMASS2)
			CON = CON1 * CON2 * CON3
	
                if(.not.cxwf) then  ! all of these arrays had 1/R factors
                 WF1C= FFR4(R1*RIN,FORML(1,KNF),MINT)
!                 WF2C= FFR4(R2*RIN,FORML(1,KNT),MINT)
                 WF1CD= FFR4(R1*RIN,FORMLFD(1,INF),MINT)
!                 WF2CD= FFR4(R2*RIN,FORMLTD(1,INF),MINT)
                 else  ! cxwf
                 WF1C= FFC4(R1*RIN,FORMC(1,KNF),MINT)
!                 WF2C= FFC4(R2*RIN,FORMC(1,KNT),MINT) 
                 WF1CD= FFC4(R1*RIN,FORMCFD(1,INF),MINT)
!                 WF2CD= FFC4(R2*RIN,FORMCTD(1,INF),MINT)
                 endif ! ~ cxwf
C
           DO 49 IL=1,KNL
              LT = LVAL(IABS(CHNO(IL+NL0,1)))
              LF = LVAL(IABS(CHNO(IL+NL0,2)))
           if(LISTCC>6) write(KO,*) ' IN,IL,NM =',IN,IL,NM(IN,IL),':',LT,LF
     x  , 'wf2,wf2d =',real(WF2C(INF)),real(WF2CD(INF)),INF,KNT,LN
           if(LISTCC>6)  call flush(KO) 
               DO 48 IM=1,NM(IN,IL)
                if(Abs(MCG(3,IM,IN,IL))<1d-20) go to 48
                  M   = NINT(MCG(0,IM,IN,IL))
                  MV  = NINT(MCG(1,IM,IN,IL))
                  MVP = NINT(MCG(2,IM,IN,IL))
                  MP = M + MVP - MV
		           SPH1 =  0d0 
		  if(M==0) SPH1 =  YLMCF(LF,M) * AL(LF+1,IABS(M)+1) 
		  SPH2 =  YLMCF(LNP,MVP) * ALNP(LNP+1,IABS(MVP)+1)
		S1=0d0; S2=0d0; S3=0d0; S4=0d0; S5=0d0; S=0d0
	     if(CHNO(IL+NL0,5)==0)  then  ! wf

		if(M==0) then

		 S1 = WF2CD(INF) * SPH1 * SPH2 * WF1C
           if(LISTCC>6) write(KO,35) ' S1 =',WF2CD(INF),SPH1,SPH2,WF1C
35	   format(10x,a,9g12.5)
		
		 S2A = 0
		 do lam=-1,1
		 if(lnp-1>=0 .and.abs(mvp-lam)<=lnp-1) 
     x	 S2A  = S2A + cleba(lnp,mvp,lam)  ! cleb6(lnp-1d0,mvp-lam+Z,1d0,lam+Z,lnp+Z,mvp+Z)   
     x                    * YLMCF(1,lam) * ALN(1+1,IABS(lam)+1) 
     x         * YLMCF(lnp-1,mvp-lam) * ALNP(lnp-1+1,IABS(mvp-lam)+1) 
		  enddo
		 S2B = sqrt(4d0*PI*lnp*(2*lnp+1)/3d0) * WF1C 
		 S2C = - WF2C(INF) * SPH1 * SP/R1 
		 S2 = S2C * S2B * S2A
		 
		 S3A =  WF1CD - LNP*WF1C/R1
		 S3B = COST*COSF + SINT*SINF
		 S3C = - WF2C(INF) * SPH1 * SPH2 * SP
		 S3 =  S3C * S3B * S3A
!	   if(LISTCC>4) write(KO,35) ' S3 =',real(WF2C(INF)),S3C,S3B,S3A,S3
		endif ! M==0	 		 
	
		 S4A = 0
		 do lam=-1,1
		 if(LF-1>=0 .and.abs(M-lam)==0) then  ! M=lam as hat(R)=hat(z)
		   if(LISTCC>6) write(KO,35) ' S4A =',real(lam),
     X                      clebb(LF,M,lam) ! cleb6(LF-1d0,M-lam+Z,1d0,lam+Z,LF+Z,M+Z) 
     x                    , YLMCF(1,lam) , ALN(1+1,IABS(lam)+1) 
     x  , YLMCF(LF-1,M-lam) , AL(LF-1+1,IABS(M-lam)+1),LF+Z,M-lam+Z 
       	 S4A  = S4A + clebb(LF,M,lam) ! cleb6(LF-1d0,M-lam+Z,1d0,lam+Z,LF+Z,M+Z) 
     x                    * YLMCF(1,lam) * ALN(1+1,IABS(lam)+1) 
     x                    * YLMCF(LF-1,M-lam) * AL(LF-1+1,IABS(M-lam)+1) 
     			endif
		  enddo ! lam
		 S4B = sqrt(4d0*PI*LF*(2*LF+1)/3d0)
		 S4C = - WF2C(INF) * WF1C * SPH2 * SPP/RF 
		 S4 = S4C * S4B * S4A
              if(LISTCC>6) write(KO,35) ' S4 =',S4C,S4B,S4A

	     else if(M==0) then  ! DERIV term  (CHNO(IL+NL0,5)==1) (only when M==0)
	
		S5A = - SPP * COST
		S5B = WF2C(INF) * WF1C
		S5C = SPH1 * SPH2
		S5 = S5A * S5B * S5C
	
	     endif
		 S = S1 + S2 + S3 + S4 + S5
!		 S = S1

 43            CG =  CON * MCG(3,IM,IN,IL) * S 
     X                * ALT(LT+1,IABS(MP)+1)
     X                * ALN(LN+1,IABS(MV)+1)
       F = 1!sd6
     	 if(LISTCC>4 .and.abs(S)>1d-100) 
     x    write(KO,45) RT,RF,TH*180/pi,INF,IL,IM,M,CHNO(IL+NL0,5),
     X   CON,MCG(3,IM,IN,IL),R1,
     X   real(F*S1),real(F*S2),real(F*S3),real(F*S4),real(F*S5)
!     X                , ALT(LT+1,IABS(MP)+1)
!     X                , ALN(LN+1,IABS(MV)+1)
     x   ,real(S),real(CG)
45 	 format(' RT,RF,TH =',2f8.3,f8.1,' @',4i3,i2,' : ',
     x          12g13.5,g13.5)
               IF(FFREAL) THEN ! sum over all M values
                RUM(IL) = RUM(IL) +  CG
               ELSE
                SUM(IL) = SUM(IL) +  CG
               ENDIF
 48            CONTINUE
 49            CONTINUE
C
 50    CONTINUE
       CALL SPLINT( DNL/(HNL*NLT) + NLC    ,NLO,IJ,NP,COFFIN)
       DO 80 IL=1,KNL
               IF(FFREAL) THEN
                  NLOC(J)=NLOC(J) + RUM(IL)**2
               ELSE
                  NLOC(J)=NLOC(J) + ABS(SUM(IL))**2
               ENDIF
        DO 70 M=1,NP
            T  =  COFFIN(M) / real(NLT)
            JJ = IJ + M - 1
         IF(FFREAL) THEN
            FNL(I,JJ,IL) = FNL(I,JJ,IL) + RUM(IL) * T
         ELSE
            FNC(I,JJ,IL) = FNC(I,JJ,IL) + SUM(IL) * T
         ENDIF
70      CONTINUE
80    CONTINUE

85    CONTINUE
100   CONTINUE
	deallocate(YLMCF)
      IF(DRY) RETURN
C
      DO 150 IL=1,KNL
      IF(FFREAL) THEN
             WRITE(12) ((FNL(I,J,IL),I=1,NLL),J=1,NLO)
      ELSE
             WRITE(12) ((FNC(I,J,IL),I=1,NLL),J=1,NLO)
      ENDIF
	scalr=0
 	do I=1,NLL
         RT = (I-1) * RINTO
	   PR = mod(I-1,5)==0 .and. RT<2*RHO.and.PPR
         RFMAX = (RHO-B*RT)/A
         RFMIN = (RHO+B*RT)/A
	   if(PR)write(400+IL+NL0,101) RT,RFMIN,RFMAX
101	   format('# RT,RFMIN,RFMAX =',f8.3,' :',2f8.3)
	do J=1,NLO
         DNL = (J - NLM/2 - 1) * HNL  +  CENTRE
         RT0 = RT*HFHT
         RF  = RT0    + DNL
       if(FFREAL) T = abs(FNL(I,J,IL))
	 if(.not.FFREAL) T = abs(FNC(I,J,IL))
	 if(RF>0.and.RF<2*RHO.and.RF<RFMAX+RINTO*2.and.PR) 
     x      write(400+IL+NL0,102) RF,FNL(I,J,IL)+(I-1)*0.1
102	 format(f8.3,2f12.4)
	   if(T>scalr) then
  	     scalr = T
	     Ipos = I;  Jpos = J
	   endif
	 enddo
	 	if(PR) write(400+IL+NL0,*) '&'
	 enddo
	 	if(PPR) write(400+IL+NL0,*) '#   scalr=',real(scalr)

	if(scalr>5.0) then   ! strangely
         RT = (Ipos-1) * RINTO
         DNL = (Jpos - NLM/2 - 1) * HNL  +  CENTRE
         RT0 = RT*HFHT
         RF  = RT0    + DNL
	 WRITE(KO,105) CHNO(IL+NL0,1),CHNO(IL+NL0,2),IL,CP,scalr
     x       ,RT,RF
 105  FORMAT(' NL m-interaction v-c(',I4,') c(',I4,') of IL,CP=',
     x       2I3,' has max ',f13.4,' at RT,RF =',2f8.3)
	endif
C
      IF(NLPLOT.LE.0) GO TO 150
      WRITE(KO,110) CHNO(IL+NL0,1),CHNO(IL+NL0,2),IL,CP
 110  FORMAT(' NL m-interaction v-c(',I4,') c(',I4,') of IL,CP=',
     x       2I3,' is'/)
      IF(FFREAL) CALL DISPLY(FNL(1,1,IL),NLL,NLO,NLL,SCALE)
      IF(.NOT.FFREAL) THEN
         CALL DISPLR(FNC(1,1,IL),NLL,NLO,NLL,SCALR)
           SCALI = SCALR * RERR * 1E3
         CALL DISPLI(FNC(1,1,IL),NLL,NLO,NLL,SCALI)
         SCALE = MAX(SCALR,SCALI)
         ENDIF
      IF(SCALE.lt.1e-10) GO TO 150
      NLPLOT = NLPLOT - 1
C       IF(FFREAL) THEN
C         DO 95 I=2,NLL
C         FNLD(I,1) = 0
C         DO 95 J=1,NLO
C95       FNLD(I,1) = FNLD(I,1) + FNL(I,J,IL) * HNL *MLT
C         WRITE(KO,98) (FNLD(I,1),I=2,NLL)
C98       FORMAT(1X,18F7.2)
C       ENDIF
150   CONTINUE
C
      CALL CCTIME(IPTIME)
       PT = IPTIME * 0.01
      IF(LISTCC.GE.4) WRITE(KO,992) PT,(THM(I)*180/PI,I=2,NLL)
992   FORMAT(/' Theta - maxima (after',F6.2,' secs) are'/(1X,20F6.1))
      RETURN
      END
