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
*****INFORM*************************************************************
      SUBROUTINE INFORM(FORML,FORMC,NSP,FORMF,NF,PTYPE,TR,QVAL,BPHASE,
     &   QNF,D0,BEE,NCHAN,NEX,ENEX,JEX,BAND,MASS,AFRAC,CCFRAC,ITC,COPY,
     &   N,NLN,MINT,NNN,RNN,RMI,HCM,RIN,MR,NNU,ERANGE,DK,PI,NAME,
     &   NBINS,EMID,DELTAE,KMINX,NKBIN,NORBIN,BSMAT,RMAS)
	use parameters
	use io
	use drier
	use searchdata
	use searchpar
	use fresco1, only: FATAL,cxwf,sumccbins,
     & 			pluto,npluto,ppots,nextpot
	use trace, only: cdcc
	use parallel, only: iame
      IMPLICIT REAL*8(A-H,O-Z)
C				mp should be defined here the same as in fread
	parameter (mp=200)
      REAL*8 FORMR(MAXM),FORML(MAXNLR,MSP,2),FFR4,KMINX(2,MSP)
      REAL*8 EMID(MSP),DELTAE(MSP),autowid
      COMPLEX*16 FORMF(MAXM,MLOC),TC,FORMC(MAXNLC,MSP,2),FRMC(MAXM),
     x           FFC4,BSMAT(NCHBINMAX,NCHBINMAX,max(1,NKMAX),MSP)
      INTEGER QNF(19,MSP),BAND(2,MXP,MXX),NEX(MXP),TR,TNT(4,mp),
     &        PTYPE(12,MLOC),PARITY,POTK,ITC(MXP,MXX),COPY(2,MXP,MXX,2),
     & 	      NKBIN(2,MSP),NORBIN(MSP)
      REAL*8 D0(MSP,2),JEX(6,MXP,MXX),MASS(4,MXP),WL(20),WLD(20),
     &       J,JN,JNMAX,JNMIN,JCORE,JCOM,KCORE,KCOM,MA,MB,
     &       JNA,JNB,JCOREA,JCOREB,KCOREA,KCOREB,K,QVAL(MXP),RMAS(MSP)
      REAL*8 ENEX(2,MXP,MXX),AFRAC(MXPEX,MXPEX,2,MSP),BEE(MSP,5),
     &       TRITON(8),COEF(mp),BPHASE(2,NKMAX,MSP),FORMIN(2*MAXM),HCM,
     &       FORMII(2*MAXM)
      REAL*8 CCFRAC(MXPEX,2,MXPEX)
      real*8,allocatable:: PSI(:,:),CC(:,:,:,:)
      complex*16,allocatable:: PSIC(:,:)
      CHARACTER*3 VAR(31),ADJ,VARK*4,FFKIND*7,CHME*4,cbin*3
      CHARACTER*8 NAME(2,MXP),TFACTS,TBIN*6,NORMS*4,KFA*2
      CHARACTER CH1,CH2,PSIGN(3)
      data PSIGN / '-','?','+' /
      CHARACTER*140 COMMENT
      LOGICAL FAIL3,FRAC,LAST1,EIGEN,TRES,TDEL,op,keep,cmb,rlb
      LOGICAL TKNRM,PWF,TKMAT,REN(MLOC),FUSED(1000),LIAP,VPOT
     	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN
	namelist/overlap/ kn1,kn2,ic1,ic2,in,kind,ch1,nn,l,lmax,sn,
     &         ia,j,ib,kbpot,krpot,be,isc,ipc,nfl,nam,ampl,keep,
     & 	       dm,nk,er,e,rsmin,rsmax,ppower,nlag,phase,autowid
	namelist/dalitz/ triton
	namelist/twont/ tnt,coef
      DATA VAR /'BE','VR','WR','VSO','USO','VTR','UTR','T2R','DEF',
     &          'DEF','DF0','DF1','DF2','DF3','DF4','DF5','DF6','ALL',
     &          '','','SSC','USR','','','LVD','LDD','LVC',
     &          'LDC', '','','Lsq'/
      DATA TRITON /  .224, 0.50, 1.00, 1.38, 5.5, 0.0, 1.38, 4.2405 /
      FRAC(X) = ABS(X-NINT(X)).GT.1E-5
      FAIL3(X,Y,Z) = FRAC(X+Y+Z) .OR. X.GT.Y+Z .OR. X.LT.ABS(Y-Z)
C
C     the 'QNF' array gives the channel & quantum numbers of
C                the form factors (both 1 & 2 particle bound states)
C         for each KN=KN1,KN2   QNF(i,KN) gives
C
C     QNF(1 : KN1 = no. of a parent form factor for which
C                   fractional parentages etc are defined
C         2 : ICR =  core partition
C         3 : IA  =  core excitation pair (or zero, if not specified)
C         4 : ICP =  compound nucleus partition
C         5 : IB  =  compound nucleus excitation pair (or zero)
C         6 : IN = 1 for projectile state, = 2 for target state
C
C         7 : KIND of bound state (0-3 = 1N  and  6-9 = 2N)
C                  = 0 :  ln,sn; jn                (any IA,IB)
C                  = 1 :  ln, (sn,Jcore)jn; Jcom   (fixed IA,IB)
C                  = 2 :  ln,sn; jn part of sn/K/parity deformed state
C                  = 3 : (ln,sn)jn, Jcore; Jcom    (fixed IA,IB)
C                  = 4 :  (Dalitz-Thacker)
C                  = 5 :
C                  = 6 :  ln, (l,s)sn; jn            & (.5,5)T
C                  = 7 : (ln,l)jn, (s,Jcore)sn; Jcom & (.5.5)T,Tcor;Tcom
C                  = 8 :
C                  = 9 : (ln,(l,s)sn)jn,Jcore; Jcom  & (.5.5)T,Tcor;Tcom
C
C         8 : NN = no. of nodes (incl origin & not infinity, so NN.ge.1)
C         9 : ln = L values relative to core nucleus
C        10 : 2*sn= 2 * total spin of bound cluster  (excl. ln & Jcore)
C                               (but for KIND 7, excluding l too)
C        11 : 2*jn= 2 * (ln + sn as vectors)
C                               (but for KIND 7, jn = ln + l = 'lambda',
C                                and for KIND 1, jn = sn + Jcore = 'S')
C        12 : n  =  0 for single-particle states
C                =  index of nucleon-nucleon separation for 2N states
C        13 : l  =  l value between two nucleons        for 2N states
C        14 : 2*S=  2 * combined NN spin  = 2*(0. or 1.)for 2N states
C                  so (vec) sn = l + S                  for 2N states
C        15 : 2*T=  2 * combined NN isospin=2*(0 or 1)  for 2N states
C                  so l + S + T is odd                  for 2N states
C        16 : Type = character identifier for otherwise-structerless
C                    bound clusters.  A-M for + parity, O-Z for -
C        17 : Number of coupled channels (1, for single channel)
C        18 : IL = Incoming wave for multichannel continuum bins
C        19 : KN1 = no. of parent formfactor for IL=1
C
C        where  Jcore = spin of core nucleus = JEX(IN,ICR,IA)
C               Jcom  = spin of compound nuc = JEX(IN,ICP,IB)
C                  so (vec) Jcom = Jcore + jn
C               Kcore = K-value of core nuc  = JEX(IN+2,ICR,IA)
C               Kcom  = K-value of compound  = JEX(IN+2,ICP,IB)
C                  and      K = Kcom - Kcore
C               Tcore = isospin of core nuc  = JEX(IN+4,ICR,IA)
C               Tcom  = isospin of compound  = JEX(IN+4,ICP,IB)
C                  so (vec) Tcom = Tcore + T
C          and  Parity = sign(1, BAND(IN,ICP,IB) * BAND(IN,ICR,IA))
C
C     the 'BEE' array gives the energies and norms of
C                the form factors (both 1 & 2 particle bound states)
C         for each KN=KN1,KN2   BEE(KN,i) gives
C
C      BEE(KN,1 : BE  = binding energy of state (negative if unbound)
C           ,2 : NORM= root-mean-square norm of wavefunction
C           ,3 : BEE(KN,3) = rms radius of state
C           ,4 : BEE(KN,4) = ANC of state
C           ,5 : ETAP=2*k*ETA      "     "   "    "
C  but BEE(KN,1 : k**2= BE*CONV  (during call to EIGCC)
C
      EPSCON= 1E-4
      PWF = .false. 	! may need in the future
      PWF = .true. 	! may need in the future
      SMALL = 1D0/SQRT(FPMAX)
      MAXC = 40
      Z = 0.0
      IF(sumccbins)THEN
       CCFRAC(1:MXPEX,1:2,1:MXPEX) = 0.0
      ELSE
       AFRAC(1:MXPEX,1:MXPEX,1:2,1:MSP) = 0.0
      ENDIF
C
C
C    single-particle form factors and their parameters
C    -------------------------------------------------
      DO 20 KN=1,MSP
      DO 20 I=1,17
20    QNF(I,KN) = 0
      NSP = 0
      NBINS = 0
      LAST1 = .FALSE.
      FUSED(:) = .false.
	FORMR(:) = 0.
C
	inquire(iolength=iol) FORMR
	inquire(8,opened=op) 
      if(.not.op) then
      IF(MACH==8.or.MACH==3.or.MACH==4) THEN
           OPEN(8,RECL=4*iol,FILE=TMPN(1:lnbl(TMPN))//'08',
     X     ACCESS='DIRECT',FORM='UNFORMATTED')
        ELSE
           OPEN(8,ACCESS='DIRECT',RECL=4*iol,STATUS='SCRATCH',
     X     FORM='UNFORMATTED')   ! 4 rather than 2, to allow for complex bins
        ENDIF
	endif
!	write(6,*) ' FILE 8 opened, recl =',2*iol
	FFKIND='   real'
	if(cxwf) FFKIND='complex'
21    FORMAT(//' ',132('*'),//,
     &' The following ',a7,' SINGLE-PARTICLE FORM FACTORS are ',
     &  ' constructed:'/)
22    FORMAT(/' No.   P1 P2 IN KIND T N  L  S1 IA J/S  IB ',
     &  'BIND XFER  BE    SC PC FIL AFRAC Adjust  to',
     &  '  Z  Mass   K     Norm    rms      D0     D   ANC/Gsp')

      DO 900 KNP=1,10000
!      READ(KI,852) KN1,KN2,IC1,IC2,INI,KIND,CH1,NN,L,LMAX,SN,IAK,J,IB,
!     &         KBPOT,KRPOT,E,ISC,IPC,NFL,NAM,AMPL
	KN1=0;IC1=0;IC2=0;IN=0; BE=0.; IA=0;BE=0.0;IA=0;NK=0
	read(KI,nml=overlap)
852   FORMAT(2I3,4I2,1X,A1,3I2,F4.1,I2,F4.1,I2,2I3,F8.4,4I3,3F8.4)
	IAK=IA
	E = BE
      IF(TR.GE.4)
     &WRITE(KO,852) KN1,KN2,IC1,IC2,IN,KIND,CH1,NN,L,LMAX,SN,IAK,J,IB,
     &         KBPOT,KRPOT,E,ISC,IPC,NFL,NAM,AMPL,rsmin,rsmax
      KN2 = MAX(KN2,KN1)
      KN2I= KN2
      INI = ABS(IN)
      IN = ABS(IN)
      IF(KN1*IC1*IC2*IN.LT.1)GO TO 1000
      CALL CHECK(KN2,MSP,3)
C     if(knp.eq.1)  write(KO,8858) 8,n,'complex*8'
C8858       FORMAT(/' File ',I3,' needs RL =',I3,' ',A9,' numbers')
      IF(KNP.EQ.1) WRITE(KO,21) FFKIND
      DZ = MASS(2+IN,IC1) - MASS(2+IN,IC2)
      DMM = MASS(IN,IC1) - MASS(IN,IC2)
      ICR = IC2
      ICP = IC1
      IF(DMM.LT.0) ICR = IC1
      IF(DMM.LT.0) ICP = IC2
      IF(IC1.GT.NCHAN.OR.IC2.GT.NCHAN.OR.(IN.NE.1.AND.IN.NE.2))THEN
            WRITE(KO,23) KN1,KN2,IC1,IC2,INI
23    FORMAT(/' UNUSABLE COMBINATION OF CHANNEL AND/OR QUANTUM NUMBERS :
     & '   /,' ',2I4,':',3I3,1X,'?'/)
            GO TO 890
            ENDIF
C     IF(IN.EQ.INH) HDN = HP(ICR)
C     IF(IN.NE.INH) HDN = HP(ICP)
	isearch=0
C             STEP SIZES HERE ARE INDEPENDENT OF PARTITION!
      HDN = HCM
C     z1z2 = mass(2+in,icr) * sign(dz,dm)
      ZC = MASS(2+IN,ICR)
      AC = MASS(IN,ICR)
      ETAS = - (AC-2*ZC)/AC  * (DM-2*ZC)/DM
      NKHERE = 0
      ERHERE = ERANGE
	if(nk/=0) NKHERE = nk
	if(abs(er)>1e-9) ERHERE = er
       AMPL = AMPL * SQRT(REAL(max(NAM,0)))
      RM = DM * MASS(IN,ICR) / (DM + MASS(IN,ICR))
      CONV = FMSCAL * RM
      IF(TR.GE.4) WRITE(KO,*) NF,ICR,ICP,DZ,ZC*ABS(DZ),DM,RM,CONV
      IAI = MIN(MAX(IAK,1),NEX(ICR))
!      IB  = MIN(IB,        NEX(ICP))
      IF(IN==2)THEN
       IB  = MIN(IB,        NEX(ICP)+NEX(ICR))
      ELSE
       IB  = MIN(IB,        NEX(ICP))
      ENDIF
         AB = 1.0
         IF(L.LT.0) AB = -1.0
         L = MAX(L,0)
      IAMIN = 1
      IAMAX = NEX(ICR)
      LMAX = MAX(LMAX,L)
      IF(KIND.LE.2 .OR. KIND.EQ.6) THEN
         IAMIN= IAI
         IAMAX= IAI
         ENDIF
	LMX1 = 0
C
      IF(KIND.LE.5) THEN
      IF(KNP.EQ.1.OR..NOT.LAST1) WRITE(KO,22)
C                          one-nucleon bound state(s)
C                          --------------------------
      LMIN = 0
      JNMIN = 0.0
      JNMAX = J
      JCOM = JEX(IN,ICP,MAX(IB,1))
C
      IF(KIND.EQ.0) THEN
C                          (ln,sn) jn  coupling order
         LMIN = L
         LMAX = L
         JNMIN= J
         JNMAX= J
C
      ELSE IF(KIND.EQ.1) THEN
C                             (ln, (sn,Jcore)s; Jcom) coupling order
         IF(IB.LE.1) IB = 1
         TCOM = JEX(IN,ICP,IB)
C+++++++++++++++++++++++++++++++ PATCH MARCH 1993
         JNMIN = J
C
      ELSE IF(KIND.EQ.2) THEN
C                             (ln,sn) jn part of a sn/k/parity state
C                             in a deformed potential.
         IB = MAX(IB,1)
         KCOM = JEX(IN+2,ICP,IB)
         KCORE= JEX(IN+2,ICR,IAI)
         K = KCOM - KCORE
         JNMIN = ABS(K)
         JNMAX = LMAX + SN
C
      ELSE IF(KIND.EQ.3) THEN
C                             (ln,sn)jn, jcore; jcom coupling order
         IB = MAX(IB,1)
         JNMAX = LMAX + SN
      ELSE IF(KIND.EQ.4) THEN
C                              construct form for dalitz-thacker triton
!         IF(NN.GE.2) READ(KI,*) TRITON
         IF(NN.GE.2) READ(KI,nml=dalitz)
         IF(IPC.GE.1) WRITE(KO,44) TRITON
44       FORMAT(' Dalitz-Thacker parameters =',8F9.4)
         BEE(KN1,1) = TRITON(8)  * CONV
         NU = TRITON(2)/HDN + 1.1
	 allocate(PSI(MINT,1))
         DO 45 I=1,NU
45       PSI(I,1) = 0.0
         DO 46 I=NU+1,MINT
            X = (I-1)*HDN
            XR= X - TRITON(2)
46     PSI(I,1) = (EXP(-TRITON(1)*XR) + (TRITON(4)-1)*EXP(-TRITON(3)*XR)
     &             - TRITON(7)*EXP(-TRITON(5)*XR) )  * SQRT(X)
         M = 1
         KN2 = KN1
         VARY = 0.0
         ADJ = 'DZT'
         EIGEN = .TRUE.
          QNF(1,KN1) = KN1
          QNF(2,KN1) = ICR
          QNF(4,KN1) = ICP
          QNF(6,KN1) = IN
          QNF(7,KN1) = KIND
          QNF(8,KN1) = 1
          QNF(9,KN1) = 0
          QNF(17,KN1) = 1
         IF(MOD(IPC,2).EQ.1) WRITE(KO) 826,HDN*MR,(PSI(I,1),I=1,MINT,MR)
826       FORMAT(' DT wavefunction at intervals of',F6.3,' from 0 is',
     &           /(1X,12F10.5))
         GO TO 50
      ENDIF ! KINDS
C
      IF(E.EQ.0.0) E = 10.0
      EIGEN = E.GT.0.0
C                      if eigen then bound state, else continuum bin
         IF(ISC.EQ.0.AND.EIGEN) THEN
            THETA = CONV
            VARY = E
         ELSE
            THETA = 0.0
            VARY = 1.0
         ENDIF
      IF(.NOT.EIGEN) THEN
C                      continuum bin 
	    if(ISC==0)  ISC=2
            I10 = mod(abs(ISC),10)
            TRES  = I10.EQ.3 .OR. I10.EQ.4 .OR. I10.EQ.7 .OR. I10.EQ.8
            TDEL  = I10.EQ.1 .OR. I10.EQ.2
            TKMAT = I10.GE.5 .AND. I10.LE.8
            TKNRM = ISC.ge.10
            TFACTS = ' without'
            IF(TRES) TFACTS = '    WITH'
            IF(TDEL) TFACTS = 'phase of'
            IF(TKMAT)TFACTS = ' (1-iK) '
            IF(TKMAT.and.TRES)TFACTS = 's(1-iK) '
            TBIN =           '  Real'
!           if(M>1) TBIN =   '~~Real'     ! M only known below
            if(cxwf) TBIN = 'Complx'
         ENDIF  ! not EIGEN
      KN = KN1-1
      PARITY = 0
       IK = 0
      DO 30 IA=IAMIN,IAMAX
c      DO 30 IA=IAMAX,IAMIN,-1
         JCORE = JEX(IN,ICR,IA)
         KCORE = JEX(IN+2,ICR,IA)
C        tcore = jex(in+4,icr,ia)
         IF( COPY(IN,ICR,IA,1).NE.0 ) GO TO 30
         IF(MOD(KIND,2).EQ.1)
     &      PARITY = SIGN(1,BAND(IN,ICP,IB)*BAND(IN,ICR,IA))
         IF(IPC.GE.5) WRITE(KO,*) 'TRY IA =',IA,' SO PARITY =',PARITY
      DO 29 LN=LMIN,LMAX
c      DO 29 LN=LMAX,LMIN,-1
         IF(IPC.GE.5) WRITE(KO,*) 'TRY LN =',LN
          IF(MOD(KIND,2).EQ.1 .AND. (-1)**LN.NE.PARITY  .OR.
     &       MOD(KIND,2).EQ.0 .AND. (-1)**(L+LN).NE.1)  GO TO 29
!      DO 28 JN=JNMIN,JNMAX,0.5
      NJN=NINT(JNMAX-JNMIN)*2
      DO 28 IJN=0,NJN
c      DO 28 IJN=NJN,0,-1
      JN=JNMIN+IJN*0.5
C
      IF(KIND.NE.1 .AND. FAIL3(LN+Z,JN,SN)) GO TO 28
         IF(IPC.GE.5) WRITE(KO,*) 'TRY JN =',DBLE(JN)
      IF(KIND.EQ.1 .AND.(FAIL3(SN,JCORE,JN) .OR.
     &                   FAIL3(LN+Z,JN,JCOM))) GO TO 28
      IF(KIND.EQ.3 .AND. FAIL3(JN,JCORE,JCOM)) GO TO 28
      KN = KN + 1
         IF(KN.GT.KN2) GO TO 28
      QNF(1,KN)  = KN1
      QNF(2,KN)  = ICR
      QNF(3,KN)  = IA
      QNF(4,KN)  = ICP
      QNF(5,KN)  = IB
      QNF(6,KN)  = IN
      QNF(7,KN)  = KIND
      QNF(9,KN)  = LN
      QNF(10,KN) = NINT(2.*SN)
      QNF(11,KN) = NINT(2.*JN)
         DO 25 I=12,15
25       QNF(I,KN) = 0
      QNF(16,KN) = ICHAR(CH1)
      IF(IA.EQ.IAI .AND. L.EQ.LN .AND. ABS(JN-J).LT..1) IK=KN
	LMX1 = max(LMX1,LN+1)
      BEE(KN,1) = E * CONV
      IF(KIND.EQ.1 .OR. KIND.EQ.3)
!     &   BEE(KN,1) = (E + ENEX(IN,ICR,IA) - ENEX(IN,ICR,IAI)) * CONV
     &   BEE(KN,1) = (E + ENEX(IN,ICR,IA)) * CONV
      BEE(KN,1) = BEE(KN,1) - THETA * VARY
      BEE(KN,5) = ZC*ABS(DZ) * ETACNS * 2 * RM  * SQRT(FMSCAL) * AB
      IF(sumccbins .AND. IB.GE.1)THEN
       CCFRAC(ITC(ICP,IB),IN,IB) = AMPL
      ELSEIF(IA.GE.1 .AND. IB.GE.1)THEN
       AFRAC(ITC(ICP,IB),ITC(ICR,IA),IN,KN) = AMPL
      ENDIF
      IF(IPC.GE.5) WRITE(KO,27) KN,(QNF(I,KN),I=1,16),BEE(KN,1)/CONV
27    FORMAT(' BOUND PARTIAL WAVE AT',I3,' IS',16I4,', ASYMP BE =',F8.4)
28    CONTINUE
29    CONTINUE
30    CONTINUE
      IF(KN.GT.KN2) THEN
           WRITE(KO,31) KN-KN2
31         FORMAT(/' NO ROOM FOR',I4,' MORE CHANNELS: INCREASE KN2')
           KN = KN2
        ENDIF
      KN2 = KN
      M = KN2 - KN1 + 1
      QNF(17,KN1:KN2) = M
         IF(M.EQ.0) THEN
          WRITE(KO,*) '  NO SUITABLE CHANNELS FOUND ????  INPUT DATA WER
     XE:'
        WRITE(KO,852) KN1,KN2,IC1,IC2,INI,KIND,CH1,NN,L,LMAX,SN,IAK,J,IB
     &        ,KBPOT,KRPOT,E,ISC,IPC,NFL,NAM,AMPL
               GO TO 890
             ENDIF
      IL= IK  - KN1 + 1
        IF(IL.le.0.or.IL.gt.M) IL=1
      QNF(18,KN1:KN2) = IL
      QNF(19,KN1:KN2) = QNF(1,KN1)
      IF(IL>1)THEN
       DO KNA=KN1-1,1,-1   ! search back for parent with IL=1
        LIAP=.TRUE.
        DO IQNF=1,11     ! parent will have identical quantum numbers
         IF(IQNF==3 .OR. IQNF>=9)THEN
          IF(QNF(IQNF,KNA)/=QNF(IQNF,KN1) .OR.
     &       abs(JEX(QNF(6,KNA),QNF(4,KNA),QNF(5,KNA))
     &           -JEX(QNF(6,KN1),QNF(4,KN1),QNF(5,KN1)))>0.001 .OR.
     &       abs(BEE(KNA,1)-BEE(KN1,1)/CONV)>0.001 .or.
     &       QNF(18,KNA)/=1 ) LIAP=.FALSE.
         ENDIF
        ENDDO
        IF(LIAP)THEN
         QNF(19,KN1:KN2) = QNF(1,KNA)
         EXIT
        ENDIF
       ENDDO
      ENDIF
      IF(IPC.GE.5)
     &WRITE(KO,*)' IEXT=',QNF(5,KN1),' ALPHA=',QNF(5,QNF(19,KN1))
C
      IF(NFL.GT.0) THEN
C                       Just read in wave functions from file NFL
        ADJ = 'FIL'
!        VARY = 0.0
	 rlb = EIGEN 
	 cmb = .not.rlb
        GO TO 50
       ENDIF
C                                          couplings
C                                          ---------
	allocate(CC(M,M,NF,4))
      DO 40 INA=1,M
         KNA = KN1 + INA - 1
         LNA = QNF(9,KNA)
         JNA = QNF(11,KNA)*0.5
         JCOREA = JEX(IN,ICR,QNF(3,KNA))
         KCOREA = JEX(IN+2,ICR,QNF(3,KNA))
      DO 40 INB=1,M
         KNB = KN1 + INB - 1
         LNB = QNF(9,KNB)
         JNB = QNF(11,KNB)*0.5
         JCOREB = JEX(IN,ICR,QNF(3,KNB))
         KCOREB = JEX(IN+2,ICR,QNF(3,KNB))
C
      REN(:) = .false.
      DO 35 JF=1,NF
         KP = PTYPE(1,JF)
         POTK = PTYPE(2,JF)
         DO 32 I=1,4
32       CC(INA,INB,JF,I) = 0.0
          IF(PTYPE(4,JF).LT.0) GO TO 35
          IF(PTYPE(3,JF).LT.0.OR.(KP.NE.KBPOT.AND.KP.NE.KRPOT)) GO TO 35
          IF(PTYPE(2,JF).GE.10 .AND.
     X       PTYPE(3,JF).EQ.0  .AND. INA.NE.INB) GO TO 35
      IF(KIND.EQ.0)
     &   T = TENS0 (POTK,LNA,SN,JNA,JCOREA,LNB,SN,JNB,JCOREB,
     &              PTYPE(3,JF),ABS(DZ),ZC)
      IF(KIND.EQ.1)
     &   T = TENSLS(POTK,LNA,SN,JNA,JCOREA,JCOM,LNB,SN,JNB,JCOREB,
     &           PTYPE(3,JF),SN,KCOREA,SN,KCOREB,ABS(DZ),ZC,PTYPE(5,JF))
      IF(KIND.EQ.2)
     &   T = TENDEF(POTK,LNA,SN,JNA,LNB,JNB,PTYPE(3,JF),K,ABS(DZ),ZC)
      IF(KIND.EQ.3)
     &   T = TENSOR(POTK,LNA,SN,JNA,JCOREA,JCOM,LNB,SN,JNB,JCOREB,
     &           PTYPE(3,JF),SN,KCOREA,SN,KCOREB,ABS(DZ),ZC,PTYPE(5,JF),
     &           ETAS)
      IF(ABS(T).LT.SMALL) GO TO 35
      if(POTK/=9) T = - T * CONV
     X      * (-1)**NINT((JCOREA-JCOREB - ABS(JCOREA-JCOREB))/2.)
C           The above phase factors with JCORE etc
C           are there only because of definition of M(Ek) matrix element

         IF(KP.EQ.KBPOT) THEN
	    IF(POTK.eq.9) THEN   ! effective mass correction
                I = 4
            ELSE IF(ISC.EQ.0 .OR. POTK.EQ.0 .OR. .NOT.EIGEN) THEN ! bin potls not adjusted yet
                I = 1
            ELSE IF(abs(ISC).EQ.POTK) THEN
                I = 2
            ELSE IF((abs(ISC).EQ.8.OR.abs(ISC).EQ.9) .AND.
     & 		    PTYPE(4,JF).GT.0) THEN
                I = 2
            ELSE IF(abs(ISC).GE.10.AND.PTYPE(4,JF).GT.0
     &                       .AND.PTYPE(3,JF).EQ.abs(ISC)-10) THEN
                I = 2
            ELSE
                I = 1
            ENDIF
C                               i=1 is the fixed  & i=2 the varied part
C                               i=3 is the reference potl & i=4 is the effective mass correction
         CC(INA,INB,JF,I) = T
	 REN(JF) = REN(JF) .or. I==2
         ENDIF
         IF(KP.EQ.KRPOT)     CC(INA,INB,JF,3) = T
      IF(IPC.GE.5) WRITE(KO,34) KNA,LNA,JNA,JCOREA,KNB,LNB,JNB,JCOREB,
     &      JF,(PTYPE(I,JF),I=1,5),(-CC(INA,INB,JF,I)/CONV,I=1,3)
34    FORMAT(' Between',I4,' (',I3,2F4.1,') &',I4,' (',I3,2F4.1,') ',
     &       ' for potl at',I3,' (',5I3,') --- get',3F10.4)
35       CONTINUE
40     CONTINUE

	 IF(MOD(KIND,2).NE.1)  then  ! provide JCOM and PARITY in uncoupled cases
	   PARITY= (-1)** QNF(9,KN1) !  else keep spatial parity
	   JCOM = QNF(11,KN1)*0.5       ! just l+s=j
	   endif
C-------------------------------------------------------------------
	rlb = EIGEN 
	cmb = .not.rlb
	FFKIND='   real'
	if(cmb) FFKIND='complex'
	
	if(nlag>0) then
	
	 if(NFL==0)  NFL = -33
	  I = ichar('0')
          CHME = char(I+mod(iame/100,10))//char(I+mod(iame/10,10))//
     X       char(I+mod(iame,10))//'.'
         !if(iame==0) CHME
	 nbin=abs(NFL)
 	  cbin = char(I+mod(nbin/100,10))//char(I+mod(nbin/10,10))//
     x         char(I+mod(nbin,10))

	  write(KO,*) 
!	  write(KO,*) ' ************ '
	  write(KO,400) FFKIND,nlag,'fort.'//trim(CHME)//cbin
400	format(' Find ',a7,' overlap with ',i4,' Lagrange mesh basis',
     x 	    ', with file ',a)
	   op = .false.
	  do i=1,npluto
	  op = op .or. pluto(i)==KBPOT
	  enddo
	  if(.not.op) then
	   write(6,*) ' Potential KP=',KBPOT,' not in pluto(:)! Stop'
	   stop
	  endif
	  VPOT = (EIGEN.and.isc/=0) .or.(phase/=0. .and..not.EIGEN)
	  if(VPOT.and.EIGEN) then 
           if(phase/=0.) write(KO,401) phase
401	  format('   Fix potential giving phase shift of',f8.3,' deg')
           if(abs(autowid)>1e-5) write(KO,4011) autowid
4011      format('    and set bin width =',f8.3,'*Gamma')
	  endif
          open(nbin,file='fort.'//trim(CHME)//cbin,form='formatted')

	  open (900,file='pluto.'//CHME//'pin',form='formatted')
	  do i=1,nextpot(KBPOT)-1
	  write(900,'(A)') ppots(i,KBPOT)
	  enddo
402	  continue
	IA = IAMIN ! guess
	 NA = IAMAX-IAMIN+1
	 write(900,403)  'particle', DM,abs(DZ),SN,
     x           NAME(IN,ICR),MASS(IN,ICR),ZC,abs(NFL),
     X           NA,(JEX(IN,ICR,IA),IA=IAMIN,IAMAX) 
	 write(900,4034) (ENEX(IN,ICR,IA),IA=IAMIN,IAMAX)
403	 format('  name1=''',a8,'''  mass1=',f8.4,' z1=',f8.3,
     x           ' spin1=',f5.1,/
     x           '  name2=''',a8,'''  mass2=',f8.4,' z2=',f8.3,/
     x           '  nbin =',i4/
     x           '  ncore=',i1,' icore=',5f9.1)
4034     format( '        ',1x,'encore=',5f9.5)
 	 write(900,404) nlag,0,isc,cmb,iame
404	 format('  nlag=',i3,' njt=',i2,'  plot(:)=0, bin=',i3,
     x          ' nearest=T complexbin=',L1,' bloch=T iame=',i4)

      endif ! nlag>0

      IF(EIGEN) THEN ! bound statess
C                                              **** BOUND STATES
        if(nlag==0)then
	 allocate(PSI(MINT,M))
         ADJ = VAR(abs(ISC)+1)
        IF(IL.le.0.or.IL.gt.M) IL=1
      CALL EIGCC(PSI,FORMF,CC,NF,QNF(1,KN1),BEE(KN1,5),AB,MR,
     &           BEE(KN1,1),THETA,VARY,MAXM,IL,NN,MINT,HDN,M,IPC,
     &           2*M+1,MAXC,EPSCON,IFAIL,LMX1,MINT)
       if(IFAIL>0) then
        if(number_calls>5) then
	 penalty = penalty + fine
	 write(6,41) number_calls,KN1,IFAIL,penalty
41	 format('  At call ',i5,', sp state ',i3,': IFAIL =',i3,:,
     & 		' so  penalty =',1p,e10.1)
	else
	 write(6,41) number_calls,KN1,IFAIL
	 if(FATAL)  stop 'EIGCC FAILURE'
	endif
       endif
       else !nlag>0
	 if(VPOT) then ! vary potential for E
	 isearch=2
	 if(NA>1) isearch=3  ! vary deformed potential, not original WS
	 write(900,4051) isearch,1,NN,-E,IL,JCOM,JCOM,PARITY,PARITY,NN
4051	 format('  search=',i2,' sjt=',i1,' enodes=',i2,
     x     ' eigen=',f10.6,' echan =',i2,
     x     '  jtot=',2f6.1,' parity=',2i3,' pnodes=',i2)

	  else         ! find E in fixed potential  (ignore input E)
 	 write(900,4057) NN,IL,-E,JCOM,PARITY
4057	 format('   pnodes =',i3,' echan =',i2,' eigen=',f10.6,
     x          ' jtot=',f6.1,' parity=',i2,' emaxlist = 0')
	 endif
       endif !nlag

       ELSE ! not eigen
C                                              **** CONTINUUM BINS
	 allocate(PSIC(MINT,M))
         FK = SQRT(FMSCAL * (-E) * RM)
        IF(ERHERE.GT.0) THEN
C                            E RATIO = ERANGE
          IF(ERHERE.LT.1.) ERHERE = 1./ERHERE
          T = ERHERE ** 0.25
          MB = FK * T
          MA = FK / T
        ELSE
C                   E DIFFERENCE = ERANGE
         EMIN = -E - ABS(ERHERE)*0.5
         EMIN = MAX(EMIN,1D-5)
         EMAX = EMIN + ABS(ERHERE)
         MB = SQRT(FMSCAL * EMAX * RM)
         MA = SQRT(FMSCAL * EMIN * RM)
        ENDIF   ! ERHERE
          NK = MAX(10, NKHERE)
	  if(abs(DK)>1e-9) NK = MAX( INT((MB-MA)/DK), NK)
	  NK = min(NK,NKMAX)
	  NBINS = NBINS+1
	  NKBIN(1,NBINS) = KN1
	  NKBIN(2,NBINS) = NK
	  NORBIN(NBINS) = ISC
	  EMID(NBINS) = -E
          DELTAE(NBINS) = EMAX - EMIN
	  KMINX(1,NBINS) = MA
	  KMINX(2,NBINS) = MB
         NORMS = '  no'
         IF(MOD(ISC,2).EQ.1) NORMS = 'unit'
         KFA = '  '
         IF(TKNRM) KFA = '*k'
         T = 1.0 / (FMSCAL * RM)
!        IF(IL.le.0.or.IL.gt.M) IL=1         ! moved earlier
      WRITE(KO,49) TBIN,MA**2 * T, MB**2 * T, NK, TFACTS, KFA,NORMS
49    FORMAT(/' ',A6,' Bin from',F11.6,' MeV to',F11.6,' MeV (cm) in',
     X  I4, ' steps, ',A8,' T-matrix (or D0) factors',
     X A2,' and ',A4,' normalising:')
      if(M.gt.1) WRITE(KO,491) IL
491   FORMAT(' Incoming waves in channel ',I3)
         I=1; if(cdcc>1) I=KN1
         ADJ = 'BIN'
         QNF(8,KN) = -1

      if(nlag==0) then  ! numerov bins

      CALL BINCC(PSIC,FORMF,CC,NF,M,BEE(KN1,1),IL,CONV,BPHASE(1,1,KN1),
     &  ISC,MA,MB,NK,QNF(1,KN1),BEE(KN1,5),MINT-1,HDN,IPC,TRES,TDEL,
     &  TKMAT,TKNRM,LMX1,MAXM,ANC,BSMAT(1,1,1,I))

        else 		!  nlag>0 bins

	 de= max(abs(ER)/NK,1d-10)              ! precision of eminscat printing!
	 ! E = -BEE(KN1,1)
	if(VPOT) then 				! vary potential to phase
	 isearch=2
	 write(900,406) isearch,1,NN,IL,-E,phase,autowid,
     x                  JCOM,JCOM,PARITY,PARITY
406	 format('  search=',i2,' sjt=',i1,' enodes=',i2,' echan =',i2,
     x     ' eigen=',f10.6,' bloch=T phase=',f6.2,' autowid=',f9.3,/,
     x     '  jtot=',2f6.1,' parity=',2i3)
	else					! simply get phases from given potential
	 isearch=0
	 write(900,407) JCOM,PARITY,1,IL
407	format('  jtot=',f6.1,' parity=',i3,' plot(4)=',i1,' echan=',i3)
        endif ! VPOT

	 ex = 0.0
	 KN = KN1+IL-1
         ICR = QNF(2,KN) ; IA  = QNF(3,KN) ; IN  = QNF(6,KN)
         IF(KIND.EQ.1 .OR. KIND.EQ.3) ex = ENEX(IN,ICR,IA)
	write(900,408) -E+ex-abs(er)*0.5,-E+ex+abs(er)*0.5,de
408    format('   eminscat=',f15.10,' emaxscat=',f15.10,' de=',1p,e10.2)

      endif ! nlag

      ENDIF ! EIGEN

      if(nlag>0)then  ! Lagrange mesh, eigen or bin
	 write(900,*) '&'
	 close(900)
       call system('pluto < pluto.'//CHME//'pin > pluto.'//CHME//'pout')
	 Adj='PL'
	 if(VPOT) Adj='PV'
	 if(EIGEN.and.ISC==0) Adj='PE'
	endif ! nlag>0
C-------------------------------------------------------------------
50       WNORM = 0.0
         WRMS = 0.0
         WD0 =0.0
         WD = 0.0
         IF(NFL.NE.0) then
	    call openif(ABS(NFL))
	    FUSED(abs(NFL)) = .true.
	    endif
        I = 0
          if(IB>0) I = SIGN(1,BAND(IN,ICP,IB))
          IF(M.GE.2) WRITE(KO,51) M,KN1,JCOM,PSIGN(I+2)
51     FORMAT(/' The following',I3,' components are in a group labelled'
     &        ,I5,' for ',f5.1,a1)
	   if(cxwf.and..not.allocated(PSIC)) allocate(PSIC(MINT,M))
	   if(.not.cxwf.and..not.allocated(PSI)) allocate(PSI(MINT,M))
      KM = KN1 - 1
      DO 860 KN=KN1,KN2
         IC = KN - KN1 + 1
         L = QNF(9,KN)
         BEE(KN,2) = 0.0
         BEE(KN,3) = 0.0
         BEE(KN,4) = 0.0
      IF(NFL.GT.0.or.nlag>0) THEN
!OLD         READ(NFL,*) (FORMR(I),I=1,MINT)
         READ(abs(NFL),'(a)',end=72)  COMMENT
!	if(nlag==0.or.isearch==0) then
         READ(abs(NFL),*) NPOINTS,RSTEP,RFIRST
!        else
!         READ(abs(NFL),*) NPOINTS,RSTEP,RFIRST,VARY  ! get also potential scaling for eigensolution
!        endif
	 go to 74
72	 write(KO,73) IC,KN
73	 format(' EOF reached when looking for channel #',I3,' @',i4)
	 stop 'EOF NFL'

74       if(IC==1.or.IPC>2)  then
	 if(rlb) VARK='real'
	 if(cmb) VARK='comp'
	 WRITE(KO,75) COMMENT,VARK,NPOINTS,RSTEP,RFIRST
75      FORMAT('  Input:',a140,/'  Reading ',a4,1x,I4,'*2 wf+vertex ',
     X   'points at ',F8.4,' intervals, starting from r =',F8.4)
	 endif
	 if(COMMENT(1:8)=='Pluto EV') then
 	  read(COMMENT,'(8x,f10.5,f9.5)') VARYE,VARYV
	  if(ISC==0) VARY = VARYE
	  if(ISC/=0)  VARY = VARYV
	 endif
	if(NPOINTS>2*MAXM) then
	  write(KO,*) ' Only room to read in ',2*MAXM,' potential',
     x     ' points, not ',NPOINTS
	  stop
	  endif

          RSTEPI = 1./RSTEP
	if(rlb) then
          READ(abs(NFL),*) (FORMIN(I),I=1,NPOINTS)
	IF(IPC>5.and.final)  then
	  do i=1,NPOINTS
          write(258,*) (i-1)*RSTEP+RFIRST,FORMIN(I)
	  enddo
	  write(258,*) '&'
	 endif
          DO I=1,MINT
            R = (I-1)*HDN
            FORMR(I) = FF4(R-RFIRST,RSTEPI,FORMIN,NPOINTS)
	  enddo
	    if(cxwf) FRMC(:) = FORMR(:)
	 else ! complex
          READ(abs(NFL),*) (FORMIN(I),FORMII(I),I=1,NPOINTS)
	IF(IPC>5.and.final)  then
	  do i=1,NPOINTS
          write(258,*) (i-1)*RSTEP+RFIRST,FORMIN(I),FORMII(I)
	  enddo
	  write(258,*) '&'
	 endif
	   if(.not.allocated(PSIC)) allocate(PSIC(MINT,M))
          DO I=1,MINT
            R = (I-1)*HDN
            FRMC(I) = cmplx(FF4(R-RFIRST,RSTEPI,FORMIN,NPOINTS),
     x                      FF4(R-RFIRST,RSTEPI,FORMII,NPOINTS))
	  enddo
	    if(.not.cxwf) FORMR(:) = real(FRMC(:))
	 endif ! rlb or cmb

        DO I=1,MINT-1
         R = (I-1)*HDN
	if(cxwf) PSIC(I,IC) = FRMC(I)*R
	if(.not.cxwf) PSI(I,IC) = FORMR(I)*R
	enddo
	ELSE ! eigcc or bincc

!			copy all channels, not just IC, as need for vertex fn calculation
 	 if(cxwf.and.rlb) PSIC(:,:) = PSI(:,:)
 	 if(.not.cxwf.and.cmb) then
 	    PSI(:,:) = PSIC(:,:)
	    T = maxval(abs(imag(PSIC(:,:))))
	    if(IC==1.and.T>1e-10) write(KO,*) 
     x         ' *** IMAGINARY part of wf dropped !! ***'
     x        ,'    ( max value =',real(T),')'
	   endif

        ENDIF  ! reading

	  if(final) call openif(46)

         BEE(KN,1) = BEE(KN,1)/CONV 
         if(.not.EIGEN) BEE(KN,1) = BEE(KN,1) - ENEX(IN,ICR,IAI)   ! agree with TKJ = K2(IL) + KP-K2(J) in BINCC
         IF(ISC.EQ.0.AND.EIGEN) BEE(KN,1) = BEE(KN,1) + VARY
         FK = SQRT(CONV * ABS(BEE(KN,1)))

        IF(mod(IPC,2).eq.1.and.EIGEN.and.final)  then
          WRITE(46,58) KN,QNF(9,KN),QNF(11,KN)*0.5,QNF(3,KN)
          written(46) = .true.
	 endif
  58   FORMAT('# For Channel #',I3,': l =',I2,', j =',F5.1,', c#',I3)
	  if(PWF) then
	  inquire(55,opened=op)
	  if(.not.op) then
	      call rewop(55); call rewop(59)

        IF(IC==1.and.IPC.GE.1) WRITE(KO,22)

	  endif
	  endif
	if(RSMIN>HDN .or. RSMAX<(MINT-2)*HDN) then ! trim radial wfs

        DO I=1,MINT
         R = (I-1)*HDN
	 if(R<RSMIN .or. R>RSMAX) then
	   if(.not.cxwf) PSI(I,IC) = 0.0
	   if(cxwf)     PSIC(I,IC) = 0.0
	 endif
	 ENDDO
	 if(RSMIN>HDN) WRITE(KO,591) 'below',RSMIN
	 if(RSMAX<(MINT-2)*HDN) WRITE(KO,591) 'above',RSMAX
591	 format(/' *** The next WFN is set to zero ',a5,f8.4,' fm ***')
	endif

        IF(mod(IPC,2)>0.and.final)  then
		call openif(58)
	 if(PWF) WRITE(55,154) KN1,IC,M,QNF(9,KN),QNF(11,KN)*.5,
     x            -real(BEE(KN,1)),cmb,rlb,cxwf
154	 format(4i5,f5.1,f10.5,3l2,' : KN,in-set,total-set,l,j,E')
         if(PWF) WRITE(55,*) MINT,real(HDN),0.,cxwf,' : N,H,R0,cmp'
          WRITE(58,155) KN1,-real(BEE(KN,1)),QNF(9,KN),QNF(11,KN)*0.5,
     X		 QNF(11,KN)+1,MINT,cmb,rlb,cxwf
	if(PWF)
     x    WRITE(59,155) KN1,-real(BEE(KN,1)),QNF(9,KN),QNF(11,KN)*0.5,
     X		 QNF(11,KN)+1,MINT,cmb,rlb,cxwf
155       format('# KN1,E,l,j,2j+1,N,cmb =',i4,f10.5,i4,f5.1,2i4,3l2)
          if(PWF) written(55) = .true.
          written(58) = .true.
          if(PWF) written(59) = .true.
	 endif
              ETA  = 0.5*BEE(KN,5)/FK
        DO 60 I=1,MINT-1
         R = MAX(I-1.,0.01)*HDN
	   if(.not.cxwf) then
!          IF(NFL.LE.0.and.nlag==0) FORMR(I) = PSI(I,IC) / R
            FORMR(I) = PSI(I,IC) / R
	    if(mod(IPC,2)>0.and.final) write(58,156) real(R),PSI(I,IC)
           WWW = FORMR(I)
	   else
!           IF(NFL.LE.0.and.nlag==0) FRMC(I) = PSIC(I,IC) / R
            FRMC(I) = PSIC(I,IC) / R
	    if(mod(IPC,2)>0.and.final) 
     x          write(58,156) real(R),PSIC(I,IC)
!           WWW = abs(FRMC(I))  ! abs only
           WWW = sign(abs(FRMC(I)),real(FRMC(I)))  ! abs only, with sign of real part
	   endif
156		format(f8.3,2g15.5)
           W2 = WWW**2 * R*R * HDN
         BEE(KN,2) = BEE(KN,2) + W2
         BEE(KN,3) = BEE(KN,3) + W2 * R*R
         IF(EIGEN.and.((mod(IPC,2).eq.1.and.I.ge.10)
     x             .or.I.eq.MINT-1)) THEN
             	IE = 0
       	    WWW = WWW*R
             	CALL WHIT(ETA,R,FK,E,L,WL,WLD,IE)
           	ANC = WWW*exp(dble(IE))/WL(L+1)
	       if(final) then
		call openif(46)
             	write(46,59) R,ANC,WWW,WL(L+1)
		written(46) = .true.
59       	format(f9.3,3g15.6)
	      endif
         ENDIF
60    CONTINUE
	    if(mod(IPC,2)>0.and.final) write(58,*) '&'
	if(mod(IPC,2)>0.and.NFL.le.0.and.PWF) then
	   if(rlb)  write(55,64) (PSI(I,IC),I=1,MINT)
	   if(cmb) write(55,64) (PSIC(I,IC),I=1,MINT)
	   call flush(55)
	   endif
        BEE(KN,4) = ANC
         IF(mod(IPC,2).eq.1.and.EIGEN.and.final) then
            write(KO,61) R,KN,ANC
61           format(/'    At R=',F9.3,', ANC in ch',I3,' is',f12.5)
           write(46,*)  '&'
         endif
      IF(.not.cxwf.and.NFL.LE.0) FORMR(1) = 2. * FORMR(2) - FORMR(3)
      IF(     cxwf.and.NFL.LE.0) FRMC(1) = 2. * FRMC(2) - FRMC(3)
C     if(abs(be(kn,2)).lt.1e-20) go to 860
      KM = KM + 1
      IF(KM.NE.KN) THEN
         DO 62 I=1,17
62       QNF(I,KM) = QNF(I,KN)
         DO 63 I=1,4
63       BEE(KM,I) = BEE(KN,I)
         IF(sumccbins)THEN
          CCFRAC(ITC(ICP,IB),IN,QNF(5,KN)) = 0.0
         ELSE
          AFRAC(ITC(ICP,IB),ITC(ICR,QNF(3,KN)),IN,KN) = 0.0
         ENDIF
       ENDIF
       IF(NFL.LT.0.and.nlag==0) THEN
          WRITE(ABS(NFL),635) KM,IC,M,NAME(IN,ICR),NAME(IN,ICP),
     x       qnf(9,km),qnf(10,km)*0.5,qnf(11,km)*.5,
     x       JEX(QNF(6,KM),QNF(2,KM),QNF(3,KM)),JCOM,-BEE(KM,1)
635	  format(' Overlap wf',i3,'#',i3,'/',i3,' of <',a8,'|',a8,'>:',
     x       ' lsj,I;J =',i3,2f5.1,',',f5.1,';',f5.1,' @',f9.5,' MeV')
          WRITE(ABS(NFL),*) MINT,real(HDN),0,cxwf
          if(.not.cxwf) WRITE(ABS(NFL),64) (FORMR(I),I=1,MINT)
          if(cxwf) WRITE(ABS(NFL),64) (FRMC(I),I=1,MINT)
	  written(abs(NFL)) = .true.
         ENDIF
64     FORMAT(1P,6E12.4)
         if(.not.cxwf) WRITE(8,REC=KM) (FORMR(I),I=1,N)
         if(cxwf) WRITE(8,REC=KM) (FRMC(I),I=1,N)
         WNORM = WNORM + BEE(KM,2)
         WRMS  = WRMS  + BEE(KM,3)
         BEE(KM,3)= SQRT(BEE(KM,3)/(BEE(KM,2) + SMALL))
         BEE(KM,2) = SQRT(BEE(KM,2))
         IF(EIGEN) then
  	   if(cxwf) BEE(KM,2) = BEE(KM,2) * SIGN(1D0,real(FRMC(6)))
  	   if(.not.cxwf) BEE(KM,2) = BEE(KM,2) * SIGN(1D0, FORMR(6))
           endif
         QNF(8,KM) = 1
       DO 65 I=1,MINT
       if(.not.cxwf) FORMR(I) = FORMR(I)/(MAX(I-1.,0.01)*HDN)**QNF(9,KM)
       if(cxwf) FRMC(I) = FRMC(I)/(MAX(I-1.,0.01)*HDN)**QNF(9,KM)
       IF(I.GE.3.and..not.cxwf) then
          if(FORMR(I)*FORMR(I-1).LT.0..and.I<MINT)  then
	  	QNF(8,KM) = QNF(8,KM) + 1
		endif
        endif
       IF(I.GE.3.and.cxwf) then
          if(real(FRMC(I))*real(FRMC(I-1)).LT.0..and.I<MINT)  then
	  	QNF(8,KM) = QNF(8,KM) + 1
		endif
        endif
65       CONTINUE

        if(.not.cxwf) FORMR(1) = 2.0 * FORMR(2) - FORMR(3)
        if(cxwf) FRMC(1) = 2.0 * FRMC(2) - FRMC(3)
c
	if(.not.cxwf) then   !  put into FRMC in any case
	 FRMC(:) = FORMR(:)
	endif
       DO 70 I=1,NLN
          R = (I-1)/RIN / HDN
	TC = FFC4(R,FRMC,MINT)
       if(.not.cxwf) FORML(I,KM,1) = TC   ! real part only
70     if(     cxwf) FORMC(I,KM,1) = TC
C
       FORMR(:) = 0.0
       FRMC(:) = 0.0
      IF(NFL.LE.0 .AND. KIND.NE.4.and.nlag==0) THEN  ! get VERTEX functions
!	write(6,*) 'Vertex fn from wfs'
      DO 86 JF=1,NF
      DO 86 IP=1,M  !!		all channels available, in this case !
      IF(KRPOT.GT.0) THEN
         T = - CC(IC,IP,JF,3)/CONV
       ELSE
         T = - (CC(IC,IP,JF,1) + VARY*CC(IC,IP,JF,2)) / CONV
      ENDIF
!      IF(ABS(T).LT.SMALL) GO TO 86
         DO 80 I=2,MINT
         R = (I-1)*HDN
         if(.not.cxwf) FORMR(I) = FORMR(I)+PSI(I,IP)*FORMF(I,JF)*T/R
         if(cxwf) FRMC(I)  = FRMC(I)+PSIC(I,IP)*FORMF(I,JF)*T/R
80       CONTINUE
!	  I = 3.0/HDN+1
!	  write(98,*) 'i,j,jf,F,C,V(3)=',IC,IP,JF,FORMF(I,JF),T,FORMR(I)*3.0
86    CONTINUE
       ELSE IF(NFL.LE.0 .AND. KIND.EQ.4.and.nlag==0) THEN
         DO 87 I=2,MINT-1
         R = (I-1)*HDN
         FORMR(I) =-(PSI(I+1,1) - 2*PSI(I,1) + PSI(I-1,1)) /(HDN*HDN *R)
 87      CONTINUE
       ELSE
C             NFL > 0  .or. nlag>0
	 if(rlb) then
          READ(abs(NFL),*) (FORMIN(I),I=1,NPOINTS)
	IF(IPC>5.and.final)  then
	  do i=1,NPOINTS
          write(259,*) (i-1)*RSTEP+RFIRST,FORMIN(I)
	  enddo
	  write(259,*) '&'
	endif
          DO I=1,MINT
            R = (I-1)*HDN
            FORMR(I) = FF4(R-RFIRST,RSTEPI,FORMIN,NPOINTS)
	  enddo
	   if(cxwf) FRMC(:) = FORMR(:)
         else  ! complex
          READ(abs(NFL),*) (FORMIN(I),FORMII(I),I=1,NPOINTS)
	IF(IPC>5.and.final)  then
	  do i=1,NPOINTS
          write(259,*) (i-1)*RSTEP+RFIRST,FORMIN(I),FORMII(I)
	  enddo
	  write(259,*) '&'
	 endif
          DO I=1,MINT
            R = (I-1)*HDN
            FRMC(I) = cmplx(FF4(R-RFIRST,RSTEPI,FORMIN,NPOINTS),
     x                      FF4(R-RFIRST,RSTEPI,FORMII,NPOINTS))
	  enddo
	   if(.not.cxwf) FORMR(:) = real(FRMC(:))
	 endif
       ENDIF

       IF(NFL.LT.0.and.nlag==0) THEN  
           if(.not.cxwf) then
              FORMR(1) = 2.0 * FORMR(2) - FORMR(3)
	      WRITE(ABS(NFL),64) (FORMR(I),I=1,MINT)
	      endif
           if(cxwf) then
              FRMC(1) = 2.0 * FRMC(2) - FRMC(3)
	      WRITE(ABS(NFL),64) (FRMC(I),I=1,MINT)
	      endif
	   written(abs(NFL)) = .true.
        ENDIF
         D0(KM,1) = 0.0
         D0(KM,2) = 0.0
      DO 88 I=1,MINT
         R = MAX(I-1.,1E-5)*HDN
         R1 = R ** QNF(9,KM)
         if(.not.cxwf) TC = FORMR(I)
         if(cxwf) TC = FRMC(I)
	    if(mod(IPC,2)>0.and.PWF.and.rlb) write(59,156) R,FORMR(I)*R
	    if(mod(IPC,2)>0.and.PWF.and.cmb) write(59,156) R,FRMC(I)*R
         D0(KM,1) = D0(KM,1) + R * TC *         R  * R1
         D0(KM,2) = D0(KM,2) + R * TC * SINH(FK*R) * R1 / FK
         IF(.not.cxwf.and.R1.NE.0.0) FORMR(I) = FORMR(I) / R1
         IF(cxwf.and.R1.NE.0.0) FRMC(I) = FRMC(I) / R1
88	continue
	    if(mod(IPC,2)>0.and.PWF) write(59,*) '&'
c      
      if(QNF(9,KM).GT.0) then
         IF(.not.cxwf) FORMR(1) = 2.0*FORMR(2)-FORMR(3)
         IF(cxwf) FRMC(1) = 2.0 * FRMC(2) - FRMC(3)
	endif
c
         R = HDN * SQRT(4*PI)
            DO 89 I=1,QNF(9,KM)
89          R = R / ((2.*I+1.)*(2.*I))
         D0(KM,1) = D0(KM,1) * R
         D0(KM,2) = D0(KM,2) * R
            WD0 = WD0 + D0(KM,1)
            WD  = WD  + D0(KM,2)
       DO 90 I=1,NLN
          R = (I-1)/RIN / HDN
          if(.not.cxwf) TC = FFR4(R,FORMR,MINT)
          if(cxwf) TC = FFC4(R,FRMC,MINT)
       if(.not.cxwf) FORML(I,KM,2) = TC
       if(     cxwf) FORMC(I,KM,2) = TC
            R = MAX(I-1.,1E-5)/RIN
            R1 = R ** (QNF(9,KM)+1)
      IF(IPC.GE.6) THEN
          if(I==1) WRITE(132,*) 'Vertex function @',KM,QNF(9,KM)
       if(.not.cxwf) write(132,'(f8.3,4g12.4)') R,FORML(I,KM,:)*R1
     x           , FORML(I,KM,2)/FORML(I,KM,1),FORML(I,KM,1)
       if(     cxwf) write(132,'(f8.3,8g12.4)') R,FORMC(I,KM,:)*R1
     x           , FORMC(I,KM,2)/FORMC(I,KM,1),FORMC(I,KM,1)
      endif
90     continue
      IF(IPC.GE.6) THEN
          WRITE(KO,992) KM,QNF(9,KM)
          if(.not.cxwf) WRITE(KO,997) ((I-1)/RIN,
     x         ((I-1)/RIN)**(QNF(9,KM)+1)*FORML(I,KM,2),I=1,NLN)
          if(     cxwf) WRITE(KO,998) ((I-1)/RIN,
     x         ((I-1)/RIN)**(QNF(9,KM)+1)*FORMC(I,KM,2),I=1,NLN)
         ENDIF
! 992 FORMAT(/' V.Psi/R**L @ channel #',I3,' with L =',I2,
!    &       //'  R   V.Psi/R**L')
  992 FORMAT(/' V.Psi @ channel #',I3,' with L =',I2,
     &       //'  R   V.Psi')
  997 FORMAT(4(0p,F6.2 ,1p,E12.4,2X))
  998 FORMAT(3(0p,F6.2 ,1p,2E12.4,2X))
C
	CH2 = ' '
	if(M>1.and.IC==IL) CH2 = '*'
        IF(NFL.GT.0) VARY=0.
	RMAS(KM) = RM
      if(abs(BEE(KM,4))<99) then
      WRITE(KO,858) KM,ICR,ICP,IN,KIND,CH1,CH2,QNF(8,KM),QNF(9,KM),SN,
     &    QNF(3,KM),QNF(11,KM)*.5,QNF(5,KM),KBPOT,KRPOT,
     &    BEE(KM,1),ISC,IPC,NFL,AMPL,  ADJ,VARY,
     &    ABS(DZ),DM,FK,BEE(KM,2),BEE(KM,3),D0(KM,1),D0(KM,2),BEE(KM,4)
	else
      WRITE(KO,8581) KM,ICR,ICP,IN,KIND,CH1,CH2,QNF(8,KM),QNF(9,KM),SN,
     &    QNF(3,KM),QNF(11,KM)*.5,QNF(5,KM),KBPOT,KRPOT,BEE(KM,1),
     &    ISC,IPC,NFL,AMPL,  ADJ,VARY,
     &    ABS(DZ),DM,FK,BEE(KM,2),BEE(KM,3),D0(KM,1),D0(KM,2),BEE(KM,4)
	endif
858   FORMAT(/1x,I4,':',4I3,2X,2A1,i2,I3,F4.1,I3,F4.1,I4,I4,I4,F9.4,2I3,
     &  I3,F7.4,1X,A3,F7.4,';',F3.0,F6.3,F7.4,F7.4,F7.3,F7.1,F7.1,F8.4)
8581  FORMAT(/1x,I4,':',4I3,2X,2A1,i2,I3,F4.1,I3,F4.1,I4,I4,I4,F9.4,2I3,
     &  I3,F7.4,1X,A3,F7.4,';',F3.0,F6.3,F7.4,F7.4,F7.3,F7.1,F7.1,
     &  1p,E10.2)
      NSP = MAX(NSP,KM)
860   CONTINUE
         DO 861 KN=KM+1,KN2
         DO 861 I=1,16
861      QNF(I,KN) = 0
      IF(KN2.GT.KN1) WRITE(KO,859) WNORM,sqrt(WNORM),
     &                 SQRT(WRMS/MAX(WNORM,Z+SMALL)),WD0,WD
859   FORMAT(/' Overall sq,rms norms & radius =',2F8.4,',',f8.4,
     &     ', & Overall D0 & D =',2F9.2)
 	 if(allocated(PSI)) deallocate(PSI)
 	 if(allocated(PSIC)) deallocate(PSIC)

         IF(KM.GT.KN1) WRITE(KO,*) ('-',I=1,132)
      LAST1 = .TRUE.
		do id=1,datasets
		if(data_type(id)==5) then
		 do ip=1,datalen(id)
 		  if(KN1==bs_no(ip,id)) then
		   T = VARY
	 	   EE = (T-datavals(ip,id))/dataerr(ip,id)
                   data_chisq(id) = data_chisq(id) + EE**2
                   theoryvals(ip,id) = T
		   endif
	          enddo
		  endif
	         enddo
      ELSE !if(kind.gt.5)
C
C    two-particle form factors
C    -------------------------
       if(NNN.lt.5) then
        write(6,*) ' Only ',NNN,' points for two-particle transfers!'
        stop 'RNN'
       endif
      IT = KBPOT
      KNZR = KRPOT
      EP2 = E
      LMIN = L
      SMAX = 1.0
      SMIN = SN
      JN = J
      IF(KIND.EQ.9) THEN
         JN = 0.0
         IB = MAX(IB,1)
         ENDIF
      NK = NN
      IF(ABS(EP2).LT.1E-5) EP2 = 1.0
      IF(KNP.EQ.1.OR.LAST1) WRITE(KO,8625)
      WRITE(KO,8621) KN1,KN2,IC1,IC2,INI,KIND,CH1,J,IT,NK,KNZR,LMIN,LMAX
     &       ,SMIN,NFL,ISC,IPC,EP2,AMPL
8621  FORMAT(/1X,6I4,3X,A1,F5.1,',',I1,3I5,I7  ,F5.1,3I4,F6.2,'%',f8.4)
8625  FORMAT( /' The following TWO-PARTICLE FORM FACTORS are constr'
     & ,'ucted :'//'  KN1 KN2  P1  P2  IN KIND TYP J12,T  PAIRS KNZR ',
     &  'LMIN LMAX SMIN FILE ISC IPC  EP2%  AMPL')
      TNT(1,1)=0
      FK= 0.
      IF(NK.GT.0) THEN
!            READ(KI,*) ((TNT(I,JJ),I=1,4),COEF(JJ),JJ=1,NK)
            READ(KI,nml=twont)
8626  FORMAT(3(4I3,F8.4))
	EGS = 0.0
      IF(TNT(1,1).GE.1) THEN
      EGS = BEE(TNT(1,1),1) + BEE(TNT(2,1),1)
         IF(KNZR.NE.0) EGS = EGS - BEE(KNZR,1)
      ENDIF
           WRITE(KO,8627) ((TNT(I,JJ),I=1,4),COEF(JJ),JJ=1,NK)
8627       FORMAT(/5(4I4,F8.4,','))
        ENDIF
	if(cxwf) then
	do 8629 KN=1,MSP
	do 8629 ip=1,2
	do 8629 I=1,NLN
8629	FORML(I,KN,ip) = FORMC(I,KN,ip)   ! to get real single-nucleon states
	endif
	if(IPC>0.and.TNT(1,1)>0) then
	KN = TNT(1,1)
	write(6,*) 'SP state :',KN
	write(6,'(10f9.4)') FORML(1:10,KN,1)
	endif

      CALL TWONN(KN1,KN2,HDN,MR,RNN,ISC,IPC,NLN,NNN,QNF,EGS,
     &       KIND,IN, IAK,IAMIN,IAMAX,IB,JEX,BAND,MXP,MXX,ICR,ICP,COPY,
     &       KNZR,LMIN,LMAX,SMIN,SMAX,JN,J,IT,D0,RIN,DM,MASS(IN,ICR),
     &   NK,TNT,COEF,NFL,cxwf,FORML,MAXNLN,MSP,FORMR,RMI,FMSCAL,
     &   NNN*2,NNU,BEE,FK,N,EP2,NAME(IN,ICP),MAX(NK,1))
      NSP = MAX(NSP,KN2)
      DO 867 KN=KN1,KN2
            IF(QNF(12,KN).EQ.1) WRITE(KO,865)
      WRITE(KO,865) KN,QNF(8,KN),QNF(9,KN),QNF(12,KN),QNF(13,KN),
     &      QNF(14,KN)*.5,QNF(10,KN)*.5,QNF(11,KN)*.5,(BEE(KN,I),I=1,3),
     &      (D0(KN,I),I=1,2)
865   FORMAT(27X,I4,':',' NN LL,(n l s)j; J12 =',2I3,',(',2I2,F4.1,')',
     &  F4.1,';',F4.1,'  BE,Norm =',F8.3,F8.4,',',F6.3,F8.1,F7.1)
      QNF(2,KN) = ICR
      QNF(4,KN) = ICP
      QNF(5,KN) = IB
      QNF(6,KN) = IN
      QNF(7,KN) = KIND
      IA = QNF(3,KN)
      IF(sumccbins)THEN ! shouldn't make it here, 2 particle kinds
       CCFRAC(ITC(ICP,IB),IN,IB) = AMPL ! cc formfactors not implemeted
       stop
      ELSEIF(IA.GE.1 .AND. IB.GE.1)THEN
       AFRAC(ITC(ICP,IB),ITC(ICR,IA),IN,KN) = AMPL
      ENDIF
      QNF(16,KN) = ICHAR(CH1)
      QNF(17,KN) = KN2-KN1+1
	if(cxwf) then
	do 866 ip=1,2
	do 866 I=1,NLN
866	FORMC(I,KN,ip) = FORML(I,KN,ip)   ! to get back two-nucleon states
	endif
867   CONTINUE
      IF(IPC.GT.0) WRITE(KO,871)
871   FORMAT(' ',132('-'))
      LAST1 = .FALSE.
      ENDIF !if(kind<=5 or kind>5)
890	if(allocated(PSI)) deallocate (PSI)
	if(allocated(CC)) deallocate (CC)
	if(EIGEN.and.ISC<0) THEN ! renormalise all varied potentials permanently!
         DO JF=1,NF
	 if(REN(JF)) then
	   FORMF(:,JF) = FORMF(:,JF) * VARY
           ADJ = VAR(abs(ISC)+1)
	   write(KO,895) ADJ,JF,VARY
895	   format('  ### ',a3,' potential at ',i3,
     x		  ' permanently scaled by ',f8.4)
	   endif
	 ENDDO
	 ENDIF

      if(NFL/=0) close(abs(NFL))
900   CONTINUE
C
1000      inquire(55,opened=op)
	  if(op) close(55) 
	  inquire(46,opened=op)
	  if(op) close(46)
	  do i=1,1000
	  if(FUSED(i)) rewind i
	  enddo
      RETURN
      END
*****EIGCC**************************************************************
      SUBROUTINE EIGCC(PSI,FORMF,CCF,NFM,QNF,ETAP,AB,MR,KAP2,THETA,P,
     &  MAXN,MC,NODES,NP,H,M,PCON,MM2,MAXC,EPS,IFAIL,LMX1,MINT)
	use io
	use factorials
	use drier
	! NP here is MINT
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 PSI(MINT,M),F(NP,M,M),CCF(M,M,NFM,4),
     &       KAP2(M),ETAP(M),MAT(MM2,MM2),EMASS(MAXN),
     &       RM(MAXN),DRM(MAXN),D2RM(MAXN)
      COMPLEX*16 FORMF(MAXN,NFM)
      INTEGER QNF(19,M),PCON,COUNT,BC,CCP,CC
      REAL*8 CENT(M),ZI(M,M),ZM(M,M),ZP(M,M)
     &             ,COUT(M),COUTP(M),COUPL(M,M)
      REAL*8 ETA,K,WL(LMX1),WLD(LMX1)
      LOGICAL SING,EM
C
      COUNT = 1
      CCP = 1
      BC = 0
      PP=P
      N = NP-1
      RN = (NP-1)*H
      HP = H*H
      R12= 1./12.
      HP12=HP/12.
      NR = 2*M+1
      IFAIL = 0
      SMALL = 1D0 / FPMAX
      SMALLR = SQRT(SMALL)
      SMALLQ = SQRT(SMALLR)
C     IF(NP.GT.MAXN.OR.NR.GT.MM2) STOP 101
      IF(NP.GT.MAXN.OR.NR.GT.MM2) CALL ABEND(8)
      IF(PCON.GE.3) WRITE(KO,207) NODES,MC
  207 FORMAT(/'    Parameters      Mismatch   Nodes ->',I2,' in ch',i3/)
	EMASS(:) = 1.0
	EM = .false.
      DO JF=1,NFM   ! just use effective mass defined for channel MC
!	write(48,'(a,i3,10f8.4)') 'EMass',JF,CCF(MC,MC,JF,4),
!     x          (DBLE(FORMF(I*10+1,JF)),I=0,5)
       EM = EM .or. abs(CCF(MC,MC,JF,4))>1e-10
       EMASS(:) = EMASS(:) - DBLE(FORMF(:,JF)) * CCF(MC,MC,JF,4) 
      ENDDO
	RM(:) = 1d0/EMASS(:)      ! 1/M(r)
      call DERIV(RM,DRM,H,NP)   ! first derivative of 1/M
      call DERIV(DRM,D2RM,H,NP) ! second derivative of 1/M
      if(EM.and.PCON>4) then
	do I=1,NP
	write(148,'(f8.3,4f10.5)') (I-1)*H,EMASS(I),RM(I),DRM(I),D2RM(I)
	enddo
	write(148,*) '&'
	endif
102   DO 21 J=1,M

      CENT(J) = QNF(9,J) * (QNF(9,J)+1)
      K = SQRT(ABS(KAP2(J) + THETA*P)) * AB
      ETA  = 0.5*ETAP(J)/K
      L = QNF(9,J)
      IE = 0
      CALL WHIT(ETA,RN+H,K,E,L,WL,WLD,IE)
      COUTP(J) = WL(L+1)
      CALL WHIT(ETA,RN,K,E,L,WL,WLD,IE)
21    COUT(J) = WL(L+1)
!      if(COUNT==1.and.abs(THETA)>0) write(149,*) '&',P,THETA,KAP2(MC)
!      if(COUNT==1.and.abs(THETA)>0) write(151,*) '&',P,THETA,KAP2(MC)

      I   =(NP + 1) *3/4
         IMAX = NP/2
         DEL = -SQRT(FPMAX)
      MAM = 0
22    I=I-1
      UDIAG = -KAP2(MC) - THETA*P
      DO 23 JF=1,NFM
23    UDIAG = UDIAG + DBLE(FORMF(I,JF)) * 
     &           ( CCF(MC,MC,JF,1) + P * CCF(MC,MC,JF,2))
!      V = UDIAG - (-KAP2(MC) - THETA*P)
!      T = UDIAG
      UDIAG = EMASS(I) * (UDIAG - DRM(I)**2/RM(I)/4d0 + D2RM(I)/2d0)
!      V2 = EMASS(I) * (V - DRM(I)**2/RM(I)/4d0 + D2RM(I)/2d0)
!      if(COUNT==1.and.abs(THETA)>0)
!     &  write(149,'(f8.3,4f10.5)') (I-1)*H,-(UDIAG+KAP2(MC))/THETA - P,
!     &   -(T+KAP2(MC))/THETA - P,-V/THETA,-V2/THETA
!      if(COUNT==1.and.abs(THETA)>0)
!     &  write(151,'(f8.3,8f12.5)') (I-1)*H,-V2/THETA,
!     &    - EMASS(I) * (V - DRM(I)**2/RM(I)/4d0 + D2RM(I)/2d0)/THETA,
!     &    - EMASS(I) * (V )/THETA,
!     &    - EMASS(I) * (  - DRM(I)**2/RM(I)/4d0 )/THETA,
!     &    - EMASS(I) * ( D2RM(I)/2d0)/THETA,
!     &       - V/THETA
     
      DEN = -CENT(MC)/((I-1)*H)**2 + UDIAG
         IF(DEL.LT.DEN) THEN
            DEL = DEN
            IMAX = I
            ENDIF
      IF(DEN.LT.0 .AND. I.GT.10) GO TO 22
      if(MAM==0) then
        MAM = I
        IF(I.EQ.10) MAM = IMAX
        MAP = MAM+1
      endif
      IF(COUNT.EQ.1 .AND. PCON.GE.7)
     &WRITE(KO,*)MC,QNF(9,MC),K,ETA,RN+H,COUTP(MC),RN,COUT(MC),MAM,MAP
     &             ,CENT(MC)
      if(I>5) go to 22
103   DO 30 J=1,M
      DO 25 IT=1,M
      ZM(IT,J) = 0.0
25    ZI(IT,J) = 0.0
      ZM(J,J) = COUTP(J)
30    ZI(J,J) = COUT(J)
C
C      outer integration,  zi from np to map
      NT = -1
      NF = NP
      NO = NT * (MAP-NP) + 1
40    DO 60 III=1,NO
         I = NF + (III-1)*NT
         II= (I + MIN0(NT,0)) + 1
         RRI= 1.0/(II-1.)**2
      DO 42 IT=1,M
      DO 42 J=1,M
42    F(I,J,IT) = ZI(IT,J) * (1. + CENT(J)*RRI*R12 )
       DO 45 J=1,M
       DO 45 L=1,M
         C = 0.0
         DO 425 JF=1,NFM
            T = CCF(L,J,JF,1)  + P * CCF(L,J,JF,2)
            IF(T.EQ.0.0) GO TO 425
         C = C + T * DBLE(FORMF(II,JF)) 
425      CONTINUE
         IF(L.EQ.J) C = C - KAP2(J) - THETA * P
       if(EM) then
         if(L==J) C = C - DRM(II)**2/RM(II)*.25d0 + .5d0*D2RM(II)
         C = C * EMASS(II)
       endif
       COUPL(L,J) = C 
       IF(C.EQ.0) GO TO 44
       C = C * HP12
          DO 43 IT=1,M
43         F(I,L,IT) = F(I,L,IT) - C * ZI(IT,J)
44    CONTINUE
45    CONTINUE
      DO 49 IT=1,M
         DO 49 L=1,M
49       MAT(IT,L) = 0.0
      DO 54 L=1,M
      DO 54 J=1,M
      C = COUPL(L,J)
      IF(C.EQ.0.0) GO TO 54
      C = C * HP
      DO 53 IT=1,M
53    MAT(IT,L) = MAT(IT,L) + C * F(I,J,IT)
54    CONTINUE
      DO 55 J=1,M
      DO 55 IT=1,M
      ZP(IT,J) = 2*ZI(IT,J) - ZM(IT,J) - MAT(IT,J)
     &                      + F(I,J,IT) * CENT(J) * RRI
      ZM(IT,J) = ZI(IT,J)
55    ZI(IT,J) = ZP(IT,J)
C                              now check for incipient overflows:
      DO 59 IT=1,M
      C = 0.
      DO 57 J=1,M
57    C = MAX(C,ABS(F(I,J,IT)))
      IF(C .LT. FPMAX) GO TO 59
      C = SMALLR
      DO 58 J=1,M
         ZI(IT,J) = ZI(IT,J) * C
         ZM(IT,J) = ZM(IT,J) * C
         DO 58 L=1,III
            I = NF + (L-1)*NT
58       F(I,J,IT) = F(I,J,IT) * C
59    CONTINUE
60    CONTINUE
C
      NT = -NT
      IF(NT.LE.0) GO TO 3
      DO 65 J=1,M
      DO 64 IT=1,M
      ZM(IT,J) = 0.0
64    ZI(IT,J) = 0.0
65    ZI(J,J) = 1E-10 * H**(QNF(9,J)+1) / EXP(0.5 * FACT(QNF(9,J)+1))
C      inner integration,  zi from 1 to mam
      NF = 1
      NO = NT * (MAM-1) + 1
      GO TO 40
3     COUNT = COUNT + 1
C   now calc. derivatives at matching pt. inner(zm) & outer(zp)
      DO 80 J=1,M
         MAT(J,NR) = 0.0
         MAT(J+M,NR) = 0.0
      DO 75 IT=1,M
      DEL =147.0*F(MAM,J,IT)-360.0*F(MAM-1,J,IT)+450.0*F(MAM-2,J,IT)
     1 -400.0*F(MAM-3,J,IT)+225.0*F(MAM-4,J,IT)
     2 -72*F(MAM-5,J,IT)+10.0*F(MAM-6,J,IT)
      ZM(IT,J) = 1  * DEL / (60.0 * H)
      DEL =147.0*F(MAP,J,IT)-360.0*F(MAP+1,J,IT)+450.0*F(MAP+2,J,IT)
     1 -400.0*F(MAP+3,J,IT)+225.0*F(MAP+4,J,IT)
     2 -72*F(MAP+5,J,IT)+10.0*F(MAP+6,J,IT)
      ZP(IT,J) =(-1)* DEL / (60.0 * H)
C     WRITE(KO,*) j,it,zm(it,j)/f(mam,j,it)-zp(it,j)/f(map,j,it)
      MAT(J,IT) = F(MAM,J,IT)
      MAT(J,IT+M) = -F(MAP,J,IT)
      MAT(J+M,IT) = ZM(IT,J)
75    MAT(J+M,IT+M) = - ZP(IT,J)
80    CONTINUE
      DO 78 IT=1,MM2
78    MAT(M+MC,IT) = 0.0
      MAT(M+MC,M+MC) = 1.0
      MAT(M+MC,NR)  = 1.0
c     do 81 j=1,nr
c     WRITE(KO,*) j,(mat(j,it),it=1,nr)
c81        continue
       CALL GAUSSR(2*M,MM2,MAT,SING,DET,SMALL,.false.)
       IF(SING) GO TO 200
       DO 88 L=1,M
          DO 85 I=1,NP
85        PSI(I,L) = 0.0
       DO 88 IT=1,M
       DO 87 I=2,MAM
87     PSI(I,L) = PSI(I,L) + F(I-1,L,IT) * MAT(IT,NR)
       DO 88 I=MAP,NP
88     PSI(I,L) = PSI(I,L) + F(I,L,IT) * MAT(IT+M,NR)
       DERIN = 0.0
       DEROUT= 0.0
       DO 90 IT=1,M
       DERIN = DERIN + ZM(IT,MC) * MAT(IT,NR)
       DEROUT= DEROUT+ ZP(IT,MC) * MAT(IT+M,NR)
90     CONTINUE
       NCO = 1
       DO 91 I=3,NP
C 91     IF(PSI(I,MC)*PSI(I-1,MC).LT.0.) NCO = NCO + 1
       IF(PSI(I,MC).GT.0 .AND. PSI(I-1,MC).LT.0.) NCO = NCO + 1
91     IF(PSI(I,MC).LT.0 .AND. PSI(I-1,MC).GT.0.) NCO = NCO + 1
       DEN = 0.0
       DO  96 J=1,M
       DO  96 L=1,M
         DO 94 JF=1,NFM
         T = CCF(J,L,JF,2)
         IF(ABS(T).LT.SMALL) GO TO 94
            DO 92 IMAX=NP,1,-1
C92          IF(ABS(DBLE(FORMF(IMAX,JF))).GT.SMALL) GO TO 925
92          IF(ABS(DBLE(FORMF(IMAX,JF))).GT.SMALLR .AND.
     X         ABS(PSI(IMAX,L)).GT.SMALLQ) GO TO 925
925         DO  93 I=2,IMAX
93          DEN = DEN + PSI(I,J) * T * DBLE(FORMF(I,JF)) * PSI(I,L) 
     x                   * EMASS(I)
94       CONTINUE
         IF(THETA.EQ.0.0 .OR. J.NE.L) GO TO 96
         DO 95 I=2,NP
95       DEN = DEN - PSI(I,J) * THETA * PSI(I,L)
96    CONTINUE
       DEL =(DEROUT - DERIN) / PSI(MAP,MC)
       IF(PCON.GE.3) WRITE(KO,217) P,DEL,NCO,MAM
C 217    format(e14.6,e15.3,i6)
  217    FORMAT(1X,F13.6,F15.8,2I6)
       POLD = P
       IF(COUNT.GT.MAXC) GO TO 210
       IF(abs(P).lt.1d-4) GO TO 210
       IF(NCO.NE.NODES) GO TO 6
       IF(ABS(DEL).LT.EPS) GO TO 7
         IF(DEN.EQ.0.0) GO TO 220
       Q =-DEL * PSI(MAP,MC)**2 / (DEN * H)
       IF(P.LT.Q) Q=P
       IF(P.LT.-Q) Q = -0.5*P
       P = P+Q
       BC = 1
       Q = 2*Q
99     IF(THETA.NE.0.0) GO TO 102
       if(P>POLD*1.8) GO TO 102  ! recalculate matching point MAM for big potl changes!
       if(P<POLD*0.6) GO TO 102  ! recalculate matching point MAM for big potl changes!
       GO TO 103
6     IF(BC.NE.0) GO TO 61
      CC = 1
      IF(NODES.LT.NCO) CC = -1
      IF(DEN.LT.0.0) CC = -CC
      IF(CC.LT.0) CCP = -IABS(CCP)
      IF(CCP.GT.0) PP = P
      IF(CCP.LT.0) PP = 0.5*PP
      P = PP*CC + P
      GO TO 99
61    Q = 0.5*Q
      P = P - 0.5*Q
      GO TO 99
7     Q = 0.0
      DO 700 J=1,M
      DO 700 I=1,NP
      PSI(I,J) = PSI(I,J)*sqrt(EMASS(I))
700   Q = Q + PSI(I,J)**2
         PIC = 1.0
      Q = SIGN(PIC,PSI(3,MC)) / SQRT(Q * H)
      DO 710 J=1,M
      DO 710 I=1,NP
710   PSI(I,J) = PSI(I,J) * Q
      IF(PCON.EQ.0) RETURN
      IF(MOD(PCON,2).EQ.0)      RETURN
      DO 996 J=1,M
      WRITE(KO,992) J,QNF(9,J),QNF(11,J)*0.5,QNF(3,J)
  992 FORMAT(/' For Channel #',I3,' with l =',I2,',  j =',F5.1,
     &  ', around core state #',I3,             //,'   R    Y(R)')
  996 WRITE(KO,998) ((I-1)*H,PSI(I,J),I=1,NP,MR)
  998 FORMAT(6(1X,0P,F6.2,2X,1P,E11.4,', '))
      RETURN
200   IFAIL = 1
      WRITE(KO,202) DET
202   FORMAT(' ***** FAIL IN EIGCC : DETERMINANT =',1P,2E12.4)
      RETURN
210   IFAIL = 2
      WRITE(KO,212) DEL,MAXC
212   FORMAT(' ***** FAIL IN EIGCC : DISCREPANCY =',1P,E12.4,
     &       ' EVEN AFTER',I3,' ITERATIONS'/)
      RETURN
220   IFAIL = 3
      WRITE(KO,222)
222   FORMAT(' ***** FAIL IN EIGCC : NO VARIABLE PART FOUND !!!!!'/)
      RETURN
      END
      SUBROUTINE GAUSSR(N,NR,A,SING,DET,EPS,SHOW)
C
C    solve by Gaussian elimination sum(j): A(i,j).P(j) = A(i,N+1)
C             where A(j) is left in A(j,N+1)
C
	use io
       IMPLICIT REAL*8(A-H,O-Z)
       PARAMETER(M = 1)
       REAL*8 A(NR,N+M),DET,RA
       LOGICAL SING,SHOW
       SING  = .FALSE.
      NPM = N + M
C      DO 201 I=1,N
C201   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
      DET = 1
      DO 9 K = 1, N
      DET = DET * A(K,K)
      IF (ABS (A(K,K)) .GT. EPS ) GO TO 5
         SING = .TRUE.
         WRITE(KO,3) K,DET
    3    FORMAT(//' THE MATRIX IS SINGULAR AT',I3,'  determinant is ',
     &  E16.8/)
      DO 401 I=1,N
401   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
402   FORMAT( 1X,20F6.3/(1X,20F6.3))
         RETURN
    5 KP1 = K + 1
         RA = 1.0/A(K,K)
      DO 6 J = KP1, NPM
    6 A(K,J) = A(K,J) * RA
      A(K,K) = 1
      DO 9 I = 1, N
      IF (I .EQ. K  .OR. ABS (A(I,K)) .EQ. 0) GO TO 9
Cdir$ ivdep
         DO 8 J = KP1, NPM
    8    A(I,J) = A(I,J) - A(I,K)*A(K,J)
         A(I,K) = 0
    9 CONTINUE
      IF(SHOW) WRITE(KO,15 ) DET
15    FORMAT(/' The determinant is ',E16.8)
      RETURN
      END
C     ================================================================
C     Five point differentiation formula (real functions)
C
      subroutine DERIV(dpaFCT,dpaDF,dpStep,iRadialMax)

      implicit integer(i)
      implicit double precision(d)

      dimension dpaA(5,5),dpaFCT(iRadialMax),dpaDF(iRadialMax)
      data dpaA(1,1),dpaA(1,2),dpaA(1,3),dpaA(1,4),dpaA(1,5)
     &     /-50.d0,96.d0,-72.d0,32.d0,-6.d0/
      data dpaA(2,1),dpaA(2,2),dpaA(2,3),dpaA(2,4),dpaA(2,5)
     &     /-6.d0,-20.d0,36.d0,-12.d0,2.d0/
      data dpaA(3,1),dpaA(3,2),dpaA(3,3),dpaA(3,4),dpaA(3,5)
     &     /2.d0,-16.d0,0.d0,16.d0,-2.d0/
      data dpaA(4,1),dpaA(4,2),dpaA(4,3),dpaA(4,4),dpaA(4,5)
     &     /-2.d0,12.d0,-36.d0,20.d0,6.d0/
      data dpaA(5,1),dpaA(5,2),dpaA(5,3),dpaA(5,4),dpaA(5,5)
     &     /6.d0,-32.d0,72.d0,-96.d0,50.d0/
      data dpEMFact/24.d0/

      do iJ=1,iRadialMax
       iK=3
       if(iJ.lt.3)iK=iJ
       if(iJ.gt.iRadialMax-2)iK=iJ-iRadialMax+5
       dpSUM=0.d0
       do i=1,5
        iJJ=iJ+i-iK
        dpSUM=dpSUM+dpaA(iK,i)*dpaFCT(iJJ)
       end do
       dpaDF(iJ)=dpSUM/(dpStep*dpEMFact)
      end do

      return
      end

C     ================================================================
C     Five point differentiation formula (real functions)
C
      subroutine DERIVC(dpaFCT,dpaDF,Step,iRadialMax)

      implicit integer(i)
      implicit double complex(d)
      real*8 Step

      dimension dpaA(5,5),dpaFCT(iRadialMax),dpaDF(iRadialMax)
      data dpaA(1,1),dpaA(1,2),dpaA(1,3),dpaA(1,4),dpaA(1,5)
     &     /-50.d0,96.d0,-72.d0,32.d0,-6.d0/
      data dpaA(2,1),dpaA(2,2),dpaA(2,3),dpaA(2,4),dpaA(2,5)
     &     /-6.d0,-20.d0,36.d0,-12.d0,2.d0/
      data dpaA(3,1),dpaA(3,2),dpaA(3,3),dpaA(3,4),dpaA(3,5)
     &     /2.d0,-16.d0,0.d0,16.d0,-2.d0/
      data dpaA(4,1),dpaA(4,2),dpaA(4,3),dpaA(4,4),dpaA(4,5)
     &     /-2.d0,12.d0,-36.d0,20.d0,6.d0/
      data dpaA(5,1),dpaA(5,2),dpaA(5,3),dpaA(5,4),dpaA(5,5)
     &     /6.d0,-32.d0,72.d0,-96.d0,50.d0/
      data dpEMFact/24.d0/

      do iJ=1,iRadialMax
       iK=3
       if(iJ.lt.3)iK=iJ
       if(iJ.gt.iRadialMax-2)iK=iJ-iRadialMax+5
       dpSUM=0.d0
       do i=1,5
        iJJ=iJ+i-iK
        dpSUM=dpSUM+dpaA(iK,i)*dpaFCT(iJJ)
       end do
       dpaDF(iJ)=dpSUM/(Step*dpEMFact)
      end do

      return
      end
      
      
!	 N = nint(10./H) + 1
!	write(301,816) KP/conv,W(N,1),YFAC*INTP/DK,DELTAI
!816     format(f10.5,2f10.5,2f8.4,',',f8.4)
