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
*****BINCC*************************************************************
      SUBROUTINE BINCC(Y,FORMF,CC,NF,M,K2,IL,CONV,BPHASE,ISC,KMIN,KMAX,
     &                 NK,QNF,ETAP,NMAX,H,PCON,TRES,TDEL,TKMAT,TKNRM,
     &                 LMX1,MAXN,ANC,BSMAT)
	use factorials
	use io
	use drier
	use parameters, only: nchbinmax
	use trace, only: cdcc
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 KMIN,KMAX,K,KP,INTP,K2(M),CC(M,M,NF,3),BPHASE(2,NK)
      INTEGER PCON,QNF(19,M),IOP(M)
      LOGICAL SING,TRES,TKNRM,TDEL,TKMAT,TRA
      COMPLEX*16 F(NMAX+1,M,M),CH,CHD,TMAT,TMATI,YFAC,
     &           ZI(M,M),ZM(M,M),ZP(M,M),
     &           MAT(2*M,3*M),C,T,COUPL(M,M),CI,
     x           BSMAT(NCHBINMAX,NCHBINMAX,NK)
      COMPLEX*16 FORMF(MAXN,NF),Y(NMAX+1,M),W(NMAX+1,M),WY
      REAL*8 ETA,KJ,CF(LMX1),CG(LMX1),CFP(1),CGP(1),KEIG,ETAP(M)
      parameter(nhsp=4)
      real radhsp(nhsp),phres(nhsp)
      data radhsp / 4., 6., 8., 10. /
C
C     pcon    trace iterations     final iteration    list wave function
C      0          no                    no               no
C      1          no                    yes              yes
C      2          no                    yes              no
C      3          yes                   yes              yes
C      4          yes                   yes              no
C
C  If TDEL, multiply wave functions by phase of T before integrating
C  If TRES, then multiply wave functions by abs(T) before integrating.
C  If TKNRM, then multiply wave functions by k before integrating.
C
      NMAXP=NMAX+1
	MAXNR = 2*M
	MAXNR1 = 3*M
      SMALL = 1D0/FPMAX
         CALL CHECK(NMAXP,MAXN,7)
	TRA = PCON>4 .and. TKMAT
      HP=H*H
      R12 = 1./12.
      CI = (0.,1.)
      SQFPI = SQRT(4.*PI)
      RADEG = 180.0/PI
      HP12 = HP*R12
       II = MIN(NK/10,1)
      MATCH = NMAX-10
      DK = (KMAX - KMIN)/(NK-1)
C     INTP = DK/SQRT(KMAX-KMIN) * SQRT(2.0/PI)
      INTP = DK * SQRT(2.0/PI)
      NP = 2*M + 1
       DO 12 L=1,M
         DO 10 N=1,NMAXP
10       Y(N,L) = 0.0
12     CONTINUE
        YINT = 0.0
        IF(PCON.GE.5) then
       WRITE(KO,998)
       write(42,*) DK,NK,KMIN
       write(42,*) H,NMAXP,0
       write(42,*) M,(QNF(9,L),l=1,M)
       endif
C
	ANC = 1e+6
	DELTAL = 0.0
	deltalast = 0.
	deltaiadd = 0.
	ELAST = 0.0
      DO 90 IK=1,NK
         K = KMIN + (IK-1)*DK
         KP= K*K
         NOP = 0
         do 20 J=1,M
         TKJ = K2(IL) + KP-K2(J)
         if(TKJ>0.) then
            NOP = NOP+1
            IOP(NOP) = J
            endif
20	  continue
C
      DO 30 L=1,M
          X = (K2(IL) + KP-K2(L))/CONV
 	  if(Abs(X).lt.0.002) then
 		write(KO,15) KP/CONV,L,X
 15		format(' Skipping energy',f8.4,' as ch',i3,' is ',
     X    ' at energy',f8.4,': too close to threshold!')
    		go to 90
 		endif
      DO 28 J=1,M
      ZI(J,L) = 0.0
28    ZM(J,L) = 0.0
30    ZI(L,L) = H**(QNF(9,L)+1) / EXP(0.5 * FACT(QNF(9,L)+1))
C
      F(1,:,:) = 0d0
      DO 60 I=2,NMAXP
         RRI= 1.0/(I-1.)**2
      DO 42 IT=1,M
      DO 42 J=1,M
42    F(I,J,IT) = ZI(IT,J) * (1. + QNF(9,J)*(QNF(9,J)+1)*RRI*R12 )
       DO 45 J=1,M
       DO 45 L=1,M
         C = 0.0
         DO 425 JF=1,NF
            T = CC(L,J,JF,1)
            IF(ABS(T).LT.SMALL) GO TO 425
         C = C + T * FORMF(I,JF)
425      CONTINUE
         IF(L.EQ.J) C = C + K2(IL) + KP-K2(J)
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
     &                      + F(I,J,IT) * QNF(9,J)*(QNF(9,J)+1) * RRI
      ZM(IT,J) = ZI(IT,J)
55    ZI(IT,J) = ZP(IT,J)
60    CONTINUE
      X = (NMAXP-1) * H
      DO 65 J=1,MAXNR1
      DO 65 L=1,MAXNR
65       MAT(L,J) = 0.0
      DO 75 J=1,M
         DO 70 IT=1,M
            MAT(J,IT) = F(NMAXP,J,IT)
70          MAT(J+M,IT)=F(NMAX ,J,IT)
         TKJ = K2(IL) + KP-K2(J)
         KJ = SQRT(ABS(TKJ))
         L= QNF(9,J)
         XL= L
         ETA = ETAP(J) * 0.5 / KJ
         IF(TKJ.GT.0.0) THEN
            CALL COULFG(KJ*X,ETA,0D0,XL,CF,CG,CFP,CGP,2,0,I,M1)
             CH  = CMPLX(CG(L+1),CF(L+1)) * (0.,.5)
            CALL COULFG(KJ*(X-H),ETA,0D0,XL,CF,CG,CFP,CGP,2,0,I,M1)
             CHD = CMPLX(CG(L+1),CF(L+1)) * (.0,.5)
         ELSE
         IE = 0
           CALL WHIT(ETA,X,KJ,E,L,CF,CG,IE)
           CH = CF(L+1) * (0.,.5)
           CALL WHIT(ETA,X-H,KJ,E,L,CF,CG,IE)
           CHD = CF(L+1) * (0.,.5)
         ENDIF
      IF(PCON.GE.5) WRITE(142,74) J,TKJ/CONV,KJ,X,L,ETA,CH,CHD
     x       			,MAT(J,J),MAT(J+M,J)
74	format(' Bincc #',i3,3f8.3,i3,f8.3,2f10.4,',',2f10.5,1p,4e12.4)
          MAT(J,J+M) = CH
          MAT(J+M,J+M)=CHD
!          IF(J.EQ.IL) THEN
!             MAT(J,NP) = - CONJG(CH)
!             MAT(J+M,NP)=- CONJG(CHD)
!             ENDIF
	   MAT(J,2*M+J) = - CONJG(CH)
	   MAT(J+M,2*M+J)=- CONJG(CHD)
75      CONTINUE
C
       CALL GAUSS5(2*M,MAXNR,MAT,M,SING,T,SMALL,PCON>=5)
C         IF(SING) STOP 'SINGULAR BIN'
          IF(SING) THEN
             WRITE(KO,*) 'SINGULAR ENERGY SOLUTION OMITTED FROM BIN,'
             WRITE(KO,77) ((K2(IL) + KP-K2(J))/CONV,J=1,MIN(M,5))
77           FORMAT(' at channel energies =',5F10.4)
             GO TO 90
           ENDIF
C
C**********************************************************************
C*************Inserted code for core ex********************************
C**********************************************************************
C
      W(:,:) = 0.
      IOUT = IL
      IN = IOP(IOUT)
      NP = 2*M + IN

      DO 810 J=1,M		! so W = scattering wf at 1 energy
      DO 810 IT=1,M
      DO 810 N=1,NMAXP
      W(N,J) = W(N,J) + F(N,J,IT) * MAT(IT,NP)  
810    continue    

	if(cdcc>1.and.M<=NCHBINMAX) 
     x	BSMAT(1:M,1:M,IK) = MAT(M+1:M+M,2*M+1:2*M+M)
!	 if(PCON>=3.and.M>1.and.cdcc>1) then
!	  write(KO,*) 
!	  write(KO,*) ' S matrix at Ein=',KP/CONV,' K=',K
!	    do I=1,M
!	 	write(KO,'(i3,10(2f8.4,:,'',''))') I,BSMAT(I,1:M,IK)
!	    enddo
!	 endif

      TMATI = (MAT(IN+M,NP)-1.)/(0.,2.)
      DELTAI = 0.0
       X = 0.0
           IF(abs(MAT(IN+M,NP)).gt.1e-20) then
             X = (0.,-.5) * LOG(MAT(IN+M,NP))
             DELTAI = RADEG * X
            endif
      WY = DELTAI
	 if(IK>1) then
           if(deltai<deltalast-90) deltaiadd=deltaiadd+180
           if(deltai>deltalast+90) deltaiadd=deltaiadd-180
	 endif
          deltalast=deltai
          deltai = deltai + deltaiadd
	
	BPHASE(1,IK) = DELTAI/RADEG
	BPHASE(2,IK) = K
      IF(PCON.GE.3)  WRITE(KO,999) IK,IN,K,KP/conv,
     X   MAT(il+M,NP), DELTAI,-sqfpi*tmati / conv

      YFAC = 1.0
        IF(TDEL) YFAC = EXP(-CI*BPHASE(1,IK))
        IF(.not.TDEL) BPHASE(1,IK)=0.0  ! no phase for CDCC processing
        IF(TRES) YFAC = CONJG(TMATI)
        IF(TKNRM) YFAC = YFAC*K
      YINT = YINT + ABS(YFAC)**2 * DK
C
      DO 815 J=1,M
      DO 815 N=1,NMAXP
      Y(N,J) = Y(N,J) + W(N,J) * YFAC * INTP
815    continue

       DERIV = 1e-6
       if(IK>1) DERIV = (DELTAI-DELTAL)/(RADEG * (KP/conv - ELAST))
	 DELTAL = DELTAI
	 ELAST =  KP/conv
	 ANC = min(ANC,2d0/DERIV)
      IF(PCON.GE.3)  then
	 WRITE(43,9997) K,WY,KP/conv,-tan(DELTAI/RADEG)/K
         ETA = ETAP(1) * 0.5 / KJ
         if(IK<2) WRITE(44,9998) KP/conv,DELTAI,0d0,ETA
         if(IK>1) WRITE(44,9998) KP/conv,DELTAI,2d0/DERIV,ETA
	 written(43) = .true.
	 written(44) = .true.
	   do N=1,nhsp
       	  L= QNF(9,IL); XL= L;
       	  KJ = sqrt(KP)
       	  ETA = ETAP(IL) * 0.5 / KJ
           CALL COULFG(KJ*radhsp(N),ETA,0D0,XL,CF,CG,CFP,CGP,2,0,I,M1)
             HSP = -atan2(CF(L+1),CG(L+1))*RADEG
          PHRES(N) = DELTAI-HSP
!	    STR = SIN(PHRES(N)/RADEG)
	    enddo
	    write(143,817) KP/conv,(PHRES(N),N=1,nhsp)
	    write(144,817) KP/conv,(sin(PHRES(N)/radeg),N=1,nhsp)
817	    format(f8.3,4(1x,f10.3))	 
	 endif
C
      IF(PCON.GE.5) then			! for information only
      DO 80 J=1,M
         BE = (K2(IL) + KP-K2(J))/CONV
         KJ = SQRT(ABS(BE)*CONV)
       IN = 0
       IF(J.eq.IL) IN =1
      TMAT = (MAT(J+M,NP)-IN)/(0.,2.)
      DELTA = 0.0
      IF(abs(MAT(J+M,NP)).gt.1e-20)
     .DELTA = RADEG * (0.,-.5) * LOG(MAT(J+M,NP))
      IF(DELTA.LT.0.0) DELTA = DELTA + 180.0
      C =-SQFPI * TMAT / CONV
      IF(IK/II*II.EQ.IK) THEN
       if(J.eq.IL) then
          IF(BE.GT.0) WRITE(KO,999) IK,J,K,BE,MAT(J+M,NP),DELTA,C
          IF(BE.LE.0) WRITE(KO,9999) IK,J,K,BE,MAT(J+M,NP)
       else
          IF(BE.GT.0) WRITE(KO,999) IK,J,K,BE,MAT(J+M,NP)
          IF(BE.LE.0) WRITE(KO,9999) IK,J,K,BE,MAT(J+M,NP)
       endif
      ENDIF
 80   continue
      endif

90    CONTINUE
      IF(PCON.ge.0) then
        IF(PCON.ge.3) write(43,*) '&'
        IF(PCON.ge.3) write(44,*) '&'
      endif
      YINT = 1.0/SQRT(YINT)
C
      WOLD = -1.
      IF(MOD(ISC,2).EQ.1) THEN
       X = 0.0
C                        NORMALISE WAVE FUNCTION TO UNITY EXACTLY!!
       DO 92 J=1,M
       DO 92 N=2,NMAX
  92   X = X + abs(Y(N,J))**2
       C = 1.0 / SQRT(X*H)
       WOLD = YINT/ABS(C)
      ELSE
C                   Normalise according to YINT
C                    (if not TRES, then YINT = KMAX - KMIN)
       C = YINT
      ENDIF
       DO 95 N=1,NMAX
       DO 93 J=1,M
93     Y(N,J) = Y(N,J) * C
95     CONTINUE
C
      K = (KMIN+KMAX) * 0.5
      WNORM = 0.0
      RMS   = 0.0
      CD0   = 0.0
      CD1   = 0.0
      DO 100 J=1,M
	if(QNF(9,J)>0) Y(1,J) = 0.0
      DO 100 N=2,NMAX
      X=(N-1)*H
      WNORM = WNORM + abs(Y(N,J))**2
      RMS   = RMS   + abs(Y(N,J))**2 * X*X
         WY = 0.0
         DO 98 L=1,M
         DO 98 JF=1,NF
98       WY = WY + CC(L,J,JF,1) * FORMF(N,JF) * Y(N,L)
      CD0 = CD0 + X**(QNF(9,J)+1) * WY *H
      CD1 = CD1 + X**(QNF(9,J)+1) * Y(N,J)*H
      IF(N.EQ.MATCH.AND.J.EQ.IL) CWN = Y(N,J) * EXP(X*K) / SQFPI
100   CONTINUE
      WNORM = SQRT(WNORM*H)
      RMS   = SQRT(RMS*H)/WNORM
      CD0 =-CD0 * SQFPI      /CONV
      CD1 =-CD1 * SQFPI * K*K/CONV
      D1 = CD1
      D0 = CD0
      Q =CWN * 4*PI / CONV
      IF(WOLD.LT.0.) WOLD=WNORM
      !IF(PCON.GT.0) WRITE(KO,991) WNORM,WOLD,RMS,CD0,CD1
      RMS = WNORM
      Q = D1
      IF(PCON.EQ.0)RETURN
C
  991 FORMAT(  / '0 Norm =',F8.4,
     &  ' from =',F8.4,', RMS =',F8.4,' & D0 =',F12.2,' or',1p,e12.2/)
C
      IF(MOD(PCON,2).EQ.0)      RETURN
      DO 996 J=1,M
      WRITE(KO,992) J,QNF(9,J)
  992 FORMAT(/' For channel #',I3,' with L =',I2//,'   N    Y(R)')
!      if(ISC>0)  WRITE(KO,997) (N,Y(N,J),N=1,NMAXP,5)
!      if(ISC<0) WRITE(KO,998) (N,Y(N,J),N=1,NMAXP,5)
      WRITE(KO,998) (N,Y(N,J),N=1,NMAXP,5)
  996 continue
  997 FORMAT(4(I4,2X,E11.4,3X))
  998 FORMAT(3(I4,2X,2E11.4,3X))
 9997 FORMAT(F11.7,2F9.4,1x,F11.6,f9.4)
 9998 FORMAT(F11.7,F9.4,1x,F11.6,f9.4)
 9999 FORMAT(' For #',i6,i2,', K,BE =',2F11.7,': S =',2F12.2)
  999 FORMAT(' For #',i6,i2,', K,BE =',2F11.7,': S =',2F9.4,:,
     & ' & Del =',F8.3,',  D0 =',2F9.2)
      RETURN
      END
      SUBROUTINE GAUSS5(N,NR,A,M,SING,DET,EPS,SHOW)
C
C    SOLVE BY GAUSSIAN ELIMINATION SUM(J): A(I,J).P(J) = A(I,N+1)
C             WHERE P(J) IS LEFT IN MAT(J,N+1)
C
       IMPLICIT REAL*8(A-H,O-Z)
	parameter(KO=142)
       COMPLEX*16 A(NR,N+M),DET,RA
       LOGICAL SING,SHOW
       SING  = .FALSE.
      NPM = N + M
      DO 201 I=1,N
201   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
      DET = 0.
      DO 9 K = 1, N
      IF (abs(A(K,K)) .NE. 0.0 ) GO TO 5
         SING = .TRUE.
         WRITE(KO,3) K,DET/LOG(10.)
    3  FORMAT(//' THE MATRIX IS SINGULAR AT',I3,', Log10 determinant is
     &  ',2E16.8/)
      DO 401 I=1,N
401   IF(SHOW) WRITE(KO,402) (            A(I,J)  ,J=1,NPM)
C402   FORMAT( 1X,20F6.1/(1X,20F6.3))
402   FORMAT( 1X,1P,14E9.1/(1X,24E9.1))
         RETURN
    5 KP1 = K + 1
      DET = DET + LOG(A(K,K))
         RA = 1.0/A(K,K)
      DO 6 J = KP1, NPM
    6 A(K,J) = A(K,J) * RA
      A(K,K) = 1
      DO 9 I = 1, N
      IF (I .EQ. K  .OR. ABS (A(I,K)) .EQ. 0) GO TO 9
CDIR$ IVDEP
         DO 8 J = KP1, NPM
    8    A(I,J) = A(I,J) - A(I,K)*A(K,J)
         A(I,K) = 0
    9 CONTINUE
      IF(SHOW) WRITE(KO,15 ) DET/LOG(10.)
15    FORMAT(/' Log10 determinant is ',2F10.5)
      IF(SHOW) WRITE(KO,402) (A(I,N+1),I=1,N)
      RETURN
      END
