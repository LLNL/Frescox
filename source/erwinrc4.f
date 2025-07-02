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
****ERWINRC4**************************************************************
      SUBROUTINE ERWINRC4(ECM,COEF,FORMF,NF,FORMFR,CH,REPEAT,NBL,
     $  EL,SMAT,L,JVAL,CORESP,LL1,NEQS,N,H,M,MD,SCL,RENORM,MINTL,
     $  AL,RE,ICUTC,SHOW,SHSMAT,SMATEL,CUTOFF,PTYPE,FIM,FIMD,LOCFIL,
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
      COMPLEX*16 FORMF(MAXM,NF),FORMFR(NF),CLIST(MAXCH,MAXCH,MCLIST),C,
     &    S,CZ,SE,CH(2,MAXCH),SMAT(MAXCH),ONEC,ZI2,CI,ZDL,ZHPLUS,ZHMINUS
      REAL*8     FIMD(NEQS,NEQS),FIM(NEQS,NEQS),ZP(NEQS),F8
      REAL*8, allocatable:: COUPL(:,:,:),DIAG(:,:),ZI(:,:,:),ZM(:,:,:),
     &    	FI(:,:,:)
      COMPLEX*16, allocatable:: MAT(:,:)
      REAL*8 COEF(MAXCH),JVAL(MAXCH),CORESP(MAXCH),H2C(MAXCH),MAGN(MAXN)
     &          ,ECM(MAXCH),H(MAXCH),LL1(MAXCH),
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),CRCRAT(MAXCH)
      INTEGER EL,L(MAXCH),CUTOFF,SHOW,
     X        CUTVAL(MAXCH),NFLIST(MAXCH,MAXCH,MCLIST),
     X        NCLIST(MAXCH,MAXCH),PTYPE(12,NF),IPVT(2*NEQS)
      LOGICAL SING,SHSMAT,SMATEL,FCWFN,FJSWTCH,BLAS,REPEAT,LOCFIL
	real*4 TM
      TM(I) = real(SECOND() - TIME0)

       call system_clock(ic,icr,icm)
       if(icr.ne.0) SEC0 = real(ic)/real(icr)
	

	time0=0.0; I=0
	UNC = 1d-100  ! Threshold for negligible couplings
C
        numthread = 1
!$      numthread = OMP_GET_MAX_THREADS()
	if(NEQS>0.and.(SHSMAT.and.SMATEL).and.final)
     X   write(6,*) 'ERWINRC4 solutions =',NEQS,numthread,NBL
     X   	,SHSMAT,SMATEL
      write(48,*) 'ERWINRC4 @',NEQS,numthread,NBL,' @ ',TM(I),N,M
	allocate(DIAG(NEQS,numthread))
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
       if(SHOW.ge.3) then
	write(6,*) 'ERWINRC4: CUTOFF,IS =',CUTOFF,IS
        write(6,*) 'ERWINRC4: CUTVAL =',(CUTVAL(I),I=1,NEQS)
	endif
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
4     FORMAT(' ERWINRC4 step sizes are',10F8.5)
5     FORMAT(' ERWINRC4 given',4I6,L4)
      CALL CHECK(NEQS,MAXB,15)
	if(REPEAT) go to 61
C
	T = NBL*8e-9*NEQS*numthread * 3
	write(48,*) ' Allocate ERWINRC4 arrays in ',real(T),' GB'
	allocate (FI(NBL,NEQS,numthread),ZI(NBL,NEQS,numthread),
     x            ZM(NBL,NEQS,numthread))

	MAGN(:) = 0.0
C
	do K=1,NEQS
	    DIAG(K,:) = 0d0
          H2C(K) = H(K)**2 / COEF(K)

         FI(:,K,:) = 0; ZM(:,K,:)  = 0; ZI(:,K,:)  = 0
	   if(MINTL>1) FIM(:,K) = 0; 
	   FIMD(:,K)= 0
	 do J=1,NEQS
	 do nc=1,NCLIST(J,K)
	   IF(mod(PTYPE(6,NFLIST(J,K,NC)),2)/=0) then
		write(6,*) ' Some couplings erroneously include in CC'
		stop
		endif
	 enddo
	 enddo
	enddo
C
        DO 141 K=1,NEQS
 141     IF(SHOW.GE.9) WRITE(KO,142) K,H(K),COEF(K),(CH(II,K),II=1,2)
 142     FORMAT(' FOR CH.',I4,' H,COEF,CH(1,2)=',10F10.5)
C-----------------------------------------------------------------------

!$OMP  PARALLEL DO  
!$OMP&  PRIVATE (ITHR,I,K,J,IT,I1,I2,II,JF,R,T,F8,RI2,C,NC)

	DO 60 ITHR=1,numthread
	I1 = (ITHR-1)*NBL + 1
	I2 = min(ITHR*NBL,NEQS)
	II = I2-I1+1
		write(48,*) 'Thread ',ITHR,' range: ',I1,I2,'>',II
	  if(LOCFIL) then
	  rewind 19
	  do I=1,IS-1
	   read(19) 
	  enddo
	  endif
	  
         DO 55 I=IS,NM1
          RI2 = 1D0 / DBLE(I-1)**2
	    if(LOCFIL) read(19) FORMFR
       call system_clock(ic,icr,icm)
       if(icr.ne.0) SEC = real(ic)/real(icr)

	write(48,12) ITHR,TM(I),SEC-SEC0,(I-1)*H(1)
12	format(' Thread ',i3,': cpu,elap =',2f9.2,' @ R=',f6.2)
!	write(48,*) ' OpenMP thread ',ITHR,' at time ',TM(I),SEC-SEC0,
!     x           ', R =',real((I-1)*H(1))
	TT0 = second()

	  
          DO 13 K=1,NEQS
             C = 0d0
	     if(LOCFIL) then
             DO NC=1,NCLIST(K,K)
		   C=C+CLIST(K,K,NC) * FORMFR(NFLIST(K,K,NC))
	       ENDDO
	     else
             DO NC=1,NCLIST(K,K)
	 	   C=C+CLIST(K,K,NC) * FORMF(I,NFLIST(K,K,NC))
	       ENDDO
	     endif
!             C = C * (0.,1.)**(L(J)-L(K)) * H2C(K)  ! but diagonal:
             C = C * H2C(K)
 
           DIAG(K,ITHR) = -LL1(K)*RI2 - ECM(K) * H2C(K) + C
	     
           IT = K
              IF(I==CUTVAL(K).and.I1<=IT.and.IT<=I2) then
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
              ZI(IT-I1+1,K,ITHR) = SCL * EXP(T)
!		  IF(SHOW.GE.3) write(KO,1302) I,IT,K,ZI(IT-I1+1,K,ITHR),0.0
		  endif
 13        CONTINUE

		IF(SHOW.GE.3) then
		do 14 K=I1,I2
		IT = K
 14		 if(I==CUTVAL(K)) write(KO,1302) I,IT,K,ZI(IT-I1+1,K,ITHR)
 1302         FORMAT(' AT I=',I3,' ZI(',I5,',',I5,')=',1P,2E12.2)
		endif

         DO 158 K=1,NEQS
	     T = (1d0 - DIAG(K,ITHR) * 
     x                  (R12 - DIAG(K,ITHR)*(ENA2 + DIAG(K,ITHR)*ENA3)))
 158       FI(1:II,K,ITHR) = ZI(1:II,K,ITHR) * T

	   IF(SHOW.GE.7.and.I1==1) THEN
          DO K=1,NEQS
		WRITE(KO,165) I,K,H(K)*(I-1),L(K),
     &    DIAG(K,ITHR)/H2C(K)+ECM(K),0.0,FI(1,K,1),ZI(1,K,1)!*(0.,1.)**(L(EL)-L(K))
!     &               ,FI(K,K,1),ZI(K,K,1) !, COUPL(K,1,ITHR)
          ENDDO
 165      FORMAT(' AT I,K,R =',I6,I4,F7.3,
     &          ' L,PE,PSI =',I4,2F9.3,1P,6E10.1)
          ENDIF
	write(48,*) ' ERWINRC4: 158 done: ',real(SECOND()-TT0),ITHR
	    
			TTT0 = SECOND()
	    DO 24 K=1,NEQS
	     DO 24 J=1,NEQS
            if(K==J.or.NCLIST(K,J)==0) go to 24
	 	C = 0d0
           if(LOCFIL) then
            DO 18 NC=1,NCLIST(K,J)
 18          C=C+CLIST(K,J,NC) * FORMFR(NFLIST(K,J,NC))
           else
            DO 181 NC=1,NCLIST(K,J)
 181          C=C+CLIST(K,J,NC) * FORMF(I,NFLIST(K,J,NC))
           endif
	      if(abs(C)>UNC) then
             C = C * R12 *(0.,1.)**(L(K)-L(J)) * H2C(K)
!		write(49,19) 'ITHR,I,K,J,C = ',ITHR,I,K,J,real(C)
!!     x		,' @ ',real(SECOND()-TT0),real(SECOND()-TTT0)
 19		format(a,4i6,2g12.3,a,f8.4,f10.6)
			TTT0 = SECOND()
 		FI(1:II,K,ITHR) = FI(1:II,K,ITHR) -C*ZI(1:II,J,ITHR)
!	 	call daxpy(II,C,ZI(1,J,ITHR),1,FI(1,K,ITHR),1)
		endif
 24       CONTINUE
	write(48,*) ' ERWINRC4: 24 done: ',real(SECOND()-TT0),ITHR
	    
         DO 43 K=1,NEQS
	   
	    ZP(I1:I2) = 2d0*ZI(1:II,K,ITHR)-ZM(1:II,K,ITHR) 
     x                     - DIAG(K,ITHR)*FI(1:II,K,ITHR)
	    
           DO 35 J=1,NEQS
            if(K==J.or.NCLIST(K,J)==0) go to 35
	 	C = 0d0
           if(LOCFIL) then
            DO 28 NC=1,NCLIST(K,J)
 28          C=C+CLIST(K,J,NC) * FORMFR(NFLIST(K,J,NC))
           else
            DO 29 NC=1,NCLIST(K,J)
 29           C=C+CLIST(K,J,NC) * FORMF(I,NFLIST(K,J,NC))
           endif
	      
	      if(abs(C)>UNC) then
	 	C = C * (0.,1.)**(L(K)-L(J)) * H2C(K)
		ZP(I1:I2) = ZP(I1:I2) - C * FI(1:II,J,ITHR)
		endif
 35	    continue

            ZM(1:II,K,ITHR) = ZI(1:II,K,ITHR)
            ZI(1:II,K,ITHR) = ZP(I1:I2)
 43       CONTINUE
	write(48,*) ' ERWINRC4: 43 done: ',real(SECOND()-TT0),ITHR
	    
       IF(I.EQ.M) THEN
         if(MINTL>1) FIM(I1:I2,:) = FI(1:II,:,ITHR)
!	   exit at this point, so use FI rather than save in FIM	
	write(48,*) 'Found M =',M,MINTL
       ELSE IF(I.EQ.MD) THEN
        FIMD(I1:I2,:) = FI(1:II,:,ITHR)
	write(48,*) 'Found MD =',MD
       ENDIF

            T = 0.0
         DO 47 IT=1,II
            DO 47 K=1,NEQS
            F8 = FI(IT,K,ITHR)
            T = MAX(T, abs(F8))
   47       IF(T .GT. BIG ) GO TO 48
	 if(numthread>1) go to 51 ! skip renorms for now
            MAGN(I) = T
         IF(I.NE.NM1) GO TO 51
 48      IF(T.EQ.0.0) GO TO 51
            MAGN(I) = T
         T = RENORM/T
         IF(SHOW.GE.1) WRITE(KO,49) T,I,RENORM,ITHR
 49      FORMAT(' Renormalising by',E12.4,' at',I4,' to max. WF.',
     x         E12.4,I3)
         if(t.eq.0d0) then
          write(KO,*) 'RENORMALISING TO 0! at I=',I,' as lWF=',MAGN(I)
          write(KO,*) ' FI:',FI
          T = SMALL
          endif
         FMIN = small/T
            ZI(1:II,:,ITHR) = ZI(1:II,:,ITHR) * T
            ZM(1:II,:,ITHR) = ZM(1:II,:,ITHR) * T
            if(MINTL>1)  FIM(I1:I2,:) = FIM(I1:I2,:) * T
            if(MINTL==1)  FI(I1:I2,:,ITHR) = FI(I1:I2,:,ITHR) * T
            FIMD(I1:I2,:) = FIMD(I1:I2,:) * T
C
 51      CONTINUE
	write(48,*) ' ERWINRC4: 51 done: ',real(SECOND()-TT0),ITHR
	 if(I.eq.M) go to 60 ! so use FI not FIM

 55       CONTINUE
 60       CONTINUE
!$OMP END PARALLEL DO
	deallocate (DIAG,ZM,ZI)
	if(MINTL>1) deallocate(FI)
C-----------------------------------------------------------------------
c
 61     write(48,*) 'Allocate MAT'
	allocate(MAT(2*NEQS,2*NEQS+1))
       DO 65 I=1,NR
 65    MAT(I,1:NP) = 0d0
         DO 70 K=1,NEQS
            MAT(K+NEQS,1:NEQS) = FIMD(1:NEQS,K)
	    if(MINTL>1) then
             MAT(K,1:NEQS)   = FIM(1:NEQS,K)
            else
		DO 68 ITHR=1,numthread
		I1 = (ITHR-1)*NBL + 1
		I2 = min(ITHR*NBL,NEQS)
		II = I2-I1+1
 68            	MAT(K,I1:I2)   = FI(1:II,K,ITHR)
	    endif
 70          CONTINUE
	if(MINTL==1) deallocate(FI)
        write(48,*) 'Allocated MAT'
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
	if(.true.) then
            CALL ZGETRF (NR,NR,MAT,NR,IPVT,IER)
            IF ( IER .NE. 0 ) THEN
               WRITE (IOUT,1000) NR,IER
 1000 FORMAT (' ERWINRC4 ',i6,' : ERROR RETURN FROM ZGETRF, IERR = ',I6)
		WRITE(IOUT,1001) (MAT(I,I),I=1,NR)
 1001 FORMAT (1p,10e12.4)
               GO TO 600
            ENDIF
            CALL ZGETRS('N',NR,1,MAT,NR,IPVT,MAT(1,NP),NR,IER)
            IF ( IER .NE. 0 ) THEN
               WRITE (IOUT,1010) IER
 1010 FORMAT (' ERWINRC4 : ERROR RETURN FROM ZGETRS, IERR = ',I6)
               GO TO 600
            ENDIF
	  S = 0d0
	  do K=1,NR
	   T = 1.0
	   if(IPVT(K)/=K) T = -1.0
	   S = S + log(MAT(K,K)*T)
	  enddo
	else

       CALL GAUSS(NR,NR,MAT,SING,S,SMALL,.FALSE.)
          IF(SING) GO TO 600
	endif
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
1157  FORMAT(' ERWINRC4: LOG10(DETERMINANT) =',2F11.3,
     X  ', ACCURACY LOSS =', 1P,E10.2)
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
205   FORMAT('0****** WARNING : ACCURACY LOSS IN ERWINRC4 =',F5.1,
     X ' DIGITS, SO EXPECT ERRORS OF',F9.4,' % OF UNITARITY'/)
      IF(AL*RE.GT..03) WRITE(KO,*) AL,ALUN,RE
	endif
	
      deallocate(MAT)
      RETURN
600   continue
      DO 605 I=IS,NM1
605   IF(MAGN(I).NE.0.0) MAGN(I) = LOG10(MAGN(I))
      WRITE(KO,610) (MAGN(I),I=IS,NM1)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      CALL ABEND(4)
      END
	
****ERWINCC4**************************************************************
      SUBROUTINE ERWINCC4(ECM,COEF,FORMF,NF,FORMFR,CH,REPEAT,NBL,
     $  EL,SMAT,L,JVAL,CORESP,LL1,NEQS,N,H,M,MD,SCL,RENORM,
     $  AL,RE,ICUTC,SHOW,SHSMAT,SMATEL,CUTOFF,PTYPE,FIM,FIMD,LOCFIL,
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
C        SOLVE 'NEQ' COUPLED SCHROEDINGERS EQUATIONS BY EXACT CC, COMPLEX COUPLINGS
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
      COMPLEX*16 FIMD(NEQS,NEQS),FIM(NEQS,NEQS),C,ZP(NEQS),F8
      COMPLEX*16, allocatable:: COUPL(:,:,:),DIAG(:,:),ZI(:,:,:),
     &    	ZM(:,:,:),FI(:,:,:)
      REAL*8 COEF(MAXCH),JVAL(MAXCH),CORESP(MAXCH),H2C(MAXCH),MAGN(MAXN)
     &          ,ECM(MAXCH),H(MAXCH),LL1(MAXCH),
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),CRCRAT(MAXCH)
      INTEGER EL,L(MAXCH),CUTOFF,SHOW,
     X        CUTVAL(MAXCH),NFLIST(MAXCH,MAXCH,MCLIST),
     X        NCLIST(MAXCH,MAXCH),PTYPE(12,NF),IPVT(NEQS*2)
      LOGICAL SING,SHSMAT,SMATEL,FCWFN,FJSWTCH,BLAS,REPEAT,LOCFIL
      AMD1(C) = ABS(C) 
C
        numthread = 1
!$      numthread = OMP_GET_MAX_THREADS()
	call flush(6)
	if(NEQS>0.and.(SHSMAT.and.SMATEL).and.final)
     X   write(48,*) 'ERWINCC4 solutions =',NEQS,numthread,NBL
     X   , NEQS,SHSMAT,SMATEL
	call flush(6)
	call flush(48)
      write(48,*) 'ERWINCC4 solutions =',NEQS,numthread,NBL
	call flush(48)
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
       if(SHOW.ge.3) then
	write(6,*) 'ERWINCC4: CUTOFF,IS =',CUTOFF,IS
        write(6,*) 'ERWINCC4: CUTVAL =',(CUTVAL(I),I=1,NEQS)
	endif
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
4     FORMAT(' ERWINCC4 step sizes are',10F8.5)
5     FORMAT(' ERWINCC4 given',4I6,L4)
      CALL CHECK(NEQS,MAXB,15)
	if(REPEAT) go to 61
C
	allocate (FI(NBL,NEQS,numthread),ZI(NBL,NEQS,numthread),
     x            ZM(NBL,NEQS,numthread))

	MAGN(:) = 0.0
C
	do K=1,NEQS
	    DIAG(K,:) = 0d0
          H2C(K) = H(K)**2 / COEF(K)

         FI(:,K,:) = 0;   ZM(:,K,:)  = 0;   ZI(:,K,:)  = 0;  
	 FIM(:,K) = 0;   FIMD(:,K)= 0; 
	enddo
C
        DO 141 K=1,NEQS
 141     IF(SHOW.GE.9) WRITE(KO,142) K,H(K),COEF(K),(CH(II,K),II=1,2)
 142     FORMAT(' FOR CH.',I4,' H,COEF,CH(1,2)=',10F10.5)
C-----------------------------------------------------------------------

!$OMP  PARALLEL DO  
!$OMP&  PRIVATE (ITHR,I,K,J,IT,I1,I2,II,JF,R,T,F8,RI2,C,NC)

	DO 60 ITHR=1,numthread
	I1 = (ITHR-1)*NBL + 1
	I2 = min(ITHR*NBL,NEQS)
	II = I2-I1+1
		write(48,*) 'Thread ',ITHR,' range: ',I1,I2,'>',II
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
	   IF(mod(PTYPE(6,NFLIST(K,K,NC)),2)==0) 
     x		   C=C+CLIST(K,K,NC) * FORMFR(NFLIST(K,K,NC))
	       ENDDO
	     else
             DO NC=1,NCLIST(K,K)
	   IF(mod(PTYPE(6,NFLIST(K,K,NC)),2)==0) 
     x	 	   C=C+CLIST(K,K,NC) * FORMF(I,NFLIST(K,K,NC))
	       ENDDO
	     endif
           DIAG(K,ITHR) = -LL1(K)*RI2 - ECM(K) * H2C(K) + C*H2C(K)
	     
           IT = K
              IF(I==CUTVAL(K).and.I1<=IT.and.IT<=I2) then
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
              ZI(IT-I1+1,K,ITHR) = SCL * EXP(T)
!		  IF(SHOW.GE.3) write(KO,1302) I,IT,K,ZI(IT-I1+1,K,ITHR),0.0
		  endif
 13        CONTINUE

		IF(SHOW.GE.3) then
		do 14 K=I1,I2
		IT = K
 14		 if(I==CUTVAL(K)) write(KO,1302) I,IT,K,ZI(IT-I1+1,K,ITHR)
 1302         FORMAT(' AT I=',I3,' ZI(',I5,',',I5,')=',1P,2E12.2)
		endif

         DO 20 K=1,NEQS
           COUPL(K,:,ITHR) = 0.0
           IF(I.LT.ICUTC) GO TO 20
           DO 19 J=1,NEQS
            C = 0d0
	    if(K==J) go to 19
	   if(LOCFIL) then
            DO 18 NC=1,NCLIST(K,J)
 18 	   IF(mod(PTYPE(6,NFLIST(K,J,NC)),2)==0) 
     x       C=C+CLIST(K,J,NC) * FORMFR(NFLIST(K,J,NC))
	   else
            DO 181 NC=1,NCLIST(K,J)
 181	   IF(mod(PTYPE(6,NFLIST(K,J,NC)),2)==0) 
     x        C=C+CLIST(K,J,NC) * FORMF(I,NFLIST(K,J,NC))
	   endif
 19         COUPL(K,J,ITHR) = C  * H2C(K)
 20	 CONTINUE
C THE DIAGONAL PART COUPL(K,K,ITHR) IS ZERO

         DO 158 K=1,NEQS
	     S = (1d0 - DIAG(K,ITHR) * 
     x                  (R12 - DIAG(K,ITHR)*(ENA2 + DIAG(K,ITHR)*ENA3)))
 158       FI(1:II,K,ITHR) = ZI(1:II,K,ITHR) * S

	   IF(SHOW.GE.7.and.I1==1) THEN
          DO K=1,NEQS
		WRITE(KO,165) I,K,H(K)*(I-1),L(K),
!    &    DIAG(K,ITHR)/H2C(K)+ECM(K),FI(1,K,1),ZI(1,K,1)
     &    DIAG(K,ITHR)/H2C(K)+ECM(K),FI(1,K,1)
          ENDDO
 165      FORMAT(' AT I,K,R =',I6,I4,F7.3,
!    &          ' L,PE,PSI =',I4,6F9.3,1P,4E10.1)
!    &          ' L,PE,PSI =',I4,1p,10E11.2)
     &          ' L,PE,PSI =',I4,2F9.3,1P,20E10.1)
          ENDIF
	    
	    DO 24 K=1,NEQS
	     DO 24 J=1,NEQS
             C = COUPL(K,J,ITHR) * R12
	      if(C/=0d0) FI(1:II,K,ITHR) = FI(1:II,K,ITHR)
     &                                    -C*ZI(1:II,J,ITHR)
 24       CONTINUE
	    
         DO 43 K=1,NEQS
	   
	    ZP(I1:I2) = 2d0*ZI(1:II,K,ITHR)-ZM(1:II,K,ITHR) 
     x                     - DIAG(K,ITHR)*FI(1:II,K,ITHR)
	    
           DO 34 J=1,NEQS
	      C = COUPL(K,J,ITHR)
 34	      if(C/=0d0) ZP(I1:I2) = ZP(I1:I2) - C * FI(1:II,J,ITHR)

            ZM(1:II,K,ITHR) = ZI(1:II,K,ITHR)
            ZI(1:II,K,ITHR) = ZP(I1:I2)
 43       CONTINUE
	    
       IF(I.EQ.M) THEN
        FIM(I1:I2,:) = FI(1:II,:,ITHR)
       ELSE IF(I.EQ.MD) THEN
        FIMD(I1:I2,:) = FI(1:II,:,ITHR)
       ENDIF

            T = 0.0
         DO 47 IT=1,II
            DO 47 K=1,NEQS
            F8 = FI(IT,K,ITHR)
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
            ZI(1:II,:,ITHR) = ZI(1:II,:,ITHR) * T
            ZM(1:II,:,ITHR) = ZM(1:II,:,ITHR) * T
            FIM(I1:I2,:) = FIM(I1:I2,:) * T
            FIMD(I1:I2,:) = FIMD(I1:I2,:) * T
C
 51      CONTINUE

 55       CONTINUE
 60       CONTINUE
!$OMP END PARALLEL DO
	deallocate (COUPL,DIAG,FI,ZM,ZI)
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
	if(.false.) then
            CALL ZGETRF (NR,NR,MAT,NR,IPVT,IER)
            IF ( IER .NE. 0 ) THEN
               WRITE (IOUT,1000) IER
 1000 FORMAT (' ERWINRC4 : ERROR RETURN FROM ZGETRF, IERR = ',I6)
               STOP
            ENDIF
            CALL ZGETRS('N',NR,1,MAT,NR,IPVT,MAT(1,NP),NR,IER)
            IF ( IER .NE. 0 ) THEN
               WRITE (IOUT,1010) IER
 1010 FORMAT (' ERWINRC4 : ERROR RETURN FROM ZGETRS, IERR = ',I6)
               STOP
            ENDIF
	  S = 0d0
	  do K=1,NR
	   T = 1.0
	   if(IPVT(K)/=K) T = -1.0
	   S = S + log(MAT(K,K)*T)
	  enddo
	else
       CALL GAUSS(NR,NR,MAT,SING,S,SMALL,.FALSE.)
	endif
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
	    if(ECM(K)>0.) then
             AL = MAX(AL, H2C(K))
             ALUN = MAX(ALUN, H2C(K)*ABS(SMAT(K)))
	    endif
!	     write(106,*) ALUN,K,WVD,H2C(k),smat(K),AL
      ENDIF
150      CONTINUE
      IF(SHOW.ge.1) WRITE(KO,1157) S,AL
1157  FORMAT(' ERWINCC4: LOG10(DETERMINANT) =',2F11.3,
     X  ', ACCURACY LOSS =', 1P,E10.2)
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
205   FORMAT('0****** WARNING : ACCURACY LOSS IN ERWINCC4 =',F5.1,
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
	
