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
****ERWIMR1**************************************************************
! Use BCAST for getting support data to the helpers
! Use scalapack for simeqs for boundary conditions.

      SUBROUTINE ERWIMR1(ECM,COEF,FORMF,NF,CH,REPEAT,NBL,
     $  EL,SMAT,L,JVAL,CORESP,LL1,NEQS,N,H,M,MD,SCL,RENORM,MINTL,
     $  AL,RE,IC,SHOW,SHSMAT,SMATEL,CUTOFF,PTYPE,FIM,FIMD,LOCFIL,
     $  CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIS,NCLIST,NFLIS,
     X  PFLIS,MFLI)
	use io
	use factorials
	use parameters
	use drier
	use searchpar, only: final
	use fresco1, only: TIME0
	use parallel
!$      use omp_lib
      IMPLICIT REAL*8(A-H,O-Z)
C
C        SOLVE 'NEQ' COUPLED SCHROEDINGERS EQUATIONS BY EXACT CC, REAL COUPLINGS
C        SOLVE  
C             (COEF(K). D2/DR2 + EN(R,K)).W(R,K)
C                  + SUM(J): COUPL(R,K,J).W(R,J) + INHOMG(R,K) = 0
C           WHERE COUPL(R,K,J) = SUM(NC from 1 to NCLIST(K,J)): 
C                                FORMF(R,JF)*CLIS(PFLIS(K,J)+NC)
C           WHERE    JF = NFLIS(PFLIS(K,J)+NC)
C
C     ASSUMED FOR ARRAYS HERE THAT N<=MAXN AND NEQS <= MAXCH
C
C  Using 'Enhanced Numerov' of Thorlacius & Cooper (JCP 72(1987) 70)
C   with 5 terms in cosh(sqrt(12T)) expansion, but only diagonal potl.
C
      COMPLEX*16 FORMF(MAXM,NF),CLIS(MFLI),C,
     &    S,CZ,SE,CH(2,MAXCH),SMAT(MAXCH),ONEC,ZI2,CI,ZDL,ZHPLUS,ZHMINUS
!      REAL*8     FIMD(NEQS,NEQS),FIM(NEQS,NEQS)
	REAL*8 FIMD(NBL,NEQS,numthread),FIM(NBL,NEQS,numthread)
	REAL*8, allocatable:: FI(:,:,:),DIAG(:,:),ZI(:,:,:),ZM(:,:,:)   	
      COMPLEX*16, allocatable:: MAT(:,:),MATL(:,:),RHSL(:)
      integer DESCM(9),DESCR(9)
      REAL*8 COEF(MAXCH),JVAL(MAXCH),CORESP(MAXCH),H2C(MAXCH),MAGN(MAXN)
     &          ,ECM(MAXCH),H(MAXCH),LL1(MAXCH),ZP(NEQS),F8,
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),CRCRAT(MAXCH)
      INTEGER EL,L(MAXCH),CUTOFF,SHOW,ihelper,MODEQ,
     X        CUTVAL(MAXCH),NFLIS(MFLI),PFLIS(MAXCH,MAXCH),
     X        NCLIST(MAXCH,MAXCH),PTYPE(12,NF),IPVT(2*NEQS),iv(12)
      INTEGER IAM,NPROCS,ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER NBROW,NBCOL,LDA,LDB,ICN
      LOGICAL SING,SHSMAT,SMATEL,FCWFN,FJSWTCH,BLAS,REPEAT,LOCFIL,OTHERS
	real*4 TM
      TM(I) = real(SECOND() - TIME0)
	I=0
	IOUT = KO
C
	if(NEQS>0.and.(SHSMAT.and.SMATEL))
     X   write(6,*) 'ERWIMR1-v4 solutions =',NEQS,numthread,NBL
     X   , NEQS,SHSMAT,SMATEL
      write(48,*) 'ERWIMR1 @',NEQS,numthread,NBL,' @ ',TM(I),N,M
	allocate(DIAG(NEQS,numthread))
!      call flush(48)
      NM1 = N-1
      NR = 2*NEQS
      NP = NR + 1

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
	write(6,*) 'ERWIMR1: CUTOFF,IS =',CUTOFF,IS
        write(6,*) 'ERWIMR1: CUTVAL =',(CUTVAL(I),I=1,NEQS)
	endif
      NN = NM1 - IS + 1

        SMALL = 1.0/FPMAX
        EPS = SQRT(SMALL)
        BIG = 1./EPS
C
      IF (FCWFN) THEN
C                           print CRCWFN
           IF(SHOW.ge.4) write(6,*) ' CG2,CF2(M) CG1,CF1(MD) ='
       DO 720 K=1,NEQS
       DO 720 K2=1,NEQS
           IF(SHOW.ge.4.and.ABS(CGMAT(K2,K,2)).gt.0.0) 
     *       WRITE(6,716) K,K2,CGMAT(K2,K,2),CFMAT(K2,K,2),
     *           CGMAT(K2,K,1),CFMAT(K2,K,1)
716            FORMAT(1X,2I3,': CRCWFN= ',4(D15.7,3X))
720	 CONTINUE
	ENDIF

      IF (.not. FJSWTCH)  THEN  ! solve radial equations

C
      IF(SHOW.GE.2) WRITE(KO, 5) NF,M,NEQS,NR,REPEAT
      IF(SHOW.GE.2) WRITE(KO, 4) (H(K),K=1,NEQS)
4     FORMAT(' ERWIMR1=v4 step sizes are',10F8.5)
5     FORMAT(' ERWIMR1-v4 given',4I6,L4)
      CALL CHECK(NEQS,MAXB,15)
	if(REPEAT) go to 61
C
	T = NBL*8e-9*NEQS*numthread * 3
	write(500+iame,*) ' Allocate ERWIMR1 arr',real(T),' GB @',TM(I)
	allocate (FI(NBL,NEQS,numthread))

	MAGN(:) = 0.0
C
	do K=1,NEQS
          H2C(K) = H(K)**2 / COEF(K)

         FI(:,K,:) = 0; 
         if(MINTL>1) FIM(:,K,:) = 0; 
	   FIMD(:,K,:)= 0
	 do J=1,NEQS
	 do nc=1,NCLIST(J,K)
		 JF = NFLIS(PFLIS(J,K)+NC)
		IF(mod(PTYPE(6,JF),2)/=0) then
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
	write(500+iame,*) 'begin',mpihelp
      LAST = mpihelp
#ifdef MPI
	DO 40 ITHR=2,mpihelp    ! call helpers
	 ihelper = ithr-1
	 write(500+iame,*) ' Node ',iame,' to send ITHR=',ITHR,
     x                   ' to ',ihelper,' (node =',ihelper+iame,')'
      call MPI_send(ITHR,1,MPI_INTEGER,ihelper,1,
     x                                 commgroup,ierr)
40	continue
! 					! broadcast supporting data
      call MPI_BCAST(NEQS,1,MPI_INTEGER,0,commgroup,ierr)
     
	 iv(1)=KO; iv(2)=NBL; iv(3)=IS; iv(4)=NM1; iv(5)=NF
	 iv(6)=SHOW; iv(7)=MINTL; iv(8)=MAXCH; iv(9)=MAXM
	 iv(10)=M; iv(11)=MD; iv(12)=MFLI
	 
	call MPI_BCAST(iv,12,MPI_INTEGER,0,commgroup,ierr)      

      call MPI_BCAST(LOCFIL,1,MPI_LOGICAL,0,commgroup,ierr)
      call MPI_BCAST(SCL,1,MPI_DOUBLE_PRECISION,0,commgroup,ierr)
      call MPI_BCAST(FORMF,MAXM*NF,MPI_DOUBLE_COMPLEX,0,commgroup,ierr)
      call MPI_BCAST(CLIS,MFLI,MPI_DOUBLE_COMPLEX,0,commgroup,ierr)
      call MPI_BCAST(NFLIS,MFLI,MPI_INTEGER,0,commgroup,ierr)
      call MPI_BCAST(H2C,NEQS,MPI_DOUBLE_PRECISION,0,commgroup,ierr)
      call MPI_BCAST(ECM,NEQS,MPI_DOUBLE_PRECISION,0,commgroup,ierr)
      call MPI_BCAST(H,NEQS,MPI_DOUBLE_PRECISION,0,commgroup,ierr)
      call MPI_BCAST(L,NEQS,MPI_INTEGER,0,commgroup,ierr)
      call MPI_BCAST(CUTVAL,NEQS,MPI_INTEGER,0,commgroup,ierr)
      call MPI_BCAST(NCLIST,MAXCH**2,MPI_INTEGER,0,commgroup,ierr)      
      call MPI_BCAST(PFLIS,MAXCH**2,MPI_INTEGER,0,commgroup,ierr)
! 40       CONTINUE
 	LAST = 1
#endif /* MPI */
 
 !  Do my own work in the meantime ! (or all, if no MPI!)
 
  	  DO ITHR = 1,LAST
  	  write(500+iame,*) ' Master does ',ITHR,', MINTL =',MINTL
 	CALL ERWICMR1(ITHR,NBL,NEQS,LOCFIL,IS,NM1,M,MD,NF,SHOW,MINTL,
     x    H2C,FORMF,CLIS,NFLIS,NCLIST,PFLIS,ECM,H,L,SCL,CUTVAL,MAGN,
     X    FI(1,1,ITHR),FIM(1,1,ITHR),FIMD(1,1,ITHR),
     X    MAXCH,MAXM,MFLI,FACT,FPMAX,KO,RENORM,mpihelp,mfact,iame)
 	  enddo 
 	  
#ifdef MPI
 ! collect distant results !
 	DO 60 ITHR=2,mpihelp    ! from helpers
	 ihelper = ithr-1
	 write(500+iame,*) ' Node ',iame,' to receive stuff from ',
     x    ihelper, ' in commgroup: node =',ihelper+iame

!            write(500+iame,*) 'FIMD to receive as ',NBL*NEQS,' real'      
 	call MPI_recv(FIMD(1,1,ITHR),NBL*NEQS,MPI_DOUBLE_PRECISION,
     x                ihelper,100,commgroup,status,ierr)
!	 write(500+iame,*) ' Node ',iame,'  received FIMD from ',ihelper
      if(MINTL>1) then
	call MPI_recv(FIM(1,1,ITHR),NBL*NEQS,MPI_DOUBLE_PRECISION,
     x                            ihelper,101,commgroup,status,ierr)
	 write(500+iame,*) ' Node ',iame,'  received FIM from ',ihelper
	 else
      call MPI_recv(FI(1,1,ITHR),NBL*NEQS,MPI_DOUBLE_PRECISION,ihelper,
     x                                 102,commgroup,status,ierr)
	 write(500+iame,*) ' Node ',iame,'  received FI from ',ihelper
	 endif
      call MPI_recv(MAGN,NM1,MPI_DOUBLE_PRECISION,ihelper,103,
     x                                 commgroup,status,ierr)
      call MPI_recv(RENORM,1,MPI_DOUBLE_PRECISION,ihelper,104,
     x                                 commgroup,status,ierr)
	 write(500+iame,*) ' Node ',iame,'  received all from ',ihelper
 60	continue
#endif /* MPI */
	if(MINTL>1) deallocate(FI)
C-----------------------------------------------------------------------
c
 61   continue  

	
#ifdef SCP
	ENDIF       ! FJSWTCH
 	write(500+iame,*) 'Finished Numerov @',TM(I)
	MODEQ = 1
#else /* not SCP */	

	T = 2*NEQS*(2*NEQS+1)*8e-9
 	write(500+iame,*) 'Allocate MAT @',TM(I),' with ',real(T),' GB'
	allocate(MAT(2*NEQS,2*NEQS+1))
       DO 65 I=1,NR
 65    MAT(I,1:NP) = 0d0
	 DO 70 ITHR=1,numthread
		I1 = (ITHR-1)*NBL + 1
		I2 = min(ITHR*NBL,NEQS)
		II = I2-I1+1
		
           DO 67 K=1,NEQS
             MAT(K+NEQS,I1:I2) = FIMD(1:II,K,ITHR)
	     if(MINTL>1) then
             MAT(K,I1:I2)   = FIM(1:II,K,ITHR)
            else
             MAT(K,I1:I2)   = FI(1:II,K,ITHR)
	      endif
 67	    continue	    
 70          CONTINUE
!	if(MINTL==1) deallocate(FI)
        write(500+iame,*) 'Allocated MAT @',TM(I)
C
      ELSE  ! FJSWTCH true
c                              skip nmrov , match CRCWFN to zero
  	 MAT(:,:) = 0.0
         DO 780 IT=1,NEQS
               MAT(IT+NEQS,IT) = ONEC
               MAT(IT,IT)   = ONEC*CRCRAT(IT)
 780       CONTINUE
	ENDIF       ! FJSWTCH

C                    MATCH TO REQUIRED ASYMPTOTIC CONDITIONS
C
      IF (FCWFN) THEN
C                           match Numerov to CRCWFN
c
       DO 82 K=1,NEQS
       DO 82 K2=1,NEQS
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
100        IF (SHOW.GE.3) WRITE(KO,105) I,(MAT(I,J),J=1,NP)
105      FORMAT(' MAT(',I2,',*) =', /(6(1P,E10.2,E9.1,1X)))
C
 	if(SHOW>10) then
 		DO I=1,NR
 		DO J=1,NR
 		write(465,106) I,J,MAT(I,J),ZMAT(I,J),MAT(I,J)-ZMAT(I,J)
 		enddo
 		enddo
106		format(' MAT(',I4,',',i4,') =',4G12.3,' -? ',2G12.3)
		endif
 	MODEQ = 2
#endif /* SCP */	

	if(MODEQ==1) then  ! Scalapack
#ifdef SCP
	! my scalapack context is ICNTXT(master+1)=ICN
	ICN = ICNTXT(master+1)
	! send out MAT array to various helpers
	NBCOL = NBL ! same blocking factor as for the wfs
	NBROW = NBCOL ! scalapack needs square blocks
	 write(500+iame,*) ' START PZGETR',NBCOL,' @ ',TM(I)
       CALL BLACS_GRIDINFO(ICNTXT(master+1),NPROW,NPCOL,MYROW,MYCOL)
        write(500+iame,*) ' BLACS_GRIDINFO: R,C =',MYROW, MYCOL,
     x   ' in ',NPROW,NPCOL
      LDA = NUMROC(NR,NBROW,MYROW,0,NPROW)
      LDB = NUMROC(NR,NBCOL,MYCOL,0,NPCOL)
         write(500+iame,*) ' MATL local size =',LDA,LDB,' from ',NBL
   

	T = 8e-9*LDA*LDB
	write(500+iame,*) ' Allocate MATL arr',real(T),' GB @',TM(I)
	allocate (MATL(LDA,LDB),RHSL(LDA))  ! local storage
	call DESCINIT(DESCM,NR,NR,NBROW,NBCOL,0,0,ICN,LDA,info)
	 write(500+iame,*) ' Done DESCM ',info	
	 write(500+iame,'(9i6)') DESCM	
	call DESCINIT(DESCR,NR,1 ,NBROW,NBCOL,0,0,ICN,LDA,info) ! RHS here only
	 write(500+iame,*) ' Done DESCR ',info
	 write(500+iame,'(9i6)') DESCR	
	
! distribute MAT array to the MATL
!	call pzladist0(MAT,NR,NR,NR, MATL,DESCM, 0,0)
	call pzladist0F(ZMAT,NR,NR, MATL,DESCM, 0,0)

	 write(500+iame,*) ' Done distribute MAT ',' @ ',TM(I)
	
	 write(500+iame,*) ' Call PZGETRF with square',NR
      CALL PZGETRF(NR,NR,MATL,1,1,DESCM,IPVT,INFO)
	 write(500+iame,*) ' PZGETRF done',INFO,' @ ',TM(I)
	 
!	call pzladist0(MAT(1,NP),NR,1,NR, RHSL,DESCR, 0,0)
	call pzladist0F(ZRHS,NR,1, RHSL,DESCR, 0,0)
	 write(500+iame,*) ' Done distribute RHS ',' @ ',TM(I)

	call PZGETRS('N',NR,1,MATL,1,1,DESCM,IPVT,RHSL,1,1,DESCR,INFO)
	 write(500+iame,*) ' PZGETRS done',INFO,' @ ',TM(I)
!      	write(500+iame,'(10f10.5)') RHSL

	S = 0  ! FIX RENORMING LATER
!	 call flush(500+iame)
	
! Retrieve solution by each making part
	allocate(MAT(NR,1))  ! rump vector, just for the solution
	NP = 1

	MAT(:,NP) = 0.
      DO I = 1, NR
	CALL PZELGET( 'A', ' ', MAT(I,NP), RHSL, I, 1, DESCR )
      END DO
	 write(500+iame,*) ' Done gather RHS '	,' @ ',TM(I)
!      	write(500+iame,'(10f10.5)') MAT(1:NR,NP)
	
#endif /* SCP */

	else if(MODEQ==2) then  ! Lapack
            CALL ZGETRF (NR,NR,MAT,NR,IPVT,IER)
            IF ( IER .NE. 0 ) THEN
               WRITE (IOUT,1000) NR,IER
 1000 FORMAT (' ERWIMR1 ',i6,' : ERROR RETURN FROM ZGETRF, IERR = ',I6)
		WRITE(IOUT,1001) (MAT(I,I),I=1,NR)
 1001 FORMAT (1p,10e12.4)
               GO TO 600
            ENDIF
            CALL ZGETRS('N',NR,1,MAT,NR,IPVT,MAT(1,NP),NR,IER)
            IF ( IER .NE. 0 ) THEN
               WRITE (IOUT,1010) IER
 1010 FORMAT (' ERWIMR1 : ERROR RETURN FROM ZGETRS, IERR = ',I6)
               GO TO 600
            ENDIF
	  S = 0d0
	  do K=1,NR
	   T = 1.0
	   if(IPVT(K)/=K) T = -1.0
	   S = S + log(MAT(K,K)*T)
	  enddo
	else ! MODEQ=3 : plain

       CALL GAUSS(NR,NR,MAT,SING,S,SMALL,.FALSE.)
          IF(SING) GO TO 600
	endif
        write(500+iame,*) 'Solved MAT @',TM(I)
	if(MINTL==1) deallocate(FI)
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
	 DO 149 ITHR=1,numthread
		I1 = (ITHR-1)*NBL + 1
		I2 = min(ITHR*NBL,NEQS)
		II = I2-I1+1
!	 DO 149 IT=1,NEQS
	 DO 149 I=1,II
            IT = I1+I-1
         IF(MAT(IT,NP).EQ.CZ) GO TO 149
            WVD = WVD + MAT(IT,NP)*FIMD(I,K,ITHR)
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
1157  FORMAT(' ERWIMR1: LOG10(DETERMINANT) =',2F11.3,
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
205   FORMAT('0****** WARNING : ACCURACY LOSS IN ERWIMR1 =',F5.1,
     X ' DIGITS, SO EXPECT ERRORS OF',F9.4,' % OF UNITARITY'/)
      IF(AL*RE.GT..03) WRITE(KO,*) AL,ALUN,RE
	endif
	 write(500+iame,*) ' Node ',iame,'  done set'
	
      deallocate(MAT)
      RETURN
600   continue
      DO 605 I=IS,NM1
605   IF(MAGN(I).NE.0.0) MAGN(I) = LOG10(MAGN(I))
      WRITE(KO,610) (MAGN(I),I=IS,NM1)
610   FORMAT(' Maximum magnitudes were (Log 10) ',/(1X,20F6.1))
      CALL ABEND(4)
****************************************************** MATRIX FUNCTION	
	CONTAINS
	FUNCTION ZMAT(I,J)
	COMPLEX*16 ZMAT
	ZMAT = 0d0
	if(FJSWTCH) then ! no interior integration
	  IT = J
	  if(I==IT+NEQS) ZMAT = 1d0
	  if(I==IT)      ZMAT = CRCRAT(IT)
	else ! not FJSWTCH
! from FIMD etc
	K  = I - NEQS
	K2 = J
	 if(K>0.and.K2>0.and.K<=NEQS.and.K2<=NEQS) then !get from FIMD
	 ITHR = (K2-1)/NBL + 1
	 III = K2 - (ITHR-1)*NBL
	 ZMAT = FIMD(III,K,ITHR)	
	 endif
	K = I 
	 if(K>0.and.K2>0.and.K<=NEQS.and.K2<=NEQS) then !get from FIMD
	 ITHR = (K2-1)/NBL + 1
	 III = K2 - (ITHR-1)*NBL
	 
	 if(MINTL>1) then
	   ZMAT = FIM(III,K,ITHR)	
	   else
	   ZMAT = FI(III,K,ITHR)
	   endif
	 endif	
	endif
	
	if(FCWFN) then ! coupled asymptotic wf
	  K = J - NEQS
	  K2= I
	  if(K>0.and.K2>0.and.K<=NEQS.and.K2<=NEQS) then
         ZDL = CI**(L(K)-L(K2))
         ZHPLUS = CGMAT(K2,K,2)+CI*CFMAT(K2,K,2)	! M	  
	   ZMAT   = ZI2*ZDL*ZHPLUS
	   endif
	  K2 = I - NEQS
	  if(K>0.and.K2>0.and.K<=NEQS.and.K2<=NEQS) then
         ZDL = CI**(L(K)-L(K2))
         ZHPLUS = CGMAT(K2,K,1)+CI*CFMAT(K2,K,1)   	! MD	  
	   ZMAT   = ZI2*ZDL*ZHPLUS
	   endif	
	
	else  ! uncoupled asymptotic wf
	  K = J - NEQS
	  K2 = I
	  if(K==K2.and. K>0.and.K<=NEQS) then	  
	   ZMAT   = CH(1,K) ! M	  
	   endif
	  K2 = I - NEQS
	  if(K==K2.and. K>0.and.K<=NEQS) then	    
	   ZMAT   = CH(2,K)  ! MD
	   endif		
	endif
	
	END FUNCTION ZMAT
	
	
	FUNCTION ZRHS(I,J)
	COMPLEX*16 ZRHS
	ZRHS = 0d0
	if(J>1) return
	if(FCWFN) then ! coupled asymptotic wf
	  K = I 
	  if(K>0.and.K<=NEQS) then
         ZDL = CI**(L(EL)-L(K))
         ZHMINUS = CGMAT(K,EL,2)-CI*CFMAT(K,EL,2)	! M	  
	   ZRHS   = ZI2*ZDL*ZHMINUS
	   endif
	  K = I - NEQS
	  if(K>0.and.K<=NEQS) then
         ZDL = CI**(L(EL)-L(K))
         ZHMINUS = CGMAT(K,EL,1)-CI*CFMAT(K,EL,1)   	! MD	  
	   ZRHS   = ZI2*ZDL*ZHMINUS
	   endif	
	
	else  ! uncoupled asymptotic wf
	  if(I==EL     ) ZRHS = - CONJG(CH(1,EL))  ! M	  
	  if(I==EL+NEQS) ZRHS = - CONJG(CH(2,EL))  ! MD
	
	endif
	
	END FUNCTION ZRHS
      END
#ifdef MPI
****ERWINHELPER**************************************************************
	SUBROUTINE ERWINHELPER
	use parallel
	use parameters
	use factorials
	use fresco1
	use drier
        use mpi
	implicit none

! INPUT & OUTPUT
	integer NBL,NEQS,IS,NM1,NF,SHOW,MINTL,KO
	integer ITHR,MFLI,iv(12)
	real*8 SCL
	complex*16, allocatable:: FORMF(:,:),CLIS(:)
      REAL*8, allocatable::  H2C(:),ECM(:),H(:)
      INTEGER, allocatable:: L(:),CUTVAL(:),
     X        NCLIST(:,:),PFLIS(:,:),NFLIS(:)

! INPUT & OUTPUT
	real*8 RENORM

! OUTPUT
      real*8, allocatable:: FIMD(:,:),FIM(:,:),FI(:,:),MAGN(:)
      
! SCALAPACK WORKING
      INTEGER, allocatable:: IPVT(:)
      COMPLEX*16, allocatable:: MATL(:,:),RHSL(:),RHS(:)
      integer DESCM(9),DESCR(9),NBC,IC,NP,NR,I1,I2,II,info,ICN,I,J
	INTEGER NPROW,NPCOL, MYROW,MYCOL, NBROW,NBCOL, LDA,LDB,NUMROC      

      write(500+iame,*) 'BEFORE'   ,iame   
!	write(500+iame,*) NBL,NEQS,IS,NM1,NF,SHOW,MINTL,MAXCH,MAXM,KO,M,MD
!	write(500+iame,*)  ITHR,mpihelp,MFLI
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

1	write(500+iame,'(/80(''#''))')
	write(500+iame,*) 'Helper awaiting from commgroup 0 ',master	
      call MPI_recv(ITHR,1,MPI_INTEGER,0,1,
     x                                 commgroup,status,ierr)

	if(ITHR==0) then
           call MPI_finalize(ierr)
           write(500+iame,*) ' Stopping as requested'
           stop  ! end of Fresco run: stop
	   endif
     	write(500+iame,*)  'Request for ',ITHR,' of ',mpihelp
!	write(500+iame,*) 'OK:',MAXCH,LOCFIL,MAXM,M,MD,mpihelp
      call MPI_BCAST(NEQS,1,MPI_INTEGER,0,commgroup,ierr)
      if(NEQS==0) go to 1 ! stop  now, and look for another set

	call MPI_BCAST(iv,12,MPI_INTEGER,0,commgroup,ierr)      
 
       KO=iv(1); NBL=iv(2); IS=iv(3); NM1=iv(4); NF=iv(5)
	 SHOW=iv(6); MINTL=iv(7); MAXCH=iv(8); MAXM=iv(9)
	 M=iv(10); MD=iv(11); MFLI=iv(12)	 
     
      call MPI_BCAST(LOCFIL,1,MPI_LOGICAL,0,commgroup,ierr)
      call MPI_BCAST(SCL,1,MPI_DOUBLE_PRECISION,0,commgroup,ierr)

! future: remove MAXCH,LOCFIL,MAXM,M,MD

     	allocate(FORMF(MAXM,NF),CLIS(MFLI))
      allocate(H2C(NEQS),ECM(NEQS),H(NEQS),L(NEQS),CUTVAL(NEQS))
      allocate(NCLIST(MAXCH,MAXCH),PFLIS(MAXCH,MAXCH),NFLIS(MFLI))     
     
      call MPI_BCAST(FORMF,MAXM*NF,MPI_DOUBLE_COMPLEX,0,commgroup,ierr)
      call MPI_BCAST(CLIS,MFLI,MPI_DOUBLE_COMPLEX,0,commgroup,ierr)
      call MPI_BCAST(NFLIS,MFLI,MPI_INTEGER,0,commgroup,ierr)
      call MPI_BCAST(H2C,NEQS,MPI_DOUBLE_PRECISION,0,commgroup,ierr)
      call MPI_BCAST(ECM,NEQS,MPI_DOUBLE_PRECISION,0,commgroup,ierr)
      call MPI_BCAST(H,NEQS,MPI_DOUBLE_PRECISION,0,commgroup,ierr)
      call MPI_BCAST(L,NEQS,MPI_INTEGER,0,commgroup,ierr)
      call MPI_BCAST(CUTVAL,NEQS,MPI_INTEGER,0,commgroup,ierr)
      call MPI_BCAST(NCLIST,MAXCH**2,MPI_INTEGER,0,commgroup,ierr)      
      call MPI_BCAST(PFLIS,MAXCH**2,MPI_INTEGER,0,commgroup,ierr)
      
      write(500+iame,*) 'AFTER'      
!	write(500+iame,*) NBL,NEQS,IS,NM1,NF,SHOW,MINTL,MAXCH,MAXM,KO,M,MD
!	write(500+iame,*)  ITHR,MFLI,SCL
	
!!	 KO = 500+iame ! debugging
	
      allocate(FIMD(NBL,NEQS),FIM(NBL,NEQS),FI(NBL,NEQS),MAGN(NM1))
	
	call ERWICMR1(ITHR,NBL,NEQS,LOCFIL,IS,NM1,M,MD,NF,SHOW,MINTL,
     x    H2C,FORMF,CLIS,NFLIS,NCLIST,PFLIS,ECM,H,L,SCL,CUTVAL,MAGN,
     x    FI,FIM,FIMD,MAXCH,MAXM,MFLI,FACT,FPMAX,KO,RENORM,
     X    mpihelp,mfact,iame)      
      write(500+iame,*) 'ERWICMR1 has finished'      
      
            write(500+iame,*) 'FIMD to send as ',NBL*NEQS,' complex'      

	call MPI_send(FIMD,NBL*NEQS,MPI_DOUBLE_PRECISION,0,100,
     x                                 commgroup,ierr)
      write(500+iame,*) 'ERWICMR1 FIMD sent to',0,MINTL  
      if(MINTL>1) then
	call MPI_send(FIM,NBL*NEQS,MPI_DOUBLE_PRECISION,0,101,
     x                                 commgroup,ierr)
      else
      call MPI_send(FI,NBL*NEQS,MPI_DOUBLE_PRECISION,0,102,
     x                                 commgroup,ierr)
      endif
      call MPI_send(MAGN,NM1,MPI_DOUBLE_PRECISION,0,103,  ! CHANGE THIS
     x                                 commgroup,ierr)
      call MPI_send(RENORM,1,MPI_DOUBLE_PRECISION,0,104,  ! CHANGE THIS!
     x                                 commgroup,ierr)
      write(500+iame,*) 'ERWICMR1 results sent to',master      
     
	deallocate(FORMF,CLIS,H2C,ECM,H,L,CUTVAL,NCLIST,PFLIS,NFLIS)
	deallocate(FIMD,FIM,FI,MAGN)
#ifdef SCP	
! Now do scalapack work for the simultaneous equations:
	NR=NEQS*2; NP=NR+1
	! my scalapack context is ICNTXT(master+1)=ICN
	ICN = ICNTXT(master+1)
	! send out MAT array to various helpers
	NBCOL = NBL ! same blocking factor as for the wfs
	NBROW = NBCOL ! scalapack needs square blocks
	 write(500+iame,*) ' START PZGETR',NBCOL
       CALL BLACS_GRIDINFO(ICNTXT(master+1),NPROW,NPCOL,MYROW,MYCOL)
        write(500+iame,*) ' BLACS_GRIDINFO: R,C =',MYROW, MYCOL,
     x   ' in ',NPROW,NPCOL
      LDA = NUMROC(NR,NBROW,MYROW,0,NPROW)
      LDB = NUMROC(NR,NBCOL,MYCOL,0,NPCOL)
         write(500+iame,*) ' MATL local size =',LDA,LDB
    

	call DESCINIT(DESCM,NR,NR,NBROW,NBCOL,0,0,ICN,LDA,info)
	 write(500+iame,*) ' Done DESCM ',info	
	 write(500+iame,'(9i6)') DESCM	
	call DESCINIT(DESCR,NR,1 ,NBROW,NBCOL,0,0,ICN,LDA,info) ! RHS here only
	 write(500+iame,*) ' Done DESCR ',info
	 write(500+iame,'(9i6)') DESCR	
	allocate (MATL(LDA,LDB),RHSL(LDA),IPVT(NR))  ! local storage

	call pzladist1F(NR,NR, MATL,DESCM, 0,0)

	 write(500+iame,*) ' Done distribute MAT '
	 write(500+iame,*) ' Call PZGETRF with square',NR
      CALL PZGETRF(NR,NR,MATL,1,1,DESCM,IPVT,INFO)
	 write(500+iame,*) ' PZGETRF done',INFO

	call pzladist1F(NR,1, RHSL,DESCR, 0,0)

	 write(500+iame,*) ' Done distribute RHS '
           
	call PZGETRS('N',NR,1,MATL,1,1,DESCM,IPVT,RHSL,1,1,DESCR,INFO)
	 write(500+iame,*) ' PZGETRS done',INFO
!	 call flush(500+iame)
!      	write(500+iame,'(10f10.5)') RHSL
      
! Retrieve solution by each making part
 	allocate(RHS(NR))
	RHS(:) = 0.
      DO I = 1, NR
	CALL PZELGET( 'A', ' ', RHS(I), RHSL, I, 1, DESCR )
      END DO
	 write(500+iame,*) ' Done gather RHS '
!	 call flush(500+iame)
!      	write(500+iame,'(10f10.5)') RHS(1:NR)
	 
	deallocate(MATL,IPVT,RHSL,RHS)
#endif /* SCP */	
	go to 1
      end
#endif /* MPI */
****ERWICMR1**************************************************************
      SUBROUTINE ERWICMR1(ITHR,NBL,NEQS,LOCFIL,IS,NM1,M,MD,NF,SHOW,MINTL
     x   ,H2C,FORMF,CLIS,NFLIS,NCLIST,PFLIS,ECM,H,L,SCL,CUTVAL,MAGN,
     x    FI,FIM,FIMD,MAXCH,MAXM,MFLI,FACT,FPMAX,KO,RENORM,
     X    mpihelp,mfact,iame)
	implicit none
! INPUT
	integer NBL,NEQS,IS,NM1,NF,SHOW,MINTL,MAXCH,MAXM,KO,M,MD
	integer ITHR,mpihelp,MFLI,mfact,iame
	logical LOCFIL
	real*8 SCL
	complex*16 FORMF(MAXM,NF),CLIS(MFLI)
      REAL*8  H2C(MAXCH),ECM(MAXCH),H(MAXCH),FACT(mfact)
      INTEGER L(MAXCH),CUTVAL(MAXCH),
     X        NCLIST(MAXCH,MAXCH),PFLIS(MAXCH,MAXCH),NFLIS(MFLI)

! INPUT & OUTPUT
	real*8 RENORM

! OUTPUT
      REAL*8 FIMD(NBL,NEQS),FIM(NBL,NEQS),FI(NBL,NEQS),MAGN(NM1)
	     
! LOCAL
      REAL*8 ZP(NEQS),FORMFR(NF)
      REAL*8 DIAG(NEQS),ZI(NBL,NEQS),ZM(NBL,NEQS),MPI_WTIME
      complex*16 C
	integer I1,I2,II,IT,I,J,K,NC,ic,icr,icm,LL
	real*8 RI2,SEC,R12,ENA2,ENA3,TT0,TTT0,FMIN,F8,SMALL,R,T,
     x       UNC,FPMAX,EPS,BIG,SECOND,TMIN,TMAX,SEC0
	logical tr

       call system_clock(ic,icr,icm)
       if(icr.ne.0) SEC0 = real(ic)/real(icr)
	tr = SHOW>0
!	tr = .true.
       
        SMALL = 1.0/FPMAX
        EPS = SQRT(SMALL)
        BIG = 1./EPS
        UNC = 1d-100  ! Threshold for negligible couplings
	
      TMAX = 20.
      TMIN = -125.
      R12 = 1D0/12D0
      ENA2 = 2D0/5D0 * R12**2
      ENA3 = - 4D0/35D0 * R12**3
	
	ZM(:,:) = 0d0; ZI(:,:) = 0d0
	
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
	  
!        write(470,12) -ITHR,MPI_WTIME()

         DO 55 I=IS,NM1
          RI2 = 1D0 / DBLE(I-1)**2
	    if(LOCFIL) read(19) FORMFR
       call system_clock(ic,icr,icm)
       if(icr.ne.0) SEC = real(ic)/real(icr)

	if(tr) write(500+iame,12) ITHR,SEC-SEC0,(I-1)*H(1)
12	format(' Thread ',i3,': elap =',f9.2,' @ R=',f6.2)
!	write(500+iame,*) ' OpenMP thread ',ITHR,' at time ',SEC-SEC0,
!     x           ', R =',real((I-1)*H(1))
	TT0 = second()
	  
          DO 13 K=1,NEQS
             C = 0d0
	     if(LOCFIL) then
             DO NC=1,NCLIST(K,K)
               LL = PFLIS(K,K)+NC
		   C=C+CLIS(LL) * FORMFR(NFLIS(LL))
	       ENDDO
	     else
             DO NC=1,NCLIST(K,K)
               LL = PFLIS(K,K)+NC
	 	   C=C+CLIS(LL) * FORMF(I,NFLIS(LL))
	       ENDDO
	     endif
!             C = C * (0.,1.)**(L(J)-L(K)) * H2C(K)  ! but diagonal:
             C = C * H2C(K)
 
           DIAG(K) = -L(K)*(L(K)+1d0)*RI2 - ECM(K) * H2C(K) + C
	     
           IT = K
              IF(I==CUTVAL(K).and.I1<=IT.and.IT<=I2) then
                 J = L(K) + 1
                 R = (I-1)*H(K)
                 T = LOG(R)*J - FACT(J)*0.5
                 T = MAX(TMIN, MIN(TMAX, T ))
              ZI(IT-I1+1,K) = SCL * EXP(T)
		  endif
 13        CONTINUE

		IF(SHOW.GE.3) then
		do 14 K=I1,I2
		IT = K
 14		 if(I==CUTVAL(K)) write(KO,1302) I,IT,K,ZI(IT-I1+1,K)
 1302         FORMAT(' AT I=',I3,' ZI(',I5,',',I5,')=',1P,2E12.2)
		endif

         DO 158 K=1,NEQS
	     T = (1d0 - DIAG(K) * 
     x                  (R12 - DIAG(K)*(ENA2 + DIAG(K)*ENA3)))
 	       FI(1:II,K) = ZI(1:II,K) * T
158	    continue

	   IF(SHOW.GE.7.and.I1==1) THEN
          DO K=1,NEQS
		WRITE(KO,165) I,K,H(K)*(I-1),L(K),
     &    DIAG(K)/H2C(K)+ECM(K),0.0,FI(1,K),ZI(1,K)!*(0.,1.)**(L(EL)-L(K))
!     &               ,FI(K,K,1),ZI(K,K,1) !, COUPL(K,1)
          ENDDO
 165      FORMAT(' AT I,K,R =',I6,I4,F7.3,
     &          ' L,PE,PSI =',I4,2F9.3,1P,6E10.1)
          ENDIF
	if(tr) write(500+iame,*) ' ERWIMR1: 158 done: ',real(SECOND()-TT0),ITHR
	    
			TTT0 = SECOND()
	    DO 24 K=1,NEQS
	     DO 24 J=1,NEQS
            if(K==J.or.NCLIST(K,J)==0) go to 24
	 	C = 0d0
           if(LOCFIL) then
            DO 18 NC=1,NCLIST(K,J)
               LL = PFLIS(K,J)+NC
 18          C=C+CLIS(LL) * FORMFR(NFLIS(LL))
           else
            DO 181 NC=1,NCLIST(K,J)
               LL = PFLIS(K,J)+NC
 181          C=C+CLIS(LL) * FORMF(I,NFLIS(LL))
           endif
	      if(abs(C)>UNC) then
             C = C * R12 *(0.,1.)**(L(K)-L(J)) * H2C(K)
!		write(49,19) 'ITHR,I,K,J,C = ',ITHR,I,K,J,real(C)
!!     x		,' @ ',real(SECOND()-TT0),real(SECOND()-TTT0)
 19		format(a,4i6,2g12.3,a,f8.4,f10.6)
			TTT0 = SECOND()
 		FI(1:II,K) = FI(1:II,K) -C*ZI(1:II,J)
!	 	call daxpy(II,C,ZI(1,J),1,FI(1,K))
		endif
 24       CONTINUE
	if(tr) write(500+iame,*) ' ERWIMR1: 24 : ',real(SECOND()-TT0),ITHR
	    
         DO 43 K=1,NEQS
	   
	    ZP(I1:I2) = 2d0*ZI(1:II,K)-ZM(1:II,K) 
     x                     - DIAG(K)*FI(1:II,K)
	    
           DO 35 J=1,NEQS
            if(K==J.or.NCLIST(K,J)==0) go to 35
	 	C = 0d0
           if(LOCFIL) then
            DO 28 NC=1,NCLIST(K,J)
               LL = PFLIS(K,J)+NC
 28          C=C+CLIS(LL) * FORMFR(NFLIS(LL))
           else
            DO 29 NC=1,NCLIST(K,J)
               LL = PFLIS(K,J)+NC
 29           C=C+CLIS(LL) * FORMF(I,NFLIS(LL))
           endif
	      
	      if(abs(C)>UNC) then
	 	C = C * (0.,1.)**(L(K)-L(J)) * H2C(K)
		ZP(I1:I2) = ZP(I1:I2) - C * FI(1:II,J)
		endif
 35	    continue

            ZM(1:II,K) = ZI(1:II,K)
            ZI(1:II,K) = ZP(I1:I2)
 43       CONTINUE
	if(tr) write(500+iame,*) ' ERWIMR1: 43 : ',real(SECOND()-TT0),ITHR
	    
       IF(I.EQ.M) THEN
         if(MINTL>1) FIM(1:II,:) = FI(1:II,:)
!	   exit at this point, so use FI rather than save in FIM	
	write(500+iame,*) 'Found M =',M,MINTL
       ELSE IF(I.EQ.MD) THEN
        FIMD(1:II,:) = FI(1:II,:)
	write(48,*) 'Found MD =',MD
       ENDIF

	 if(mpihelp>1) go to 51 ! skip renorms for now
            T = 0.0
         DO 47 IT=1,II
            DO 47 K=1,NEQS
            F8 = FI(IT,K)
            T = MAX(T, abs(F8))
   47       IF(T .GT. BIG ) GO TO 48
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
            ZI(1:II,:) = ZI(1:II,:) * T
            ZM(1:II,:) = ZM(1:II,:) * T
            if(MINTL>1)  FIM(1:II,:) = FIM(1:II,:) * T
            if(MINTL==1)  FI(1:II,:) = FI(1:II,:) * T
            FIMD(1:II,:) = FIMD(1:II,:) * T
C
 51      CONTINUE
	if(mod(I,M/4)==0) 
     x    write(500+iame,*) ' ERWIMR1: I=',I,' @',real(SECOND()-TT0)

	if(tr) write(500+iame,*) ' ERWIMR1: 51 : ',real(SECOND()-TT0)
	 if(I.eq.M) go to 60 ! so use FI not FIM

 55       CONTINUE
  60        continue
       call system_clock(ic,icr,icm)
       if(icr.ne.0) SEC = real(ic)/real(icr)

        !write(470,12) -ITHR,MPI_WTIME()
        !write(500+iame,12) ITHR,SEC-SEC0
        !if(iame>0) write(470,12) ITHR,SEC-SEC0
	return
 	 END
#ifdef SCP
!	include 'pzladist.f'
	include 'pzladistf.f'
#endif /* SCP */
