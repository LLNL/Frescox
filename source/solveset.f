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
****SOLVESET**************************************************************
!	SUBROUTINE SOLVESET
!      	implicit real*8(a-h,o-z)
C
C    SOLVING THE SET OF COUPLED CHANNELS FOR EACH INCOMING CHANNEL
C    -------------------------------------------------------------
      REPEAT = .FALSE.
	SKIP(:,1) = .false.
      DO 700 JIN=1,MINTL
         EL = INITL(JIN)
         JPSET = JPSET + 1
         IF(JSET.GT.0 .AND. JPSET.GT.JSET) GO TO 700
         DO 430 C=1,NCH
          SOMA(C,1) = 0.0
          SIMPLE(C) = .TRUE.
          SKIP(C,2) = .FALSE.
C                      SKIP(C,1) = SOURCE SAME AS LAST TIME
C                      SKIP(C,2) = FED(C) = SOURCE TERM NON ZERO
         SRC(:,C) = 0.0
         IF(DISCIT) write(8,rec=NCH+C) SRC(1:N,C)
         IF(DISCIT) written(8) = .true.
        IF(C<=NICH) then
         PSI(:,C) = 0.0
         IF(DISC8)write(8,rec=C) PSI(1:N,C)
	endif
  430    CONTINUE
C
      FAIL = .FALSE.
      AL = 0.0
      BEST = 1E5
      NBEST = 0
       NSOLI = NSOL
       if(FJSWTCH) NSOLI = 1
       if(FJSWTCH) NBEST = 1
      DO 520 ITNL=1,NSOLI
       call flush(6)
      IF(.NOT.REPEAT) THEN
C
         DO 438 C=1,NCH
          IT = ITC(PART(C,1),EXCIT(C,1))
         L = LVAL(C)
           IF(ISOCEN.eq.1) L = JTOTAL + 0.1
           IF(ISOCEN.eq.2) L = LVAL(EL)
           LL1(C) = L*(L+1d0)
         DO 437 I=1,2
  437    CH(I,C) = CHL(L+1,IT,I) * CI * HALF
         IF(CDETR.GE.2) WRITE(KO,1300) C,(CH(I,C),I=1,2)
 1300      FORMAT('0',I6, ' MATCH TO',1P,2E13.4, ' AND',2E13.4)
  438    CONTINUE
      ENDIF
C
      WREQ = ITNL.LT.NSOLI .OR. WOUT
      SHSMAT = SMATL.GE.5.AND.PADE.EQ.0.OR.SMATS.GE.6.OR.SMATL.GE.3.AND.
     X ITNL.EQ.NSOLI
      SMATEL = SMATS.GE.2.AND.ITNL.EQ.1
	RERR=ACC8
	NQ = 4*NCH
	NCH2 = 2*NCH
C
      ! if(IEX.lt.NCH.and.FALLOC.and.IBLOCK.lt.99) then
      ! if(IBLOCK.lt.99) then
      ! if(IEX.lt.NCH.or.FALLOC) then
      ! if(IEX.lt.NCH.or.WOUT.or.ITER>0) then

	GVAL = 0.
      if(FALLOC) then  
	
      CALL ERWIN(PSI,ECM(1,1),ECM(1,2),IEX,FORMF,NF,FORMF,SRC,CH,
     X EL,SMAT,LVAL,JVAL,JTARG,LL1,NCH,NICH,N,ECM(1,3),M,MD,TSTD,RENORM,
     X REPEAT,WREQ,BLOCKD,AL,RERR,ICUTC,CDETR ,SHSMAT, SMATEL,
     X CUTOFF,SKIP(1,1),SKIP(1,2),PTYPE,ferw,FIM,FIMD,LOCFIL,
     X CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)

      else if(.not.CCREAL) then
c		Erwin version optimised for pure CC S-matrix solutions:

	if(.true.) then
      CALL ERWINCC4(ECM(1,1),ECM(1,2),FORMF,NF,FORMF,CH,REPEAT,NBL,EL,
     X SMAT,LVAL,JVAL,JTARG,LL1,NCH,N,ECM(1,3),M,MD,TSTD,RENORM,
     X AL,RERR,ICUTC,CDETR ,SHSMAT, SMATEL,CUTOFF,PTYPE,FIM,FIMD,LOCFIL,
     X CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)
 	else
      CALL ERWINCC(ECM(1,1),ECM(1,2),FORMF,NF,FORMF,CH,REPEAT,EL,
     X SMAT,LVAL,JVAL,JTARG,LL1,NCH,N,ECM(1,3),M,MD,TSTD,RENORM,
     X AL,RERR,ICUTC,CDETR ,SHSMAT, SMATEL,CUTOFF,PTYPE,FIM,FIMD,LOCFIL,
     X CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,CLIST,NCLIST,NFLIST)
	endif


	else
c		Erwin version for pure CC S-matrix, REAL couplings!
	deallocate(CLIST)
      CALL ERWIMR1(ECM(1,1),ECM(1,2),FORMF,NF,CH,REPEAT,NBL,EL,
     X SMAT,LVAL,JVAL,JTARG,LL1,NCH,N,ECM(1,3),M,MD,TSTD,RENORM,MINTL,
     X AL,RERR,ICUTC,CDETR ,SHSMAT, SMATEL,CUTOFF,PTYPE,FIMR,FIMDR,
     X LOCFIL,CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,CUTVAL,
     X CLIS,NCLIST,NFLIS,PFLIS,MFLI)
	if(numthread>1) 
     x 	 write(500+iame,*) ' Node ',iame,'  done cc set 1'
	allocate(CLIST(MAXCH,MAXCH,MCLIST))

	endif

C            READ IN INITIAL WFS, if requested
      IF(INITWF.NE.0.and..not.INITWFE) THEN
	  NRF = abs(INITWF)
	  PSI(:,:) = 0d0     ! all wfs=0 EXCEPT those read in!
!          call openif(NRF)
          NA = (N-1)/2*2 + 1 ! read these
      if(ITNL==1) then ! first iteration
      IF(INITWF>0) READ(NRF,657,END=1692) NAI,HPI,ENLABI,
     x                                    JTOTALI,PARITYI
      IF(INITWF<0) READ(NRF,END=1692) NAI,HPI,ENLABI,JTOTALI,PARITYI
	if(NA/=NAI) then
 	  write(0,*) 'STOP: read wfs NA ',NA,' neq ',NAI; stop
	 endif
	if(abs(JTOTAL-JTOTALI)>0.1.or.PARITY/=PARITYI) then
 	  write(0,*) 'STOP: read J/pi ',JTOTALI,PARITYI
	  stop ' J/pi in wrong sequence'
	 endif
!657   FORMAT(I4,2F8.4,F8.1,I3,2f12.6,2f8.3)
      DO 1690 C2=1,NCH+1     ! look for at most NCH+1 channel wfs !
        IF(INITWF.GT.0)
     X     READ(NRF,660) ITI,LVALI,JVALI,JTOTALI,LVALE,JVALE,SMATI
        IF(INITWF.LT.0)
     X     READ(NRF) ITI,LVALI,JVALI,JTOTALI,LVALE,JVALE,SMATI
	  if(ITI<0) go to 1695
            IF(LVAL(EL).NE.LVALE .or. abs(JVAL(EL)-JVALE)>0.1) then
 	    write(0,*) 'STOP: wrong elastic I,J: ',LVALE,JVALE
	    stop  ' Elastic channels in wrong sequence'
	    endif
          DO 1670 C=1,NICH  !  look for this channel in current set
            IF(ITC(PART(C,1),EXCIT(C,1)).NE.ITI) GOTO 1670
            IF(LVAL(C).NE.LVALI.or.abs(JVAL(C)-JVALI)>0.1) GOTO 1670
	    write(ko,1696) 1,C,SMATI
	     SMAT(C) = SMATI
	     READWF(C) = .true.
           IF(INITWF>0) then
             READ(NRF,'(6e12.4)') (PSI(I,C),I=1,NA)
	    ELSE
             DO 1665 I=2,NA
 1665        READ(NRF) PSI(I,C)
            ENDIF
 1670      CONTINUE
! 681      FORMAT(1P,6E13.6)
 1690  CONTINUE
       go to 1695
 1692	INITWFE = .true. ! reached EOF
 1695  CONTINUE
       else !  ITNL > 1, on later iterations, just read from file 8
	!write(6,*) 'READWF alloc2',allocated(READWF); call flush(6)
	 do C=1,NICH
         if(READWF(c)) then 
 	   read(8,rec=C) PSI(1:N,C)
  	   SMAT(C) = SOMA(C,1) 
	    write(ko,1696) ITNL,C,SMAT(C)
 1696    format(' Iteration #',i2,': fixed wf for  channel ',i4,
     x          ', so S =',2f12.8)
	   endif
	 enddo  ! C
       endif    ! ITNL==1
       ENDIF    ! INITWF/=0
C
      IF(NPOST.AND.ITNL.LT.NSOLI.AND.ITNL.GT.1)
     X    CALL RENO(PSI,PHI,N,ECM(1,3),NCH,CUTOFF,NLN,NLO,MR,MLM,0,
     X         NF,  0,NL,EMPTY,MLT,NLC,SMATL,CLIST,NCLIST,NFLIST,CHNO,
     X              NPOST,NPRIOR,SRC,ECM(1,2),SIMPLE,IEX,
     X              FORMF,FORMF,LOCFIL,ECM(1,1),LVAL,PTYPE,ICUTC,NICH)
C
      IF(FLIP.AND.ITNL.GT.1) CALL XCH(PSI,N,NCH,NICH,SMAT,JTOTAL,LVAL,
     X      JVAL,JPROJ,JTARG,PART,EXCIT,COPY,MXP,MXX,MAXN,SMATL,SAME,
     X      .false.,EXCH,MAXCH)
C
C      PADE ACCELERATION
C      -----------------
        DO 470 C=1,NCH
470     SRATIO(C) = SMAT(C)
      IF(PADE.GE.1)    CALL PAD(PADE,SMAT,NCH,ITNL,SPAD,MAXCH,MAXIT+2)
        DO 480 C=1,NCH
        SRATIO(C) = SMAT(C) / (SRATIO(C) + (1E-20,0.0))
 480    IF(.NOT.PSIREN) SRATIO(C) = (1.0,0.0)
C
C               THAT IS THE PADE APPROXIMANT BY THE EPSILON ALGORITHM
C               -----------------------------------------------------
      FRACT = 0.0
      DO 500 C=1,NCH
         T = ABS(SMAT(C))
         EMPTY(C) = T.LE.1E-12
         SAME(C) = .TRUE.
         DIFF = ABS(SMAT(C) - SOMA(C,1)) * 100.
     X            * (2.*JTOTAL+5.)/(2.*JTMAX+5.)
         SAME(C) = DIFF.LE.1E-6
         IF(ECM(C,1).LT.0.) SMAT(C)=0.  ! closed channels give no xs, but can still vary
         IF(WREQ.AND.DISC8.and.C<=NICH)  write(8,rec=C) PSI(1:N,C)
         IF(FRACT.GT.DIFF.OR.EXTRA(2,PART(C,1),EXCIT(C,1))) GO TO 490
            FRACT = DIFF
            C2 = C
  490    SOMA(C,1) = SMAT(C)
  500  CONTINUE
       IF(FRACT.LT.BEST.AND.ITNL-1.GE.ITMIN.and.
     X		(MOD(ITNL,2)==1.or.FRACT<abs(IPS))) THEN
          BEST = FRACT
          NBEST = ITNL
C   Store best SMAT in SOMA(*,2) and wfns in file 18
          IF(WOUT) THEN
           IF(.NOT.OPEN18) THEN
            IF(MACH==6.or.MACH==7.or.MACH==8.or.MACH==3.or.MACH==4)THEN
               OPEN(18,ACCESS='SEQUENTIAL',
     X            FORM='UNFORMATTED',FILE=trim(TMPN)//'18')
            ELSE
               OPEN(18,ACCESS='SEQUENTIAL',
     X            FORM='UNFORMATTED',STATUS='SCRATCH')
            ENDIF
           ENDIF
           OPEN18 = .TRUE.
           REWIND 18
          ENDIF
          DO 501 C=1,NCH
          SOMA(C,2) = SMAT(C)
          IF(WOUT.and.C<=NICH) WRITE(18) (PSI(I,C),I=1,N)
          IF(WOUT            ) WRITE(18) (SRC(I,C),I=1,N)
  501     continue
	  if(WOUT) written(18) = .true.
       ENDIF
      DO 505 I=1,N
  505 IF(ITNL.eq.1) VPOT(I,3) = PSI(I,EL) / SRATIO(EL)
      IF(PSIREN.AND.WOUT.AND.ABS(SRATIO(EL)-1.).GT.0.1.AND.ITNL.GE.3)
     X  WRITE(KO,1395) SRATIO(EL),FRACT
 1395 FORMAT(' Psi(EL) renorm. by ',1P,2E9.2,', Conv =',0P,F9.5,' %')
       IF(NSOLI.GT.4.AND.ITNL-1.GE.1.AND.
     X     (SMATL.GE.3.OR.SMATL.GE.2.AND.IT0.GT.0)) THEN

          FUSIO = 0.0
          DO C=1,NCH
             IC = PART(C,1)
             IA = EXCIT(C,1)
             IT = ITC(IC,IA)
          IF(ECMC(IC,IA) .GT. 0.0) THEN
          AMDSQS = ABS(SMAT(C))**2
          IF(C.EQ.EL) AMDSQS = 1 - AMDSQS
!@@
             IF(RMASS(IC).lt.1e-5) then
                XSC = JCOEF* AMDSQS * 
     X                ( RMASS(PEL)*amu/ (HBC*K(PEL,EXL)) )
             ELSE IF(RMASS(IC).lt.1e-5) then
                XSC = JCOEF* AMDSQS / 
     X                ( RMASS(IC)*amu/ (HBC*K(IC,IA)) )
             ELSE
                XSC = JCOEF* AMDSQS *
     X                K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
             ENDIF
!            see justification on frxx3.f
!
!            XSC = JCOEF* AMDSQS *
!    X             K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))
!@@
          IF(C.EQ.EL) THEN
             FUSIO = FUSIO + XSC
            ELSE
             FUSIO = FUSIO - XSC
            ENDIF
          endif
          ENDDO
              IF(PADE.GT.0) WRITE(KO,1400) ITNL-1,FRACT,C2,SMAT(C2),
     X                           FUSIO,SPAD(C2,MIN(ITNL,MAXIT+2),2)
              IF(PADE.EQ.0) WRITE(KO,1400) ITNL-1,FRACT,C2,SMAT(C2),
     X				 FUSIO
        ENDIF
 1400 FORMAT(' So max change at iter. #',I3,' =',F10.4,' %, @ CH',I3,
     X ' with',2F10.5,', Fus=',F10.5,:,' prev S=',1p,2e12.4)
      IF(SMATS.GE.5.AND.PADE.GE.1) WRITE(KO,1430) (SMAT(C),C=1,NCH)
      IF(ITNL-1.GE.ITMIN.AND.FRACT.LT.ABS(IPS) .OR. DRY) GO TO 530
      WERR = 3.* AL*RERR
      IF(FRACT.LT.WERR*0.1.AND.ITNL-1.GT.ITMIN.AND.IPS.GT.0.0) THEN
        WRITE(KO,1401) ITNL-1,FRACT,WERR
1401      FORMAT('0STOPPING AFTER',I4,' ITERATIONS, as change',F9.3,
     X         ' % is less than 0.1 of ERWIN accuracy loss ',F9.3,' %'/)
         GO TO 530
         ENDIF
      IF(FRACT .GT. 100.*BEST .AND. NBEST-1.GT.ITMIN.AND.IPS.GT.0.) THEN
C               Give up, as appears to be diverging
         WRITE(KO,1402) ITNL-1,FRACT,BEST,NBEST-1
1402     FORMAT('0STOPPING AFTER',I4,' ITERATIONS, as change of ',F9.3,
     x ' % is worse than 100 times best change of',F9.3,'% at IT =',I3/)
         GO TO 530
         ENDIF
      IF(ITNL.EQ.NSOLI.AND.VEFF.EQ.0) GO TO 520
C
      CALL SOURCE(SRC,PSI,N,ECM(1,3),NCH,IEX,FORMF,NF,FORMF,CUTOFF,
     X  ICUTC,SIMPLE,.TRUE.,NLN,NLO,MR,NL,EMPTY,SAME,SHFL,
     X  WAVES,LVAL, MLT,SKIP,SKIP(1,2),NLC,CHNO,
     X  MLM,NICH,MMXICH,PTYPE,LOCFIL,CLIST,NCLIST,NFLIST)
C
!     IF(NPRIOR.AND.ITNL.GT.1)
      IF(NPRIOR)
     X CALL RENO(SRC,PHI,N,ECM(1,3),NCH,CUTOFF,NLN,NLO,MR,MLM,1,
     X      NF,  NCH,NL,EMPTY,MLT,NLC,SMATL,CLIST,NCLIST,NFLIST,CHNO,
!    X    NPOST.AND.ITNL+1.LT.NSOLI,NPRIOR,SRC,ECM(1,2),SAME,IEX,
     X    NPOST,NPRIOR,SRC,ECM(1,2),SAME,IEX,
     X        FORMF,FORMF,LOCFIL,ECM(1,1),LVAL,PTYPE,ICUTC,NICH)
C
      DO 510 C=1,NCH
  510  IF(DISCIT) write(8,rec=NCH+C) SRC(1:N,C)
  520 REPEAT = .TRUE.
      ITNL = NSOLI
C
c      IF(IPS.NE.0..AND.ITER.GT.4) THEN
      IF(IPS.NE.0..AND.NSOLI.GT.4) THEN
        WRITE(KO,1410) ITER
 1410   FORMAT(' FAILED TO CONVERGE AFTER',I3,' ITERATIONS !!!!!!!',/)
        IF(FATAL) WRITE(KO,*) ' Set ITER negative to allow continuation'
        FAIL = .TRUE.
        penalty = penalty + fine
       if(number_calls>5) then
	 write(6,41) number_calls,JTOTAL,PSIGN(parity+2),penalty
41	format('  At call ',i5,', cc set ',f6.1,a1,
     & 		' failed to converge, so  penalty =',1p,e10.1)
	else
C       IF(FATAL) STOP
        IF(FATAL) CALL ABEND(2)
       endif
      ENDIF
C
C    NOW HAVE CONVERGED RESULT in itnl iterations
C    --------------------------------------------
C 530 IF(ITNL.NE.NBEST.AND.NBEST.NE.0.OR.FAIL) THEN
  530 IF(IPS.NE.0.0.AND.ITNL.NE.NBEST.AND.ITER.GT.4.OR.FAIL) THEN
C               Find S-matrix elements and wfns from ITNL = NBEST
C               Get SMAT from SOMA(*,2), wfns etc. from file 18.
          FRACT = BEST
          IF(WOUT) REWIND 18
          DO 531 C=1,NCH
          SMAT(C) = SOMA(C,2)
          IF(WOUT.and.C<=NICH) READ(18) (PSI(I,C),I=1,N)
  531     IF(WOUT) READ(18) (SRC(I,C),I=1,N)
        ENDIF
C
      IF(PADE.GE.1.OR.ITER.GT.4) THEN
      ITNL = ITNL - 1
      IF(SMATL.GE.2) WRITE(KO,1420) JTOTAL,PSIGN(PARITY+2),
     X                          NBEST-1,ITNL,FRACT,JTMAX,SMAT(EL)
      IF(SMATL.GE.3) WRITE(KO,1430) (SMAT(C),C=1,NCH)
 1420 FORMAT(' Final S-matrices (',F7.1,A1,') at',I3,' after',I3,
     X       ' accurate to',
     X F9.4,' % of J=',F7.1,' unitarity:',F10.5,' +i*',F8.5,',')
 1422 FORMAT(F10.1,2F12.8,'i: elastic S-matrix  @@',
     x       f10.2,
!$   x       '/',f8.2,
!    x       i3,i4,l2,i3,2i5)
     x       i3,i4,l2,i3,2i3,l3,1p,2e11.3)
 1430 FORMAT(5(1X,F11.5,' +i*',F9.5,',') )
!1431 FORMAT(5(1X,1p,e11.3,' +i*',e9.1,',') )
      ENDIF
      MMXIT = MAX(MMXIT, ITNL)
C
	if(WOUT) then
        DO 550 C=1,NICH
        DO 550 I=1,N
        PSI(I,C) = PSI(I,C) * SRATIO(C)
550     SRC(I,C) = SRC(I,C) * SRATIO(c)
	endif
C
      IF(SMATL.GE.2) then
          WRITE(KO,1422) JTOTAL,SMAT(EL),TM(I)-TIME0J,
!$   x                  OMP_GET_WTIME()-TIME0W,
     X			DROPPED,NSTEPD,FJSWTCH,iams ,
     x                  SMATL,SMATS,final !,DONES,DONE
!     x                 ,(SMAT(EL)-1.0)/(0.,2.)
	ENDIF
	phasin(JIN) = log(SMAT(EL))*(0.,-0.5)*180./pi
	linel(JIN) = EL

      NLN1 = NLN
      include 'usescatwf.f'
C
!      IF(CDETR.GT.1) CDETR = CDETR - 1
C                next JIN    :
  700 CONTINUE

	if(numthread>1) 
     x    write(500+iame,*) ' Node ',iame,'  done cc set 2'

!        RETURN
!	END SUBROUTINE SOLVESET
