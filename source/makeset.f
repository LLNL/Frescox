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

!	SUBROUTINE MAKESET
!      	implicit real*8(a-h,o-z)
***	integer I1,I2,IX1,IX2,LAM1,LAM2
!      	implicit none

        !write(6,*) 'ENLAB, elpmax:',ENLAB,elpmax
C
C      IF(LISTCC.GE.5) WRITE(KO,*) JTOTAL,PSIGN(PARITY+2)
C
C    PARTIAL WAVE SETS FOR EACH J-TOTAL/PARITY,  AND THEIR PARAMETERS
C    ----------------------------------------------------------------
       KINTL = -1
!       KINTL = 1000000
       IF(JSET.GT.0) KINTL = JSET-JPSET
	NSMALL = 5
!	if(SMALLJ.and.ITER==0.and.ITCM>1.and..not.CCREAL)  then	! Change from  CC  to DWBA
!	   IEX = 0
!	   ITER = 1
!	   write(KO,*) '  CHANGING TO DWBA!!!!!'
!	   ENDIF
	  

!				Find number of channels:
	CALL NUMCC(NCH,IEX,MINTL,MEL,JTOTAL,PARITY,JTMIN,KINTL,
     X        NEX,NCHAN,GIVEXS,PEL,EXL,LMAX,JEX,COPY,NFUSCH,LPMAX,
     X        RMASS,MAL1,ITC,IBLOCK,BAND,SMALLS,NSMALL,CHPRES,
     X        ENLAB,elpmax)
        IF(NCH.EQ.0) go to 350 
        IF(MINTL==0.and.melfil==0) go to 350 
!        NICH = max(IEX,MEL)
        NICH = MEL
        if(ITER>0)  NICH = max(NICH,IEX)
        if(ITER>1.or.TWOW.or.IBLOCK<0)  NICH = NCH
        if(IBLOCK<0)  NICH = MEL
	 NICH = max(NICH,NFUSCH)
	if(WOUT0) NICH = NCH
        MMXCH=MAX(MMXCH,NCH)
        NBL = (NCH-1)/max(1,numthread) + 1

       JCCSET = JCCSET + 1
	if(mpisets>1.and.MOD(JCCSET-1,mpisets).eq.iams) then
         T = TM(I)
!        write(51,*) 'Node ',iams,' finds ',nint(JTOTAL),
!    X                    PSIGN(PARITY+2),' has ',NCH,MINTL,real(T)
	 write(51,10) iams,JTOTAL,PSIGN(PARITY+2),NCH,MINTL,T
10	format(' Node',i3,' finds J=',f7.1,a1,' has ',2i3,' chs ',
     x		' at ',f9.2,' secs')
	 else
!         write(51,*) 'Node ',iams,' rejects ',nint(JTOTAL)
	endif
         call flush(51)
      IF(.NOT.SCHONP) JUMP(JBLOCK,2) = JUMP(JBLOCK,2) + 1
       SCHONP = .TRUE.
      JUMP(JBLOCK,3) = JUMP(JBLOCK,3) + 1
C
C     IN CONCURRENT MACHINES, JUST DO EVERY NUMNODE CC SET
C     ----------------------------------------------------
      IF(MPIC) THEN
CCCCC   Included for iPSC/860 & other concurrent
       IF(MOD(JCCSET-1,mpisets).ne.iams) THEN
          T = TM(I)
         if(say)write(48,*) 'Node ',iams,' not to do ',nint(JTOTAL),
     X                    PSIGN(PARITY+2),JCCSET,real(T)
	call flush(48)
          go to 350 
          ENDIF
        T = TM(I)
       if(say)write(48,*) 'Node ',iams,' to do ',nint(JTOTAL),
     X                    PSIGN(PARITY+2),JCCSET,real(T)
       call flush(48)
       ENDIF
!#	write(48,*)'NCH>MAXCH.or.NICH>MAXICH.or.IEX+2>MAXB'
!#	write(48,*) NCH,MAXCH,NICH,MAXICH,IEX,MAXB
	FALLOC = WOUT.or.ITER>0.or.IEX<NCH
!	if(.not.FALLOC) NFDEC=1
	if(.not.FALLOC) NFDEC=N*2*NICH
	if(NCH>MAXCH.or.NICH>MAXICH.or.IEX+2>MAXB) then
!					Deallocate arrays:
	if(MAXCH>0) then
	  deallocate (CLIST,NCLIST,NFLIST,CH,ECM,XS,INCOME,LAMERR,
     X    CRCRAT,CUTVAL,JVAL,JPROJ,JTARG,CHNO,EXCH,CFMAT,CGMAT,READWF,
     X    LVAL,PART,EXCIT,INITL,BLOCKD,LL1,LAM,GDER,SMAT,SAME,PSI)

	  if(.not.rterms) then
            deallocate(SRC,SRATIO,SOMA,SPAD,SIMPLE,EMPTY,SKIP,ferw)
          if(CCREAL.and..not.FALLOC) then    ! only for ERWINRC4, not yet ERWINCC4
	      deallocate(FIMR,FIMDR)
	    else
	      deallocate(FIM ,FIMD)
	    endif
           endif  ! not rterms

	  if(cfalloc) deallocate(CFF)
	endif
!				Allocate arrays:
	MAXCH = NCH
	if(VEFF.ne.0) MAXCH = max(MAXCH,NSA)
	MAXICH = NICH
	MAXXCH = max(MAXXCH,MAXCH)
	MAXCH = MAXXCH
	MFNL=NCH*NCH
	MAXB = IEX+2
!       NFDEC = 2 * N *(NICH*2+IEX*NICH)  + 4*MAXB*NCH
!       NFDEC = 2 * N *MAXB*NICH  + 4*MAXB*NCH
        NFDEC =  N *MAXB*NICH  !complex
!#	write(48,*) 'Makeset: NFDEC = ',NFDEC,N,MAXB,NICH,' set'
	if(NOSOL) then
	   NFDEC=1; MAXICH=max(IEX,MEL); MAXB=0; MAXIT=0
	   endif
	if(say.and.final) write(KO,*) 'Allocate arrays for ',MAXCH,
     X	 ' channels, of which ',MAXICH,' need wfs.'
	if(NFDEC>10 000 000) write(KO,*) 
     x		'Allocate NFDEC =',NFDEC,' in ',nint(NFDEC*16E-6),' MB'
      memuse = MAXCH*16e-9*MAXCH*MCLIST
      if(say.and.final)write(48,'(a,3i7,a,f7.3,a)')
     &'Allocate CLIST(',MAXCH,MAXCH,MCLIST,') in ',memuse,' GB' 
	T = (MAXCH*2)**2 * 16e-09
       if(T>.1) write(48,'(a,f8.3,a)') ' (if not SCP) Allocate MAT in '
     x     ,T,' GB'
      
      allocate(ECM(MAXCH,3),XS(MAXCH),INCOME(MAXCH),GDER(2*MAXCH))
      allocate(CUTVAL(MAXCH),BLOCKD(MAXCH),LL1(MAXCH),LAM(2*MAXCH))
      allocate(LVAL(MAXCH),PART(MAXCH,3),EXCIT(MAXCH,3),INITL(MAXCH))
      allocate(JVAL(MAXCH),JPROJ(MAXCH),JTARG(MAXCH),READWF(MAXCH))
       T = MAXCH**2 * 1e-9 ; X = T*4*(MCLIST+1)
      if(say.and.final)write(48,*) ' Allocate NCLIST in ',real(X),' GB'
      allocate(CLIST(MAXCH,MAXCH,MCLIST),CRCRAT(MAXCH))
      allocate(NCLIST(MAXCH,MAXCH),NFLIST(MAXCH,MAXCH,MCLIST))
      allocate(SMAT(MAXCH),CHNO(MFNL,6),CH(2,MAXCH),SAME(MAXCH))
      allocate(LAMERR(2*MAXCH))

      if(.not.rterms) then
       T = MAXN * 1e-9 ; X = T*(MAXCH*16+MAXICH*16)
	if(say)write(48,*) ' Allocate SRC/PSI etc in ',real(X),' GB'
      	 allocate(SRC(MAXN,MAXCH),PSI(MAXN,MAXICH),SRATIO(MAXCH))
	 allocate(SOMA(MAXCH,2),SPAD(MAXCH,MAXIT+2,MAXIT+2))
      	 allocate(ferw(NFDEC))

	 if(CCREAL.and..not.FALLOC) then    ! only for ERWINRC4, not yet ERWINCC4
	    IF = MAXCH
	    if(MINTL==1) IF = 1   ! reduced FIM  if ERWINs called only once
            allocate(FIMR(NBL,IF,numthread),FIMDR(NBL,MAXCH,numthread))
       	     T = NBL*numthread * (IF+MAXCH)
	   else  		  ! normal
            allocate(FIM(MAXB,MAXCH),FIMD(MAXB,MAXCH))
       	     T = MAXB* MAXCH*2 * 2
	   endif
        if(say)write(48,*) ' Allocate FIM(D) in ',real(T*8e-9),' GB ',T

      	 allocate(SIMPLE(MAXCH),EMPTY(MAXCH),SKIP(MAXCH,2))
	else
	 allocate(PSI(1,1))
	endif  ! not rterms
        
	if(flip) then
	       T = MAXCH*MAXCH*8e-9
		if(say)write(48,*) ' Allocate EXCH in ',real(T),' GB'
	   allocate(EXCH(MAXCH,MAXCH))
	  else
	   allocate(EXCH(1,1))
	  endif
!       write(6,*) ' Alloc ',allocated(CFMAT),allocated(CGMAT),
!    x			allocated(CFF),':',FCWFN,cfalloc
        if(FCWFN) allocate(CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2))
        if(.not.FCWFN) allocate(CFMAT(1,1,2),CGMAT(1,1,2))
	if(cfalloc) then
	   T = (MAXCH*2)**2 * 16e-09
	   if(T>0.10) write(KO,'(a,f8.3,a)') ' Allocate CFF in ',T,' GB'
	   allocate(CFF(MAXCH,MAXCH,IPMAX))
	   endif

!	else
!	write(KO,*) 'Have arrays for ',MAXCH,' channels, of which ',
!     X	 MAXICH,' need wfs.'
	endif   ! need to reallocate
!#	write(48,*) 'Makeset: NFDEC = ',NFDEC,N,MAXB,NICH,' proceed'

      CALL CCSET(JTOTAL,PARITY,ETOTAL,JTMIN,KINTL,LPMAX,
     X        NEX,NCHAN,GIVEXS,QVAL,ENEX,PEL,EXL,LMAX,JEX,ECM,ECMC,
     X        LVAL,JVAL,PART,EXCIT,COPY,JPROJ,JTARG,CUTL,CUTR,CUTOFF,
     X        HP,RMASS,INCOME,BLOCKD,LUSED,MAL1,LJMAX,GAM,
     X        MINTL,INITL,ITC,IBLOCK,IEX,NCH,NCHPART,CHBASE,
     X        K,ETA,CHL,NAME,BAND,N,ISOCEN,RTURN,SMALLS,NSMALL,DROPPED,
     x        ENLAB,elpmax)
!		write(6,*)  ' CCSET: IEX,MAXB =',IEX,MAXB
!         WRITE(48,*) JTOTAL,PSIGN(PARITY+2),'CCSET:',NCH,MINTL,IEX
C
C    READ IN ANY J/P-DEPENDENT POTENTIALS
C    ----------------------------------
      IB = INT(JTOTAL) + 1
      DO 240 JF=1,NF0
        JFT = PTYPE(2,JF)
          NFL = PTYPE(5,JF)
        IF(JFT.LE.7 .AND. NFL.GT.0.and.NFL<30) THEN
	  OPEN(NFL,ACCESS='DIRECT',FORM='FORMATTED',RECL=72,
     X			STATUS='OLD',file='potent25')
          IF(NFL.GE.24) IB = INT(JTOTAL) + 1
          IF(NFL.LT.24) IB = IPARIT-1
         NR = (N-1)/3 + 1
         IR0 = (IB-1)*NR
       WRITE(KO,*) 'Read in form at ',JF,' for block ',
     X               IB, ' with records from card ',IR0+1
         DO 230 IR=1,NR
           I0 = (IR-1)*3
230       READ(NFL,REC=IR+IR0,FMT=232)
     X                  (FORMF(I+I0,JF),I=1,MIN(3,N-I0))
232    FORMAT(6E12.4)
        IF(TRENEG.GT.0) WRITE(KO,235) ((I-1)*HCM,FORMF(I,JF),I=1,M,MR)
235         FORMAT(5(1X,F5.1,':',2F9.4,1X))
	close(NFL)
       ENDIF
240   CONTINUE
C
!     CLIST(:,:,:) = 0.0   ! zero-init not necessary (and time-consuming for big cases!)
      NCLIST(:,:) = 0	   ! this DOES need to be set to zero at start!
      NL = 0
C
C    OPTICAL POTENTIALS FOR EACH CHANNEL & TENSORS BETWEEN CHANNELS
C    --------------------------------------------------------------
!	EPS=-1		! DEBUG ONLY
      ICH = 0
           IF(LISTCC.GE.1) WRITE(48,*) JTOTAL,PSIGN(PARITY+2)
!           WRITE(480,*) 
!           WRITE(480,*) ENLAB,JTOTAL,PSIGN(PARITY+2)
      DO 325 IC=1,MXP
	NEXK = NCHPART(IC)
	ETATG = (MASS(2,IC)-2*MASS(4,IC))/MASS(2,IC)
	ETAPR = (MASS(1,IC)-2*MASS(3,IC))/MASS(1,IC)
	ETAS = - ETATG * ETAPR
	I0 = CHBASE(IC)
	 repeat = .false.
         DO 302 C=I0+1,I0+NEXK
         IA = EXCIT(C,1)
         KPA = CPOT(IC,IA)
	 WPA = 1.0

	   if(J40P) then
	   do JF=1,NF0
	   if (PTYPE(1,JF)==KPA) then
	   if (PTYPE(5,JF)==40) then
	     IB = IPARIT-1			! 1 for IPARIT=1 (+), 2 for IPARIT=3 (-1) as PARITY=(-1)**IPARIT
	     KPA = nint(FORMDEF(IB,JF))   ! parity-dependent choice of KPA
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C,' uses parity-dep potential: ',KPA
	   endif
	   if (PTYPE(5,JF)==41) then
	     L = LVAL(C)
	     KPA = nint(FORMDEF(min(L,6)+1,JF))   ! L-dependent choice of KPA
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C,' uses L-dep potential: ',KPA
	   endif
	   if (PTYPE(5,JF)==42) then
	     L = int(JTOTAL)
	     KPA = nint(FORMDEF(min(L,6)+1,JF))   ! J-dependent choice of KPA
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C,' uses J-dep potential: ',KPA
	   endif
	   if (PTYPE(5,JF)==43) then
	     L = LVAL(C)
	     IB=MOD(L,2)+1 		 ! 1 for even L, 2 for odd L
	     KPA = nint(FORMDEF(IB,JF))   ! L-parity-dependent choice of KPA
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C,' uses L-parity-dep potential: ',KPA
	   endif
           if (PTYPE(5,JF)==45) then !            E-dependent choice of KPA
	    if(ENLAB<FORMDEF(2,JF)) then
	     KPA=nint(FORMDEF(1,JF)); KPA2=0; WPA=1.; !WPA2=0.
	     endif
	    J = 0
 	    do I=4,12,2
     	      if(FORMDEF(I-3,JF)>0.and.FORMDEF(I-1,JF)>0) then
	  	J = I
	    if(ENLAB<FORMDEF(I,JF).and.ENLAB>=FORMDEF(I-2,JF)) then
	      KPA = nint(FORMDEF(I-3,JF)); KPA2 = nint(FORMDEF(I-1,JF))
	      T = (ENLAB-FORMDEF(I-2,JF))
     x              /(FORMDEF(I,JF)-FORMDEF(I-2,JF))
	      WPA  = 1. - T   !  weight for lower energy potential
	     endif
	     endif  ! KPs>0
	    enddo
	    if(J>1) then
	    if(ENLAB>=FORMDEF(J,JF)) then  ! beyond last energy listed with KP>0
	     KPA=nint(FORMDEF(J-1,JF)); KPA2=0; WPA=1.; !WPA2=0.
	     endif
	     endif ! J>1
           IF(LISTCC.GE.1) then
           write(48,*) ' Channel ',C,' uses E-dep potentials: ',
     x           KPA,KPA2,' with weights ',real(WPA),real(1.-WPA)
           write(6,*) ' Channel ',C,' uses E-dep potential: ',KPA,
     x                 ' with weight ',real(WPA), ' A'
	   endif
           endif


	   endif   ! KPA
	   enddo   ! JF
	   endif ! J40P

!@         C2I = I0+1; C2F = I0+NEXK			!@ fixed March 2014: new DO 302 and MAXLAM0 condition
!@         if(MAXLAM0==0) then  ! monopoles only !!
!@           C2I = C; C2F = C
!@           endif
!	 write(48,*) ' For C=',C,' the C2 goes from ',C2I,' to ',C2F

!@        DO 302 C2=C2I,C2F
         DO 302 C2=I0+1,I0+NEXK
	  if(abs(LVAL(C)-LVAL(C2))>MAXLAM0) go to 302

            IB = EXCIT(C2,1)
            KPB = CPOT(IC,IB)
	    WPB = 1.

	   if(J40P) then
	   do JF=1,NF0
	   if (PTYPE(1,JF)==KPB) then
!	   IF(LISTCC.GE.0) WRITE(KO,*) C,'@@',JF,':',(PTYPE(II,JF),II=1,6)
	   if (PTYPE(5,JF)==40) then
	     KPB = nint(FORMDEF(IPARIT-1,JF))   ! parity-dependent choice of KPB
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C2,' uses parity-dep potential: ',KPB
	   endif
	   if (PTYPE(5,JF)==41) then
	     L = LVAL(C2)
	     KPB = nint(FORMDEF(min(L,6)+1,JF))   ! L-dependent choice of KPB
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C2,' uses L-dep potential: ',KPB
	   endif
	   if (PTYPE(5,JF)==42) then
	     L = int(JTOTAL)
	     KPB = nint(FORMDEF(min(L,6)+1,JF))   ! J-dependent choice of KPB
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C2,' uses J-dep potential: ',KPB
	   endif
	   if (PTYPE(5,JF)==43) then
	     L = LVAL(C2)
	     KPB = nint(FORMDEF(mod(L,2)+1,JF))   ! L-parity-dependent choice of KPB
           IF(LISTCC.GE.1)
     x	   write(48,*) ' Channel ',C2,' uses L-parity-dep potential: ',KPB
	   endif
           if (PTYPE(5,JF)==45) then !            ! E-dependent choice of KPB
            if(ENLAB<FORMDEF(2,JF)) then
             KPB=nint(FORMDEF(1,JF)); KPB2=0; WPB =1.; ! WPB2=0
             endif
            do I=4,12,2
     	      if(FORMDEF(I-3,JF)>0.and.FORMDEF(I-1,JF)>0) then
	      J = I
            if(ENLAB<FORMDEF(I,JF).and.ENLAB>=FORMDEF(I-2,JF)) then
              KPB2 = nint(FORMDEF(I-3,JF)); KPB = nint(FORMDEF(I-1,JF))
              WPB = (ENLAB-FORMDEF(I-2,JF))
     x              /(FORMDEF(I,JF)-FORMDEF(I-2,JF))  ! weight for upper-energy potential
              !WPB2  = 1. - WPB
             endif
             endif ! KPs>0
            enddo
	    if(J>1) then ! found some last energy

            if(ENLAB>=FORMDEF(J,JF)) then  ! beyond last energy listed with KP>0
             KPB=nint(FORMDEF(J-1,JF)); KPB2=0; WPB=1.; !WPB2=0.
             endif
             endif ! J>1
           IF(LISTCC.GE.1) then
           write(48,*) ' Channel ',C,' uses E-dep potentials: ',
     x                 KPB,KPB2,' with weights ',real(1-WPB),real(WPB)
           write(6,*) ' Channel ',C,' uses E-dep potential: ',KPB,
     x                 ' with weight ',real(WPB), ' B'
	   endif
           endif

	   endif
	   enddo
	   endif ! J40P

!         if(KPA/=KPB) go to 302     !! ALLOW COUPLINGS BETWEEN DIFFERENT KP VALUES!! THEY ADD
!					 This feature is required for the E-dependent linear interpolation

         DO 300 JF=1,NF0
!       IF(LISTCC.GE.4) WRITE(KO,*) C,C2,'@@@',JF,':',
!    x   		(PTYPE(II,JF),II=1,6),',',KPA,KPB
            IF(PTYPE(3,JF).LT.0 .OR. PTYPE(4,JF).LT.0)    GO TO 300

           OP = .true.
           if(KPA /= KPB) then          ! how to deal with this case:
            if(MIXPOT(IC)==0) then      !  no couplings 
              OP =  .false.

	    else if(MIXPOT(IC)==1) then ! couplings only to or from the first state in that partition
              OP = IA==1 .or. IB==1

	    else if(MIXPOT(IC)==2) then ! couplings only according to KP value of the destination level
              OP = PTYPE(1,JF) == KPA

	    else if(MIXPOT(IC)==3) then ! all couplings produced by either KP value
              OP = .true.
            endif
           endif
            if(listcc>6 .or. listcc>4.and.OP) write(6,4801) 
     #       C,C2,MIXPOT(IC),PTYPE(1,JF),KPA,KPB,OP
4801	    format('C,C2=',2i4,' MIXPOT =',I2,' potl',I3,' KPA,B =',2i3,
     #          ' => OP:',L1)
 
           if( .not. OP) go to 300

! OLD coding (ignore!)
!           IF(.not.(PTYPE(1,JF)==KPA .or. 
!    x               PTYPE(1,JF)==KPB.or.abs(WPB)>1e-8) )  GO TO 300
!            IF(PTYPE(1,JF).NE.KPA.and.abs(WPB)<1e-8) go to 300   !! WHEN KPA/=KPB, potential effective that is listed for destination
!            IF(PTYPE(1,JF).NE.KPA) go to 300   !! WHEN KPA/=KPB, potential effective that is listed for destination

	    IF(PTYPE(2,JF).le.3 .and. C/=C2) then
      		IF(LISTCC.GE.7) write(6,*) 'Skip as ',PTYPE(2,JF),' <= 3'
		go to 300 ! scalar => no off-diagonal for central vol, deriv & projectile vso
		endif
	    WPOT = 0.
	    if(PTYPE(1,JF)==KPA) WPOT = WPA
	    if(PTYPE(1,JF)==KPB) WPOT = WPB
!            write(480,481) JF,PTYPE(1,JF),KPA, WPA,KPB,WPB
!481	 format(2i4,":",i4,f8.3,i4,f8.3)
            if(abs(WPOT)<1e-20) go to 300 

            JFT = PTYPE(2,JF)
	    JFTT = JFT
	    NFD = PTYPE(4,JF)
!       IF(LISTCC.GE.5) WRITE(KO,*) 'Try 1 form @',JF,':',JFT
	   SCALEL = 1d0
	   if(LDEP(JF)>0) then
	     T = JTOTAL
	     if(LDEP(JF)==2) T = (LVAL(C)+LVAL(C2))/2
	     RH = (T- VARYL(1,JF))/VARYL(2,JF)
	     EE = DDEXP(-RH)
	     if(LSHAPE(JF)==0) then
	 	SCALEL = 1d0/(1d0 + 1d0/EE)
	     else if(LSHAPE(JF)==1) then
	 	SCALEL = 1d0/(1d0 + 1d0/EE)**2
	     else if(LSHAPE(JF)==2) then
	 	SCALEL = DDEXP(-RH*RH)
	     endif
	    if(LISTCC>50) write(KO,*) 'JL-dependence for form ',JF,
     X		': LDEP,VARY,T,RH,EE,SCALEL=',
     X             LDEP(JF),VARYL(1:2,JF),T,RH,EE,SCALEL
	   endif
      IF(EXCIT(C,3)/=EXCIT(C2,3) .AND.(JFT==10.OR.JFT==12.OR.JFT==14)
     X .OR.
     X   EXCIT(C,2)/=EXCIT(C2,2) .AND.(JFT==11.OR.JFT==13))
     X  GO TO 300
C**********************************************************************
      IF(JFT.GE.10 .AND. JFT.LE.15 .AND.
     X   (EXCIT(C,2).NE.EXCIT(C2,2) .AND. EXCIT(C,3).NE.EXCIT(C2,3)))
     X  GO TO 300
      IF(JFT.LT.10 .AND.
     X   (EXCIT(C,2).NE.EXCIT(C2,2) .OR. EXCIT(C,3).NE.EXCIT(C2,3)))
     X  GO TO 300

       T = 1.0
	 if(JFT==0) then
           if(IC.eq.PEL.and.IA.eq.EXL) then
            if(mod(PLANE,2)==1) T = 0.  ! elastic channel
           else
            if(    PLANE/2 ==1) T = 0.  ! nonelastic channels
           endif
	 endif

       KK = PTYPE(3,JF)
       IF(KK.eq.7) THEN
C                       Take K=7 as an inelastic monopole
          IF(IA.eq.IB) GO TO 300 ! exclude diagonal monopole
          KK = 0
       ELSE if(KK==0) then
C			Normal KK=0 should exclude off-diagonal monopoles
C                       BUT: Allow off-diagonal couplings for spin.spin monopoles
          if(C/=C2.and.JFT/=8.and.JFT/=4) GO TO 300
            if(JFT==12) T = sqrt(2*JPROJ(C)+1.)   ! projectile couplings
            if(JFT==13) T = sqrt(2*JTARG(C)+1.)   ! target couplings
       ENDIF
C
        KP = (-1)**KK   ! allow parity change for odd multipoles
        if (JFTT.eq.10.and.BAND(1,ic,ia).ne.BAND(1,ic,ib)*KP .or.
     x      JFTT.eq.11.and.BAND(2,ic,ia).ne.BAND(2,ic,ib)*KP ) go to 300

C
      IF(JFT>=12.and.JFT<=13.and.PTYPE(3,JF)>0) THEN
C                           LOOK UP TABLE OF ALLOWED COUPLINGS for KK>0 or off-diagonal monopole (7)
         T = 0.0
         DO 294 IX=1,NIX
           IF(MATRIX(1,IX).NE.EXCIT(C,JFT-10) .OR.
     X        MATRIX(2,IX).NE.EXCIT(C2,JFT-10) .OR.
     X        MATRIX(3,IX).NE.PTYPE(3,JF) .OR.
     X        MATRIX(4,IX).NE.JF) GO TO 294
         T = T + MEK(IX)
294      CONTINUE
C
      ELSE IF(JFT>=14.and.JFT<=16) THEN		! Deal with later
         GO TO 302

      ELSE IF(JFT.EQ.17) THEN
	 if(PTYPE(3,JF)/=IC) go to 300
C		This last condition allows only the *first* JF slot.
C                           LOOK UP TABLE OF ALLOWED COUPLINGS :
         T = 0.0
         DO 298 IX=1,NIX
           IF(MATRIX(1,IX).NE.EXCIT(C,1) .OR.
     X        MATRIX(2,IX).NE.EXCIT(C2,1).OR.
     X        MATRIX(5,IX).NE.NFD) GO TO 298
         T = T + MEK(IX)
         KK = MATRIX(3,IX)
	 if(KK==7) KK=0
	 JFTT = 0
         IF(EXCIT(C,3).NE.EXCIT(C2,3)) JFTT=13
         IF(EXCIT(C,2).NE.EXCIT(C2,2)) JFTT=12
	 if(JFTT==0) JFTT = MATRIX(6,IX)+11
298      ENDDO
	if(.not.repeat) then
	 allocate (OO(NEXK,NEXK),EVEC(NEXK,NEXK))
	 OO(:,:) = 0d0
	endif
         JF0 = JF
	 repeat = .true.
       ENDIF
C
          S =  TENSOR(JFTT,LVAL(C),JPROJ(C),JVAL(C),JTARG(C),
     X                JTOTAL,LVAL(C2),JPROJ(C2),JVAL(C2),JTARG(C2),
     X                KK,JEX(3,IC,IA),JEX(4,IC,IA),
     X                 JEX(3,IC,IB),JEX(4,IC,IB),
     X           MASS(3,IC),MASS(4,IC),PTYPE(5,JF),ETAS)
	  S = S
C    X      * CI**NINT(JPROJ(C2)+JTARG(C2) - JPROJ(C) - JTARG(C))
     X      * CI**NINT( - ABS(JPROJ(C2)-JPROJ(C))
     X                   +     JPROJ(C2)-JPROJ(C)
     X                   - ABS(JTARG(C2)-JTARG(C))
     X                   +     JTARG(C2)-JTARG(C))
C           The above phase factors with JPROJ & JTARG etc.,
C           are there only because of definition of M(Ek) matrix element
!	if(JFT<=11.or.JFT>=18) 
	if(JFT/=17) S = S * CI**(LVAL(C2)-LVAL(C)) * SCALEL *  WPOT
	if(repeat.and.JFT==17) then
C					Use couplings as deformation lengths
	 OO(C-I0,C2-I0) = OO(C-I0,C2-I0) + T*REAL(S)
!	 OO(C2-I0,C-I0) = OO(C2-I0,C-I0) + T*REAL(S)
           IF(LISTCC.GT.1.and.abs(T)>EPS) WRITE(KO,1330) C,C2,JF,
     X	(PTYPE(II,JF),II=1,6),KK,LDEP(JF),JFTT,T,S,SCALEL
!         IF(LISTCC>2.and.abs(T)>EPS) write(KO,*) OO(C-I0,C2-I0),C,C2,I0
!          IF(LISTCC>2.and.abs(T)>EPS) write(KO,*) OO
	endif
	   NC = NCLIST(C,C2)+1
	if(abs(S*T)>EPS.and.JFT/=17) then
	   if(NC>MCLIST) then
             write(KO,*) 'For channels ',C,C2
             write(KO,*) 'Need NC=',NC,' > MCLIST=',MCLIST,' !!!'
             write(KO,*)'So far use forms ',NFLIST(C,C2,:),' need ',JF
             call check(NC,MCLIST,30)
             stop
             endif
	   CLIST(C,C2,NC) = S*T
	   NFLIST(C,C2,NC) = JF
	   NCLIST(C,C2) =  NC
           IF(C.NE.C2.AND.S*T.NE.0.0) ICH = MAX(ICH,C2)
           IF(LISTCC.GT.1) WRITE(KO,1330) C,C2,JF,
     X    (PTYPE(II,JF),II=1,6),KK,LDEP(JF),NC,T,S,SCALEL,WPOT
 1330      FORMAT(' For',I5,' from',I5,' by form',I5,' (',6I3,')',2i3,
     X	    i10,' get',F8.4,2F15.6,3F8.4)
	 endif
  300    CONTINUE
  302  CONTINUE
!            write(6,*) 'Finished TENSOR 302 @ ',real(TM(I))

! 					Deal with 14,15,16 now
	if(J14T16) then
         DO 308 C=I0+1,I0+NEXK
         IA = EXCIT(C,1)
         KP = CPOT(IC,IA)
         DO 308 C2=I0+1,I0+NEXK
            IB = EXCIT(C2,1)
            KPB = CPOT(IC,IB)
         DO 307 JF=1,NF0
!       IF(LISTCC.GE.2) WRITE(KO,*) 'To ',C,' from  ',C2,', ',
!     X	' TRY ',JF,' @ ',(PTYPE(II,JF),II=1,6)
            IF((PTYPE(1,JF).NE.KP .and. PTYPE(1,JF).ne.KPB) .OR.
     X         PTYPE(3,JF).LT.0 .OR. PTYPE(4,JF).LT.0)    GO TO 307
            JFT = PTYPE(2,JF)
       IF(JFT<14.or.JFT>17) GO TO 307		! Deal with 14,15,16 only
	    NFD = PTYPE(4,JF)
	   SCALEL = 1d0
	   if(LDEP(JF)>0) then
	     T = JTOTAL
	     if(LDEP(JF)==2) T = (LVAL(C)+LVAL(C2))/2
	     RH = (T- VARYL(1,JF))/VARYL(2,JF)
	     EE = DDEXP(-RH)
	     if(LSHAPE(JF)==0) then
	 	SCALEL = 1d0/(1d0 + 1d0/EE)
	     else if(LSHAPE(JF)==1) then
	 	SCALEL = 1d0/(1d0 + 1d0/EE)**2
	     else if(LSHAPE(JF)==2) then
	 	SCALEL = DDEXP(-RH*RH)
	     endif
	    if(LISTCC>50) write(KO,*) 'JL-dependence for form ',JF,
     X		': LDEP,VARY,T,RH,EE,SCALEL=',
     X             LDEP(JF),VARYL(1:2,JF),T,RH,EE,SCALEL
	   endif
!      IF(LISTCC.GE.2) WRITE(KO,*) ' TRY2 ',JF,' for ',
!     X   EXCIT(C,3),EXCIT(C2,3),EXCIT(C,2),EXCIT(C2,2)
      IF((EXCIT(C,3)/=EXCIT(C2,3) .AND.JFT==14)  .OR.
     X   (EXCIT(C,2)/=EXCIT(C2,2) .AND.JFT==15)) GO TO 307
!      IF(LISTCC.GE.2) WRITE(KO,*) ' TRY3 ',JF,' * ',SCALEL
C**********************************************************************
!?      IF(JFT==16 .AND.
!?     X   (EXCIT(C,2).eq.EXCIT(C2,2) .or. EXCIT(C,3).eq.EXCIT(C2,3)))
!?     X  GO TO 307

       KK = PTYPE(3,JF)
       if(KK==7) KK=0
C	 I1 and I2 are 1(proj) or 2(targ) for the two transitions
       IF(JFT>=14.and.JFT<=15) THEN		
	 I1 = JFT-13			! Both = 1 for 14, 2 for 15
	 I2 = JFT-13
       ELSE IF(JFT==16) THEN
	 I1=1				! Projectile
	 I1=2				! Target
       ENDIF
!       IF(LISTCC.GE.2) WRITE(KO,*) ' TRY4 ',JF,KK,I1,I2,NFD
C                           LOOK UP TABLE OF ALLOWED COUPLINGS :
         DO 305 IX1=1,NIX
         DO 305 C1=I0+1,I0+NEXK
         DO 305 IX2=1,NIX
           LAM1 = +MATRIX(3,IX1)
           LAM2 = -MATRIX(3,IX2)
!       IF(LISTCC.GE.2) WRITE(KO,*) ' TRY5 ',IX1,IX2,' @ ',LAM1,LAM2,C1
           IF(JFT>=14.and.JFT<=15) LAM1=abs(LAM1)

           IF(MATRIX(1,IX1).NE.EXCIT(C,I1+1) .OR.
     X        MATRIX(2,IX1).NE.EXCIT(C1,I1+1) .OR.
     X        LAM1<0           .OR.
     X        MATRIX(5,IX1).NE.NFD) GO TO 305
         TH1 = MEK(IX1)

           LAM2 = -MATRIX(3,IX2)
           IF(JFT>=14.and.JFT<=15) LAM2=abs(LAM2)
           IF(MATRIX(1,IX2).NE.EXCIT(C1,I2+1) .OR.
     X        MATRIX(2,IX2).NE.EXCIT(C2,I2+1) .OR.
     X        LAM2<0           .OR.
     X        MATRIX(5,IX2).NE.NFD) GO TO 305
         TH2 = MEK(IX2)

	 IF(KK<abs(LAM1-LAM2) .or. KK>LAM1+LAM2) go to 305
!       IF(LISTCC.GE.2) WRITE(KO,*) ' Found ',IX1,IX2,' @ ',LAM1,LAM2,C1
	 T = TH1*TH2 * 0.5 			!0.5 from 1/2!
C
          S =  TENSOR2(JFT,LVAL(C),JPROJ(C),JVAL(C),JTARG(C),
     X                JTOTAL,LVAL(C2),JPROJ(C2),JVAL(C2),JTARG(C2),
     X                KK,LAM1,LAM2,JEX(3,IC,IA),JEX(4,IC,IA),
     X                 JEX(3,IC,IB),JEX(4,IC,IB),
     X                 MASS(3,IC),MASS(4,IC),PTYPE(5,JF))
	  S = S
     X      * CI**NINT( - ABS(JPROJ(C2)-JPROJ(C))
     X                   +     JPROJ(C2)-JPROJ(C)
     X                   - ABS(JTARG(C2)-JTARG(C))
     X                   +     JTARG(C2)-JTARG(C))
C           The above phase factors with JPROJ & JTARG etc.,
C           are there only because of definition of M(Ek) matrix element
	S = S * CI**(LVAL(C2)-LVAL(C)) * SCALEL

	   NC = NCLIST(C,C2)+1
	if(abs(S*T)>EPS*0.) then
	   if(NC>MCLIST) then
		write(KO,*) 'For channels ',C,C2
		write(KO,*) 'Need NC=',NC,' > MCLIST=',MCLIST,' !!!'
		write(KO,*) 'So far forms ',NFLIST(C,C2,:),'; need ',JF
		call check(NC,MCLIST,30)
		stop
		endif
	   CLIST(C,C2,NC) = S*T
	   NFLIST(C,C2,NC) = JF
	   NCLIST(C,C2) =  NC
           IF(C.NE.C2.AND.S*T.NE.0.0) ICH = MAX(ICH,C2)
           IF(LISTCC.GT.1) WRITE(KO,1331) C,C1,C2,JF,
     X    (PTYPE(II,JF),II=1,6),KK,LDEP(JF),NC,TH1,TH2,S,SCALEL
 1331      FORMAT(' For',I5,' via',I5,' from',I5,' by form',I5,
     X		' (',6I3,')',3i3,' get',2F8.4,2F15.6,F8.4)
	 endif
  305    CONTINUE
  307    CONTINUE
  308  CONTINUE
	endif ! J14T16

	if(repeat) then  
           IF(LISTCC.GT.2) then
               write(KO,*) ' OO matrix:'
              do i=1,NEXK
               write(KO,3101) i,(OO(i,j) , j=1,NEXK)
	      enddo
 3101          format(1x,i3,12f10.3,:,/(4x,12f10.3))
	     endif
C
C    ALL-ORDER COUPLING POTENTIALS BETWEN CHANNELS
C    ---------------------------------------------
	call HDIAG(OO,NEXK,NEXK,0,EVEC,I)

           IF(LISTCC.GT.1) WRITE(KO,*) ' Diagonalised matrix',NEXK,I
           IF(LISTCC.GT.2) then
	  	write(KO,311)(OO(IE,IE),IE=1,NEXK)
311		format('  Eigenvalues are:'/(1x,10f8.4))
	 	endif
	IF = JF0
        do 322 ia=1,NEXK
	do 322 ib=ia,NEXK
	  c = ia+I0
	  c2= ib+I0
!	if(LISTCC>1) write(KO,312) c,c2,IF,NFD
!312	format(30X,'  All-order coupling ',i3,' <->',i3,' is at ',i3,
!    X		', deforming potential at ',i3)
               CALL CHECK(IF,MLOC,24)
C					The monopole is replaced completely!
!               PTYPE(4,NFD) = -JF0
	 FORMF(1:N,IF) = 0d0
	  do 320 I=1,N
	   R1 = (I-1)*HP(IC)
	   do 320 IE=1,NEXK
	   X = OO(IE,IE)*R4PI
	   R = R1 - X
	   S = FFC(R/HP(IC),FORMF(1,NFD),N)
            FORMF(I,IF)=FORMF(I,IF) + S * EVEC(IE,IA)*EVEC(IE,IB)
320	   continue
	L = 1
	if(C/=C2) L=2
C				Do forward (IE=1) and reverse (IE=2) couplings
	Do 321 IE=1,L
	if(IE==2) then
	 I = C
	 C = C2
	 C2 = I
	 endif
C						Put in phase factors:
	  S =     CI**(LVAL(C2)-LVAL(C))
	  
	   NC = NCLIST(C,C2)+1
	   if(NC>MCLIST) then
		write(KO,*) 'For channels ',C,C2
		write(KO,*) 'Need NC=',NC,'>MCLIST=',MCLIST,' !!!'
		write(KO,*) 'So far use forms ',NFLIST(C,C2,:),' need',IF
		call check(NC,MCLIST,30)
		stop
		endif
	   CLIST(C,C2,NC) = S
	   NFLIST(C,C2,NC) = IF
	   NCLIST(C,C2) =  NC
           IF(C.NE.C2) ICH = MAX(ICH,C2)
           IF(LISTCC.GT.1) 
     X	WRITE(KO,1330) C,C2,IF,(PTYPE(II,IF),II=1,6),NC,0,0,1.,S
321	continue
      IF(TRENEG*LISTCC.GE.1) THEN
            WRITE(KO,110) IF
110         FORMAT(/' Potential Form at',I4,' is')
            WRITE(KO,120) ((I-1)*HP(IC),FORMF(I,IF),I=1,N,MR)
120         FORMAT(5(1X,F5.1,':',2F9.4,1X))
           ENDIF
	 IF = IF+1
322	 continue
	deallocate(OO,EVEC)
	endif
	 
  325  CONTINUE
!	write(48,*) ' lanecoup =',lanecoup
	if(lanecoup) then
           IF(LISTCC.GT.1) write(6,*) 'LANE COUPLINGS: '
      DO 1325 IC1=1,MXP
	NEXK = NCHPART(IC1)
	ETATG = (MASS(2,IC1)-2*MASS(4,IC1))/MASS(2,IC1)
	ETAPR = (MASS(1,IC1)-2*MASS(3,IC1))/MASS(1,IC1)
	ETAS = - ETATG * ETAPR
	 T = 2*sqrt(abs(ETAS)/MASS(2,IC1))   ! coefficient for charge-exchange coupling
	I0 = CHBASE(IC1)
	 if(say)write(48,*) ' Partition ',IC1,' has chs ',I0+1,I0+NEXK
         DO 1322 C=I0+1,I0+NEXK
         IA = EXCIT(C,1)
         KPA = CPOT(IC1,IA)
	 if(say)write(48,*) ' Look for Lane couplings IC,IA =',IC1,IA, 
     x 				' with KP =',KPA
	   do JF=1,NF0
	   if (PTYPE(1,JF)==KPA) then
            JFTT = PTYPE(2,JF)
	   if (JFTT==26.or.JFTT==27) then   ! Vol or Surface Lane form factors
                KK = PTYPE(3,JF)
       		if(KK==7) KK=0
	 	if(say)write(48,*) ' Found Lane couplings IC,IA =',IC1,IA, 
     x 				' with KP =',KPA,', T,K =',JFTT,KK
	DO 1320 IC2=1,MXP   
	 IB = IA 		!  assume identical IA values coupled by the Lane coupling
	 KPA2= CPOT(IC2,IA)   
	 if(KPA==KPA2 .and.IC2.ne.IC1) then    ! only off-diagonal couplings here: others already done
	 	if(say)write(48,*) ' Match Lane couplings IC2,IB =',
     x                             IC2,IB, ' cf ',IC1
	  I2 = CHBASE(IC2)
	  do 1312 C2=I2+1,I2+NCHPART(IC2)
	  if(ABS(JPROJ(C2)-JPROJ(C))>.1 .or. ABS(JTARG(C2)-JTARG(C))>.1
     x   .or.LVAL(C2).ne.LVAL(C)) go to 1312 			! no change of projectile or target spins!

          S =  TENSOR(JFTT,LVAL(C),JPROJ(C),JVAL(C),JTARG(C),
     X                JTOTAL,LVAL(C2),JPROJ(C2),JVAL(C2),JTARG(C2),
     X                KK,JEX(3,IC1,IA),JEX(4,IC1,IA),
     X                 JEX(3,IC2,IB),JEX(4,IC2,IB),
     X          MASS(3,IC1),MASS(4,IC1),PTYPE(5,JF),T)
          S = S * CI**NINT( - ABS(JPROJ(C2)-JPROJ(C))
     X                   +     JPROJ(C2)-JPROJ(C)
     X                   - ABS(JTARG(C2)-JTARG(C))
     X                   +     JTARG(C2)-JTARG(C))
C           The above phase factors with JPROJ & JTARG etc.,
C           are there only because of definition of M(Ek) matrix element
        S = S * CI**(LVAL(C2)-LVAL(C))   !* SCALEL

          NC = NCLIST(C,C2)+1
           if(NC>MCLIST) then
                write(KO,*) 'For channels ',C,C2
                write(KO,*) 'Need NC=',NC,'>MCLIST=',MCLIST,' !!!'
                write(KO,*) 'So far use forms ',NFLIST(C,C2,:),' need',IF
                call check(NC,MCLIST,30)
                stop
                endif
           CLIST(C,C2,NC) = S
           NFLIST(C,C2,NC) = JF
           NCLIST(C,C2) =  NC
           IF(C.NE.C2) ICH = MAX(ICH,C2)
           IF(LISTCC.GT.1)
     X  WRITE(KO,1330) C,C2,JF,(PTYPE(II,JF),II=1,6),NC,0,0,1.,S

1312 	CONTINUE   ! C2

	 endif     ! KPA==KPA2
1320 	CONTINUE   ! IC2
	endif
	endif
	enddo      ! JF
1322 	CONTINUE   ! C
1325 	CONTINUE   ! IC1
	endif ! lanecoup
!	    write(6,*) 'Finished TENSOR @ ',real(TM(I))

C
C    COUPLING COEFFICIENTS BETWEEN PAIRS OF PARTIAL WAVES
C    ---------------------------------------------------
      if(OPEN12) REWIND 12
      IF(REW14) REWIND 14
      DO 290 CP=1,NCP
            IF(.NOT.COUPLE(CP).OR.JTOTAL.GT.JMAX(CP)) GO TO 290
            IF(.NOT.LOCAL(CP).AND.ITER==0.and..not.rterms) GO TO 290
	    IF(JTMIN.LT.Z .AND. JTOTAL.LE.-JTMIN-.1) GO TO 290
         IC1 = ICTO(CP)
         IC2 = ICFROM(CP)
            IN1 = MOD(ABS(KIND(CP))-1,2)+1
         IF(LISTCC.GT.2) WRITE(KO,1329) CP,ICTO(CP),ICFROM(CP),KIND(CP)
     x                  ,NLL(CP),NLN
 1329  FORMAT(/' Coupling #',I4,' to',I4,' from',I4,' of KIND=',i2,2i4)

      CALL CPAIR(CP,IC1,IC2,KIND(CP),QQ(CP),KPCORE(CP),REV(CP),LOCAL(CP)
     X          ,ICOM,ICOR,LOCF(CP),NKP(1,CP),GPT(1,1,CP),XA(CP),XB(CP),
     X          XP(CP),XQ(CP),IREM(CP),FILE(CP),CHNO,ffreal(CP),KLT(CP),
     X    NCH,LVAL,JVAL,JPROJ,JTARG,JTOTAL,PART,EXCIT,JEX,COPY,BAND,
     X    QNF,AFRAC,CCFRAC,CUTOFF,ALOSSM,ACC8,FPT(1,1,CP),NLL(CP),
     X    NKP(1,CP),BE(1,1),CLIST,NFLIST,NCLIST,DNL,EPS,NFORML1,NFORML2,
     X    HP,ITC,NEX(ICOR(CP,IN1)),ICFROM(CP),ICTO(CP),MATRIX,MEK,
     X    CPSO(CP),SPINTR,K,POTCAP(1,1,CP),LTRANS(1,CP),NLOC(1,CP),
     X    MLCALL(1,CP),MLCALL(2,CP),FORML,VPOT,VPOT(1,2),NIB(CP),MASS,
     x    FORMC,INFILE(CP),BETAR(CP),BETAI(CP),PARITY,RMAS,
     x    NPWCOUP,PWCOUP,JPWCOUP,IP4(CP),CPOT,NF0,FORMF,PTYPE,NEX)

  290 CONTINUE
      IF(NL.GT.0) JTEST = JTEST + 1
      X = ACC8*ALOSSM
      IF(JTEST.EQ.4 .AND.X.GT.1E-4.and..not.DRY) then
         WRITE(KO,291) LOG10(MAX(ALOSSM,1D0)),X*100.
  291 FORMAT('0***** WARNING : SOME LARGE L-TRANSFER GIVES ACCURACY LOSS
     X IN KERNELS OF',F5.1,' DIGITS,',/,
     X '0          SO MUST EXPECT ITS ERROR TO BE',   F10.5,' %',
     X '    MTMIN could be decreased'/)
         if(X > 0.1) then
                write(6,'(//'' THIS ERROR IS TOO LARGE. STOPPING'')')
                stop
                endif
         endif
!	    write(6,*) 'Finished CPAIR @ ',real(TM(I))

      NFNL = MAX(NFNL,NL)
      NL = MIN(NL,MFNL)


C
C    CHECK COULOMB MATCHING!  (at M, not N, as read-in Coulomb potentials do not go that far!)

       I = 0
       DO 335 C=1,NCH
          IC = PART(C,1)
          IA = EXCIT(C,1)
          R0 = (M-1)*HP(IC)
          T = MASS(3,IC) * MASS(4,IC) * COULCN / R0
          if(JFT==0) then
           if(IC.eq.PEL.and.IA.eq.EXL) then
            if(mod(PLANE,2)==1) T = 0.  ! elastic channel
            else
            if(    PLANE/2 ==1) T = 0.  ! nonelastic channels
           endif
          endif

          CF = 0.0
	if(LOCFIL) then
	  rewind 19
	  do IB=1,M-1
	  read(19)
	  enddo
	  read(19) FORMF(1:NF,1)
	  endif
       DO 330 NC=1,NCLIST(C,C)
	JF = NFLIST(C,C,NC)
       if(LOCFIL) then
       R = CLIST(C,C,NC) * FORMF(JF,1) 
       else
       R = CLIST(C,C,NC) * FORMF(M,JF) 
       endif
       if(LAMBDA(JF)>1.and.FCWFN.AND.PWFLAG(IC)) 
     X   R = R - CLIST(C,C,NC) * STREN(JF)/R0**(LAMBDA(JF)+1)
        IF(ABS(R).GT.1E-5 .AND. ABS(HP(IC)-HPOT(JF)).GT.1E-5)
     X      WRITE(KO,*) '***** CHANNEL ',C,' USING POTENTIAL AT',JF,
     X          ', WHICH HAS STEP SIZE ',HPOT(JF),', NOT',HP(IC),'!!'
  330  CF = CF + R
       T = T - CF
       IF(final.and.ABS(T).GT.0.020) WRITE(KO,1340) C,IC,EXCIT(C,1),T,R0
c      IF(ABS(T).GT.0.03 ) I = I + 1
 1340  FORMAT(/' ****** ERROR : CHANNEL',I4,' PARTITION',I3,' LEVEL',I3,
     X' HAS MATCHING DEFICIENCY OF',F8.4,' MEV. AT R =',F8.2)
  335  CONTINUE
C      IF(I.GT.0) STOP 'COULOMB MONOPOLES'
       IF(I.GT.0) write(6,*) ' BAD Coulomb deficiences!'
       IF(I.GT.0..and.number_calls<0) go to 350 
!      IF(I.GT.0) CALL ABEND(4)
      if(TRENEG>=1) then
      write(89,210) JTOTAL,PARITY,NCH
210	format('#',f7.1,i3,i5,' : J,pi,NCH')
      DO 215 C1=1,NCH
	write(89,212) C1,LVAL(C1),JVAL(C1),JPROJ(C1),JTARG(C1),
     x                PART(C1,1),EXCIT(C1,1)
      DO 215 C2=1,NCH
      IN = NCLIST(C1,C2)
      write(89,211) C1,C2,IN
	do I=1,IN
	write(89,213) NFLIST(C1,C2,I),CLIST(C1,C2,I)
	enddo
211	format('#',4i5)
212	format('##',i4,' LJ:',i5,3f6.1,2i4)
213	format('#',i5,1p,2e14.6)
215   continue
	write(89,210) 0.,0,0
      endif
C
 1358 FORMAT(' ',F9.3,11F11.4/(11X,11F11.4) )
 
	if(cfalloc)
     x 	CALL ASYMPTOPIA(NCH,ipmax,CI,STREN,LVAL,MASS,PART,
     X    ryev,LAMBDA,NF,CLIST,NFLIST,NCLIST,CFF)

	if(MELFIL/=0)
     X    call WRITEMEL(PEL,jtotal,psign,parity,nch,part,MINTL,
     X	  ipmax,ETOTAL,LVAL,JVAL,JPROJ,JTARG,CI,NLN,FORMF,STREN,
     X    LAMBDA,NF,MR,EXCIT,symm,nproj,ntarg,levelp,levelt,MASS,INITL,
     X    CLIST,NFLIST,NCLIST,ECM(1,1),1./(FMSCAL*RMASS(PEL)),CFF)
       If(NOSOL) go to 350 
C
! C   Only use only CRCWFN if R-turn > |RMATCH| + GAP
!         GAP = 4.5
!         if(cutr.lt.0.) GAP=abs(cutr)
C   Only use only CRCWFN if R-turn > CUTOFF*HP(1)
         GAP = 4.5
         if(CUTR.lt.-1.) GAP = max(RTURN - CUTOFF*HP(1),1d0)

      FJSWTCH = FCWFN.AND.JTOTAL.GT.AJSWTCH
     X         .and.(RTURN.gt.RMATR+GAP .or.CUTOFF.gt.MRM*5/6)
      IF(FJSWTCH .and. JSWITCH.le.0.1) JSWITCH=JTOTAL
      IF(CDETR.gt.0) WRITE(KO,*) 'At J=',JTOTAL,
     X               ', Rturn =',RTURN,':',FJSWTCH
C
!	write(KO,*) 'LAMBDA(:) =',LAMBDA(1:MLOC),' so LAMAX =',LAMAX
!	write(106,*) 'STREN(:) =',STREN(1:MLOC)

!@	if(NGAIL>0) then
!@	   DEGENY = 1d-8
!@	   APEPS  = 1d-6
!@	   APFLG(:) = (/ 0, 0, 0, 1, 0 /)
!@	   prwf = .false.
!@	   RASYM = RTURN+50.
!@
!@           CALL GAILCO (NCH,APFLG,ECM,LVAL,CFF,MAXCH,ipmax,NGAIL,
!@     x			DEGENY,GMAX)
!@!	   write(KO,*) ' Gail coefficients found '
!@            I = 0
!@      	   	write(KO,*) 
!@338      	   CALL GAILIT (CFMAT,CGMAT,NCH,MAXCH,APFLG,RASYM,ECM,
!@     X			LVAL,NGAIL,APEPS,LMAX,AERR)
!@      	   CALL WRONSK2(NCH,MAXCH,CFMAT,CGMAT,2,ECM,WERRMAX)
!@!	   write(KO,*) ' Gail expansion evaluated '
!@      	   write(KO,1020) JTOTAL,RASYM,AERR,WERRMAX
!@      	   write(156,*) jtotal,AERR,WERRMAX
!@      	   write(157,*) jtotal,sqrt(cfmat(1,2,1)**2+cfmat(1,2,2)**2)
!@           if(max(AERR,WERRMAX)>gailacc.or.RASYM<RTURN+GAP) then
!@	      RASYM = max(RASYM + 100.,RASYM*1.1)
!@	      RASYM = min(RASYMAX,RASYM)
!@	      I = I+1
!@      	   	write(156,*) '&'
!@	      if(I<10.and.RASYM<RASYMAX) goto  338
!@	      write(KO,*) ' INCREASING RASYM TO ',real(RASYM),
!@     x		 ' DID NOT REDUCE GAILITIS ERRORS SUFFICIENTLY!!'
!@	      endif
!@
!@
!@      	    IF (PRWF) THEN
!@         	WRITE (KO,1000) 'F',RASYM
!@         	CALL WRTMAT (CFMAT,NCH,NCH,MAXCH,10,KO)
!@          	WRITE (KO,1000) 'G',RASYM
!@          	CALL WRTMAT (CGMAT,NCH,NCH,MAXCH,10,KO)
!@		call flush(KO)
!@ 1000 		FORMAT (' ',a1,' Solutions at R =',F8.2)
!@ 1020 		format(' At J,R =',f7.1,f8.1,' Gailitis error =',
!@     X		1p,e9.1,' & max Wronskian error =',e9.1)
!@      	    ENDIF
!@    
!@	   endif

        IF(FCWFN.and.RTURN.ge.max(RASYM,RMATCH)*1.5) THEN
            WRITE(KO,*) 'At J=',JTOTAL,', Rturn =',RTURN
            WRITE(KO,*) 'Turning point beyond R limit,',
     X                  ' STOPPING SOON.'
!            DONES = DONES+500
             DONES = DONES+50
c            GO TO 750
         ENDIF

      NSTEPD=0
      IF (FCWFN) CALL PWISER(LVAL,ECM,K,ETA,CLIST,NCLIST,NFLIST,
     *    PART,EXCIT,CHL,ACCRCY,RASYM,CFMAT,CGMAT,CRCRAT,NCH,M,MD,
     *    PWFLAG,LAMAX,STREN,LAMBDA,DERIV,FORMF,NF,MRM,GAP,
     *    ITC,SSWITCH,FJSWTCH,CDETR,SMATL.ge.3,NGAIL>0,NSTEPD,
     *    ALLPART)
      IF (FCWFN.and..not.PWFLAG(PEL)) then
	write(6,*) ' MUST HAVE PWFLAG SET IN ELASTIC CHANNEL!'
	stop
	endif
      IF (FCWFN.and.CDETR>0) then
	if(DERIV) then
	   write(KO,*) ' Find CRCWFN derivatives at',(mrm-1)*real(hcm)
	  else
	   write(KO,*) ' Find CRCWFN at',(m-1)*real(hcm),' and ',
     *				(md-1)*real(hcm)
	  endif
	endif
C
      WRITE(KS) JTOTAL,NCH,MINTL,IPARIT,JUMPER
      WRITE(KS) (LVAL(C),JVAL(C),PART(C,1),EXCIT(C,1),C=1,NCH)
      written(KS) = .true.
c
      DO 345 C=1,NCH
        CUTVAL(C) = CUTOFF
        IF(CUTL.lt.0) then
C       			Make lower cutoff on L(C) not JTOTAL
               CUTVAL(C) = MAX(abs(CUTL)*LVAL(C), CUTR/HP(1),0d0) + 1
         IC=PART(C,1); IA=EXCIT(C,1); T = ETA(IC,IA)
!        RTURN =(T+SQRT(T**2 + LVAL(C)*(LVAL(C)+1d0)))/K(IC,IA)
               if(CUTR.lt.0.) CUTVAL(C) = max(CUTVAL(C),
     *                  int(10.**(CUTR/(JTOTAL+1))*RTURN/HP(1)) )
c        CUTVAL(C) = 5
       if(CDETR.gt.0) 
     *  write(KO,*) 'L-cutoff at ',CUTVAL(C),ICUTC,' for C,Lin,J=',
     *       		C,LVAL(C),real(JTOTAL)
        endif
345       continue
       CUTOFF = N
      DO 346 C=1,NCH
       CUTVAL(C) = max(2,CUTVAL(C))
346       CUTOFF = min(CUTOFF,CUTVAL(C))
       if(CDETR.gt.0) write(KO,*) 'min-cutoff = ',CUTOFF
c
	READWF(:)=.false.
       if(BPM/=0 .or. VEFF /=0) then
        ELCOEF(:) = 0.
	EL = INITL(1)  !   temporary fix, ok for spin-0 projectiles
        DO NC=1,NCLIST(EL,EL)
          JF = NFLIST(EL,EL,NC)
         ELCOEF(JF) = ELCOEF(JF) + CLIST(EL,EL,NC)
         enddo
        endif

	MFLI = 0
	DO 360 C=1,NCH
	DO 360 C2=1,NCH
360   MFLI = MFLI + NCLIST(C,C2)
	if(say.and.final) write(48,*) 'Allocate ',MFLI,
     x                      ' spaces in coupling list',MAXCH
	allocate(CLIS(MFLI),NFLIS(MFLI),PFLIS(MAXCH,MAXCH))
	NFLI = 0
	DO 370 C=1,NCH
	DO 370 C2=1,NCH
	PFLIS(C,C2) = NFLI
	 DO 370 NC=1,NCLIST(C,C2)
       NFLI = NFLI + 1
       CLIS(NFLI) = CLIST(C,C2,NC)
370    NFLIS(NFLI) = NFLIST(C,C2,NC)
!    	write(KO,*) 'Used ',NFLI,' spaces in coupling list'
350	continue

       
!	END SUBROUTINE MAKESET

