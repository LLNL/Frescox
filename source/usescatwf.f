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
       TIME0J = TM(I)
       if(SMATL.ge.3) then
       write(KO,1433) 
 1433       format(' With velocity factors, final S-matrix elements:')
!	  rtrace = 0.
	  XSC = 0.
          DO 540 C=1,NCH
            IC = PART(C,1)
            IA = EXCIT(C,1)
!@@
            IF(RMASS(IC).lt.1e-5) then
                T = sqrt(RMASS(PEL)*amu/ (HBC*K(PEL,EXL)))
            ELSE IF(RMASS(PEL).lt.1e-5) then
                T = 1./sqrt(RMASS(IC)*amu/ (HBC*K(IC,IA)))
            ELSE
                T = sqrt(K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL)))
            ENDIF
!	    see justification at frxx3.f
!
!           T = sqrt(K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL)))
!@@
          SOMA(C,1) = SMAT(C) * T
!	  if(C.ne.EL.and.abs(SMAT(C))>1d-50) rtrace = abs(SMAT(C))
	  T4 = abs(SOMA(C,1))
	  if(C.ne.EL.and.T4>1d-50) XSC = T4
540       CONTINUE
             WRITE(KO,1430) (SOMA(C,1),C=1,NCH)
!	  if(abs(rtrace)>1d-50) write(155+JIN,*) ECM(EL,1),rtrace
!	  if(abs(XSC)>1d-50) write(165+JIN,*) ECM(EL,1),XSC
!#        write(KO,1434) 
!#  1434       format(' With velocity factors, final T-matrix elements:')
!			S = 1 + 2iT
!# 	 do C=1,NCH
!# 	    if (C==EL) SOMA(C,1)=SOMA(C,1) - 1d0
!#             SOMA(C,1)=SOMA(C,1)/(0d0,2d0)
!# 	 enddo

!#             WRITE(KO,1430) (SOMA(C,1),C=1,NCH)

       endif
C
      if(NFUS>0) FUSL(2:1+NFUS) = 0.0
C-------------------------------------------CORE FUSION POTENTIAL
      DO 556 IT=1,NFUS
      DO 556 C1=1,NCH
	IC = PART(C1,1); IA = EXCIT(C1,1)
      if(ITC(IC,IA)==IT) then
	C = 0
        DO 555 JF=1,NF
          IF(PTYPE(1,JF).NE.KFUS) GO TO 555
          IF(PTYPE(2,JF).EQ.0)    GO TO 555
          IF(PTYPE(2,JF).GT.2)    GO TO 555
          IF(PTYPE(3,JF).NE.0)    GO TO 555
          T4 = 0.0
          DO 554 I=1,N
554       T4 = T4 + AIMAG(FORMF(I,JF)) * ABS(PSI(I,C1))**2
          T4 = T4 * HP(IC) / (K(PEL,EXL) * ECM(EL,1))
          T = -T4 * XCOEF * (2*JTOTAL+1.) * 4.*PI
          FUSL(1+IT) = FUSL(1+IT) + T

	T4 = abs(SMAT(C1))**2 * JCOEF * 
     x       K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL))

	           FUSLL(LVAL(C1)+1,2) = FUSLL(LVAL(C1)+1,2) + T
	  if(C==0) FUSLL(LVAL(C1)+1,1) = FUSLL(LVAL(C1)+1,1) + T4
	C = 1  ! add XO for first JF only

	if(SMATS>2) then
       write(93,5545) JTOTAL,psign(PARITY+2),IT,C1,lval(C1),JF,T,T4
!     x   ,SMAT(C1),JCOEF,K(IC,IA)/RMASS(IC)/(K(PEL,EXL)/RMASS(PEL)),
!     x    FUSLL(LVAL(C1)+1,:)
5545	 format('#',f5.1,a1,' IT/C =',2i4,'(L=',i3,'): JF=',i3,
     x      ' > PF =',f10.5,' XO =',f10.5)
!     x ', S =',2f10.5,' J, KF =',2f8.4,' Cumu:',2f10.5)
	if(T>1e+10.and.JF>2) then
      write(1000+C1,5545) JTOTAL,psign(PARITY+2),IT,C1,lval(C1),JF,T,T4
	  do I=1,M
	   R = (I-1)*HP(IC)
	  write(1000+C1,5546) R,PSI(I,C1)! ,AIMAG(FORMF(I,JF))
	  enddo
	  do I=M-3,M+100
	   R = (I-1)*HP(IC)
        CALL COULFG(R*K(IC,IA),ETA(IC,IA),Z,LVAL(C1)+1D0,
     X                  CFG,CFG(1,2),CFG(1,3),CFG(1,4),2,0,J,L1)	   
            L1 = LVAL(C1)+1
            T = 0.0
            IF(C1.EQ.EL) T = 1.0
               CF = CFG(L1,1)
               CG = CFG(L1,2)
            C6 = (CG + CI*CF) * (0.,0.5)
            C7 = (CG - CI*CF) * (0.,0.5)
            C7 = T* C7 - SMAT(C1)*C6   
	   
	  write(1000+C1,5546) R,C7
5546	  format(f8.3,2f12.5,f8.4)
	  enddo
	  
	  write(1000+C1,*) '&'
	endif
	endif
555       CONTINUE
      ENDIF
556   continue
C
      IF(VEFF.NE.0) THEN
C-------------------------------------------EFFECTIVE POTENTIALS
            DO 560 I=1,NLN
  560       VPOT(I,1) = PSI((I-1)*MR+1,EL)
            CALL CHECK(IEX,MAXICH,25)
	if(NL>0.and.NRBASES/=0) write(6,*) 
     x	 	' NONLOCAL COUPLINGS NOT INCLUDED IN VEFF!!'
!	write(139,*) 'N,(N-1)*HP(1) =',N,(N-1)*HP(1)
!	write(139,*) 'NLN1,(NLN1-1)*HP(1)*MR =',NLN1,(NLN1-1)*HP(1)*MR
!	write(139,*) 'NLN,(NLN-1)*HP(1)*MR =',NLN,(NLN-1)*HP(1)*MR
!	write(139,*) 'IEX =',IEX
        DO 590 C2=1,IEX
        DO 590 NC=1,NCLIST(EL,C2)
	 IF = NFLIST(EL,C2,NC)
C     IF(EL.GT.IEX .OR. PTYPE(3,IF).LE.0) GO TO 590
      IF(.NOT.(EL.LE.IEX .AND.
     X   (PTYPE(3,IF).GT.0.OR.PTYPE(3,IF).EQ.0.AND.EL.NE.C2))) GOTO 590
           C6 = CLIST(EL,C2,NC)
!	   write(139,*) 'EL,C2,NC =',EL,C2,NC,' S=',C6, ' IF =',IF
           IF(ABS(C6).LT.1E-10) GO TO 590
           DO 580 I=1,N
  580      SRC(I,1) = SRC(I,1) + C6 * FORMF(I,IF) * PSI(I,C2)
  590   CONTINUE
      REWIND 15
C   Choose T as J-dependent factor for Weighted Equivalent Potential
         T = (2*JTOTAL+1.) * (1. - ABS(SMAT(EL))**2)
         IF(ABS(VEFF).EQ.2 .AND. ABS(SMAT(EL)) .LT. 0.1) T = 0.0
         IF(ABS(VEFF).EQ.3) THEN
          T = 0.0
          DO 595 C=1,NCH
  595       IF(C.NE.EL) T = T + ABS(SMAT(C))**2
          T = T * (2.*JTOTAL + 1.)
         ENDIF
         IF(ABS(VEFF).EQ.4) T = 1.0
!	   write(139,*) 'Weight factor =',T
C
          CALL CHECK(NSA,MAXCH,5)
         DO IMA=1,NSA
           READ(15) (FNC(I,IMA),I=1,NLN1)
	 ENDDO
         DO IMA=1,NSA
           READ(15) (WNM(I,IMA,2),I=1,NLN1)
	 ENDDO
         IMA = NINT(JVAL(EL) - LVAL(EL) + JEX(1,PEL,EXL)) + 1
         DO I=1,NLN
          J = (I-1)*MR + 1
          FNC(I,IMA) = FNC(I,IMA) + CONJG(PSI(J,EL))*SRC(J,1) * T
          WNM(I,IMA,2) = WNM(I,IMA,2) +   ABS(VPOT(I,1))**2     * T
	 ENDDO
      REWIND 15
         DO IMA=1,NSA
          WRITE(15) (FNC(I,IMA),I=1,NLN1)
	 ENDDO
         DO IMA=1,NSA
          WRITE(15) (WNM(I,IMA,2),I=1,NLN1)
	 ENDDO
      ENDIF
C
      DO 655 C=1,NICH
         IC = PART(C,1)
         IA = EXCIT(C,1)
         IT = ITC(IC,IA)
      AMDSQS = ABS(SMAT(C))**2
      IF(mod(WAVES,2).ne.0) THEN
         DO 652 I=CUTVAL(C),N
            VPOT(I,1) =PSI(I,C)
            IF(abs(WAVES).GE.4) THEN
                C7 = PSI(I,EL) + 1E-30
                if(WAVES.le.0) C7 = VPOT(I,3) + 1E-30
                GO TO  651
                ENDIF
            IF(WAVES.GE.0) GO TO 652
        CALL COULFG((I-1)*ECM(C,3)*K(IC,IA),ETA(IC,IA),Z,LVAL(C)+1D0,
     X                  CFG,CFG(1,2),CFG(1,3),CFG(1,4),2,0,J,L1)
            IF(J.GT.0) THEN
               VPOT(I,1) = 0.
               GO TO 652
            ENDIF
            L1 = LVAL(C)+1
            T = 0.0
            IF(C.EQ.EL) T = 1.0
               CF = CFG(L1,1)
               CG = CFG(L1,2)
            C6 = (CG + CI*CF) * (0.,0.5)
            C7 = (CG - CI*CF) * (0.,-.5)
            C7 = -(T* C7 + SMAT(C)*C6)
C           C7 = CF
 
  651       VPOT(I,1) = VPOT(I,1)/C7
            if(C.eq.1) VPOT(I,2) = 0.0
            VPOT(I,2) = VPOT(I,2) + VPOT(I,1)
  652     CONTINUE
       ENDIF
      IF(mod(WAVES,2).NE.0.AND.ABS(AMDSQS).GT.1E-16) THEN
      IF(abs(WAVES).le.3) THEN
        WRITE(KO,1440) C,SMAT(C),AMDSQS,(VPOT(I,1),I=CUTVAL(C),M)
 1440 FORMAT(//' For channel no.',I3,', the S-matrix element is',
     X F15.5,' +i*',F10.5,   ' (square modulus =',F15.8,')',/
     X ' with the radial wave function'//
C    @   (5(F12.5,' +i*',F9.5,',') ) )
     X      1P,            (5(E12.3,' +i*',E9.2,',') ) )
        WRITE(191,1494) C, ((I-1)*HP(IC),PSI(I,C),I=1,N-3)
1494    format('# WF in channel ',i4,':'/(f8.3,2e14.5))
         WRITE(191,*) '&'

       ELSE
        IF(C.NE.EL.or.WAVES.le.0)
     X    WRITE(KO,1441) C,SMAT(C),AMDSQS,(C,ENLAB,(I-1)*HP(IC),
     X LVAL(C),ABS(VPOT(I,2))**2,VPOT(I,1),I=MR+1,M,MR)
 1441 FORMAT(//' For channel no.',I3,', the S-matrix element is',
     X F15.5,' +i*',F10.5,   ' (square modulus =',F15.8,')',/
     X ' with ratio to elastic wavefunction:'//
C    @   (5(F12.5,' +i*',F9.5,',') ) )
     &   (' MX:',I2,2F6.2,I4,2F12.5,' +i*',F9.5) )
      ENDIF
      ENDIF
  655 CONTINUE
C
      IF(WDISK.NE.0) THEN
	 NA = (N-1)/2*2 + 1 ! print these
		if(say)write(48,*) 'WDISK: M,N,NA = ',M,N,NA
      IF(WDISK>0) then
	  call openif(17)
         WRITE(17,657) NA,HP(PEL),ENLAB,JTOTAL,PARITY,MASS(:,PEL)
	else
	  call openuf(17)
         WRITE(17) NA,HP(PEL),ENLAB,JTOTAL,PARITY,MASS(:,PEL)
657   FORMAT(I4,2F8.4,F8.1,I3,2f12.6,2f8.3)
	endif
	  written(17) = .true.
      DO 690 IT=1,ITCM
          IF(ABS(MOD(WDISK,2)).EQ.1.AND.IT.NE.ITC(PEL,EXL)) GOTO 690
          DO 670 C=1,NICH
            IF(ITC(PART(C,1),EXCIT(C,1)).NE.IT) GOTO 670
            IF(ABS(MOD(WDISK,2)).EQ.1.AND.C.ne.EL) GOTO 670
            IF(WDISK.GT.0)
     X     WRITE(17,660) IT,LVAL(C),JVAL(C),JTOTAL,LVAL(EL),JVAL(EL),
     X                   SMAT(C),ETA(PART(C,1),EXCIT(C,1))
            IF(WDISK.LT.0)
     X     WRITE(17) IT,LVAL(C),JVAL(C),JTOTAL,LVAL(EL),JVAL(EL),
     X                   SMAT(C),ETA(PART(C,1),EXCIT(C,1))
  660         FORMAT(2I4,2F6.1,I4,F6.1,2F15.10,f12.8)
           IF(WDISK.GT.0.AND.WDISK.LE.2) THEN
              WRITE(17,'(1p,6e12.4)') (PSI(I,C),I=1,NA)
           ELSE IF(WDISK.GE.3) THEN
              WRITE(17,681) (PSI(I,C),I=1,NA)
           ELSE
             DO 665 I=2,NA
             C6 = PSI(I,C)
  665        WRITE(17) C6
           ENDIF
  670      CONTINUE
  681      FORMAT(1P,6E13.6)
  690  CONTINUE
        IF(WDISK.GT.0) WRITE(17,660) -1,0,0.,0.,0,0.,0.,0.,0.
        IF(WDISK.LT.0) WRITE(17) -1,0,0.,0.,0,0.,0.,0.  ,0.
       ENDIF
C
C     PRINT OUT REACTION CROSS SECTIONS
C     ---------------------------------
      DO 692 IT=0,ITCM
 692  SIGJ(IT) = 0.0
      CALL DISPX(SMAT,NCH,JIN,EL,JCOEF,XS,SIGJ,PART,EXCIT,
     X           JEX,ITC,K,RMASS,PEL,EXL,LVAL,FUSL,OUTJ,
     X           JVAL,JPROJ,JTARG,JTOTAL,PARITY,
     X           SMATS,CHANS+DONE,JTOTAL.GE.JTMAX-3.,SCALE,
     X           ITCM,IF1,IF2)
C
      WRITE(KS) EL,(SMAT(C),C=1,NCH),FUSL(2:1+NFUS1),OUTJ
C
      if(iams==0) then ! do it now!
       FUSJ = FUSJ + FUSL(1)
       TOTJ = TOTJ + OUTJ
       if(NFUS>0) CFUSJ(1:NFUS) = CFUSJ(1:NFUS) + FUSL(2:1+NFUS)
      endif
      IF(REAL(SMAT(EL))>0.9) then
      SMALLJ =  SMALLCOUP>0.
      DO 693 IC=1,NCHAN
        NA = NEX(IC)
        DO 693 IA=1,NEX(IC)
        IT = ITC(IC,IA)
      CHSIZES(IT) = SIGJ(IT) / JCOEF
      IF(IT/=ITC(PEL,EXL).and.K(IC,IA)>1e-5) then
      IF(CHSIZES(IT)<SMALLCHAN.and.CHPRES(IT)>0) 
     X		SMALLS(IT) = SMALLS(IT) + 1
      IF(CHSIZES(IT)>SMALLCOUP) SMALLJ = .false.
        T = ETA(IC,IA)
        T4 =(T+SQRT(T**2 + JTOTAL*(JTOTAL+1d0)))/K(IC,IA)
*	write(139,*) ia,ic,jtotal,real(rturn),real(rturn+gap),real(t4)
*	call flush(139)
!		Drop a channel if its  turning point is too far outside
!		the turning point of the elastic channel.
	if(T4>RTURN+GAP.and.SMALLCHAN>0.) SMALLS(IT) = SMALLS(IT) + 1
      endif
693   CONTINUE
!	write(138,'(f7.1)')  JTOTAL
!	write(138,'(1p,8e10.2)')  CHSIZES(1:ITCM)
!	call flush(138)
      ELSE
       SMALLJ =  .false.
      ENDIF

C
      IF(.not.FAIL) THEN
      IF(SMATS.EQ.3.AND.ABS(DBLE(SMAT(EL)) - 0.5).LT.0.45) SMATL = 4
      IF(SMATS.EQ.3.AND.ABS(DBLE(SMAT(EL)) - 0.5).GT.0.45) SMATL = 2
      ENDIF
C
      IF(JTOTAL.GT.max(jleast,ABS(JTMIN)).and.
     x           (IEXCH==0.or.XS(EL)>1e-5)) then
      IF(XS(EL).GT.ABSEND) then
   	   DONES = DONES - 1
	   COMPARE = '>'
	else
      	   DONES = DONES + 1
	   COMPARE = '<'
	endif
	if(say.and.final)
     x  write(48,694) XS(EL),COMPARE,ABSEND,JTOTAL,JTMIN,DONES
694	format('  XS ',g12.3,1x,a1,f9.4,' & JT',f7.1,' >',f4.1,
     x   ' so DONES now =',i4)
	endif
       DO 695 I=1,NTHRESH
  695  IF(DBLE(SMAT(EL)).gt.RESM(I).and.THRJ(I).le..1) THRJ(I)=JTOTAL
C
