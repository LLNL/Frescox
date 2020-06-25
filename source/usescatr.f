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
       I = 0; TIME0J = TM(I)
       ! if(SMATL.ge.3.and.NCHAN.ge.3) then
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
       endif
C
      if(NFUS>0) FUSL(2:1+NFUS) = 0.0
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
      IF(SMATS.EQ.3.AND.ABS(DBLE(SMAT(EL)) - 0.5).LT.0.45) SMATL = 4
      IF(SMATS.EQ.3.AND.ABS(DBLE(SMAT(EL)) - 0.5).GT.0.45) SMATL = 2
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
