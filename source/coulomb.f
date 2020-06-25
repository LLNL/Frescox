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

!	SUBROUTINE COULOMB
!     	implicit real*8(a-h,o-z)
	 
	 CHL(:,:,:) = 0.
      DO 161 IC=1,NCHAN
      NA = NEX(IC)
      DO 160 IA=1,NA
      IT = ITC(IC,IA)
!!       if(RMASS(IC).lt.1e-5) go to 160
         RMK = (M-1)*HP(IC) * K(IC,IA)
         CRho = dcmplx(RMK)
	 if(DERIV.and.rterms) RMK = (MRM-1)*HP(IC) * K(IC,IA)
         RMKD = (MD-1)*HP(IC)  * K(IC,IA)
         IF(ECMC(IC,IA).GT.0.0) THEN
           CEta = ETA(IC,IA)
         CALL PHASES(ETA(IC,IA),MAL1,CSIG(1,IT))
!@       IF (NGAIL>=1) go to 160
          T = K(IC,IA)
       IF (FCWFN.AND.PWFLAG(IC)) THEN
          CALL COULFG(RASYM*T,ETA(IC,IA),Z,XLMAX,CFG,CFG(1,2),
     *         CFG(1,3),CFG(1,4),1,0,I,M1)
           DO 149 L1=1,MAL1
             CF = CFG(L1,1)
             CG = CFG(L1,2)
             CHL(L1,IT,1) = CG+CI*CF
!				Derivatives wrt R:
             CF = CFG(L1,3) * T
             CG = CFG(L1,4) * T
             CHL(L1,IT,2) = CG+CI*CF
 149       CONTINUE
       ELSE
        CALL COULFG(RMK,ETA(IC,IA),Z,XLMAX,CFG,CFG(1,2),
     *              CFG(1,3),CFG(1,4),1,0,I,M1)
           DO 150 L1=1,MAL1
              CF = CFG(L1,1)
              CG = CFG(L1,2)
C 150         CHL(L1,IT,1) = DCMPLX(CG,CF)
  150         CHL(L1,IT,1) = CG+CI*CF
	if(DERIV) then
!				Derivatives wrt R:
           DO  L1=1,MAL1
              CF = CFG(L1,3) * T
              CG = CFG(L1,4) * T
	      CHL(L1,IT,2) = CG+CI*CF
	   ENDDO
	else
        CALL COULFG(RMKD,ETA(IC,IA),Z,XLMAX,CFG,CFG(1,2),
     *              CFG(1,3),CFG(1,4),2,0,I,M1)
           DO 151 L1=1,MAL1
              CF = CFG(L1,1)
              CG = CFG(L1,2)
C 151         CHL(L1,IT,2)= DCMPLX(CG,CF)
  151         CHL(L1,IT,2)= CG+CI*CF
	endif
       ENDIF
         R = (M-1)*HP(IC) 
      ELSE
C              CLOSED CHANNELS:
        IE = 0
         R = (M-1)*HP(IC)
        CALL WHIT(ETA(IC,IA),R,K(IC,IA),ECMC(IC,IA),MAL1,
     x            CFG,CFG(1,2),IE)
           DO L1=1,MAL1
              CSIG(L1,IT) = 0.0
              CHL(L1,IT,1) = dcmplx(0d0,CFG(L1,1))
	   ENDDO
	if(DERIV) then
           DO  L1=1,MAL1
              CHL(L1,IT,2) = dcmplx(0d0,CFG(L1,2))
	   ENDDO
	else
         R = (MD-1)*HP(IC)
        CALL WHIT(ETA(IC,IA),R,K(IC,IA),ECMC(IC,IA),MAL1,
     x            CFG,CFG(1,2),IE)
           DO L1=1,MAL1
  	    CHL(L1,IT,2)= dcmplx(0d0,CFG(L1,1))
	    enddo
	endif
          CEta = dcmplx(0d0,-ETA(IC,IA))
          CRho = dcmplx(0d0,R*K(IC,IA))
      ENDIF
C
      if(pralpha.and.DERIV.and.pcon>3) then
      do L1=1,min(MAL1,10)
 	L=L1-1
      C6 = CHL(L1,IT,2)/CHL(L1,IT,1) * R
      Shift = real(C6)
      Pen = aimag(C6)
      T=CF2(CRho,CEta,dcmplx(L),CI,1d-7,AERR,20000,1d-15)*CRho
      AL = ENLAB
      AL = ECMC(IC,IA)
      write(1400+it*100+l,1511) AL,Shift,T,R*K(IC,IA),eta(ic,ia),
     x                          CHL(L1,IT,1),CHL(L1,IT,2)*R,Shift/T
      write(2400+it*100+l,1511) AL,CSIG(L1,IT)
 1511	format(1p,12e13.5)
 	enddo
  	endif
  160 CONTINUE
  161 CONTINUE
C              Find shift values at R-matrix pole positions
      if(btype == 'A'.or.eobs) then

	do ip=1,nvars
        if(srch_kind(ip)==3) then
         it = srch_rterm(ip)
         EPOLE = srch_value(ip)
         do i = 1,nvars   ! find widths for that pole
          if(srch_kind(i)==4.and.srch_rterm(i)==it) then
           L = srch_r_lch(i)
           IC = srch_r_ic(i)
           IA = srch_r_ia(i)
           EE = EPOLE + QVAL(IC)
           X =  sqrt(FMSCAL*RMASS(IC) * ABS(EE))   ! wave number k
           RMK = (M-1)*HP(IC) * X
           ETAS = ETACNS * MASS(2+1,IC) * MASS(2+2,IC)
     X          * SQRT(RMASS(IC)/ ABS(EE)) * GAM(IC,IA)
           if(EE>0d0) then
             CRho = dcmplx(RMK)
             CEta = dcmplx(ETAS)
           else
             CEta = dcmplx(0d0,-ETAS)
             CRho = dcmplx(0d0,RMK)
           endif

           T = CF2(CRho,CEta,dcmplx(L),CI,1d-7,AERR,20000,1d-15)*CRho
           srch_r_shift(i) = T
           if(pralpha) write(60,16101) i,srch_name(i),IC,IA,L,EE,T
16101	 format(' Width #',i3,' (',a10,') with IC,IA,L,E=',3i3,f8.4,
     x          ' S=',f10.5)
          endif! kind=4
	 enddo ! i
 	endif  ! kind=3
 	enddo  ! ip
        ENDIF ! btype=A

!	END SUBROUTINE COULOMB
