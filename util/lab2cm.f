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
C	Test of lab2cm for elastic scattering:
c
	implicit real*8 (a-h,o-z)
	a1=1; a2=24
	a3=a1; a4=a2
	ecmi = 2;
	ecmf = 2

	do ith=0,180,1
	th = ith
	call lab2cm(TH,A1,A2,A3,A4,ECMI,ECMF,XCM_LAB)
	write(1,*) ith,xcm_lab
	write(2,*) ith,th-ith
	enddo
	end

	subroutine lab2cm(TH,A1,A2,A3,A4,ECMI,ECMF,XCM_LAB)
	implicit real*8 (a-h,o-z)
	pi = 4d0*atan(1d0)
		
	x=(A1*A3*ECMI/((A2*A4)*(ECMF)))**.5
	
	thetalab = TH*pi/180.d0
	coslab=cos(thetalab)
        c2=coslab**2
        s2=1d0-c2
*
c As 0 < TH < 180 deg (at most), we have always sinlab > 0, so it is 
c ok to take only the positive square root of s2.
*
        sinlab=sqrt(s2)
*
        thetacm=thetalab+asin(x*sinlab)
*
        coscm=cos(thetacm)
        thetacm=thetacm*180.d0/pi !switching to degrees
*
        a=1.+x*coscm
        b=abs(a)
	XCM_LAB=b/(1.+x*x+2.*x*coscm)**1.5
!        sigCM=sigLAB*XCM_LAB

	TH = thetacm   ! return cm degrees
	return
	end
