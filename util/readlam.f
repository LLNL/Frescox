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
	parameter(maxmam=100,maxl=100,maxens=20)
	complex A(maxmam,0:maxl,maxens),ampl(maxmam)
	real elab(maxens)
	integer lmax(maxens)

	icreq = 1
	iareq = 1

	ie=1
1	read(5,5,end=30) elab(ie)
	write(0,6) ie,elab(ie)
5 	format(47x,g12.5)
6	format(' Energy #',i3,' = ',f10.5)

9	read(5,*) ic,ia,mam
	write(0,*) 'ic,ia,mam = ',ic,ia,mam
	if(ic==0) then  ! end of this energy
	  ie = ie + 1
	  go to 1
	 endif
	
10	read(5,11) lp,ampl(1:mam*min(1,lp+1))
	write(0,*) 'lp =',lp
11	FORMAT(4x,I5,1X,1P,10E12.4 / (10X,10E12.4))
	if(lp<0) then  ! end of thise partial waves
	  lmax(ie) = lp
          go to 9
	endif
	if(ic==icreq.and.ia==iareq) then
	 a(1:mam,lp,ie) = ampl(1:mam)
	 mamreq = mam
	 lastl = lp
	endif

	go to 10
30	stop
	end
	
