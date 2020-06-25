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
	parameter (mvars=148)
	real*8 chi,pr(mvars),vars(mvars)
	integer list(mvars)
	list(:) = 0
	list(1:1) = 94
10	read(5,11,end=99) ic,nvars,chi
	read(5,12) vars
11      format(2i6,1p,e12.4,1x,0p,6f10.5)
12	format(5e12.5)
        ip = 0
	do i=1,nvars
	if(list(i)>0) then
          ip = ip+1
          pr(ip) = vars(list(i))
	endif
	enddo
        biggest = maxval(vars)
        np = ip
        write(1,20) ic,chi,biggest,pr(1:np)
20	format(i6,1p,20e12.4)
	go to 10
99	stop
	end
