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
C read fam amplitude files from Fresco,
C and plot individual m-dependent cross sections
      real ja,jb,jap,jbp
      complex f1(100)
      if = 10
      do ist=1,100
      read(5,*,end=99)ja,jb,jap,jbp,nangl
      ni1=int((2*ja+1)*(2*jb+1)*(2*jap+1)*(2*jbp+1))
c     print *,ni1
	faci = 10.0/((2.*ja+1.)*(2.*jb+1.))
      do ith=1,nangl
      read(5,*)angl
      read(5,90)(f1(i),i=1,ni1)
      if(angl>1e-5) then
      do i=1,ni1
      write(if+i,60)  angl,abs(f1(i))**2*faci
60	format(f8.3,g14.5)
      enddo
      endif
90    format( 6e12.5)
      enddo
      	if = if+ni1
      enddo
99	stop
      END
 
