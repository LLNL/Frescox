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
      parameter(ne=1000)
      real e(ne),f(ne),ef(ne),df(ne)
C
      i = 1
1      read(5,*,err=20,end=20) e(i),f(i)
       i = i + 1
      go to 1
20    ni = i - 1
      write(0,*) 'No. of points =',ni
      do 50 i=1,ni
50    ef(i) = e(i)*f(i)
      do 60 i=4,ni-4
60    df(i) = (ef(i+3) - 2*ef(i) + ef(i-3)) / (3.0*(e(i+1)  - e(i)))**2
      do 70 i=4,ni-4
70    print *,e(i),df(i),f(i)
      STOP
      end
