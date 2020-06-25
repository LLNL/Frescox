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
      SUBROUTINE E01DAF(MX,MY,X,Y,F,PX,PY,LAMDA,MU,C,WRK,IFAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER           IFAIL, MX, MY, PX, PY
      DOUBLE PRECISION  C(MX*MY), F(MX*MY), LAMDA(MX+4), MU(MY+4),
     *                  WRK((MX+6)*(MY+6)), X(MX), Y(MY)

      write(0,*) 'Nag library not available for 2D interpolation. Stop'
	stop
	end

      SUBROUTINE E02DFF(MX,MY,PX,PY,X,Y,LAMDA,MU,C,FF,WRK,LWRK,IWRK,
     *                  LIWRK,IFAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER           IFAIL, MX, MY, PX, PY
      DOUBLE PRECISION  C((PX-4)*(PY-4)), FF(MX*MY), LAMDA(PX), MU(PY),
     *                  WRK(LWRK), X(MX), Y(MY)
      INTEGER           IWRK(LIWRK)

      write(0,*) 'Nag library not available for 2D interpolation. Stop'
	stop
	end
