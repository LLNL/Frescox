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
      include 'nlmod-mpi2.f'
C  SUBROUTINES NEEDED ON SOME MACHINES, E.G. NOT CRAY
      FUNCTION SECOND()
      use mpi
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 TARRAY(2),ETIME
      SECOND = MPI_WTIME()
      
C IBM-----------------
C     CALL CPTIME(I)
C     SECOND = I/100.0

C SUN/ALPHA---------
*     SECOND = ETIME(TARRAY)

C F90 real-time clock (NOT cpu time!)
!	call system_clock(ic,icr,icm)
!	if(icm*icr.ne.0) SECOND = real(ic)/real(icr)

      RETURN
      END
      FUNCTION ICAMAX(N,A,I)
	use io
      COMPLEX*16 A(N)
      ICAMAX = 1
      WRITE(KO,*) 'SUBROUTINE ICAMAX CALLED. ONLY AVAILABLE ON CRAY'
      CALL ABEND(64)
      END
      SUBROUTINE FLUSH(NF)
      RETURN
      END
	subroutine machine(mach)
!                0 : Unknown
!                1 : standard serial Fortran90
!                2 : standard serial Fortran90 with BLAS
!                3 : MPI with SGI SHMEM
!                4 : MPI with 
!	 	 5 : MPI with BLAS
	mach = 3
	return
	end
!				Change stdout recl on some machines
	subroutine stdout(nf,lrec)
	return
	end

	subroutine compiler(comp)
	character*30 comp
	comp = 'altix-mpi-cc'
	return
	end
