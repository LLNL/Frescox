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
      FUNCTION SECOND()
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 TARRAY(2),ETIME

      call cpu_time(etime)
      second = etime

      RETURN
      END
C  SUBROUTINES NEEDED ON SOME MACHINES, E.G. NOT CRAY
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
!                       0 : Unknown
!                       1-7 (incl.) standard serial Fortran90 
!                       8 : Intel iPSC/860 hypercube
!                       9 : CRAY T3D (planned)
!                       10 : other MIMD computer
	mach = 1 ! =2 when linked to BLAS
	return
	end
!				Change stdout recl on some machines
	subroutine stdout(nf,lrec)
	return
	end

	subroutine compiler(comp)
	character*30 comp
	comp = 'gfortran-openmp'
	call system('echo Running on `hostname`')
	return
	end
