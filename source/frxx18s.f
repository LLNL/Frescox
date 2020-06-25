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
!  Interfaces (eg. on Cray) between double-precision program and
!	single-precision libraries

      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )

      CALL CGETRF( M, N, A, LDA, IPIV, INFO )
      RETURN
      END

      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
      CALL CGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      RETURN
      END
      subroutine zaxpy(n,za,zx,incx,zy,incy)
      complex*16 zx(*),zy(*),za
      integer i,incx,incy,ix,iy,n
      call caxpy(n,za,zx,incx,zy,incy)
      return
      end
      subroutine dgemv(tr,m,n,alph,a,lda,x,incx,bet,y,incy)
      character*1 tr
      integer*4 m,n,lda,incx,incy
      real alph,bet,a,x,y
      dimension a(lda,lda),x(lda),y(lda)
      call sgemv(tr,m,n,alph,a,lda,x,incx,bet,y,incy)
      return
      end
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      REAL*8         	A( LDA, * )

      CALL RGETRF( M, N, A, LDA, IPIV, INFO )
      RETURN
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      REAL*8         A( LDA, * ), B( LDB, * )
      CALL RGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      RETURN
      END
