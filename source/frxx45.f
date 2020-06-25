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
C  INDICATOR SUBROUTINE
      logical function parall()
      parall = .true.
      return
      end
C		MPI-- SEND and RECEIVE TEXT FILES:
C
      SUBROUTINE FILESEND(KIN,NODE,MSG)
	use parallel
        use mpi
      CHARACTER*200 ALINE
!      WRITE(KIN,713)
!  713    FORMAT('EOF')
      REWIND KIN
  777 READ(KIN,'(a)',END=778) ALINE
!      IF(ALINE(1:3).eq.'EOF') GO TO 778
      do 1 L=200,1,-1
      if(ALINE(L:L).ne.' ') go to 5
1     continue
      L=0
5     continue
      call MPI_send(L,1,MPI_INTEGER,
     >               NODE,MSG,MPI_COMM_WORLD,ierr) 
      if(L.ne.0)
     >call MPI_send(ALINE,L,MPI_CHARACTER,
     >               NODE,MSG,MPI_COMM_WORLD,ierr)
      GO TO 777
  778 L = -1 ! end marker
      call MPI_send(L,1,MPI_INTEGER,
     >               NODE,MSG,MPI_COMM_WORLD,ierr) 
      REWIND KIN
      end
C      
      SUBROUTINE FILERECV(NODE,MSG,KOUTI)
	use parallel
        use mpi
      CHARACTER*200 ALINE
      KOUT = KOUTI
5     continue
      call MPI_recv(L,1,MPI_INTEGER,
     >         NODE,MSG,MPI_COMM_WORLD,status,ierr)
	if(KOUT<0.and.MSG>100) KOUT = MSG-100   ! get KOUT from MSGTYPE
      if(L<0) then
	call flush(KOUT)
	return
	endif
      if(L.ne.0) then
      call MPI_recv(ALINE,L,MPI_CHARACTER,
     >         NODE,MSG,MPI_COMM_WORLD,status,ierr) 
        WRITE(KOUT,'(a)') ALINE(1:L)
      else
        WRITE(KOUT,*) 
      endif
      GO TO 5
      end
