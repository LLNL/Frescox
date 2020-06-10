!Copyright (c) 1992-2011 The University of Tennessee.  All rights reserved.
!                        Univ. of California Berkeley,
!                        Courant Institute, Argonne National Lab, and Rice University
!
!$COPYRIGHT$
!
!Additional copyrights may follow
!
!$HEADER$
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are
!met:
!
!- Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!
!- Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer listed
!  in this license in the documentation and/or other materials
!  provided with the distribution.
!
!- Neither the name of the copyright holders nor the names of its
!  contributors may be used to endorse or promote products derived from
!  this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!
!***********************************************************************
      SUBROUTINE PZLADIST0F(ZMAT,M,N, A, DESCA, IRREAD, ICREAD)
*
*     .. Scalar Arguments ..
      INTEGER            ICREAD, IRREAD
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * ), ZMAT
      COMPLEX*16, allocatable:: WORK(:)
*     ..
*
*  Purpose
*  =======
*
*  PZLADIST0/1 distributes an array of size M*N to the process grid,
*   with the array defined by function ZMAT(I,J)
*
*  PZLADIST0 is called by {IRREAD, ICREAD} 
*  PZLADIST1 is called by all other processes 
*
*  Only the process of coordinates {IRREAD, ICREAD} accesses MAT.
*
*  WORK must be of size >= MB_ = DESCA( MB_ ).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NIN
      PARAMETER          ( NIN = 11 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            H, I, IB, ICTXT, ICURCOL, ICURROW, II, J, JB,
     $                   JJ, K, LDA, M, MYCOL, MYROW, N, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, ZGERV2D, ZGESD2D, ZMAT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( MYROW.NE.IRREAD .OR. MYCOL.NE.ICREAD ) THEN
            WRITE( *, FMT = * ) 'PZLADIST0: not called by Source node'
            WRITE( *, FMT = * ) 'Abort ...'
         CALL BLACS_ABORT( ICTXT, 0 )
	ENDIF


      IF( M.GT.DESCA( M_ ).OR. N.GT.DESCA( N_ ) ) THEN
            WRITE( *, FMT = * ) 'PZLADIST: Matrix too big to fit in'
            WRITE( *, FMT = * ) 'Abort ...'
         CALL BLACS_ABORT( ICTXT, 0 )
      END IF
*
      II = 1
      JJ = 1
      ICURROW = DESCA( RSRC_ )
      ICURCOL = DESCA( CSRC_ )
      LDA = DESCA( LLD_ )
      allocate(WORK(DESCA( NB_ )))
*
*     Loop over column blocks
*
      DO 50 J = 1, N, DESCA( NB_ )
         JB = MIN(  DESCA( NB_ ), N-J+1 )
         DO 40 H = 0, JB-1
*
*           Loop over block of rows
*
            DO 30 I = 1, M, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), M-I+1 )
               IF( ICURROW.EQ.IRREAD .AND. ICURCOL.EQ.ICREAD ) THEN
                     DO 10 K = 0, IB-1
                  	 A( II+K+(JJ+H-1)*LDA ) = ZMAT(I+K,J+H)
   10                CONTINUE
               ELSE
                     DO 20 K = 1, IB
                       WORK(K)= ZMAT(I+K-1,J+H)
   20                CONTINUE
                     CALL ZGESD2D( ICTXT, IB, 1, WORK, DESCA( MB_ ), 
     $                             ICURROW, ICURCOL )
               END IF
               IF( MYROW.EQ.ICURROW )
     $            II = II + IB
               ICURROW = MOD( ICURROW+1, NPROW )
   30       CONTINUE
*
            II = 1
            ICURROW = DESCA( RSRC_ )
   40    CONTINUE
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + JB
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   50 CONTINUE
*
	deallocate (WORK)
      RETURN
*
*     End of PZLADIST0
*
      END
      
      SUBROUTINE PZLADIST1F(M,N, A, DESCA, IRREAD, ICREAD)
*
*     .. Scalar Arguments ..
      INTEGER            ICREAD, IRREAD
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLADIST distributes an array ZMAT of size M*N to the process grid.
*
*  Only the process of coordinates {IRREAD, ICREAD} accesses MAT.
*
*  PZLADIST0 is called by {IRREAD, ICREAD} 
*  PZLADIST1 is called by all other processes 
*
*  WORK must be of size >= MB_ = DESCA( MB_ ).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NIN
      PARAMETER          ( NIN = 11 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            H, I, IB, ICTXT, ICURCOL, ICURROW, II, J, JB,
     $                   JJ, K, LDA, M, MYCOL, MYROW, N, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, ZGERV2D, ZGESD2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*
      IF( MYROW.EQ.IRREAD .AND. MYCOL.EQ.ICREAD ) THEN
            WRITE( *, FMT = * ) 'PZLADIST1:  called by Source node'
            WRITE( *, FMT = * ) 'Abort ...'
         CALL BLACS_ABORT( ICTXT, 0 )
	ENDIF
*
      II = 1
      JJ = 1
      ICURROW = DESCA( RSRC_ )
      ICURCOL = DESCA( CSRC_ )
      LDA = DESCA( LLD_ )
*
*     Loop over column blocks
*
      DO 50 J = 1, N, DESCA( NB_ )
         JB = MIN(  DESCA( NB_ ), N-J+1 )
         DO 40 H = 0, JB-1
*
*           Loop over block of rows
*
            DO 30 I = 1, M, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), M-I+1 )
               
             IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               CALL ZGERV2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                             LDA, IRREAD, ICREAD )
             END IF
               IF( MYROW.EQ.ICURROW )
     $            II = II + IB
               ICURROW = MOD( ICURROW+1, NPROW )
   30       CONTINUE
*
            II = 1
            ICURROW = DESCA( RSRC_ )
   40    CONTINUE
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + JB
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   50 CONTINUE
*
      RETURN
*
*     End of PZLADIST1
*
      END      
      
