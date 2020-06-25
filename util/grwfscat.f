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
c  Simple program to read FRESCO scattering state wave functions (complex)
c  output when NFL<0.
c
c  If job.wf is the output file, use this program as
c	grwf < job.wf > job.plot
c	xmgr job.plot


      parameter(MAXN=10000)
      character*100 COMMENT
      dimension dr(MAXN),di(MAXN),vertr(MAXN),verti(MAXN)
C
1        READ(5,'(a)',end=999)  COMMENT
         READ(5,*) NPOINTS,STEP,RFIRST
         WRITE(0,71) COMMENT,NPOINTS,STEP,RFIRST
71      FORMAT(/'  Input:',a80,/'  Reading ',I4,' data points at '
     X   ,F8.4,' intervals, starting from r =',F8.4)
        if(NPOINTS.gt.MAXN) then
          write(0,*) ' Only room to read in ',MAXN,' potential',
     x     ' points, not ',NPOINTS
          stop
          endif
        READ(5,*) (Dr(I),Dr(I),I=1,NPOINTS)
        READ(5,*) (VERTr(I),VERTi(J),I=1,NPOINTS)
C
      totd = 0.0
      rms = 0.0
      do 15 i=1,NPOINTS
        r = (i-1)*STEP
	dr(i) = dr(i) * r
	di(i) = di(i) * r
        de = (dr(i)**2 + di(i)**2) * step
        totd = totd + de
15      rms = rms + de * r*r 
      rms = sqrt(rms/totd)
      write(0,16) totd,rms
16    format(' Wfntn volume integral =',F8.4,', rms radius =',F8.3)

      DO 200 I=1,NPOINTS
      R = (I-1)*STEP

      WRITE(6,*) real(R),Dr(I),Di(I) !,VERTr(I),VERTi(I)
200   CONTINUE
      write(6,*) '&'
	go to 1
C
999   STOP
      END
