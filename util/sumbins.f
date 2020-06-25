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
      parameter(mset=100, mlen=2000)
      real x(mlen),y(mlen,mset),ysum(mlen)
      integer incl(mset)
      character*80 card
      character first
      
       nset = 1
       nel   = 1
       n = 0
10    continue
      read(5,'(a)',err=20,end=30) card
      first = card(1:1)
      if(first.eq.'@' .or. first.eq.'#') go to 10
      read(card,*,err=20) x(nel),y(nel,nset)
      n = max(n,nel)
      nel = nel + 1
      if(nel.gt.mlen) then
         write(0,*) 'ONLY ',mlen,' values allowed'
         stop
         endif
      go to 10
20    nset = nset+1
      if(nset.gt.mset) then
         write(0,*) 'ONLY ',mset,' sets allowed'
         stop
         endif
      nel = 1
      go to 10
30    if(nel.eq.1) nset = nset -1
      write(0,*) 'Read in ',nset,' sets, each of ',n,' values'
      ninc = nset-1
         do 40 i=1,ninc
40       incl(i) = i+1
       write(0,*) 'Including:',(incl(i),i=1,ninc)
       scale = 1.0
       write(0,*) 'scaling by ',scale
       do 60 j=1,n
       ysum(j) = 0
       do 60 i=1,ninc
60      ysum(j) = ysum(j) + y(j,incl(i))
       do 70 j=1,n
70     write(6,*) x(j),ysum(j)*scale
        stop
       end
