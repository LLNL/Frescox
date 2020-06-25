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
      program sumbins_cc
      implicit none
      real sumit(5000)
      real sumbu
      real sumia(5000)
      real f1,f2,f3,f4,energy,be
      real jexp,jex(500),sumjpi(5000),jold,jmax
      real jcore(10)
      integer bandp,par(5000),set(5000),iset,bandpold
      integer parcore(10)
      integer i1,i2,i3,ia(5000),iia,nset,count,icore
      integer i,ni,k,nk,it,itcm,order,maxj,nb,kmax,nc,ic
      character text*80
      character*1 sign(-1:1)
      
      sign(-1)='-'
      sign(1) ='+'

* read in itcm= number of channels
* from xst file, fort,13
      read(13,*) text
      read(13,*) i1,i2,i3,f1,be
      read(13,*) i1,itcm
      nb=0
      iset=0
      do it=1,itcm
        read(13,'(2(f5.1,i3,f8.4),1p,e12.4)') 
     x          jexp,bandp,energy,f2,i1,f3,sumit(it)
* count number of bound states, as exluded from bu x sec
        if(energy<be)then
          nb=nb+1
          iset=iset+1
          set(it)=iset
          jex(iset)=jexp
          par(iset)=bandp
        else
          if((jexp/=jold.or.bandpold/=bandp).or.it==nb+1)then
            iset=iset+1
            jex(iset)=jexp
            par(iset)=bandp
          endif
          set(it)=iset
        endif
        jold=jexp
        bandpold=bandp
      enddo
      read(13,*) i1,nc
      do ic=1,nc
        read(13,'(2(f5.1,i3,f8.4),1p,e12.4)') 
     x          jexp,bandp,energy,f2,i1,f3,sumit(it)
        jcore(ic)=jexp
        parcore(ic)=bandp
      enddo
      write(6,*)'# states (itcm)=',itcm
      nset=iset

      call get_state(itcm,ia)

      write(6,*)'# bound states=',nb

      sumbu=0.
      sumjpi=0.
      sumia=0.
      
      do it=1,itcm
        if(it>nb)sumbu=sumbu+sumit(it)
        iset=set(it)
        if(it>1)sumjpi(iset)=sumjpi(iset)+sumit(it)
        if(it>nb)sumia(ia(it))=sumia(ia(it))+sumit(it)
      enddo

      do it=1,itcm
        write(6,'(a,i3,a,i3,a,f4.1,a1,a,i1,a,e12.4,a)')
     &       'channel = ',it,
     &       ' set = ',set(it),
     &       ' j/pi =',jex(set(it)),sign(par(set(it))),
     &       ' core = ',ia(it),
     &       '   =>   x sec = ',sumit(it),' mb'
      enddo
      write(6,*)
      write(6,'(a,g12.4,a)')'Breakup xsec  = ',sumbu
      do iset=2,set(itcm)
        write(6,'(a,i3,a,f5.1,a,a,g12.4,a)')
     &   'set = ',iset,
     &   ' j/pi = ',jex(iset),sign(par(iset)),
     &   '   =>   x sec = ',sumjpi(iset),' mb'
      enddo
      do iia=minval(ia(1:itcm)),maxval(ia(1:itcm))
        write(6,'(a,i2,a,f5.1,a,a,g12.4,a)')
     &   'core state = ',iia,
     &   ' j/pi = ',jcore(iia),sign(parcore(iia)),
     &   '   =>   x sec = ',sumia(iia),' mb'
      enddo
      end

      subroutine get_state(n,ia)
      implicit none
      integer:: ib,j,n,ia(n),ib1
      character:: line*10
      ia(:)=0
      ib=1
      j=0
      do while(j<2)
       read(1,'(a)')line
       if(line(1:10)=='          ')j=j+1
      enddo
      do while(ib<=n)
       if(ib<100)then
        read(1,'(26x,i2,4x,i2,45x)')ia(ib),ib1
        if(ib1/=ib.and.ib1/=0)write(0,*)'ERROR: IB MISMATCH',ib,ib1
       else
        read(1,'(26x,i2,51x)')ia(ib)
       endif
       if(ia(ib)==0)ia(ib)=1
       ib=ib+1
      enddo
      end
