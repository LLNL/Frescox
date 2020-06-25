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
      program sigma_j
* calculates partial cross breakup section from fresco output
* reads data from file fort.38
* extended to deal with half int j, seperate +ve and -ve parity
      implicit none
      real sigj(0:500,10000),sumit(0:500,10000),sumji(0:500,20000)
      real sumel,sumbu,sumbuj(10000),sumch(0:500),sumbuj1(20000)
      real j(10000),jk(10000),sumia(500)
      real f1,f2,f3,f4,energy,be
      real jexp,jex(0:500),sumjpi(0:500),jold
      integer bandp,par(0:500),set(0:500),iset,bandpold
      integer i1,i2,i3,ia(500),iia
      integer i,ni,k,nk,it,itcm,order,maxj,nb
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
     x          jexp,bandp,energy,f2,i1,f3,f4
* count number of bound states, as exluded from bu x sec
        if(energy<be)then
          nb=nb+1
          if(it>1)iset=iset+1
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
      write(0,*)'# states (itcm)=',itcm

      call get_state(itcm,ia)

      write(0,*)'nb=',nb
      write(0,*)'WARNING:assumes all bound states before breakup states'
      write(0,*)'WARNING:assumes only 2 fusion xsec output after bu'
      write(0,*)'        check fort.38'

      do i=1,10000
        read(38,'(f7.1,9x,9g12.4,:,/(16x,9g12.4))',end=10)
     &  j(i),(sigj(it,i),it=0,itcm),f1,f2
      enddo

   10 ni=i-1

      write(0,*)'# of partial waves=',ni

      k=1
      nk=1
      do it=0,itcm
        sumit(it,1)=sigj(it,1)
      enddo
      jk(k)=j(1)
      do i=2,ni
        if(j(i)==j(i-1))then
          nk=nk+1
          do it=0,itcm
            sumit(it,k)=sumit(it,k)+sigj(it,i)
          enddo
        else
          k=k+1
          nk=1
          jk(k)=j(i)
          do it=0,itcm
            sumit(it,k)=sigj(it,i)
          enddo
        endif
      enddo

      sumji=0.
* note: have to interpolate as not all partial waves calculated
* using 4-point lagrange interpolation to third order
      call interp(jk,sumit,itcm,k,sumji)

      sumel=0.
      sumbu=0.
      sumbuj=0.
      sumch=0.
      sumjpi=0.
      sumia=0.
      
      maxj=jk(k)-jk(1)
c      write(0,*)k,maxj
      do i=1,maxj+1
        sumel=sumel+sumji(0,i)
c        write(10,'(f7.1,3(1x,e10.4))')i-1+jk(1),sumji(0,i)
        do it=0,itcm
c        write(0,*)i,it
          sumch(it)=sumch(it)+sumji(it,i)
          if(it>nb)sumbuj1(i)=sumbuj1(i)+sumji(it,i)
          if(it>nb)sumbu=sumbu+sumji(it,i)
          iset=set(it)
          if(it>nb)sumjpi(iset)=sumjpi(iset)+sumji(it,i)
          if(it>nb)sumia(ia(it))=sumia(ia(it))+sumji(it,i)
        enddo
      enddo

c      write(10,'(f7.1,g12.4)')(i-1+jk(1),sumbuj1(i),i=1,maxj+1)
      
      do i=1,k
        do it=0,itcm
          if(it>=nb)sumbuj(i)=sumbuj(i)+sumit(it,i)
        enddo
      enddo


      write(0,'(a,g12.4,a)')'Elastic x sec = ',sumel,' mb'
      do it=0,itcm
        write(0,'(a,i3,a,i3,a,e12.4,a)')
     &       'channel = ',it,
     &       ' set = ',set(it),
     &       '   =>   x sec = ',sumch(it),' mb'
      enddo
      write(0,'(a,g12.4,a)')'Breakup xsec  = ',sumbu
      do iset=1,set(itcm)
        write(0,'(a,a,f5.1,a,g12.4,a)')
     &   'j/pi = ',sign(par(iset)),jex(iset),
     &   '   =>   x sec = ',sumjpi(iset),' mb'
      enddo
      do iia=minval(ia(1:itcm)),maxval(ia(1:itcm))
        write(0,'(a,i2,a,g12.4,a)')
     &   'core state = ',iia,
     &   '   =>   x sec = ',sumia(iia),' mb'
      enddo


* write out bu xsec to stout at each j calculated in fresco
* helps to see if step sizes used are good
      do i=1,k
        if(sumbuj(i)==0)stop
        write(6,'(f7.1,3(1x,e10.4))')jk(i),sumbuj(i)
      enddo
* write out interpolated bu xsec for all j values
      write(6,*)'&'
      do i=0,maxj
        write(6,'(f7.1,3(1x,e10.4))')i+jk(1),sumbuj1(i+1)
      enddo
* write out elastic xsec for each j calculated in fresco
      write(6,*)'&'
      do i=1,k
        write(6,'(f7.1,3(1x,e10.4))')jk(i),sumit(0,i)
      enddo

      end
      
      
      subroutine interp(x,y,kmax,n,yi)
* 4-point Lagrange interpolation to 3rd order
* x(1:n) is array of x values can be non regular
* y(k,1:n) is array of y values at given x(:) points
*         for each set k
* kmax is maximum number of sets k
* n is number of array points
* yi(k,1:x(n)+1) is final interpolated set for all x values upto x(n)
      implicit none
      real x(10000),y(0:500,10000),d(10000),yi(0:500,20000),xi
      integer max,n,order,k,kmax,m,ix,i,j,i1,m1
      logical interpolate
      m=3
      max=x(n)-x(1)
      do ix=1,max+1
        xi=ix-1+x(1)
        interpolate=.true.
        d=1.
        do i1=1,n
          if(x(i1)<xi)then
            cycle
          elseif(x(i1)==xi)then
            m1=i1
            interpolate=.false.
            exit
          else
            if(i1<3)then
              m1=1
            elseif(i1==n)then
              m1=n-3
            else
              m1=i1-2
            endif
            exit
          endif
        enddo
        if(interpolate)then
          do i=m1,m1+m
            do j=m1,m1+m
              if(i/=j) d(i)=d(i)*(xi-x(j))/(x(i)-x(j))
            enddo
          enddo
          do k=0,kmax
            yi(k,ix)=0.
            do i=m1,m1+m
              yi(k,ix)=yi(k,ix)+d(i)*y(k,i)
            enddo
          enddo
        else
          do k=0,kmax
            yi(k,ix)=y(k,m1)
          enddo
        endif
      enddo
      end


      subroutine get_state(n,ia)
      implicit none
      integer:: ib,j,n,ia(n),ib1
      character:: line*10
      ia(:)=0
      ib=1
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
