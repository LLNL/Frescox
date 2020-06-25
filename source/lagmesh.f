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
	subroutine setup_mesh(M,MM,N,a,pcon)
	use meshes, f=>ybas
	implicit none
	integer M,MM,N,i,j,k,pcon,kot
	real*8 x,hm,hmm,a,b,c,d,f1,Pd,PN,ascale
!---------------------------------------------------------------------
!     calculating zeros of legendre/jacobi polynomial
!---------------------------------------------------------------------
	kot = 63
	hm = 1d0/real(M)
	hmm = 1d0/real(MM-1)
	ascale = 1.0/sqrt(a)
	if(.not.allocated(xi)) allocate(xi(N),rad(N))
!	write(6,*) 'leg,d1: ',allocated(leg),allocated(legd)
	if(allocated(leg)) deallocate(leg,legd)
	allocate(leg(MM),legd(MM))
      do j=1,MM
      x=(j-1)*hmm
      call legendre(x,N,leg(j),legd(j))
!	  write(kot,*) x,leg(j),legd(j)
      end do

      k=0
      do j=1,MM-1
      x=(j-1)*hmm
      if(leg(j).gt.0.d0.and.leg(j+1).le.0.d0.or.
     x   leg(j).lt.0.d0.and.leg(j+1).ge.0.d0)  then
      k=k+1
      call newraph(x,leg(j),N,xi(k))
!	  if(PCON>3) write(48,*) 'Found node near ',x,' at ',xi(k),': ',k
      end if
      end do
		if(N/=k) stop 'Not enough nodes found in lagmesh'
      if(PCON>0) write(kot,*)'subroutine newraph done'
!	  close(48)
	  deallocate(leg,legd)

!----------------------------------------------------------------
!     calculating Lagrange functions, f, over range 0<x<1
!----------------------------------------------------------------
	  hm = 1d0/real(M)
!	  if(.not.allocated(f))  then
	  write(48,*) allocated(f),allocated(fddi),allocated(lam)
	  allocate(f(M+1,N),lam(N),fddi(N,N))
!	  endif
	  allocate(leg(M+1),legd(M+1))
      do j=1,M+1
      x=(j-1)*hm
      call legendre(x,N,leg(j),legd(j))
	  if(PCON>3) write(48,*) x,leg(j),legd(j)
      end do
	  if(PCON>3) write(48,*) '&'

      do i=1,N
		rad(i) = xi(i)*a
      do j=1,M+1
      x=(j-1)*hm
      b=(-1)**i
      c=dsqrt((1.d0-xi(i))/xi(i))
      if (dabs(x-xi(i)).le.0d0) then

      d=x*2.d0*legd(j)  ! used L'Hopital's rule here
      else
      d=x*leg(j)/(x-xi(i))
      end if
      f(j,i)=b*c*d * ascale
!		if(PCON>0) write(161,*) x*a,f(j,i)
      end do
!		if(PCON>0) write(161,*) '&'
      end do
!		if(PCON>0) close(161)

!----------------------------------------------------------------
!     calculating Lambda(i)
!----------------------------------------------------------------
!      write(*,*)'calculating Lambda(i) has started'
      do i=1,N
      call legendre(xi(i),N,PN,Pd)
      b=(-1)**i
      c=dsqrt((1.d0-xi(i))/xi(i))
      d=xi(i)*2*Pd    ! used L'Hopital's rule here
      f1=b*c*d
      lam(i)=1.d0/f1**2
      end do
      if(PCON>0) write(kot,*)'calc. of Lambda done (used sub. legendre)'
!----------------------------------------------------------------------
!     calculating fi''(xi) and fi''(xj), to be used when calculating
!     <fi|T|fj>
!----------------------------------------------------------------------
      do i=1,N
      do j=1,N

      if(i.eq.j) then
      fddi(i,j)=-lam(i)**(-0.5d0)*(N*(N+1)*xi(i)*(1.d0-xi(i)) 
     x         -3.d0*xi(i)+1.d0)/(3.d0*xi(i)**2*(1.d0-xi(i))**2)
      else

      fddi(i,j)=(-1)**(i+j)*lam(j)**(-0.5d0)*((2.d0*xi(j)**2-xi(i)
     & -xi(j))/(xi(i)*(xi(j)-xi(i))**2))*dsqrt(xi(i)*(1.d0-xi(i))
     & /(xi(j)*(1.d0-xi(j))**3))
      end if
      end do
      end do

!      write(kot,*)'fddi calculated!'
	  if(PCON>0) then
	  do i=1,N
	  write(62,*) xi(i)*a,lam(i)
	  enddo
	  close(62)
	  endif
	return
	end
	
      subroutine newraph(x,P,N,x0)
!     this routine finds the zeros of legendre functions
      double precision F,Fdash,z1,z2,x,x0,P
      integer N
      x0=x
      F=P
      z1=2.d0*x-1.d0
      do 20 while (dabs(F).gt.1.d-9)
      call legendre(x0,N,F,Fdash)
      z2=z1-F/Fdash
      z1=z2
      x0=0.5d0*(z1+1.d0)      
20    continue
      return
      end


      subroutine legendre(x,N,res,resd)
!     this routine calculates legendre function and its derivative
      double precision P(0:N),Pdash(0:N),z,i1,res,resd,x
      integer N,i
      z=2.d0*x-1.d0
      P(0)=1.d0
      if(N>0) then
	P(1)=z
      	Pdash(1)=1.d0
	endif
      do 10 i=2,N
      i1=i
      P(i)=((2.d0*i1-1.d0)*z*P(i-1)-(i1-1.d0)*P(i-2))/i1
      Pdash(i)=i*P(i-1)+z*Pdash(i-1)
10    continue
      res=P(N)
      resd=Pdash(N)
      return
      end
	

