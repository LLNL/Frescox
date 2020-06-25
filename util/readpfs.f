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
	parameter (mxx=100, maxj=30)
	real JEX,br(mxx,0:maxj),en(mxx,0:maxj)
	
	read(5,1) NFUS
1	format(1x,i3)
	write(0,*) '# States =',NFUS
	if(NFUS>mxx) stop 'mxx'

	jmax = -1
	slats = -1
	KE = 0
	do IA=1,NFUS
	read(5,14656) II,CORFUS,SIGR,BSPR,JEX,EEX,EOUT
14656 format(i3,3f10.5,f5.1,2f8.3)
	if(IA>1) then
	j = nint(JEX-0.1)
	jmax = max(j,jmax)
	KE = KE + 1
	if(abs(JEX-SLAST)>0.1) then
	  KE = 1 ; SLAST = JEX
	  write(0,*) ' New spin ',j,' at ',IA
	 endif
	 br(KE,j) = BSPR
	 en(KE,j) = EOUT
	 write(1,*) KE,j,br(KE,j)
	endif ! IA>1
	enddo ! IA

	write(6,20) ' Spreading fractions'
20    FORMAT( '@subtitle "',A,'"'/,'@legend ON',/,
     X       '@legend x1 0.2',/,'@legend y1 0.8'/
     X       '@xaxis label "L value"',/,
     X       '@yaxis label "Probability"')

	IS = 0
	do IE=KE,1,-1
	
	write(6,25) IS,en(IE,0)
25	format('@legend string ',i2,' "Energy(a-p) =',f8.3,'"')
	do j=0,jmax
	write(6,*) j,br(IE,j)
	enddo ! j
	write(6,*) '&'

	IS = IS+1
	enddo ! IE
	end
