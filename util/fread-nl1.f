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

	character*80 heading
	integer pset,jset,pade,buttle,pcon,eigens,vsearch,echan,enodes
	integer nrbases,nrbmin,meigs,chans,listcc,treneg,iblock
	logical bloch,pralpha
	real*8 hcm,rmatch,rintp,jtmin,jtmax,ips,rmatr,ebeta(2),phase
	real*8 elab(3)

	character*8 namep,namet
	real*8 massp,masst,zp,zt,qval
	logical pwf
	integer nex

	real*8 jp,ep,jt,et
	integer ptyp,ptyt,cpot

	integer kp,type,shape
	real*8 p(0:4)

	namelist/fresco/ hcm,rmatch,rintp,jtmin,jtmax,pset,jset,
     x		ips,it0,iter,iblock,pade,
     x		nrbases,nrbmin,buttle,pralpha,pcon,rmatr,ebeta,meigs,
     x          eigens,vsearch,echan,enodes,bloch,phase,
     x          chans,listcc,treneg,elab
	namelist/partition/ namep,massp,zp,nex,pwf,namet,masst,zt,qval
	namelist/states/ jp,ptyp,ep,cpot,jt,ptyt,et
	namelist/pot/  kp,type,shape,p
	namelist/overlap/ kn1
	namelist/coupling/ icto,icfrom,ip1,ip2,ip3
	
	read(5,'(a)') heading
	write(1,'(a)') heading
	read(5,'(a)')
	write(1,'(''NAMELIST'')')
	hcm = 0.1; rintp=0.5; rmatch=10;
	pset = 0; jset = 0; 
	ips=0.001; it0=0; iter=0; iblock=0; pade=0
	nrbases=0; nrbmin=0; buttle=0; pralpha=.false.; meigs=0
	eigens=0; vsearch=0; echan=0; enodes=2; bloch=.false.; phase = 90.
	chans=0; listcc=0; treneg=0; elab(:) = 0.

	read(5,nml=fresco)
	write(1,nml=fresco)

10	nex = 0
	pwf = .false.; qval = 0.
	read(5,nml=partition)
	if(nex>0) then
           write(1,nml=partition)
	  else
	   write(1,*) '&partition /'
	   go to 20
	  endif
	do 15 iex=1,nex
	  ptyp=0; ptyt=0; cpot=0; ep=0.; et=0.; jp=0; jt=0.
	read(5,nml=states)
	write(1,nml=states)
15	continue
	go to 10

20	kp = 0
	type =0; shape =0; p(:) = 0.
	read(5,nml=pot)
	if(kp>0) then
	  write(1,nml=pot)
	  go to 20
	 else
	  write(1,*) '&pot /'
         endif

30	kn1= 0
	read(5,nml=overlap)
	if(kn1>0) then
	  write(1,nml=overlap)
	  go to 40
	 else
	  write(1,*) '&overlap /'
	 endif

40	icto=0
	read(5,nml=coupling)
	if(kn1>0) then
	  write(1,nml=coupling)
	  go to 50
	 else
	  write(1,*) '&coupling /'
	 endif

50	stop
	end
