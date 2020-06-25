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
	implicit none
	character*80 heading
	integer pset,jset,pade,buttle,pcon,eigens,vsearch,echan,enodes
	integer nrbases,nrbmin,meigs,iblock,nearfa,iter,it0,koords
	integer chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
     x   lampl,veff,kfus,wdisk,bpm,melfil,cdcc,nfus
	integer kqmax,pp,minl,maxl,nnu,mtmin,llmax,msolve,readstates
	integer nparameters,numnode,sumform,maxcoup(3),ompform,inh
	integer pel,exl,lab,lin,lex,jump(6),jbord(7),nlab(4),plane
	logical bloch,pralpha,dry,fatal,complexbins,psiren,nosol
	real*8 hcm,rmatch,rintp,jtmin,jtmax,ips,rmatr,ebeta(2),phase
	real*8 hnl,rnl,centre,absend,cutl,cutr,cutc,epc,erange,dk
	real*8 elab(4),expand(11),hktarg,smallchan,smallcoup,
     x   unitmass,finec,weak,thmin,thmax,thinc
      real*8 hnn,rnn,rmin,rsp, rasym,accrcy,switch,ajswtch, sinjmax
	character rela,iso,ch1
      character*70 TMP,MASFIL


	character*8 namep,namet
	real*8 massp,masst,zp,zt,qval
	logical pwf
	integer nex

	real*8 jp,ep,jt,et
	integer ptyp,ptyt,cpot,copyp,copyt

	integer kp,type,shape
	real*8 p(0:7),ap,at,rc
	logical itt,nosub
	
	integer kn1,kn2,ic1,ic2,in,kind,nn,l,lmax,ia,ib,kbpot,krpot,
     x        isc,ipc,nfl,nam,nk
	real*8 sn,j,be,ampl,dm,er,e
	
	real*8 jmax,rmax,p1,p2
	integer icto,icfrom,ip1,ip2,ip3,ip4,ip5,kfrag,kcore
	integer kn,nforms
	real*8 a
	logical keep
	
	integer iex,ki

        namelist/fresco/hcm,rmatch,rintp,hnl,rnl,centre,hnn,
     X   rnn,rmin,rsp, rasym,accrcy,switch,ajswtch, sinjmax,plane,
     X   jtmin,jtmax,absend,dry,rela,nearfa,jump,jbord,pset,jset,
     x   kqmax,pp,thmin,thmax,thinc,koords,cutl,cutr,cutc,complexbins,
     x   ips,it0,iter,fatal,iblock,pade,psiren,iso,nnu,maxl,minl,mtmin,
     X   epc,erange,dk, nosol,nrbases,nrbmin,pralpha,pcon,rmatr,ebeta,
     X   meigs,buttle,weak,llmax,expand,hktarg,msolve,eigens,
     X   chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
     x   lampl,veff,kfus,wdisk,bpm,melfil,cdcc,nfus, nparameters,
     x   smallchan,smallcoup,numnode,sumform,maxcoup,ompform,
     X   inh,TMP,MASFIL,unitmass,finec,pel,exl,lab,lin,lex,elab,nlab,
     x   vsearch,echan,enodes,bloch,phase

	namelist/partition/ namep,massp,zp,nex,pwf,namet,masst,zt,qval,
     x                      readstates
	namelist/states/ jp,ptyp,ep,copyp,cpot,jt,ptyt,et,copyt
	namelist/pot/  kp,type,shape,p,itt,nosub,ap,at,rc
        namelist/overlap/ kn1,kn2,ic1,ic2,in,kind,ch1,nn,l,lmax,sn,
     &         ia,j,ib,kbpot,krpot,be,isc,ipc,nfl,nam,ampl,keep,
     &         dm,nk,er,e
     
      namelist/coupling/icto,icfrom,kind,ip1,ip2,ip3,ip4,ip5,
     X                    p1,p2,jmax,rmax,kfrag,kcore,nforms
     
      namelist/cfp/ in,ib,ia,kn,a,keep

	ki = 5
			
	read(ki,'(a)') heading
	write(1,'(a)') heading
	read(ki,'(a)')
	write(1,'(''NAMELIST'')')
	hcm = 0.1; rintp=0.5; rmatch=10;rnl=0; hnl=hcm; centre=0.
	rsp = 0; rmin=0; rasym=0; accrcy=0.; switch=0; ajswtch=0; 
	sinjmax=0; rnn=0; nnu=0; plane=0; maxl=0; minl=0; pcon=0
	jump(:) = 0; jbord(:)=0.; iso=' '; mtmin=0; epc=0.;hnn=0;
	erange= 0; nosol=.false.; llmax=-1; expand(:) = 0; msolve=0
	rmatr=0; ebeta(:)=0.; weak=1d-3; hktarg=0.2d0
	absend=-1.0; dry=.false.; rela=' '; nearfa=0; 
	pset = 0; jset = 0; dk=0.
	ips=0.001; it0=0; iter=0; iblock=0; pade=0; fatal=.true.
	nrbases=0; nrbmin=0; buttle=0; pralpha=.false.; meigs=0
	eigens=0; vsearch=0; echan=0; enodes=2; bloch=.false.; phase = 90.
	chans=0; listcc=0; treneg=0;cdetr=0; smats=0; xstabl=0; nlpl=0;
	waves=0; lampl=0; veff=0; kfus=0; wdisk=0; bpm=0; melfil=0; cdcc=0
	nfus=0; nparameters=0; smallchan=0.; smallcoup=0.; 
	numnode=0; sumform=0; maxcoup(:)=0; ompform=0; inh=0; 
	tmp=' '; masfil=' '; unitmass=0.; finec=0; 
	pel=0; exl=0; lab=0; lin=0; lex=0;
	elab(:) = 0.; nlab(:) = 0
         expand(1) = 1.2d0; expand(2) = 1.00; expand(3) = 0;
         expand(4) = 4.0; expand(5) = 0.4d0; expand(6) = 2.
         expand(7) = 1.0; expand(8) = 1.0; expand(9) = 1.0
         expand(10)= 1.0; expand(11)= 1.0	

	read(ki,nml=fresco,err=200)
	write(1,nml=fresco)

10	nex = 0
	pwf = .false.; qval = 0.; readstates=0
	read(ki,nml=partition)
	if(nex>0) then
           write(1,nml=partition)
	  else
	   write(1,*) '&partition /'
	   go to 20
	  endif
	if(readstates>0) read(readstates,*) nex
	do 15 iex=1,nex
	  ptyp=0; ptyt=0; cpot=0; ep=0.; et=0.; jp=0; jt=0.
	if(readstates==0) then
	     read(ki,nml=states)
	   else
	     read(readstates,nml=states)
	   endif
	write(1,nml=states)
15	continue
	go to 10

20	kp = 0
	type =0; shape =0; p(:) = 0.; itt=.false.; nosub=.false.
	ap=0; at=0; rc=0
	read(ki,nml=pot)
	if(kp>0) then
	  write(1,nml=pot)
	  go to 20
	 else
	  write(1,*) '&pot /'
         endif

30	keep =.false. 
	if(.not.keep) then
          kn1=0;kn2=0;ic1=0;ic2=0;in=0;kind=0;ch1=' ';nn=0;l=0;lmax=0;
          sn=0;ia=0;j=0;ib=0;dm=0;nk=0;er=erange;e=0
          kbpot=0;krpot=0;be=0;isc=0;ipc=0;nfl=0;nam=0;ampl=0
        endif

	read(ki,nml=overlap)
	if(kn1>0) then
	  write(1,nml=overlap)
	  go to 30
	 else
	  write(1,*) '&overlap /'
	 endif

40	icto=0; icfrom=0; ip1=0; ip2=0; ip3=0; ip4=0; ip5=0
	p1=0.; p2=0; jmax=0.; rmax=0.; kfrag=0; kcore=0; nforms=0
	read(ki,nml=coupling)
	if(icto>0) then
	  write(1,nml=coupling)
	   
45	  ib=0; ia=0; kn=0; a=0.0; keep=.false.
	  read(ki,cfp)
	  write(1,cfp)
	  if(in>0) go to 45
	  
	  go to 40
	 else
	  write(1,*) '&coupling /'
	 endif

50	stop
200   write(1,*) 'ERROR in fresco namelist'
	write(1,nml=fresco)
	end
