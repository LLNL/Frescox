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
	neigs = min(meigs,naa)
	allocate(val(neigs),phys(neigs))
		nlw = min(nln-2,(n-1)/mr+1)
!		write(6,*) ' Allocate gswf with ',n,nln,mr,nlw,nch
     	allocate(rvec(naa),gswf(n,nch))
C					No scattering!: just find eigenstates!
	if(meigs<naa) then	! Inverse iteration for lowest eigenstates
	 if(vsearch>0) stop '  Meigs too small'
	eps = 1e-8
	mc = 30
	E1 = ECM(INITL(1),1)
	if(tr) write(6,715) neigs,E1,EOFF
715	format(/' Inverse iteration to find ',i3,' eigenstates nearest'
     x		,F10.5,' MeV (cm).   EOFF=',f10.5)
	   call flush(6)
     	   allocate(aar(nd,nd))
		do i=1,naa
	 	aar(1:naa,i) = aa(1:naa,i)
	 	enddo
C			L.U decomposition of Hamiltonian-E matrix
	  t = second()
	  call zgetrf(naa,naa,AA,nd,ipiv,info)
	  write(138,*) 'zgetrf of ',naa,' took ',real(second()-t)
	  call flush(138)

	do 800 kv=1,neigs
c			starting vector
	do 720 i=1,naa
720	rvec(i)=1d0/sqrt(real(naa))
c	rvec(1) = 1d0
	E = 0d0
		
	do 780 ic=1,mc
	  elast = E
	   do 725 i=1,naa
725	     rhs(i,kv) = rvec(i)
        call zgetrs('N',naa,1,AA,nd,ipiv,rhs(1,kv),naa,info)

	   if(info.ne.0) then
	     write(6,*)' Singular H-E Matrix for this energy',info
	     stop 'Rmatr'
	   endif

	do 740 kvp=1,kv-1
	wn = 0d0
	if(ISNONO) then
	 do 728 i=1,naa
	 do 728 j=1,naa
728	 wn = wn + rhs(i,kv)*kk(i,j)*rhs(j,kvp)
	endif
	do 730 i=1,naa
730	wn = wn + rhs(i,kv)*rhs(i,kvp)
	do 735 i=1,naa
735	rhs(i,kv) = rhs(i,kv) - wn*rhs(i,kvp)
740	continue
	wn = 0d0
	do 745 i=1,naa
745	wn = wn + rhs(i,kv)*rvec(i)
	do 750 i=1,naa
	if(ISNONO) then
	   an = rhs(i,kv)   !   Norm matrix here is  1 + kk
	   do 748 j=1,naa
748	      an = an + kk(i,j)*rhs(j,kv)
	   rvec(i) = an
	else
	   rvec(i) = rhs(i,kv)
	endif
750	continue
	an=0d0;  tc=0d0
	do 760 i=1,naa
      	  tc = tc + rhs(i,kv)**2
760	  an = an + rhs(i,kv)*rvec(i)
	if(ISNONO) then
	   phys(kv) = an/tc
	 else
	   phys(kv) = 1d0
	 endif
	E = E1 + wn/an
	if(tr) write(6,765) ic,E,E-elast,phys(kv)
	if(tr) write(imp,765) ic,E,E-elast,phys(kv)
765	format(' Iteration ',I4,' gives ',F10.5,', change =',1P,E10.2,
     x	 0p,' (phys',f9.5,')')
	call flush(6)
	an = 1d0/sqrt(abs(an))
	do 770 i=1,naa
770	  rvec(i) = rvec(i) * an
	if(tr) write(imp,772) real(rvec(1:naa))
772	format(1x,10f7.3)
	val(kv) = E
	if(ic.gt.3 .and. abs(E-elast).lt.eps) go to 785
780	continue
	if(tr) write(6,*)  mc,kv,E
781	format(' Inverse iteration failed after',i3,
     x         ' iterations for state ',i2,' near',f10.5/)
	neigs = kv-1
	go to 801
c	stop
785	wn = 0d0
	if(ISNONO) then
	 do 788 i=1,naa
	 do 788 j=1,naa
788	 wn = wn + rhs(i,kv)*kk(i,j)*rhs(j,kv)
	endif ! ISNONO
	do 789 i=1,naa
789	wn = wn + rhs(i,kv)*rhs(i,kv)
	wn = abs(wn)
	write(6,*) ' Norms factors :',real(an),1./real(sqrt(wn))
	an = 1d0/sqrt(wn)
	do 790 i=1,naa
790	rhs(i,kv) = rhs(i,kv) * an
!	if(E>0) go to 800
!	write(6,772) real(rhs(1:naa,kv))
	an  = E1
	do i=1,naa
	do j=1,naa
	an = an + rhs(i,kv)*aar(i,j)*rhs(j,kv)
	enddo
	enddo
	if(tr) write(6,*) 'Eval at ',real(E),', <H> =',real(an)
	if(tr) write(6,*) 
800	continue
801	 deallocate(aar)
	else
     	   allocate(evec(nd,nd),aar(nd,nd))
	   E1 = ECM(INITL(1),1)
	if(ISNONO.or.vsearch>0) then
!				FIND ALL EIGENSTATES FOR NONSYMMETRIC CASE
		do i=1,naa
	         if(vsearch==0) kk(i,i) = kk(i,i) + 1d0
	         !! kk(i,i) = kk(i,i) + 1d0
		aar(1:naa,i) = aa(1:naa,i)
	 	enddo
     	   allocate(qvec(nd),vec(nd),kkt(nd,nd))
		kkt(:,:) = kk(:,:)
	   write(imp,*) 'Norm matrix:'
	   call WRRMAT(kkt,naa,naa,nd,6,imp)
C find all eigensolutions of unsymmetric matrix problem Ax=EBx
C with F02BJF
c               call f02bjf(naa,aar,nd,kkt,nd,1d-10,VAL,qvec,
c    1                      vec,.true.,evec,nd,ipiv,IT)
       write(0,*) 'Nag library not available for F02BJF. Stop'
		do i=1,naa
!		write(162,*) i,val(i),vec(i)
 		VAL(i) = VAL(i)/vec(i) 	
		if(vsearch==0) VAL(i) = VAL(i) + E1
		qvec(i) = qvec(i)/vec(i)
	 	rhs(1:naa,i) = evec(1:naa,i)
		enddo
!		write(162,*) 'IFAIL =',IT
!		write(162,805) ipiv(1:naa)
805		format(20i4)
		deallocate(vec,kkt)
	  else
!				FIND ALL EIGENSTATES FOR SYMMETRIC CASE
		do i=1,naa
	 	aar(1:naa,i) = aa(1:naa,i)
	 	enddo
             call HDIAG(aar,naa,nd,0,evec,nc)
	        phys(:)=1d0
		do i=1,naa
	 	 val(i) = aar(i,i) + E1
	 	 rhs(1:naa,i) = evec(i,1:naa)
	 	enddo
	  endif
     	   deallocate(evec,aar)
	endif

	   if(CDETR>3) then
	   write(imp,*) 'Complex eigenvectors in cols'
	   call WRCMAT(rhs,naa,neigs,nd,5,imp)
	   write(imp,*) 'For eigenvalues'
	   write(imp,'(4(G12.4,G9.2))') (VAL(i),qvec(i),i=1,naa)

	   endif
	t = 1e6; j=0
	 ac = ECM(INITL(1),1) - EOFF
	 if(vsearch>0) ac = 1.0
         do i=1,neigs
	 if(abs(val(i)-ac)<t) then
	   t = abs(val(i)-ac)
	   j = i
	  endif
	 enddo
	if(j==0) go to  950
	if(vsearch==0) write(6,810) (val(i)+EOFF,i=1,neigs)
	if(vsearch>0 ) write(6,8101) (val(i),i=1,neigs)
810	format(' ENER eigenstates + offset at ',5f10.4)
8101	format(' POTL eigenstates          at ',5f10.4)
	 eof=0.; if (vsearch>0) eof = EOFF
	write(6,8102) ac+eof,val(j)+eof,j
8102	format(' ENGY eigenstate nearest ',f10.4,' is ',f10.4,' at ',i5)
	if(allocated(qvec)) then
	if(neigs==naa.and.ISNONO.and.maxval(abs(qvec(1:neigs)))>1e-3) 
     x        then
	  write(6,811) (qvec(i),i=1,neigs)
811	  format(' Engy eigenstates imag  parts ',5f10.4)
	  deallocate (qvec)
	  endif
	  endif
	if(ISNONO.and.neigs<naa.and.CDETR>1)  
     x         write(6,812) (phys(i),i=1,neigs)
812	format(' Physical norms of Engy states',5f10.4)
	kvo = 0
	adjusted = .false.
        do kv=1,neigs
         eigen = 'Scatt'
          if(VAL(kv).lt.0) eigen = 'BOUND'
	  if(vsearch>0) eigen = 'POTL'
	  L=1
	  if(ISNONO.and.eigens>1.and.vsearch==0) L=2
!	  if(VAL(kv)>max(0d0,ECM(INITL(1),1))+2.) L=0
	  if(VAL(kv)>max(0d0,ECM(INITL(1),1))+2.or.VAL(kv)<-1d3) L=0
	  if(L>0) kvo=kvo+1
	do nb=1,L
	 ifout=200+(kvo-1)*L+nb
          write(6,820) eigen,JTOTAL,PSIGN(parity+2),VAL(kv)+EOFF
820         format(/' ',A5,' EIG',F5.1,A1,': ',F15.6)
	  nodest(:) = 1
	  gswf(:,:) = 0.0
	  pnorms(:) = 0.0
	  do 832 c=1,nch
	    phi3(c) = 0.0
	  do 830 ib=1,nbasis(c)
            ki = nbas0(c) + ib
            do i=1,n
	      gswf(i,c) = gswf(i,c) + rhs(ki,kv)*ybas(i,c,ib)
	      enddo
830	  phi3(c) = phi3(c) + rhs(ki,kv)**2
832	  pnorms(PART(c,1)) = pnorms(PART(c,1)) + phi3(c)
	  do c=1,nch
	   do i=3,n
            if(gswf(i,c)*gswf(i-1,c)<-1e-20) nodest(c) = nodest(c) +1
	   enddo
 	  enddo
!	  write(6,840) (c,phi3(c),c=1,nch)
840	 format(4(i4,F14.10,','))
!	  write(6,841) (nodest(c),c=1,nch)
841	 format(4(4x,8x,'N =',i3,',',:))
	  if(MXP>1) write(6,842) (ic,pnorms(ic),ic=1,min(5,MXP))
842	 format('  partition norms:'/,5(i4,F14.10,','))
           
	NS = 0
	if(vsearch>0.and.echan>0.and.echan<=nch.and..not.adjusted) then
 	  if(nodest(echan)==enodes.and.phi3(echan)>0.1) then
           write(6,845) vsearch,val(kv)
845         format(/' ***  Adjust potential component at',i3,
     x              ' by factor ',f10.5/)
	    adjusted = .true.; NS=1
	  endif
	endif

	if(nb==L) then
	  write(ifout,885) eigen,kv,val(kv)+eof
	  write(16,885) eigen,kv,val(kv)+eof,NS
	  write( 6,885) eigen,kv,val(kv)+eof,NS
885       format('#'/'# ',a5,' eigensolution ',i4,'  at    ',f14.6,i3,
     X      / '#    i: Part Excit   E    L   J   Jp   Jt  Blockd  Norm',
     x        '   Nodes')

	  NS = (nch-1)/10+1
	  do 896 IS=1,NS
	   nc = min(10,nch-(IS-1)*10)
	   do 887 ic=1,nc
	   I=(IS-1)*10 + ic
	   write(ifout,886) I,PART(I,1),EXCIT(I,1),
     x 		ECM(1,1)-ECM(I,1)+EOFF,LVAL(I),
     x         JVAL(I),JPROJ(I),JTARG(I),BLOCKD(I),phi3(I),nodest(I)
	   write(16,886) I,PART(I,1),EXCIT(I,1),
     x 		ECM(1,1)-ECM(I,1)+EOFF,LVAL(I),
     x         JVAL(I),JPROJ(I),JTARG(I),BLOCKD(I),phi3(I),nodest(I)
	   write( 6,886) I,PART(I,1),EXCIT(I,1),
     x 		ECM(1,1)-ECM(I,1)+EOFF,LVAL(I),
     x         JVAL(I),JPROJ(I),JTARG(I),BLOCKD(I),phi3(I),nodest(I)
886	  format('#',1x,I4,':',2I5,f6.2,i4,3f5.1,L5,f10.6,i4)
887	   continue
	  do 890 j=1,n
	    r = (j-1)*HP(PART(1,1))
          write(ifout,895) r,(gswf(j,(IS-1)*10+ic),ic=1,nc)
890	  write(16,895) r,(gswf(j,(IS-1)*10+ic),ic=1,nc)
895	  format(f7.3,1p,10e12.4)
       	write(16,*) '&'
896    	write(ifout,*) '&'
	call flush(ifout)
	endif

	if(nb==L.and.mxp>1.and.msp>1) then
	nlw = min(nln-2,(n-1)/mr+1)
	allocate (wfnxy(nlw,nlw),donell(nch,msp))
	donell(:,:) = .false.
!				Find also 3-body wave function:
	RINTP = MR*HP(1)
	write(88,'(i4,f10.5,i4,f10.5)') nlw,RINTP,-2,val(kv)+EOFF
	write(88,'(''#'')')
	write(88,'(''# Wave function'')')
	write(6,'(/'' 2D wave function:'')')
	  ac = mass(2,1)
	  an1 = mass(1,2)
	  an2 = mass(1,1) - mass(1,2)
	  am12 = an1*an2/(an1+an2)
	  am123 = ac*(an1+an2)/(ac+an1+an2)
	  rcore = 0.
	  if(abs(ac-4.)<0.01) rcore = 1.49
	 write(6,'('' Masses'',3f8.3,'' Core rms'',f8.3)') 
     X	    ac,an1,an2,rcore
		call flush(6)
	tnorm  = 0.
	rho2  = 0.
	do mc=1,100
	  qnn(:) = -1
	  changed = .false.
	do i=1,nch
	do kn=1,msp
	if(abs(afrac(itc(part(I,1),excit(I,1)),itc(2,1),1,kn))>1e-9
     X	  .and.part(i,1)==1.and..not.donell(i,kn)) then
	 l = qnf(9,kn)
	 if(qnn(1)<0) then  ! first find of new component
           qnn(1)=lval(i);qnn(2)=l;qnn(3)=qnf(10,kn);qnn(4)=qnf(11,kn)
     	   qnn(5)=nint(2*jex(1,1,2));qnn(6)=nint(2*jval(i))
             write(88,915) lval(i),l,qnf(10,kn),qnf(11,kn),
     X		1,nint(2*jex(1,1,2)),jval(i),jtotal,kn    ! new channel
             write(6,915) lval(i),l,qnf(10,kn),qnf(11,kn),
     X		1,nint(2*jex(1,1,2)),jval(i),jtotal,kn    ! new channel
915          format('#  ',6i3,2f4.1,i3,' : L,l,2*s,2*j,#I,2*I,Jp,J,kn')
	    wfnxy(:,:) = 0.
	 endif
         if(qnn(1)==lval(i).and.qnn(2)==l.and.qnn(3)==qnf(10,kn)
     X		.and.qnn(4)==qnf(11,kn).and.
     X 	   qnn(5)==nint(2*jex(1,1,2)).and.qnn(6)==nint(2*jval(i))) then
!						found another component

c                                       Plot in rectangular grid
		do ii=1,nlw
	 	r = ii*RINTP
	 	t = forml(ii+1,kn,1)*r
	 	do is=1,nlw
	 	wfnxy(ii,is) = wfnxy(ii,is) + t*gswf(is*MR+1,i)
	 	enddo
		enddo
 	    donell(i,kn)=.true.
	    changed = .true.
	    endif  ! quantum number match
	  endif  ! afrac
	  enddo  ! kn
	  enddo  ! i
	if(changed) then		! print wf found for this mc if any found
	an = 0.		
	rms(:) = 0.
	do ii=1,nlw
	 write(88,*) ii
	 write(88,942) (wfnxy(ii,is),is=1,nlw)
	 do is=1,nlw
	 t =  wfnxy(ii,is)**2 * rintp**2
	 an = an + t
		rms(1) = rms(1) + t * (ii*RINTP)**2
		rms(2) = rms(2) + t * (is*RINTP)**2
	 enddo !is
 942     format(1p,6e12.4)
	enddo !ii
	 tnorm = tnorm + an
	 rho2 = rho2 + am123*rms(2) + am12*rms(1)
	 rms(:) = sqrt(rms(:)/(an+1e-20))
	 write(88,'(3f12.6,'' = Sq.xy norm,rR rms'')') an,rms
	 write(6,'(3f12.6,'' = Sq.xy norm,rR rms'')') an,rms
	endif  ! if changed
	enddo  ! mc
	   write(88,915) -1
	   rmst = sqrt((ac*rcore**2 + rho2/(tnorm+1e-20))/(ac+an1+an2))
	   write(88,'(3f12.6,'' = Total Sq.xy norm,rho,mat rms'')') 
     X			tnorm,sqrt(rho2),rmst
	   write(6,'(3f12.6,'' = Total Sq.xy norm,rho,mat rms'')') 
     X			tnorm,sqrt(rho2),rmst
	   deallocate (wfnxy,donell)
	endif

	  if(ISNONO.and.vsearch==0) then
	  rvec(1:naa) = rhs(1:naa,kv)
	  do i=1,naa
	    an = rvec(i)
	    if(neigs==naa) an = 0.
	    do j=1,naa
	     an = an + kk(i,j) * rvec(j)
	    enddo
	    rhs(i,kv) = an
	  enddo
	  eigen = ' proj'
	  endif
	enddo
        enddo
950	if(allocated(aa)) deallocate(aa,rhs,ipiv)
	if(allocated(kk)) deallocate(kk)
	deallocate(rvec,ybas,alpha,gswf,val,phys)
