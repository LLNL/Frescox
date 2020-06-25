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
**LAGRANGE************************************************************
        module meshes
! Lagrange meshes:
        real*8,allocatable:: leg(:),legd(:),xi(:),ybas(:,:),lam(:),
     x                     fddi(:,:),rad(:)
        end module meshes

	SUBROUTINE LAGRANGE(JTOTAL,NCH,MINTL,INITL,ECM,COEF,HP,ITC,CHL,
     X	  PART,EXCIT,LVAL,JVAL,JPROJ,JTARG,N,FORMF,NF,ITCM,EOFF,
     X    CUTVAL,ISOCEN,bndx,weak,nbas,naa, KO,CDETR,WOUT,KFUS,
     X    pralpha,PCON,CH,NCHAN,KS,PARITY,PSIGN,K,RMASS,BLOCKD,
     X    CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,symm,CLIST,NCLIST,NFLIST,
     X    SIGJ,JCOEF,XS,FUSL,SMATS,SMATL,CHANS,DONE,JTMAX,JTMIN,SCALE,
     X    JEX,ABSEND,DONES,NTHRESH,RESM,THRJ,IF1,IF2,TOTJ,FUSJ,meigsi,
     X    MR,CHNO,NL,NLO,NLN1,NLC,MLT,MLM,ICUTC,ISNONO,FLIP,EXCH,GAP,
     X    ETA,FORML,AFRAC,QNF,MASS,DROPPED,NSTEPD,iams,TIME0,TIME0J,
     X    XCOEF,PTYPE,CFUSJ,SMALLCOUP,SMALLCHAN,SMALLS,CHPRES,RTURN,NEX,
     X    VEFF,FNC,WNM,NSA,WAVES,WDISK,CFG,WRITTEN,ENLAB,phasin,linel,
     x    IEXCH,FUSLL,IT0,ITER,IPS,say,ETOTAL)
	use parameters
	use searchpar
	use drier
	use fresco1, only: vsearch,echan,enodes,jleast,btype
	use meshes
      	implicit none
	integer KO,CDETR,NF,N,NCH,MINTL,ITC(MXP,MXX),ISOCEN,nbas,ITCM
	integer I,PART(MAXCH,3),LVAL(NCH),INITL(MINTL),IF,ia,JIN,nd,
     X	  EXCIT(MAXCH,3),IC,C,C2,EL,L,IT,CUTVAL(NCH),nodes,IEXCH,imp,
     X    PCON,IFAULT,ib,naa,j,ki,KS,KS1,PARITY,nrbmin,nb,meigs,meigsi,
     X    NCHAN,PEL,EXL,SMATL,SMATS,DONE,DONES,NTHRESH,CHANS,IF1,IF2,
     X    NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH),NC,neigs,
     X    MR,NL,CHNO(MFNL,6),NLO,NLC,MLT,MLM,ICUTC,mc,kv,kvp,kvo,nlw,
     X	  QNF(19,MSP),kn,qnn(6),NSTEPD,iams,nab,JF,KFUS,linel(MINTL),
     X	  PTYPE(12,MLOC),C1,NA,SMALLS(MXPEX),CHPRES(MXPEX),NEX(MXP),
     x    NSA,WAVES,VEFF,WDISK,IEX,IMA,L1,NICH,M,NLN,NLN1,nrbases,ip,
     x    DROPPED,ki1,ki2,ib2,mag,mode1,kfn,ifail,kkdim,IT0,ITER,ITNL,
     x    ik,kkk,kj,ck,ij,ncha
	logical pralpha,prrm,prkm,prtm,prbut,changed,symm,pauli,prham,
     X 		FCWFN,FJSWTCH,ISNONO,tr,BLOCKD(NCH),FLIP,prba,prmats,
     X		WOUT,FAIL,SMALLJ,WRITTEN(299),r_added,chweak(MAXCH),
     X          rm_set(mterms),drop(NCH),allorder_exch,nonel(MAXCH),red,
     x          adjusted,iterate,joined,say
	real*8 jtotal,JCOEF,ECM(MAXCH,3),JVAL(NCH),JPROJ(NCH),EOFF,
     X		COEF(NCH),HP(NCH),beta(NCH),CRCRAT(MAXCH),ABSEND,
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),Eshift,bndx,
     X		 JEX(6,MXP,MXX),TOTJ,OUTJ,FUSJ,EPS,E1,AN,E,ELAST,WN,
     X		 EXCH(MAXCH,MAXCH),PNORMS(MXP),FORML(MAXNLN,MSP,2),
     X		 AFRAC(MXPEX,MXPEX,2,MSP),RINTP,tnorm,rms(2),rmst,
     X		 MASS(4,MXP+1),ac,an1,an2,rho2,rcore,am12,am123,
     X		 TM,TIME0J,SECOND,TIME0,g,XCOEF,PI,T4,SMALLCOUP,
     X		 SMALLCHAN,CHSIZES(MXPEX),RTURN,GAP,ETA(MXP,MXX),
     X           ENLAB,FCOUL(NCH,NCH,2),FCOULP(NCH,NCH,2),
     X           rm_energy(mterms),rm_vec(nch,mterms),JTARG(NCH),
     x           rtrace,weak,rmrad1,xsc,FUSLL(LMAX1,2),IPS,ETOTAL
	REAL*8 S,T,kp,conv,K(MXP,MXX),RMASS(MXP),SIGJ(0:MXPEX),r,
     X 	 XS(NCH),JTMIN,JTMAX,THRJ(NTHRESH),RESM(NTHRESH),FUSL(1+NFUS1),
     X	 CFUSJ(NFUS),AMDSQS,CFG(MAXMUL,4),CF,CG,Z,phasin(MINTL),
     x   fc_s(lmax1),gc_s(lmax1),fcp_s(lmax1),gcp_s(lmax1),lgd,eta1,
     x   sph,cph
     
	complex*16 FORMF(MAXM,MLOC),CLIST(MAXCH,MAXCH,MCLIST)
   	COMPLEX*16 CHL(LMAX1,MXPEX,2),CH(2,MAXCH),SOMA(NCH,1),
     x     HMINp,HPLp,WFP(nch),WF,HMIN,HPL,phaze,WFN(nch),
     x     FNC(NLN1,NSA),WNM(NLN1,NSA,2),CI,C6,C7,VPOT(N,3),SRC(N,1)
	character PSIGN(3),SCALE(0:MXPEX),COMPARE
	character*5 eigen
**** Local arrays:
	real*8,allocatable:: alpha(:,:),kk(:,:),gswf(:,:),
     x	   pvec(:),qvec(:),vec(:),val(:),phys(:),evec(:,:),aar(:,:),
     x	   kkt(:,:),wfnxy(:,:),ovl(:)
	complex*16,allocatable:: aa(:,:),rhs(:,:),rvec(:),phic(:),yc(:),
     x	   rhsb(:),psi(:,:),aac(:,:,:),rch(:,:),aaa(:,:)
	integer,allocatable:: ipiv(:),ipivc(:,:)
	logical,allocatable:: donell(:,:)
	real*8 phi3(nch),EB,rhsprob(nbas*2)
	complex*16 Rmat(nch,nch),Kmat(nch,nch),Smatr(nch,nch),TC,
     X		logd(nch),Rmati,Rex0,Rmat0,SMAT(nch),Rmmat(nch,nch)
	integer info,nop,iop(nch),failed,mfails,nodest(nch),nbas0(nch),
     X          nbasis(nch),nbmax,NS,IS,ii,ifout,NN
	parameter(imp=60)
	parameter (tr=.true., prba=.false., allorder_exch=.true.)
	 TM(I) = SECOND() - TIME0

C
	pauli = .false.
         if(INH>0) stop ' Lagrange method requires INH=0'
      iterate = ITER>0 .and..not.(flip.and.allorder_exch)
     x                 .and..not.pauli .and..not.(ISNONO.or.vsearch>0)
      if(ITER>0) then
        if(iterate) WRITE(6,*) ' Iterate Lagrange mesh ',ITER,' times!'
        if(.not.iterate) WRITE(6,*) ' No Lagrange iteration allowed'
        endif

!	open(imp,form='formatted')
!	open(imp+1,form='formatted')
	call openif(imp); call openif(65)
      if(pralpha) call openif(imp+1) 
	if(PCON>=3) call openif(62)
	if(prba) call openif(70)
	prrm = pralpha.and.PCON>0 ;prkm = prrm ; prtm = pralpha
	prmats = prrm
	prbut = pralpha
	prham = pralpha .and. PCON>5
 	write(KO,*) ' PRALPHA,PCON,nbas,prrm= ',pralpha,PCON,nbas,prrm
	PI = 4d0*atan(1d0)
	Z = 0d0
	CI = (0d0,1d0)
	EL = INITL(1)
	NICH = NCH
	IEX = NCH
	M = N
	NLN = min(NLN1,(N-1)/MR+1)
	FAIL = .false.
	r_added = .false.
	rm_vec(:,:) = 0.; rm_energy(:)=0.; rm_set(:)=.false.

	rmrad1 = (n-1)*HP(1) 

        do C=1,NCH
        if(btype=='E') then
          T = sqrt(abs(bndx)*conv)
          beta(C) = sign(T,bndx)

        else if(btype=='L') then
          beta(C) = -L/rmrad1 ! so B=-L

        else if(btype=='B') then
          beta(C) = bndx/rmrad1  ! so B=bndx

        else if(btype=='k') then
          beta(C) = bndx

        else if(btype=='S') then
            if(.not.FCWFN)  then
             an2 = abs(CH(1,C))**2
             an1 = real(CH(1,C))*REAL(CH(2,C))+
     x             AIMAG(CH(1,C))*AIMAG(CH(2,C))
            else
             an2 = CGMAT(C,C,1)**2+CFMAT(C,C,1)**2                      ! F^2 + G^2
             an1 = CGMAT(C,C,1)*CGMAT(C,C,2) + CFMAT(C,C,1)*CFMAT(C,C,2)! FF' + GG'
            endif
           beta(C) = an1/an2
           if(pralpha) write(imp,*) 'C,CH:',C,CH(1:2,C),beta(C)
             !! SHIFT(C) = AC * (an1 /an2 - beta(C))
          endif
          if(pralpha) write(imp,*) ' Ch. ',C,': beta =',beta(C),' from '
     x                             ,btype,real(bndx)
         enddo
	  
       Rmat(:,:) = 0.
       if(.not.FJSWTCH) then  ! match asymptotic functions to zero here

	NN = N*5  ! that is enough to find Lagrange nodes
	call setup_mesh(N-1,NN,nbas,rmrad1,PCON)
	nbasis(:) = nbas  ! all the same now !
      do C=1,NCH
	nbas0(C) = (C-1)*nbas
         IT = ITC(PART(C,1),EXCIT(C,1))
         L = LVAL(C)
      CH(1:2,C) = CHL(L+1,IT,1:2)
!		write(60,'(a,i4,4f12.6)') 'CH2: ',C,CH(:,C)
	enddo
	nbmax = nbas

	 beta(:) = 0.
	 do C=1,NCH
	  rmrad1 = (n-1)*HP(C) 
	  conv = -1d0/coef(C)
	  if(pralpha) write(imp,*) ' Channel ',C,': beta =',beta(C)
	 enddo

	  if(pralpha.or.final) then
	  	rmrad1 = (n-1)*HP(1)  ! channel 1
	  	t = beta(1)*rmrad1
	    write(imp,191) bndx,beta(1),t,weak
191	     format('   bndx,beta(1),bp =',3f10.5,';  weak',
     x           ' if penetrability < ',1p,e10.2)
	    endif

!	if(pralpha) write(KO,*) ' H-E matrix size = ',naa
	 write(48,*) ' H-E matrix size = ',naa

	nd = naa
	meigs = meigsi
	if(meigs==0) meigs=2
	if(meigs<0) meigs=naa

	if(.not.iterate) then

     	allocate(aa(nd,nd),rch(nd,max(meigs,nch)),
     x         aaa(nbas,nbas),ipiv(nd))
     	if(ISNONO.or.vsearch>0) then  
   	   allocate(kk(nd,nd))
	   kkdim = nd
	  else
     	   allocate(kk(1,1))
	   kkdim = 1
	  endif

	do c1=1,nch
	            ki1 = nbas0(c1)
	do c2=1,nch
	            ki2 = nbas0(c2)
	call HELAGcc(c1,c2,nch,nbasis,ybas,n,beta,rmrad1,aaa,nbas,
     X	HP,FORMF,CLIST,NCLIST,NFLIST,NF,lval,coef,EL,ECM(1,1),
     X		nbas0,nbas,prba,MR,CHNO,NL,NLO,NLN,NLC,MLT,MLM,
     X		CUTVAL,ICUTC,ISNONO,kk,kkdim,vsearch,joined,PTYPE)
      aa(ki1+1:ki1+nbas,ki2+1:ki2+nbas) = aaa(1:nbas,1:nbas)
     	enddo ! c2
     	enddo ! c1

        drop(:) = .false.
!		Remix exchange terms in the Hamiltonian
       if(flip.and.allorder_exch) then
	if(prham) then
	   write(imp,*) 'H-E matrix before exchange reduction:'
	   call WRCMAT(AA,naa,naa,nd,4,imp)
	   endif
	do c1=1,nch  ! c"
	do c2=1,nch  ! c'
	T = EXCH(c1,c2)  ! exch(c",c')
	drop(c2) = drop(c2).or.abs(T)>1e-10

	  if(abs(T)>1e-10) then
	  do 26 ib=1,nbasis(C2)    ! i'
            ki1 = nbas0(C2) + ib   ! c'i'
	  do 26 ib2=1,nbasis(C1)   ! i"
            ki2= nbas0(C1) + ib2   ! c"i"
!		find overlap of basis states wn=kk(ki2,ki1) or kkt(ib2,ib)
		wn = 0.
		do i=1,n
		 wn = wn + ybas(i,ib) * ybas(i,ib2)
		 enddo
		wn = wn * HP(c2)
	  do 26 c=1,nch
	  do 26 ia=1,nbasis(C1)    ! i
            ki = nbas0(C1) + ia    ! ci
	    if(prham) then 
              WF = aa(ki,ki1) - aa(ki,ki2) * wn * T
	      write(imp,25) ki,ki1,ki2,WF,aa(ki,ki1),
     x				aa(ki,ki2),wn,T
25	    format(' EXCH: ',3i3,2f9.3,' =',2f9.3,' -',2f9.3,' *',f9.3,
     x        ' * ',f6.3)
     		endif
            aa(ki,ki1) = aa(ki,ki1) - aa(ki,ki2) * wn * T
26          continue
	  endif
            
  	enddo
  	enddo
    	DO 70 C2=1,NCH
        IF(drop(C2)) then
	  do ib=1,nbasis(C2)
            ki = nbas0(C2) + ib
    	    aa(ki,:) = 0.0
    	    aa(:,ki) = 0.0
    	    aa(ki,ki) = 1.
    	  enddo
    	  endif
70    	 CONTINUE
       endif

	if(.false.) then  !!!!! DISABLE !!!!
	if(MXP==3.and.ISNONO) then   ! TEMP FIX FOR KNOCKOUT NON-ORTHOGONALITY
	  i=0; ii=0; ib=1
	  do c=1,nch
	  if(i==0.and.part(c,1)==2)  i = ib
	  if(ii==0.and.part(c,1)==3)  ii = ib
	  ib = ib + nbasis(c)
	  enddo
	  write(KO,*) ' Assuming channels ',i,ii,' are knockout equal'
	  kk(i,ii) = 1d0
	  kk(ii,i) = 1d0
	endif
	endif


	  if(pralpha.and.ISNONO) then
     	     allocate(kkt(nd,nd),evec(nd,nd))
	     kkt(:,:) = kk(:,:)
             call HDIAG(kkt,nd,nd,0,evec,nc)
		write(65,*) ' K eigenvalues :'
		write(65,755) (kkt(i,i),i=1,nd)
755		format(10f8.4)
		do i=1,nd
		if(abs(kkt(i,i))>0.1) then
		write(65,756) kkt(i,i)
756		format(' Eigenvector for ',f10.5)
		write(65,755) evec(i,1:nd)
		endif
		enddo
	     call flush(65)
	     deallocate (kkt,evec)
	     endif

	if(pauli) then
	kvp = 0
	do c=1,nch
	if(BLOCKD(c)) kvp = kvp + nbasis(c)
	enddo
		write(KO,*) ' Allocate evec with kvp =',kvp

     	allocate(pvec(naa),qvec(naa),evec(naa,kvp),ovl(kvp))
c --------------------------------------------
C		Project Hamiltonian matrix AA off blocked states
c
c       i.e. matrix elements for QQ.H.QQ + Eshift.(1-QQ)
c	where QQ = 1 - |forbidden><forbidden|
c
C     evec(i) = overlap <SS(i)|forbidden>
	kv = 0
	do 650 c=1,nch
	if(.not.BLOCKD(c)) go to 650
	if(.not.ISNONO) then
	  write(KO,*) ' PAULI BLOCKING NEEDS NON-ORTHOGONALITY OVERLAPS'
	  go to 650
	 endif
	do 640 ib=1,nbasis(c)
	 kv = kv+1				! another state
	   if(kv>kvp) then
	   write(KO,*) ' kvp should be at least  ',kv,' nbas =',nbas
	   stop
	   endif
   	 ki = nbas0(c) + ib
	 evec(1:naa,kv) = kk(1:naa,ki)		! forbidden state
	 evec(ki,kv) = evec(ki,kv) + 1d0		!  + diagonal norm

	an1 = sum(evec(1:naa,kv)**2)
	eps = 1e-8
c			Subtract overlaps with previous basis states:
	do 630 j=1,kv-1
	ac = 0d0
	do 620 i=1,naa
620	ac = ac + evec(i,kv)*evec(i,j)
	ovl(j) = ac
630	evec(:,kv) = evec(:,kv) - ac * evec(:,j)
c
c			Find norm of orthogonalised state:
	t = sum(evec(1:naa,kv)**2)
	  if(tr) write(KO,625)  c,ib,an1,t
625	 format(' Channel ',i3,' basis #',i3,' blocked, norm =',2f10.5)
	  if(tr) write(KO,626)  (ovl(j),j=1,kv-1)
626	 format(10x,10f7.3)
	if(abs(t).lt.eps)  then
      	write(KO,*) 'Proposed PP operator omitted, as norm only',real(t)
	 kv=kv-1
	 go to 640
	 endif
c			Normalise 
	evec(1:naa,kv) = evec(1:naa,kv)/sqrt(t)
c
	tc=0d0
	do i=1,naa
	  pvec(i)=0d0
	  qvec(i)=0d0
	  do ii=1,naa
	    pvec(i)=pvec(i)+AA(i,ii)*evec(ii,kv)
	    qvec(i)=qvec(i)+AA(ii,i)*evec(ii,kv)
	  enddo
	  tc=tc+pvec(i)*evec(i,kv)
	enddo

	Eshift=1000.d0
	tc = tc + Eshift
	do i=1,naa
	do ii=1,naa
	  AA(i,ii)=AA(i,ii)-(evec(i,kv)*qvec(ii)+evec(ii,kv)*pvec(i))
     &    +tc*evec(i,kv)*evec(ii,kv)
	enddo
	enddo
640	continue	
650	continue	
     	deallocate(pvec,qvec,evec)
	endif  ! pauli

	if(prham) then
	   write(imp,*) 'H-E matrix:'
	   call WRCMAT(AA,naa,naa,nd,4,imp)
	   if(ISNONO.or.vsearch>0) then
	   write(imp,*) 'K matrix:'
	   call WRRMAT(kk,naa,naa,nd,8,imp)
 	   endif	
 	 endif	 ! prham

	if(eigens>0) then
		stop 'NO bound state with Lagrange mesh yet'	
!!		include 'boundstate.f'

		return
	  endif ! eigens

C			L.U decomposition of Hamiltonian-E matrix
	if(pralpha) write(KO,*) ' LU decomposition of matrix size ',naa
		
	  call zgetrf(naa,naa,AA,nd,ipiv,info)

	  rch(:,:) = 0d0
	  do 31 c=1,nch
	  if(.not.drop(c)) then
	  do 30 ib=1,nbasis(c)
            ki = nbas0(c) + ib
	    rch(ki,c) = ybas(n,ib)
   30	  continue
   	  endif
   31	  continue
	 	if(prham) then
 	   	write(imp,*) 'rhs source terms:'
 	   	call WRCMAT(rch,naa,nch,nd,4,imp)
  	 	endif	

          call zgetrs('N',naa,nch,AA,nd,ipiv,rch,nd,info)
	   if(info.ne.0) then
	     write(imp,*)' Error return from H-E zgetrs',info
	     stop 'zgetrs'
	   endif ! info/=0
	if(prham) then
	   write(imp,*) 'rhs solution:'
	   call WRCMAT(rch,naa,nch,nd,4,imp)
 	 endif	
	   rhsprob(:) = 0d0
	  do 40 i=1,nch
	  do 40 j=1,nch
	   tc = 0d0
	   do 35 ib=1,nbasis(i)
            ki = nbas0(i) + ib
     	   tc = tc + ybas(n,ib) * rch(ki,j)
		rhsprob(ib) = rhsprob(ib) + abs(rch(ki,j))**2
   35	   continue
	   Rmat(i,j) = - tc * coef(j) !/ rmrad1
   40	  continue
	  T = sum(rhsprob(1:nbmax))
	if(final) write(KO,42) (rhsprob(ib)/T,ib=1,nbmax)
   42	format(' Lagrange bases: uses  =',10f6.3/
     x        ('                        ',10f6.3))
     
      else ! ITERATE R-matrix solutions!!
      
         allocate(aac(nbas,nbas,nch),aaa(nbas,nbas),rhs(nbas,nch),
     x               rch(nd,nch),ipivc(nbas,nch))
     	   allocate(kk(1,1));  kkdim = 1

	do c=1,NCH
	call HELAGcc(c,c,nch,nbasis,ybas,n,beta,rmrad1,AAC(1,1,c),nbas,HP, ! diagonal
     X		FORMF,CLIST,NCLIST,NFLIST,NF,lval,coef,EL,ECM(1,1),
     X		nbas0,nbas,prba,MR,CHNO,NL,NLO,NLN,NLC,MLT,MLM,
     X		CUTVAL,ICUTC,ISNONO,kk,kkdim,vsearch,joined,PTYPE)
	if(prham) then
	   write(imp,*) 'H-E block diagonal for channel ',c
	   call WRCMAT(AAC(1,1,c),nbas,nbas,nbas,4,imp)
 	 endif	 ! prham
 	 enddo
        drop(:) = .false.

C			L.U decomposition of Hamiltonian-E matrix
	if(pralpha) write(KO,*) ' LU decompositions matrix size ',nbas**2		
	 do c=1,NCH
	  call zgetrf(nbas,nbas,AAC(1,1,c),nbas,ipivc(1,c),info)
	 enddo
	 
	 rch(:,:) = 0.0  ! start iterations	 
	 DO ITNL=0,ITER
		if(prham) write(imp,*) '======= ITERATION ',ITNL,' ======'
	 do c=1,NCH  ! iterate each channel c in turn

	 rhs(:,:) = 0.
	  do ib=1,nbasis(c)
!            ki = nbas0(c) + ib
	      rhs(ib,c)= ybas(n,ib)
	   enddo 
	
	 do ck=1,NCH  ! coupled channel k

	 if(ITNL>0.and.c.ne.ck) then
	   call HELAGcc(c,ck,nch,nbasis,ybas,n,beta,rmrad1,  ! off-diagonal
     X      aaa(1,1),nbas,HP,FORMF,CLIST,NCLIST,NFLIST,NF,lval,coef,EL,
     X	ECM(1,1),nbas0,nbas,prba,MR,CHNO,NL,NLO,NLN,NLC,MLT,MLM,
     X		CUTVAL,ICUTC,ISNONO,kk,kkdim,vsearch,joined,PTYPE)
        if(joined) then
	  if(prham) then
	   write(imp,*) 'H-E block  to channel ',c,' from ',ck,':'
	   call WRCMAT(aaa,nbas,nbas,nbas,4,imp)
 	  endif	 ! prham
 	 
	  do ib=1,nbasis(c)
            ki = nbas0(c) + ib

	    do ik=1,nbasis(ck)
	     kkk = nbas0(ck) + ik
	     do j=1,nch
	     rhs(ib,j) = rhs(ib,j) - AAA(ib,ik) * rch(kkk,j)	    
	     enddo ! j
	    enddo ! ik
	    
    	  enddo  ! ib
    	 endif ! joined
 	 endif ! ITNL > 0 & c/=ck
 	 enddo ! ck
    	  
	 	if(prham) then
 	   	write(imp,*) 'rhs source terms to c =',c,' at ',ITNL
 	   	call WRCMAT(rhs,nbas,nch,nbas,4,imp)
  	 	endif	

          call zgetrs('N',nbas,nch,AAC(1,1,c),nbas,ipivc(1,c),rhs,nbas,
     x                info)
	   if(info.ne.0) then
	     write(imp,*)' Error return from H-E zgetrs',info
	     stop 'zgetrs'
	   endif ! info/=0
	 	if(prham) then
 	   	write(imp,*) 'Solutions for c =',c,' at ',ITNL
 	   	call WRCMAT(rhs,nbas,nch,nbas,4,imp)
  	 	endif	
  	 	
	   do ij=1,nbasis(c)
	     kj = nbas0(c) + ij	  
  	    do j=1,NCH
	     rch(kj,j) = rhs(ij,j)
	    enddo ! j
	   enddo  ! ij
	
   	  enddo  !   c

	   rhsprob(:) = 0d0
	  do  i=1,nch
	  do  j=1,nch
	   tc = 0d0
	   do ib=1,nbasis(i)
            ki = nbas0(i) + ib
     	   tc = tc + ybas(n,ib) * rch(ki,j)
		rhsprob(ib) = rhsprob(ib) + abs(rch(ki,j))**2
    	   enddo
	   Rmat(i,j) = - tc * coef(j) 
   	  enddo
   	  enddo

	if(prrm) then
	   write(imp,*) ' r*R matrix at ITNL =',ITNL,':'
          do  i=1,nch
          write(imp,45) i,(Rmat(i,j) , j=1,nch)
          enddo
	  endif ! prrm

    	  ENDDO  !   ITNL
  	  
	  T = sum(rhsprob(1:nbmax))
	if(final) write(KO,42) (rhsprob(ib)/T,ib=1,nbmax)
	
      endif ! ITERATE

	if(prrm) then
	   write(imp,*) ' r*R matrix:'
          do 43 i=1,nch
   43      write(imp,45) i,(Rmat(i,j) , j=1,nch)
   45	   format(1x,i3,18f8.3,:,/(4x,18f8.3))
	   write(imp,*) ' r*R matrix diag:'
   	   write(imp,45) nch,(Rmat(j,j) , j=1,nch)
   	   write(imp,46) (Rmat(j,j)/rmrad1 , j=1,nch)
   46	   format(' :R matrix diag =',12f8.4)
	  endif ! prrm
C

	endif ! .not.FJSWTCH

!		Add in any search R-matrix term!
	do ip=1,nvars
	  if(srch_kind(ip)==3) then
	    if(abs(JTOTAL-srch_jtot(ip))<.1
     x         .and.abs(PARITY-srch_par(ip))<.1) then
              rm_energy(srch_rterm(ip)) = srch_value(ip)
              rm_set(srch_rterm(ip)) = .true.
!              write(imp,151) ki,rm_energy(ki)
 	write(imp,1501) ip,JTOTAL,srch_jtot(ip),PARITY,srch_par(ip),
     x         srch_rterm(ip),srch_value(ip)
 1501	format(' Var ',i2,' for ',2f5.1,2i3,': term #',i3,' val',f10.5)
	      endif
	  elseif(srch_kind(ip)==4.and.rm_set(max(1,srch_rterm(ip))).and.
     x            srch_r_ch(ip)<=nch) then
              rm_vec(srch_r_ch(ip),srch_rterm(ip)) = srch_value(ip)
              if(abs(srch_value(ip))>0.) r_added=.true.
 	write(imp,1501) ip,JTOTAL,srch_jtot(ip),PARITY,srch_par(ip),
     x         srch_rterm(ip),srch_value(ip)
	  endif
	enddo

!	if(prrm.or.r_added) then
	if(prrm) then
	   write(imp,*) ' r*R matrix before additional search terms:'
          do 1515 i=1,nch
 1515      write(imp,45) i,(Rmat(i,j) , j=1,nch)
	  endif ! prrm

	do ki=1,mterms
	if(rm_set(ki)) then
	do 152 i=1,nch
	do 152 j=1,nch
	 if(.not.chweak(j)) then
!	  T = sqrt(RMASS(PART(i,1))/RMASS(PART(j,1)))
	  T = sqrt(COEF(j)/COEF(i))
	  Rmat(i,j) = Rmat(i,j) +  T* (n-1)*HP(i)*
     x    rm_vec(i,ki)*rm_vec(j,ki)/(rm_energy(ki)-ETOTAL)	
	if(pralpha) 
     x  write(imp,151) ki,rm_energy(ki),i,j, rm_vec(i,ki),rm_vec(j,ki),
     x   		ECM(INITL(1),1),an1,an2,T
 151	format(' New R term',i2,' at',f8.4,:,': ',2i3,' with',2f10.5,
     x         ' at ',f8.4,' MeV',3f8.3)
!       if(i>j.and.abs(rmat(i,j))>1e-20) rtrace = Rmat(i,j)
	endif
 152	continue
 	endif
 	enddo


C****************************************************** Weak coupling limit
C                                       if requested, or needed for stability

	if(pralpha.and.r_added) then
	   write(imp,*) ' r*R matrix with additional search terms:'
          do 153 i=1,nch
  153      write(imp,45) i,(Rmat(i,j) , j=1,nch)
! 	write(155,*) ECM(INITL(1),1),abs(rtrace)
	  endif ! prrm
C****************************************************** Find open channels

	call nopen(NCH,ECM(1,1),iop,nop)
	if(pralpha.and.nop<nch) then
  	   write(KO,*) 'Only ',nop,' open channels'
	   write(imp,*) 'Only ',nop,' open channels'
	endif
	
C****************************************************** K Scattering Matrix
C
C    K =  - [G - RG~]^-1 [F - RF~] where RF~ = r R F' - r beta F, etc

	call GETK(nch,iop,nop,CFMAT,CGMAT,FCWFN,CH,MAXCH,
     X			Rmat,beta,LVAL,Kmat,imp,FCOUL,FCOULP,prmats)
C
	if(prkm) then
	   write(imp,*) ' K matrix:'
	  rtrace = 0.
          do 230 i=1,nop
	do 229 j=1,nop
        if(i>j.and.abs(Kmat(i,j))>1d-50) rtrace = Kmat(i,j)
  229   continue
  230      write(imp,231) i,(Kmat(i,j) , j=1,nop)
  	write(175,*) ECM(INITL(1),1),abs(rtrace)
  231	   format(1x,i3,8g14.4,:,/(4x,8g14.4))
	 endif
C******************************************* S Scattering Matrix
C
C    S = [1 + i K] * [1 - i K]^{-1}
C
	call GETS(Kmat,Smatr,prtm,iop,nop,nch,lval,imp)

C******************************************** S Columns for Cross sections
C
      	FUSL(:) = 0.0
	NCHA = NCH
	if(.not.WOUT) NCHA = 1
	allocate(psi(N,NCHA))
	
! done 50	SMAT(I,J) = SMAT(I,J) * (0.,1.)**(LVAL(J)-LVAL(I))
! i.e. 50	SMAT(I,EL) = SMAT(I,EL) * (0.,1.)**(LVAL(EL)-LVAL(I))

	do 600 jin=1,MINTL
	  psi(:,:) = 0.0  
	 EL = INITL(JIN)
         PEL = PART(EL,1)
         EXL = EXCIT(EL,1)
	 SMAT(:) = Smatr(:,EL) / (0.,1.)**(LVAL(EL)-LVAL(:))  ! ie. original SMAT
c
	if(WOUT.and..not.FJSWTCH) then
!		Find radial wf for channel c:
	do 250 c=1,nch
	   HMIN  = dcmplx(FCOUL (c,el,2),-FCOUL (c,el,1))
	   HMINp = dcmplx(FCOULP(c,el,2),-FCOULP(c,el,1))
	WF     = HMIN
	WFP(c) = HMINp
	 do j=1,nch  
	  HPL   = dcmplx(FCOUL (c,j,2),+FCOUL (c,j,1))
	  HPLp  = dcmplx(FCOULP(c,j,2),+FCOULP(c,j,1))
	  WF     = WF     - HPL *SMAT(j)
	  WFP(c) = WFp(c) - HPLp*SMAT(j)
	 enddo
250	WFP(c) = (0d0,.5d0)*(WFP(c) - beta(c)*WF)

	src(:,:) = 0.0  

	do 260 c=1,nch
	do 260 j=1,nch
	  rmmat(c,j) = 0.
	do ib=1,nbasis(c)
          ki = nbas0(c) + ib
	  TC  = rch(ki,j) * WFP(j) * (-coef(j))
!	  if(PRALPHA) write(imp,252) c,j,ib,TC
252	  format('TC for ',3i4,'=',2f10.5)
	if(wout) then
	do i=1,n
     	   psi(i,c) = psi(i,c) + ybas(i,ib) * TC !/ rmrad1
	enddo
	endif
	  rmmat(c,j) = rmmat(c,j) -ybas(n,ib)*rch(ki,j)*coef(j)
        enddo
  260	continue
        if(PRALPHA) write(imp,*) ' R-matrix again:'
          do 261 i=1,nch
  261      if(PRALPHA)write(imp,45) i,(rmmat(i,j) , j=1,nch)
	WFN(:) = 0.0
	do j=1,nch
	WFN(:) = WFN(:) + rmat(:,j)*WFP(j)
         if(PRALPHA) write(imp,2615) WFP(j),WFN(:)
2615    format(' Contribution from ',2f9.5,': ',20f10.5)

	enddo
        if(PRALPHA) then
	  do c=1,NCH
	 T = 0; if(C==EL) T = 1.0
	  HPL   = dcmplx(FCOUL (c,c,2),+FCOUL (c,c,1))
	  HPLp  = dcmplx(FCOULP(c,c,2),+FCOULP(c,c,1))
	   HMIN  = dcmplx(FCOUL (c,c,2),-FCOUL (c,c,1))
	   HMINp = dcmplx(FCOULP(c,c,2),-FCOULP(c,c,1))
	 Rex0 = (T*HMIN - smat(c)*HPL)/(T*HMINp - smat(c)*HPLp)
	 Rmati = WFN(c)/WFP(c)
	 WF = (T*HMIN-smat(c)*HPL)*(0.,.5)
	 C6 = (T*HMINp - smat(c)*HPLp)*(0.,0.5)
	  
!	  write(imp,262) WFP(c),WFN(c),psi(n,c)
	  write(imp,262) WFP(c),WFN(c),C6,psi(n,min(c,NCHA)),WF,Rex0,Rmati
262     format(' WFP =',2f10.5,' WF =',2f10.5,' psi-ext''='2f10.5,:,
     x     ' Psi-int,ext ',4f10.5,
     x  ';   R-ext =',2f10.5,'  R-int =',2f10.5)
	  enddo ! c
	  endif
         endif

	 DO C=1,NCH
	 SMAT(C) = Smatr(C,EL)  ! With GETS (0.,1.)**(LVAL(EL)-LVAL(C))
	 if(WOUT) psi(:,C) = psi(:,C) * (0.,1.)**(LVAL(EL)-LVAL(C))
	 ENDDO

!	C = 1         
!       ac = sqrt(-ECM(C,1)/coef(C))	! k
!       do I=1,N
!        R  = (I-1)*HP(C)
!        write(93,264) R,PSI(I,C),sin(ac*R)
!        enddo
!264	  format(f8.3,3f12.6)
!	  write(93,*) '&'

!		Remix exchange contributions.
       if(flip.and..not.allorder_exch) then
        drop(:) = .false.
	do c1=1,nch
	do c2=1,nch
	T = EXCH(c1,c2)
	SMAT(C1) = SMAT(C1) + T * SMAT(C2)
	drop(c2) = drop(c2).or.abs(T)>1e-10
        IF(SMATS.GE.5 .AND. ABS(T*SMAT(C2)).GT.1E-6)
     &                WRITE(KO,265) C1,SMAT(C1),T,C2,SMAT(C2)
265      FORMAT(' S-mat(',I3,') =',2F10.6,' after',
     &         '  adding',F9.5,' times S-mat(',I3,') =',2F10.6)
        IF(C1<=NICH.and.C2<=NICH.and.WOUT.and..not.FJSWTCH) then
          PSI(1:N,C1) = PSI(1:N,C1) + T * PSI(1:N,C2)
          endif
  	enddo
  	enddo
    	DO 270  C2=1,NCH
        IF(.NOT.DROP(C2)) GO TO 270
    	  SMAT(C2) = 0.0
    	  IF(C2<=NICH.and.WOUT.and..not.FJSWTCH) then
    	  PSI(1:N,C2) = 0.0
    	  ENDIF
270   	 CONTINUE
       endif

	phaze = log(SMAT(EL))*(0.,-0.5)*180./pi
	phasin(JIN) = phaze
	linel(JIN) = EL

      if(final) then
      IF(SMATL.GE.2)
     X     WRITE(KO,1420) JTOTAL,PSIGN(PARITY+2),EL,SMAT(EL)
      IF(SMATL.GE.3) 
     X    WRITE(KO,1430) (SMAT(C),C=1,NCH)
 1420 FORMAT(' Final S-matrices (',F7.1,A1,') Sel(',I3,') =',
     X			F10.5,' +i*',F8.5)
 1422 FORMAT(F10.1,2F12.8,'i: elastic S-matrix  @@',f10.2,i3,i4,l2,i3)
 1430 FORMAT(5(1X,F11.5,' +i*',F9.5,',') )
!1431 FORMAT(5(1X,1p,e11.3,' +i*',e9.1,',') )
      IF(SMATL.GE.2) then
         WRITE(KO,1422) JTOTAL,SMAT(EL),TM(I)-TIME0J,
     X                  DROPPED,NSTEPD,FJSWTCH,iams
         WRITE(KO,195) EL,phaze,lval(el),jval(el)
195      FORMAT( ' Elastic phase shift ',I3,' = '
     &     , 2F8.3,' deg. for the L =',I5,', J =',F7.1,' channel.')
          WRITE(45,196) ECM(EL,1),phaze,lval(el),jval(el)
196        format(f10.3,2f9.3,' for LJin =',i6,f6.1)
          written(45) = .true.
	  endif
	endif

	nrbases = 1
	include 'usescatwf.f'
	
      if(CDETR.GT.1) CDETR = CDETR - 1
C                next JIN    :
600   enddo
	if(WOUT) deallocate(psi)

	if(allocated(aa)) deallocate(aa,aaa,rch,ipiv)
	if(allocated(aac)) deallocate(aac,aaa,rhs,rch,ipivc)
	if(allocated(kk)) deallocate(kk)
	if(allocated(ybas)) deallocate(ybas,lam,fddi)
	return
	END SUBROUTINE LAGRANGE

	Subroutine HELAGcc(c1,c2,Nstate,nsturm,SS,numr,beta,rmrad,
     X		AA,nd,hcm,FORMF,CLIST,NCLIST,NFLIST,NF,LCHL,coef,
     X 		EL,ECM,nbas0,nbmax,prba,MR,CHNO,NL,NLO,NLN,
     X		NLC,MLT,MLM,CUTOFF,ICUTC,ISNONO,kk,kkdim,vsearch,
     X		joined,PTYPE)
     
!      Couplings just to c1 from c2 !  Ignore symmetries ! 
!
	use parameters
	use meshes, f=>ybas
	implicit real*8 (a-h,o-z)
	real*8 ECM(Nstate),kk(kkdim,kkdim)
	real*8 SS(numr,nbmax),beta(Nstate)
	real*8 NN,hcm(Nstate),wrad(numr),coef(Nstate),
     X	  FNL(NLN,NLO),VNL(1:numr,1:MLM),EXF(NLN),VR,wwsrch(numr)
	integer LCHL(Nstate),EL,nsturm(Nstate),nstmax,nbas0(Nstate),
     X    NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH),
     X    MR,NL,CHNO(MFNL,6),NLO,NLN,NLC,MLT,MLM,ICUTC,C,D,
     X    CUTOFF(Nstate),vsearch,c1,c2,PTYPE(12,MLOC)
	complex*16 FORMF(MAXM,MLOC),CLIST(MAXCH,MAXCH,MCLIST),VV,VV1
	complex*16 AA(nd,nd),ww(numr),SSC(nbmax),wwd(numr),
     X	  FNC(NLN,NLO),VNC(1:numr,1:MLM),ECF(NLN)
	logical pertcent,coupled,blas1,blas2,FFR,SH,NREV,ISNONO,
     X		prba,incl,forward,reverse,joined,derd
	parameter (pertcent=.false.,blas1=.true.,blas2=.false.,
     X  	   SH=.false.)
!			Adjust blas1 and blas2 for best times in your system
!	ascale = 1.0/sqrt(rmrad)
c				wrad is for integration [0,rmax], without hcm factor
	    wrad(1) = 1d0/3d0
	  g=4.d0/3.d0
	  do  i=2,numr
	    if(i.eq.numr) g = 1d0/3d0
	    wrad(i) = g
	    g=2d0-g
	   enddo

	if(prba) then
	call openif(70)
	do l1=1,nbmax
	  write(70,'(''#  Basis state '',i4,'' in ch 1'')') l1
	  do i=1,numr
	  write(70,*) (i-1)*real(hcm(1)),real(SS(i,l1))   !,wrad(i)
	  enddo
	  write(70,*) '&'
	  call flush(70)
	  !stop
	enddo
c --------------------------------------------
c       matrix elements for MM=<SS|SS> delta(lji)
c --------------------------------------------
        do l1=1,nsturm(kv1)
          do l2=1,nsturm(kv2)
            AA(l1,l2)=0d0
            do i=1,numr
              AA(l1,l2)=AA(l1,l2)+wrad(i)*hcm(kv1)*
     x                  SS(i,l1)*SS(i,l2)
            enddo
          enddo
        enddo

           write(60,*) 'MM normalisation matrix:',rmrad
           call wrcmat(AA,nd,nd,nd,4,60)
         endif

	AA(:,:) = 0d0
      N = numr
      joined = .false.
c --------------------------------------------
c       matrix elements for AA=<SS|T|SS> : kinetic energy
c --------------------------------------------
	if(c1==c2) then
	 joined = .true.
	 kv1=c1
	      L = LCHL(kv1)
	      T = -coef(kv1)/ rmrad**2
	do i=1,nsturm(kv1)
	do j=1,nsturm(kv1)
	del=0.0; if(i==j) del=1.0
	AA(i,j) = T*(-sqrt(lam(i))*fddi(j,i)+L*(L+1)*del/xi(i)**2)
	enddo
	enddo
!		write(60,*) 'T (KE) matrix:'
!	   call WRCMAT(AA,nd,nd,nd,4,60)
c --------------------------------------------
c       matrix elements for AA=<SS|L|SS> : Bloch operator
c --------------------------------------------
	kv1=c1
	      T = -coef(kv1)/ rmrad**1.5d0
	do i=1,nsturm(kv1)
	do j=1,nsturm(kv1)
      fdash=(-1)**j*sqrt((1.d0-xi(j))/xi(j))* 
     &    (-leg(N)/((1.d0-xi(j))**2)+leg(N)/(1.d0-xi(j)) 
     &   +2.d0*legd(N)/(1.d0-xi(j)))
	AA(i,j) = AA(i,j) + f(N,i)*(fdash-rmrad*beta(kv1)*f(N,j))*T
	enddo
	enddo
	
	if(.false.) then
 		   write(60,*) 'T+L matrix:'
 	   call WRCMAT(AA,nd,nd,nd,4,60)
	write(60,110) beta
110	format('Beta:',10f10.5)
	do i=1,nd
	do j=1,i-1
	if(abs(AA(i,j)-AA(j,i))>1e-5) write(60,111) i,j,AA(i,j),AA(j,i)
111	format('EA',2i5,4f15.6)
	enddo
	enddo
	endif  ! print T+L
	endif  ! c1==c2
c --------------------------------------------
c       matrix elements for AA=<SS|V2|SS>*YY2
c --------------------------------------------
	
	if(ISNONO.or.vsearch>0) kk(:,:) = 0d0
 	kv1=c1
	  ryev = -coef(kv1)
	  kv2l=Nstate
	kv2=c2
	  ww(:) = 0d0; wwd(:)=0d0
	  if(vsearch>0) wwsrch(:) = 0.
	  coupled=.false.; derd=.false.
	  do NC=1,NCLIST(kv1,kv2)
	   ifm = NFLIST(kv1,kv2,NC)
!	    incl = vsearch==0 .or. vsearch>0.and.ifm.ne.vsearch
	    incl = ifm.ne.vsearch
	   if(abs(CLIST(kv1,kv2,NC))>1e-20) then
	     coupled=.true.;  joined = .true.
            if(PTYPE(7,ifm)/=1) then ! Scalar potential
	    if(incl) then
	     if(blas1) then
	       call zaxpy(numr,CLIST(kv1,kv2,NC),FORMF(1,ifm),1,ww,1)
	     else
	       ww(:) = ww(:) + CLIST(kv1,kv2,NC)*FORMF(1:numr,ifm)
	     endif  !blas1
	    else  ! incl
	     wwsrch(:) = wwsrch(:) + CLIST(kv1,kv2,NC)*FORMF(1:numr,ifm)
	    endif  ! incl
           else   ! derivative potential
             if(incl) then
               wwd(:) = wwd(:) + CLIST(kv1,kv2,NC)*FORMF(1:numr,ifm)
               derd = .true.
             else
                write(6,*) ' Searching on derivatives not implemented'
                stop
             endif
            endif  ! PTYPE(7,IF)/=1	    
	    endif  ! CLIST
	   enddo  ! NC
	  if(coupled) then
	  do i=1,numr
	    rrad = (i-1)*hcm(kv1)
	       VV = ww(i) 
	      if(pertcent.and.kv1.eq.kv2) 
     & 		 VV = VV + ryev*(LCHL(kv1)*(LCHL(kv1)+1)
     & 		                -LCHL(EL)*(LCHL(EL)+1))/rrad**2
	      VV = VV * wrad(i) * hcm(kv1) 
     & 			* (0.,1.)**(LCHL(kv1)-LCHL(kv2))
	    do l1=1,nsturm(kv1)
	      k1=l1
	      VV1 = VV * SS(i,l1)
	      l2max = nsturm(kv2)
c					l2max so k2.le.k1 if symm.
	    if(l2max>0) then
	    if(blas2) then
!				SS is real not complex!
	        SSC(1:l2max) = SS(i,1:l2max)
	        call zaxpy(l2max,VV1,SSC,1,AA(k1,k20+1),nd)
	      else
!	    	do l2=1,l2max
!	     	AA(k1,l2) = AA(k1,l2) + SS(i,kv2,l2) * VV1
	     	AA(k1,1:l2max) = AA(k1,1:l2max) + SS(i,1:l2max) * VV1
!             	enddo
	      endif ! blas2
	      endif ! l2max>0
             enddo ! l1
	  enddo ! i
	
!                               acting on WF derivatives to right (less optimised code!)
         if(derd) then
          do i=2,numr-1
            rrad = (i-1)*hcm(kv1)
               VV = wwd(i) 
              if(pertcent) stop 'percent & derivatives not implemented'
              VV = VV * wrad(i) * hcm(kv1) / hcm(kv2) / 2.0
     &                  * (0.,1.)**(LCHL(kv1)-LCHL(kv2))
            do l1=1,nsturm(kv1)
              k1=l1
              VV1 = VV * SS(i,l1)
              l2max = nsturm(kv2)
c                                       l2max so k2.le.k1 if symm.
	     	AA(k1,1:l2max) = AA(k1,1:l2max) + 
     x            (SS(i+1,1:l2max)-SS(i-1,1:l2max)) * VV1vi fxx
             enddo ! l1
          enddo ! i
        endif ! derd	
	
!			 put search potential into kk not AA
	  if(vsearch>0) then
	   write(481,*) ' Channels ',kv1,' <- ',kv2,':  wwsrch ='
	   write(481,481) real(wwsrch(1:numr))
481	   format(1p,10e12.4)
	   sig = sign(1d0,sum(wwsrch(1:10)) )
	  do i=1,numr
!				make more sign-definite
!		if(abs(wwsrch(i))<1e-20) wwsrch(i)=sign(1d-20,sig)
	    rrad = (i-1)*hcm(kv1)
	       VV = - wwsrch(i) ! negative because V is eigenvalue
	      VV = VV * wrad(i) * hcm(kv1) 
     & 			* (0.,1.)**(LCHL(kv1)-LCHL(kv2))
	    do l1=1,nsturm(kv1)
	      k1=l1
	      VV1 = VV * SS(i,l1)
	      l2max = nsturm(kv2)
	    if(l2max>0) then
	     	kk(k1,1:l2max) = kk(k1,1:l2max) + SS(i,1:l2max) * VV1
	      endif ! l2max>0
             enddo ! l1
	    enddo ! i
	   endif ! vsearch
	  endif ! coupled

c -------------------------------------------------
c       matrix elements for AA=<SS|V(nonlocal)|SS>
c -------------------------------------------------

      HI = 1.0 / DBLE(MR)
      N = numr
      
      LASTNL = 0
       DO  INL=1,NL
         D = ABS(CHNO(INL,1))
         C = CHNO(INL,2)
    	
         NONO = MOD(CHNO(INL,3),10)
         NREV = CHNO(INL,1).LT.0 .OR. C.EQ.D
         
    	   forward = c1==D .and. c2==C 
    	   reverse = c1==c .and. c2==D .and..not.NREV
    	   if(forward.or.reverse) LASTNL = INL ! latest one needed
    	  ENDDO
      
       if(LASTNL>0) REWIND 12
      DO 50 INL=1,LASTNL
         FFR = CHNO(INL,3).LT.10
       NLL = CHNO(INL,4)

         D = ABS(CHNO(INL,1))
         C = CHNO(INL,2)
    	
         NONO = MOD(CHNO(INL,3),10)
         NREV = CHNO(INL,1).LT.0 .OR. C.EQ.D
         
    	   forward = c1==D .and. c2==C 
    	   reverse = c1==c .and. c2==D .and..not.NREV
    	   coupled = forward.or.reverse
    	   joined = joined .or. coupled
    	   if(.not.coupled) then
    	   	 read(12) ! skip this record
    	       go to 50
    	       endif
    	   
    	   IF(FFR) then
	    READ(12) ((FNL(I,J),I=1,NLL),J=1,NLO)
	    VNL(:,:) = 0d0
	   else
            READ(12) ((FNC(I,J),I=1,NLL),J=1,NLO)
	    VNC(:,:) = 0d0
	   endif
    	   
         IF(SH) WRITE(60,22) INL,D,C,NLL,.NOT.NREV,FFR,NONO
22       FORMAT(' NL coupling #',I4,' to',I4,' from Ch.',2I4,
     &         ', Reverse=',L2,', Real =',L2,', NONO =',I2)
            SQH = SQRT(hcm(C)/hcm(D))
        IF(SH.AND.FFR) CALL DISPLY(FNL,NLL,NLO,NLN,SCALE)
        IF(SH.AND..not.FFR) CALL DISPLR(FNC,NLL,NLO,NLN,SCALE)
         DO 31 J=1,MLM
            JJ = (J - NLC*MLT  - 1)
            IMIN = 1 + MAX(-JJ ,0)  +  MAX(CUTOFF(D),ICUTC)
            IMAX = N - MAX(JJ ,0)   -  5
            KMIN = MAX((IMIN-1)/MR+2,2)
            KMAX = MIN(IMAX/MR,NLL-2)
            DO 31 II=1,MR
               P = (II-1)*HI
              P1 = P - 1.
              P2 = P - 2.
              Q  = P + 1.
              X = P * P1 / 6.0  * hcm(C) * SQH
              Y = Q * P2 * 0.5  * hcm(C) * SQH
         IF(FFR) THEN
            IF(II.EQ.1) CALL EXPAND(FNL,NLN,NLL,NLO,EXF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 29 K=KMIN,KMAX
             VR = (-P2*EXF(K-1)+Q*EXF(K+2))*X + (P1*EXF(K)-P*EXF(K+1))*Y
            VNL(I,J) = VR
29           I = I + MR
          ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 30 K=KMIN,KMAX
             VV = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
            VNC(I,J) = VV
30           I = I + MR
         ENDIF
31      CONTINUE
        IF(SH.AND.FFR) CALL DISPLY(VNL,numr,MLM,numr,SCALE)
        IF(SH.AND..not.FFR) CALL DISPLR(VNC,numr,MLM,numr,SCALE)

!	A(kv1,l1;kv2,l2) = * + <SS(kv1,l1,:) | V | SS(kv2,l2,:)>
!                        = * + int(i,j) <SS(kv1,l1,i) | V(i,j) | SS(kv2,l2,j)>
	 kv1 = D
	 kv2 = C
	      k10 = 0
	      k20 = 0
	    do l1=1,nsturm(kv1)
	      k1=k10 + l1
	      l2max = nsturm(kv2)
	    do l2=1,l2max
		VV1 = 0d0
	 	VR = 0d0
		do J=1,MLM
             	JJ = (J - NLC*MLT  - 1)
             	IMIN = 1 + MAX(-JJ ,0)  +  MAX(CUTOFF(D),ICUTC)
             	IMAX = N - MAX(JJ ,0)   -  5
         	 IF(FFR) THEN
	           do I=IMIN,IMAX
	            VR = VR + SS(I,l1)*VNL(I,J)*SS(I+JJ,l2)
         	   enddo
          	 ELSE
	           do I=IMIN,IMAX
	            VV1 = VV1 + SS(I,l1)*VNC(I,J)*SS(I+JJ,l2)
         	   enddo
	         endif
	        enddo
	         if(FFR) VV1 = VR
		VV1 = VV1 * HCM(kv1)
	    IF(NONO>1) then  			 ! Do normal couplings
	     	if(forward) AA(k1,k20+l2) = AA(k1,k20+l2) + VV1
	     	if(reverse) AA(k20+l2,k1) = AA(k20+l2,k1) + VV1
	     else				 ! Do NONO overlaps
	       if(vsearch==0) then
	     	if(forward) kk(k1,k20+l2) = kk(k1,k20+l2) + VV1
	     	if(reverse) kk(k20+l2,k1) = kk(k20+l2,k1) + VV1
               endif
		if(ISNONO.and.vsearch==0) then
			stop 'NONO not implemented in Lagrange mesh'
			EK = 0 ! replace by KE operator to fix !
		 if(NONO==0) then  	! NONO post contribution to AA
	     	   if(forward) AA(k1,k20+l2) = AA(k1,k20+l2) +
     x				VV1*(EK-ECM(kv1))
	     	   if(reverse) AA(k20+l2,k1) = AA(k20+l2,k1) + 
     x				VV1*(EK-ECM(kv1))
		 else			! NONO prior contribution to AA
	     	   if(forward) AA(k1,k20+l2) = AA(k1,k20+l2) + 
     x				VV1*(EK-ECM(kv2))
	     	   if(reverse) AA(k20+l2,k1) = AA(k20+l2,k1) + 
     x				VV1*(EK-ECM(kv2))
		 endif
		endif
	     endif
            enddo
            enddo

50      CONTINUE
c ---------------------------
c -----------------------------------------------------
c       matrix elements for AA=Hamiltonian Matrix - E
c -----------------------------------------------------
	if(c1==c2) then
 	kv1=c1
	do l1=1,nsturm(kv1)
	    AA(l1,l1)=AA(l1,l1) - ECM(kv1)
          enddo
	endif ! c1==c2

	RETURN
	END
