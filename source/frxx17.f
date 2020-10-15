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
	SUBROUTINE RMATRIX(JTOTAL,NCH,MINTL,INITL,ECM,COEF,HP,ITC,CHL,
     X	  PART,EXCIT,LVAL,JVAL,JPROJ,JTARG,N,FORMF,NF,ITCM,EOFF,
     X    CUTVAL,ISOCEN,bndx,weak,nbasi,nrbmin,nd1,KO,CDETR,WOUT,KFUS,
     X    pralpha,PCON,CH,LL1,NCHAN,KS,PARITY,PSIGN,K,RMASS,BLOCKD,
     X    CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,symm,CLIST,NCLIST,NFLIST,
     X    SIGJ,JCOEF,XS,FUSL,SMATS,SMATL,CHANS,DONE,JTMAX,JTMIN,SCALE,
     X    JEX,ABSEND,DONES,NTHRESH,RESM,THRJ,IF1,IF2,TOTJ,FUSJ,meigsi,
     X    MR,CHNO,NL,NLO,NLN1,NLC,MLT,MLM,ICUTC,ISNONO,FLIP,EXCH,GAP,
     X    ETA,FORML,AFRAC,QNF,MASS,DROPPED,NSTEPD,iams,TIME0,TIME0J,
     X    XCOEF,PTYPE,CFUSJ,SMALLCOUP,SMALLCHAN,SMALLS,CHPRES,RTURN,NEX,
     X	  VEFF,FNC,WNM,NSA,WAVES,WDISK,CFG,WRITTEN,ENLAB,phasin,linel,
     x    IEXCH,FUSLL,CPSO,say,ETOTAL)
	use parameters
	use searchpar
	use drier
	use fresco1, only: vsearch,echan,enodes,jleast,btype,eobs
	use parallel, only: mpisets
      	implicit none
	integer KO,CDETR,NF,N,NCH,MINTL,ITC(MXP,MXX),ISOCEN,nbas,ITCM
	integer I,PART(MAXCH,3),LVAL(NCH),INITL(MINTL),IF,ia,JIN,nd,
     X	  EXCIT(MAXCH,3),IC,C,C2,EL,L,IT,CUTVAL(NCH),nodes,IEXCH,imp,
     X    PCON,IFAULT,ib,naa,j,ki,KS,KS1,PARITY,nrbmin,nd1,meigs,meigsi,
     X    NCHAN,PEL,EXL,SMATL,SMATS,DONE,DONES,NTHRESH,CHANS,IF1,IF2,
     X    NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH),NC,neigs,
     X    MR,NL,CHNO(MFNL,6),NLO,NLC,MLT,MLM,ICUTC,mc,kv,kvp,kvo,nlw,
     X	  QNF(19,MSP),kn,qnn(6),NSTEPD,iams,nab,JF,KFUS,linel(MINTL),
     X	  PTYPE(12,MLOC),C1,NA,SMALLS(MXPEX),CHPRES(MXPEX),NEX(MXP),
     x    NSA,WAVES,VEFF,WDISK,IEX,IMA,L1,NICH,M,NLN,NLN1,nrbases,ip,
     x    DROPPED,ki1,ki2,ib2,nbasi,mag,mode1,kfn,ifail,kkdim,ich,ip2
	logical pralpha,prrm,prkm,prtm,prbut,changed,symm,pauli,prham,
     X 		FCWFN,FJSWTCH,ISNONO,tr,BLOCKD(NCH),FLIP,prba,prmats,
     X		WOUT,FAIL,SMALLJ,WRITTEN(299),r_added,chweak(MAXCH),
     X          rm_set(mterms),drop(NCH),allorder_exch,nonel(MAXCH),red,
     x          adjusted,em,CPSO,say
	real*8 jtotal,JCOEF,ECM(MAXCH,3),JVAL(NCH),JPROJ(NCH),EOFF,
     X		COEF(NCH),HP(NCH),beta(NCH),CRCRAT(MAXCH),ABSEND,
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),Eshift,bndx,
     X		 JEX(6,MXP,MXX),TOTJ,OUTJ,FUSJ,EPS,E1,AN,E,ELAST,WN,
     X		 EXCH(MAXCH,MAXCH),PNORMS(MXP),FORML(MAXNLN,MSP,2),
     X		 AFRAC(MXPEX,MXPEX,2,MSP),RINTP,tnorm,rms(2),rmst,
     X		 MASS(4,MXP+1),ac,an1,an2,rho2,rcore,am12,am123,
     X		 TM,TIME0J,SECOND,TIME0,g,XCOEF,PI,T4,SMALLCOUP,
     X		 SMALLCHAN,CHSIZES(MXPEX),RTURN,GAP,ETA(MXP,MXX),
     X           ENLAB,FCOUL(NCH,NCH,2),FCOULP(NCH,NCH,2),fracshift,
     X           rm_vec(nch,mterms),JTARG(NCH),gameps,ETOTAL,damp,
     x           rtrace,xsc,PMAX(MAXCH),vmin(NCH),weak,rmrad1,
     x           HARDSPH(MAXCH),SHIFT(MAXCH),ggam,esh,LL1(NCH),EOF,
     x           FUSLL(LMAX1,2),rmdist(NCH),SV,JV,racah,SYMR(NCH),
     x           ET(2),penet(2),E_Brune2(0:mterms),W_Brune2(0:mterms)
        complex*16  rm_energy(mterms)
	REAL*8 S,T,kp,EIG,conv,K(MXP,MXX),RMASS(MXP),SIGJ(0:MXPEX),r,
     X 	 XS(NCH),JTMIN,JTMAX,THRJ(NTHRESH),RESM(NTHRESH),FUSL(1+NFUS1),
     X	 CFUSJ(NFUS),AMDSQS,CFG(MAXMUL,4),CF,CG,Z,phasin(MINTL),
     x   fc_s(lmax1),gc_s(lmax1),fcp_s(lmax1),gcp_s(lmax1),lgd,eta1,
     x   sph,cph,pene,rwamp,unitar(NCH),SYMS(NCH),ETACO(NCH),ENFF(NCH),
     x   RMA(nch),Pshift(nch,mterms),EobsShift,penefac,penep

	complex*16 FORMF(MAXM,MLOC),CLIST(MAXCH,MAXCH,MCLIST)
   	COMPLEX*16 CHL(LMAX1,MXPEX,2),CH(2,MAXCH),SOMA(NCH,1),
     x     HMINp,HPLp,WFP(nch),WF,HMIN,HPL,phaze,WFN(nch),
     x     FNC(NLN1,NSA),WNM(NLN1,NSA,2),CI,C6,C7,VPOT(N,3),SRC(N,1)
	character PSIGN(3),SCALE(0:MXPEX),COMPARE
	character*5 eigen,wid
**** Local arrays:
	real*8,allocatable:: ybas(:,:,:),alpha(:,:),kk(:,:),gswf(:,:),
     x	   pvec(:),qvec(:),vec(:),val(:),phys(:),evec(:,:),aar(:,:),
     x	   kkt(:,:),wfnxy(:,:),ovl(:),schro(:)
	complex*16,allocatable:: aa(:,:),rhs(:,:),rvec(:),phic(:),yc(:),
     x	   rhsb(:),psi(:,:)
	integer,allocatable:: ipiv(:)
	logical,allocatable:: donell(:,:)
	real*8 bpot(n,nch),phi3(nch),EB,rhsprob(nbasi*2)
	real*8 emass(n,nch),rm(n),drm(n),d2rm(n),vemass(n,nch)
	complex*16 Rmat(nch,nch),Kmat(nch,nch),Smatr(nch,nch),TC,
     X		logd(nch),Rmati,Rex0,Rmat0,SMAT(nch),Rmmat(nch,nch),
     x          Ssym(nch,nch)
	integer nodmin(nch),info,nop,iop(nch),failed,mfails,nodest(nch),
     X          nbasis(nch),nbas0(nch),nb0,nbmax,NS,IS,nb,ii,ifout,M1,
     x          ioi,ioj,term
	parameter(imp=60)
	parameter (fracshift=0.25, mfails=5, tr=.true., prba=.false.,
     x             allorder_exch=.true.)
	parameter(mag=1)   ! increase of basis for elastic
	 TM(I) = SECOND() - TIME0

C
	if(final) then
!	open(imp,form='formatted')
!	open(imp+1,form='formatted')
	call openif(imp); call openif(65)
	 endif
	EL = INITL(1)
      if(pralpha) write(imp,'(/''JT ='',f5.1,a1,'' E ='',f12.8/)') 
     x              jtotal,PSIGN(PARITY+2),ECM(EL,1)
      if(pralpha) call openif(imp+1) 
	if(PCON>=3) call openif(62)
	if(prba) call openif(70)
!	write(KO,*) ' PRALPHA = ',pralpha
        gameps = 1d-10 ! MeV spreading of R-matrix poles
	symm = .false.
	if(say.and.pralpha)write(48,*) 'RMATRIX, REVC: set symm=',symm

	 if(INH>0) stop ' R-matrix method requires INH=0'
	prrm = pralpha.and.PCON>0 ;prkm = prrm ; prtm = pralpha
	prmats = prrm
	prbut = pralpha
	prham = pralpha .and. PCON>4
!	write(KO,*) ' PRALPHA,prrm= ',pralpha,prrm
	PI = 4d0*atan(1d0)
	Z = 0d0
	CI = (0d0,1d0)
	NICH = NCH
	IEX = NCH
	M = N
	NLN = min(NLN1,(N-1)/MR+1)
	FAIL = .false.
	pauli = .false.
	r_added = .false.
	rm_vec(:,:) = 0.; rm_energy(:)=0.; rm_set(:)=.false.

	  rmrad1 = (n-1)*HP(1)  ! channel 1
	  conv = -1d0/coef(1)
	 beta(:) = 0.
	 do C=1,NCH
           IT = ITC(PART(C,1),EXCIT(C,1))
           L = LVAL(C)
           IF(ISOCEN.eq.1) L = JTOTAL + 0.1
           IF(ISOCEN.eq.2) L = LVAL(EL)
           CH(1:2,C) = CHL(L+1,IT,1:2)

	  rmrad1 = (n-1)*HP(C)

	if(btype=='E') then
    	  T = sqrt(abs(bndx)*conv)
    	  beta(C) = sign(T,bndx)

	else if(btype=='L') then
          beta(C) = -L/rmrad1

	else if(btype=='B') then
          beta(C) = bndx/rmrad1

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
!                Look for any B overrides
        do ip=1,nvars
          if(srch_kind(ip)==3) then
            if(abs(JTOTAL-srch_jtot(ip))<.1
     x         .and.abs(PARITY-srch_par(ip))<.1) then
              do I=1,nvars
                if(srch_rterm(I)==srch_rterm(ip).and.
     x             srch_kind(I)==4 .and.
     x             abs(srch_B(I)-nul)>1e-3) then   ! has B-override
                
              if (srch_r_ch(ip)==C) then
                beta(C) = srch_B(I)/rmrad1
              else if (PART(C,1) ==srch_r_ic(I) .and.
     x                EXCIT(C,1)==srch_r_ia(I) .and.
     x                lval(C)   ==srch_r_lch(I)) then   ! match L only, not J or S
                beta(C) = srch_B(I)/rmrad1
              endif
              endif ! kind=4 etc
              enddo ! I
            endif   ! J/pi
            endif   ! kind=3
          enddo

	if(pralpha) write(imp,*) ' Ch. ',C,': beta =',beta(C),' from '
     x                             ,btype,real(bndx),real(rmrad1)
	 enddo

	  if(pralpha.and.final) then
!	    write(imp,'(/''  Basis states for J ='',f6.1,a1,'' @'',
!    x      f8.4,'' MeV;'',f8.4,'' fm;'',i6,'' terms, buttle ='',i2/)') 
!    x       jtotal,PSIGN(PARITY+2),ECM(EL,1),rmrad1,nbasi,buttle
	  	rmrad1 = (n-1)*HP(1)  ! channel 1
	  	t = beta(1)*rmrad1
	    write(imp,191) bndx,beta(1),t,weak
191	     format('   bndx,beta(1),bp =',3f10.5,';  weak',
     x           ' if penetrability < ',1p,e10.2)
	    endif

         DO 18 C=1,NCH
	  pauli = pauli .or. BLOCKD(C)
          IT = ITC(PART(C,1),EXCIT(C,1))
         L = LVAL(C)
           IF(ISOCEN.eq.1) L = JTOTAL + 0.1
           IF(ISOCEN.eq.2) L = LVAL(EL)
           LL1(C) = L*(L+1d0)
         DO 17 I=1,2
  17     CH(I,C) = CHL(L+1,IT,I)
	    AC = (n-1)*HP(C) 
	    R = AC * K(PART(C,1),EXCIT(C,1))
	    T = R**(L+1)
	    T4 = R**(-L)
	    do i=1,L
	    T = T / (2*i+1)
	    if(i<L) T4 = T4 * (2*i+1)
	    enddo
	    PMAX(C) = R/(T**2 + T4**2)
	    HARDSPH(C) = atan(T/T4)
         if(.not.FCWFN)  then
	     an2 = abs(CH(1,C))**2
	     an1 = real(CH(1,C))*REAL(CH(2,C))+
     x             AIMAG(CH(1,C))*AIMAG(CH(2,C)) 
	    else
	     an2 = CGMAT(C,C,1)**2+CFMAT(C,C,1)**2			! F^2 + G^2
	     an1 = CGMAT(C,C,1)*CGMAT(C,C,2) + CFMAT(C,C,1)*CFMAT(C,C,2)! FF' + GG'
	    endif
	    PMAX(C) = R/an2
            SHIFT(C) = AC * (an1 /an2 - beta(C))

!         IF(CDETR.GE.2.or.final)then
!         IF(CDETR.GE.2)then
         IF(prmats)then
         if(.not.FCWFN)  then
      	    WRITE(imp,1000) C,L,(CH(I,C),I=1,2),T4,T,R,SHIFT(C)
	    else
      	    WRITE(imp,1000) C,L,
     x		(dcmplx(CGMAT(C,C,i),CFMAT(C,C,i)),I=1,2),T4,T,R,SHIFT(C)
	    endif
	  endif
 1000      FORMAT('0Channel',I3,': match with L=',I4,' to',
     X			1P,2E13.4, 2(',',2E13.4),0p,2f8.4)
  18    CONTINUE
  	if(pralpha.and.final) 
     x            write(imp,19) ECM(EL,1),PMAX(1:min(20,nch))
  19	format(f8.3,': P =',1p,20e9.1)
  	  nonel(:) = .true.
	do JIN=1,MINTL
	  nonel(INITL(JIN)) = .false.
	enddo
	do c=1,nch
	  chweak(c) = PMAX(c)<weak.and.nonel(c)
	  if(pralpha) write(imp,*) c,pmax(c),nonel(c),'weak?',chweak(c)
	enddo

	Rmat(:,:) = 0.
	if(.not.FJSWTCH) then  ! match asymptotic functions to zero here
	bpot(:,:) = 0d0
	emass(:,:) = 1.0
	em = .false.
	 do c=1,nch
	  do NC=1,NCLIST(c,c)
	    IF = NFLIST(c,c,NC)
	    if(PTYPE(2,IF)/=9) then     ! not an effective mass
	      if(PTYPE(7,IF)==0) then   ! no d/dR or L/R terms in bpot !
	   	 t = -CLIST(c,c,NC)/coef(c)  ! remember coef uses asymptotic m*
	    	 if(abs(t)>1d-20) 
     x	         bpot(1:N,c) = bpot(1:N,c) + t*REAL(FORMF(1:N,IF))
	      endif
           else
 		 emass(1:N,c) = emass(1:N,c) - REAL(FORMF(1:N,IF))
 		 em = .true.
 		endif
	  if(c<=5.and.PCON>=3) 
     x 	  write(62,'(''# channel,NC,IF,CLIST'',3i4,2f10.5)') c,NC,IF,
     X				CLIST(c,c,NC)

	  enddo 
	  if(c<=2.and.PCON>=3) then
 	  write(62,'(''# Diagonal potential in channel  '',i4,f10.5,i4)')
     X			c,coef(c),lval(c)
	  endif
	  vmin(c) = 1e10
 	  do i=5,N,5
	    r = (i-1)*HP(c)
	    t = LVAL(c)*(LVAL(c)+1)/r**2
	    vmin(c) = min(vmin(c),(bpot(i,c)+t)*(-coef(c)))
	  if(c<=5.and.PCON>=3) then
 	    write(62,'(f8.3,4f12.5)') r,bpot(i,c)*(-coef(c)),
     X	(bpot(i,c)+t)*(-coef(c)),(bpot(i,c)+t)*(-coef(c))-ECM(c,1),
     x      emass(i,c)
	   endif
	  enddo
	  if(c<=5.and.PCON>=3) then
 	    write(62,*) '&'  ,'  Vmin =',vmin(c)
	    call flush(62)
	  endif
	enddo
	if(em) symm = .false.
	
	nbas = nbasi
        if(nbasi==1) nbas = 0  ! equivalent, when buttle included
	do ip=1,nvars
	  if(srch_kind(ip)==3.and.abs(JTOTAL-srch_jtot(ip))<.1
     x       .and.abs(PARITY-srch_par(ip))<.1 .and. srch_nopot(ip)) then
	  nbas = 0
	    if(pralpha.or.final) 
     x       write(imp,'(/''  Potential disabled for J ='',f6.1,a1)') 
     x             jtotal,PSIGN(PARITY+2)
	  endif
	enddo

	allocate(ybas(n,nch,0:nbas*mag),alpha(nch,nbas*mag))

	  ybas(:,:,:) = 0.
	nb0 = 0
	red=.false.
	do c=1,nch
	  nbasis(c) = nbas
	  nbas0(c) = nb0
	  conv = -1d0/coef(c)
	  kp = ECM(C,1)*conv

c	Effective masses:
	RM(:) = 1d0/EMASS(:,c)      ! 1/M(r)
      call DERIV(RM,DRM,HP(c),n)   ! first derivative of 1/M
      call DERIV(DRM,D2RM,HP(c),n) ! second derivative of 1/M
	vemass(:,c) = - 0.25*drm(:)**2/rm(:) + 0.5*d2rm(:)
	bpot(:,c) = emass(:,c) * (bpot(:,c) - vemass(:,c))

      if(EM.and.PCON>0) then
	 do I=1,N
	   write(148,'(f8.3,6f12.5)') (I-1)*HP(c),EMASS(I,c),RM(I),
     x      DRM(I),D2RM(I),bpot(I,c)/conv,vemass(I,c)/conv
       enddo
	 write(148,*) '&'
	endif

!	  write(65,'(a,i4,8f8.3)') 'c,ecm(c,:) =',c,ECM(C,1:3),coef(c),
!     x        RMASS(PART(c,1))
	call integ(lval(c),ll1(c),n,HP(c),kp,bpot(1,c),
     x		ybas(1,c,0),emass(1,c),logd(c),imp)
	  nodest(c) = 0
	 do i=max(3,CUTVAL(c)),n-1
	  if(ybas(i,c,0)*ybas(i+1,c,0)<0d0) nodest(c) = nodest(c)+1
	 enddo
!	 write(65,*) ' c,l,ll1,nodest =',c,lval(c),ll1(c),nodest(c)
 
   !  Increase basis size if an elastic channel or nearest inelastic   NOT
	  i = nbasis(c)
	    do jin=1,MINTL
	     EL = INITL(JIN)
!	     if(c==INITL(JIN)) nbasis(c) = i*mag
         	PEL = PART(EL,1)
         	EXL = EXCIT(EL,1)
	    enddo
!	    if(PART(c,1)==PEL.and.abs(EXCIT(c,1)-EXL)<=1) 
!     x				nbasis(c)=i*mag

	  nodmin(c) = max(1,nodest(c)-nbasis(c)/2+1)
	  ib = nbasis(c)
	  if(nodmin(c)==1) 
     x		nbasis(c) = min(max(2*nodest(c),nrbmin,0),nbasis(c))
	  red = red .or. nbasis(c)<ib
 	  
	  if(pralpha) write(imp,*) ' Ch.',c,': nodest,min-node,#basis=',
     X		 nodest(c),nodmin(c),nbasis(c)

	  failed = 0
	  nodes = nodmin(c) 
          EIG = 0.
	 do ib=1,nbasis(c)
c   beta = logarithmic derivative of basis functions at Rmax
          if(ib.le.3) then
	      if(ib>1) then
		EIG = alpha(c,ib-1)
	      else if(nodest(c)<3)  then
		 EIG = ecm(c,1)
	      else if(nodest(c)>10) then
		 EIG = ecm(c,1)*(nodes/real(nodest(c)))**2
	      else 
	         EIG = vmin(c)
	      endif
	   else
              EIG = 2*alpha(c,ib-1)-alpha(c,ib-2)
	   endif

	  IFAULT=0
20	      continue
!	      write(60,*) 'EIG, vmin =',EIG,vmin(c)
	      EIG = max(EIG,vmin(c))
          call SEARCH11(nodes,lval(c),n-1,HP(c),conv,EIG,beta(C),
     X      bpot(2,c),ybas(2,c,ib),emass(2,c),PCON,IFAULT,
     x      max(1,CUTVAL(c)-1),imp+5)
	    ybas(1,c,ib) = 0d0
	  if(IFAULT.ne.0) then
	   failed = failed+1
           write(KO,201) ib,c,lval(c),nodes,nodest(c),IFAULT,EIG,beta(c)
201        format(' Basis',i3,' failure: ch.',i3,' (L =',i5,')',
     X        ' @',i3,'/',I3,' nodes.  FLT,EIG,BETA: ',i3,2f10.5)
!	   PCON = max(PCON,4)
	   if(failed>mfails.and.ib>nbas) then
	        write(KO,*) ' Keep only ',ib-1,' basis states for ch ',c
		nbasis(c) = ib-1
	        go to 207
	   endif
	   if(failed>mfails) then
	     write(KO,*) 'R-matrix basis failure: should stop'
	     write(6,*) 'R-matrix basis failure: should stop'
               if(number_calls>5) then
	 	penalty = penalty + fine
	 	write(6,41) number_calls,penalty
41	 	format('  At call ',i5,' penalty =',1p,e10.1)
	       else
	     	stop 'R-matrix basis failure'
	       endif
	     endif
	   nodes = nodes+1
	   if(ib==1) nodmin(c) = nodes
	   go to 20
	  endif
	  alpha(c,ib)=EIG

	  if(pralpha) then
	     t = - coef(1)/rmrad1
 	     if(ib.eq.1) then
		write(imp+1,*)  '&'
 	     	write(imp+1,
     x            '(''# Chan.'',2i4,2f8.2,'':  E,wf,t*wf^2,gam^2'')') 
     X				c,nodest(c),ecm(c,1),jtotal
		endif
     		ggam = 2*pmax(c)*t*ybas(n,c,ib)**2
 	     write(imp+1,'(i4,1p,4e15.5)')  
     X              nodes,alpha(c,ib),ybas(n,c,ib),t*ybas(n,c,ib)**2
     x             ,ggam
	     call flush(imp+1)
	     if(PCON>0.and.alpha(c,ib)<0.0) then
 	     	write(imp+4,'(''# Chan.'',3i4,f8.2,''  c,ib,nodes,e'')') 
     X				c,ib,nodes,alpha(c,ib)
	       do i=1,n,2
	       write(imp+4,'(f8.3,1p,e12.4)') (i-1)*HP(c),ybas(i,c,ib)
	       enddo
		write(imp+4,*)  '&'
	        call flush(imp+4)
	       endif
	  endif
	  nodes = nodes+1
	enddo
207	  nb0 = nb0+nbasis(c)
	enddo
	nbmax = maxval(nbasis(1:nch))
!!         write(0,*) ' H-E matrix size = ',nb0,nbmax
!	if(pralpha.and.PCON>0) write(ko,*) ' Lower cutoffs: ',cutval(1:NCH)

	 if(pralpha.and.nb0>0) write(KO,*) ' H-E matrix size = ',nb0
	 if(say.and.pralpha) then
         write(48,*) ' H-E matrix size = ',nb0
	if(final.and.maxval(nodest(1:nch))>0) 
     x      write(KO,21) (nodest(c),c=1,min(nch,20))
21	format(' R-matrix bases: expected  nodes',20i3)
	if(final.and.nbas>0) write(KO,22) (nodmin(c),c=1,min(nch,20))
22	format('               channel min nodes',20i3)
	if(final.and.nbmax/=nbas.or.red) 
     x          write(KO,23) (nbasis(c),c=1,min(nch,20))
23	format('               channel basis set',20i3)
	endif

	naa = nb0
	nd = naa
	meigs = meigsi
	if(meigs==0) meigs=2
	if(meigs<0) meigs=naa
	if(nb0>0) then

     	allocate(aa(nd,nd),rhs(nd,max(meigs,nch)),ipiv(nd))
     	if(ISNONO.or.vsearch>0) then  
   	   allocate(kk(nd,nd))
	   kkdim = nd
	  else
     	   allocate(kk(1,1))
	   kkdim = 1
	  endif

	call HEMATRIX(nch,alpha,nbasis,ybas,bpot,n,aa,nd,naa,HP,
     X		FORMF,CLIST,NCLIST,NFLIST,NF,lval,coef,EL,ECM(1,1),
     X		nbas0,nbmax,symm,prba,MR,CHNO,NL,NLO,NLN,NLC,MLT,MLM,
     X		CUTVAL,ICUTC,ISNONO,kk,kkdim,vsearch,emass,vemass,
     x		PTYPE,CPSO)

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
		 wn = wn + ybas(i,c2,ib) * ybas(i,c1,ib2)
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
	
		include 'boundstate.f'

		return
	  endif ! eigens

C			L.U decomposition of Hamiltonian-E matrix
	if(pralpha) write(KO,*) ' LU decomposition of matrix size ',naa
		
	  call zgetrf(naa,naa,AA,nd,ipiv,info)

	  rhs(:,:) = 0d0
	  do 31 c=1,nch
	  if(.not.drop(c)) then
	  do 30 ib=1,nbasis(c)
            ki = nbas0(c) + ib
	    rhs(ki,c) = ybas(n,c,ib)
   30	  continue
   	  endif
   31	  continue

	 if(pralpha) then  ! check convergence for Schrodinger eqn for elastic
	  allocate(schro(n))
	  c = 1
	  j = nint(2*JTOTAL)+(PARITY+1)/2
	  schro(:) = 0.
	  do 333 ib=1,nbasis(c)
            ki = nbas0(c) + ib
	   do i=1,n
	   schro(i) = schro(i) + rhs(ki,c)*ybas(i,c,ib)
	  if(mod(ib,4)==1)
     x	   write(500+j,'(f8.3,f12.5)')  HP(c)*i,schro(i) ! +(ib-1)*0.4
	   enddo
	   write(500+j,*) '&'
  333    continue
	 endif

          call zgetrs('N',naa,nch,AA,nd,ipiv,rhs,nd,info)
	   if(info.ne.0) then
	     write(imp,*)' Error return from H-E zgetrs',info
	     stop 'zgetrs'
	   endif ! info/=0
!	if(pralpha.and.CDETR>3) then
!	   write(imp,*) 'rhs solution:'
!	   call WRCMAT(rhs,naa,nch,nd,4,imp)
! 	 endif	
	   rhsprob(:) = 0d0
	  do 40 i=1,nch
	  do 40 j=1,nch
	   tc = 0d0
	   do 35 ib=1,nbasis(i)
            ki = nbas0(i) + ib
     	   tc = tc + ybas(n,i,ib) * rhs(ki,j)
		rhsprob(ib) = rhsprob(ib) + abs(rhs(ki,j))**2
   35	   continue
	   Rmat(i,j) = - tc * coef(j) 
   40	  continue
	  T = sum(rhsprob(1:nbmax))
	if(final.and.nbmax>2) write(KO,42) (rhsprob(ib)/T,ib=1,nbmax)
   42	format(' R-matrix bases: uses  =',10f6.3/
     x        ('                        ',10f6.3))

	if(prrm) then
	   write(imp,*) ' r*R matrix:'
          do 43 i=1,nch
   43      write(imp,45) i,(Rmat(i,j) , j=1,nch)
   45	   format(1x,i3,18f8.3,:,/(4x,18f8.3))
  451      format(1x,i3,18f10.5,:,/(4x,18f10.5))

	   write(imp,*) ' r*R matrix diag. (before Buttle correction):'
   	   write(imp,45) nch,(Rmat(j,j) , j=1,nch)
	   T = Rmat(1,1)/rmrad1
	   ac = -T*(ECM(1,1)-alpha(1,1))! gam^2
	   esh = -ac*SHIFT(1)		! energy shift
	   an2 = alpha(1,1) + esh	! resonance energy
	   write(178,'(f8.3,5f10.5)') ECM(1,1),T,ac,SHIFT(1),an2
	  endif ! prrm
C
	endif   ! nb0>0   (else Rmat=0 + buttle)

	if(buttle<4.and.nbas>0) then
C********************************************** `Exact' Buttle Correction
	  allocate(phic(nch),rhsb(nd))
	  if(mod(buttle,2)==0) allocate(yc(n))
	do 150 i=1,nch
	  nab = nbasis(i)
c
c       Calculate exact wave function for e.g. diagonal potential
C	at an energy EB not too close to single-channel R-matrix energy
	EB = ECM(i,1)
!	if(EB.lt.0.0) go to 150
	changed = .false.
	do 120  ib=1,nab-1
!	if(abs(alpha(i,ib)-EB).lt.abs(alpha(i,ib+1)-alpha(i,ib))*0.1) 
!    x 			then
!       EB = alpha(i,ib) + sign(0.1d0,EB-alpha(i,ib))
!    x                         *(alpha(i,ib+1)-alpha(i,ib))
!	changed = .true.
!	if(prbut) write(imp,*) ' E =',real(ECM(i,1)),' too close to ',
!    x		real(alpha(i,ib)),ib,' so move to ',real(EB) 
!       endif
         t = abs(alpha(i,ib+1)-alpha(i,ib))
         if(ib>1) t = min(t,abs(alpha(i,ib)-alpha(i,ib-1)))
	 t = min(t * fracshift,1.00d0)
	  if(nch==1) t = t*0.001d0  ! avoid coincident alpha and R-matrix poles?
	  if(buttle>=2) t=0.
         if(abs(alpha(i,ib)-EB).lt.t) then
          EB = alpha(i,ib) + sign(t,EB-alpha(i,ib))
          changed = .true.
          if(prbut) write(imp,*) ' E =',real(ECM(i,1)),' too close to ',
     x          real(alpha(i,ib)),ib,' so move to ',real(EB)
         endif
  120	continue
 
	if(mod(buttle,2)==1) then
	if(changed) call integ(lval(i),ll1(i),n,HP(i),-EB/coef(i),
     x		bpot(1,i),ybas(1,i,0),emass(1,i),logd(i),imp)
        Rex0 = 1d0/(real(logd(i)) - beta(i))
        Rmat0 = 0.0
          do 145 ib=1,nab
  145    Rmat0 = Rmat0 + ybas(n,i,ib)**2/(alpha(i,ib)-EB)*(-coef(i))
  
	else !   COMPLEX BUTTLE
	  vpot(1:N,1) = 0d0
	  do NC=1,NCLIST(i,i)
	    IF = NFLIST(i,i,NC) 
	   if(PTYPE(7,IF)==0) then
	    t = -CLIST(i,i,NC)/coef(i)
	    if(abs(t)>1d-20) vpot(1:N,1) = vpot(1:N,1) + t*FORMF(1:N,IF)
	   endif
	  enddo 
	  call integc(lval(i),ll1(i),n,HP(i),-EB/coef(i),
     x		vpot,yc,emass(1,i),logd(i),imp)
          Rex0 = 1d0/(logd(i) - beta(i))
	if(prbut) write(imp,144) i,logd(i),beta(i),Rex0
144	format('In #',i4,' logd,beta =',3f10.5,' so Rex0 =',2f10.5)
!        NOW calculate one-channel R-matrix from basis, with imag optical potl

	 if(nab>0) then
C		L.U decomposition of one-channel Hamiltonian-E matrix
	  do 146 j=1,nab
	  do 146 ki=1,nab
	  g=4.d0/3.d0
	  tc = 0d0
	  do  ii=2,n
	    if(ii.eq.n) g = 1d0/3d0
	    g=2d0-g
	    tc = tc + ybas(ii,i,j) * ybas(ii,i,ki) * g
     x          * (emass(ii,i)*(vpot(ii,1)-vemass(ii,i))-bpot(ii,i))
	   enddo
	  AA(j,ki) = tc * HP(i) * (-coef(i))
	  if(j==ki) AA(j,ki) = AA(j,ki) + alpha(i,j) - EB
  146 	  continue
!	if(prham) then
!	   write(imp,*) 'H-E matrix for channel ',i,' ::'
!	   call WRCMAT(AA,nab,nab,nd,4,imp)
! 	 endif
 	
 	 call zgetrf(nab,nab,AA,nd,ipiv,info)

	  rhsb(1:nab) = 0d0
	  do ib=1,nab
            ki = ib
	    rhsb(ki) = ybas(n,i,ib)
     	  enddo
          call zgetrs('N',nab,1,AA,nd,ipiv,rhsb,nd,info)
	   if(info.ne.0) then
	     write(imp,*)' Error return from Buttle zgetrs',info
	     stop 'zgetrs'
	   endif
	   tc = 0d0
	   do ib=1,nab
            ki = ib
     	   tc = tc + ybas(n,i,ib) * rhsb(ki)
     	   enddo
	   Rmat0 = - tc * coef(i)
	   else
	   Rmat0 = 0. 		! nab=0
	   endif

	endif
C
        phic(i)=Rex0

C 	Correct diagonal elements of full R-matrix:

	Rmati  = Rmat(i,i)
	Rmat(i,i) = Rmat(i,i) + Rex0 - Rmat0

!	if(prbut) 
!    X     write(120+i+50*(1-PARITY),'(f8.3,1p,2e12.4)') EB,Rex0-Rmat0
	if(prbut) write(imp,148) ECM(i,1),EB,i,Rmati,Rex0,Rmat0,
     X			Rmat(i,i),buttle,changed,fracshift
  148   format(' E,EB,ch =',2f8.4,i3,' R,buttle=',f8.4,f7.4,'+',
     X		f8.4,f7.4,'-',f8.4,f7.4,' =>',2f10.5,I3,L2,f5.2)
  150   continue
	endif
	if(prrm.and.buttle<4.and.nbas>0) then
!	   write(imp,*) ' Exact uncoupled r*R matrix diagonal:'
!   	   write(imp,45) nch,(phic(i) , i=1,nch)
	   write(imp,*) ' r*R matrix diagonal (with Buttle correction):'
   	   write(imp,45) nch,(Rmat(i,i) , i=1,nch)
	  endif
	  if(buttle<4.and.nbas>0) then
	     deallocate(phic,rhsb)
	     if(mod(buttle,2)==0) deallocate(yc)
	     endif

	endif ! .not.FJSWTCH

	   do c1=1,nch
	    ic = part(c1,1)
           RMA(c1)   = (n-1)*HP(C1)
	   ETACO(c1) = ETACNS * MASS(2+1,IC)*MASS(2+2,IC) *
     x                  SQRT(rmass(ic)) ! for Barker
           ENFF(c1)  = ECM(c1,1) - ETOTAL ! for Barker
	   enddo
           Pshift(:,:) = 0.0

!               Add in any search R-matrix term!
        do ip=1,nvars
          if(srch_kind(ip)==3) then
            if(abs(JTOTAL-srch_jtot(ip))<.1
     x         .and.abs(PARITY-srch_par(ip))<.1) then
              term = srch_rterm(ip)
              rm_energy(term) = srch_value(ip)  ! real part only, so far
              rm_set(term) = .true.
            endif  ! J/pi
          endif ! kind=3
         enddo !ip
        if(pralpha) write(imp,*) 'rm_set:',rm_set

        do ip=1,nvars
          ich = srch_r_ch(ip)
          term = srch_rterm(ip)
          if(srch_kind(ip)==4.and.rm_set(term).and.
     x            ich<=nch) then
          if(pralpha) write(imp,*) ip,term,' used'
             rmdist(:) = 0.0
             if(ich>=1) then ! have channel 'ich'
                rmdist(ich) = 1.0

             else ! search for channel with ic,ia,l,j values
             JV = srch_r_jch(ip)
             SV = srch_r_sch(ip)
               do c=1,nch
                if(PART(c,1) ==srch_r_ic(ip) .and.
     x             EXCIT(c,1)==srch_r_ia(ip) .and.
     x             lval(c)   ==srch_r_lch(ip)) then
                 if     (JV>-0.1) then   ! J specified
                   if(abs(jval(c)-JV)<0.1) rmdist(c)=1.0
                 else if(SV>-0.1) then   ! S specified
                    rmdist(c) = sqrt((2*SV+1.) * (2*jval(c)+1.)) *
     x            racah(lval(c)+Z,JPROJ(c),JTOTAL,JTARG(c),jval(c),SV)
      if(pralpha) write(imp,14997) lval(c)+Z,JPROJ(c),JTOTAL,JTARG(c),
     x            jval(c),SV,rmdist(c)
14997    format(' HH*RACAH(',4f5.1,';',2f5.1,') =',f8.4)
                 else
                    write(6,14998) ip,srch_r_jch(ip),srch_r_sch(ip)
14998               format(' Neither JCH nor SCH given for variable',
     x                       i3,2f5.1)
                    if(abs(srch_value(ip))>1d-30) stop
                 endif  ! JV or SV given
                C1 = c
                endif   ! part,excit, L agree
               enddo    ! c
               if(abs(sum(rmdist(:)**2)-1.) > 1e-5) then
                write(6,15001) ip,term,srch_name(ip),
     x     srch_r_ic(ip),
     x     srch_r_ia(ip),srch_r_lch(ip),srch_r_jch(ip),srch_r_sch(ip)
15001   format(' Rmatrix channel not found #',i3,' term ',i3,':',
     x           1x,a10,' ic,ia,l,J,S =',3i3,2f5.1)
         write(6,*) ip,JTOTAL,PARITY,rmdist(1:nch)

                if(abs(srch_value(ip))>1d-30) stop
               endif  ! rmdist non-unit
             endif! rmdist non trivial

             rwamp = srch_value(ip)
             if(.not.srch_rwa(ip).or.eobs) then
                ic = part(c1,1)
                E1 = rm_energy(term) + ECM(c1,1) - ETOTAL   ! using real part
                kp = -E1/coef(C1)
                wn = sqrt(abs(kp))
                S = ((n-1)*HP(C1)) * wn
                eta1 = ETACNS * MASS(2+1,IC)*MASS(2+2,IC) *
     x                  SQRT(rmass(ic)/abs(E1))
                L = srch_r_lch(ip)
                T = L
                CALL COULFG(S,eta1,0D0,T,fc_s,gc_s,fcp_s,gcp_s,1,0,J,M1)
                 T = fc_s(L+1)**2 + gc_s(L+1)**2
                pene = S / T
                Pshift(c1,term) = pene * (fc_s(L+1)*fcp_s(L+1)
     x                                 + gc_s(L+1)*gcp_s(L+1))
               do c=1,nch
                if(abs(rmdist(c))>1e-10)
     x                 Pshift(c,term) = srch_r_shift(ip)
               enddo
!           write(161,*) c1,E1,S,eta1,L,Pshift(c1,term)
               if(.not.srch_rwa(ip)) then
                rwamp = sqrt(abs(srch_value(ip)/(2*pene)))
                if(srch_value(ip)<0.0) rwamp = -rwamp
!       write(192,1501) ip,JTOTAL,srch_jtot(ip),PARITY,srch_par(ip),
!     x         term,4,srch_value(ip)
!          write(192,15004) ip,ic,e1,wn,s,eta1,l,pene,rwamp
!15004   format(2i4,4f10.5,i4,2f10.6/)
                endif
                endif

             do i=1,nch
              rm_vec(i,term) = rm_vec(i,term) + rwamp*rmdist(i)
             enddo
              if(abs(srch_value(ip))>0.) r_added=.true.
            if(pralpha) then
               wid='width'
               if(srch_rwa(ip)) wid='r.w.a'
               if(ich>0) then
             write(imp,1511) ki,wid,srch_value(ip),ich
 1511   format(' New R term',i4,' ',a5,f8.4,' in c=',i4)
               else
             write(imp,1512) ki,wid,srch_value(ip),rmdist(1:min(nch,10))
 1512   format(' New R term',i4,' ',a5,f8.4,' in chs=',10f6.3)
               endif  ! ich>0
        write(imp,1501) ip,JTOTAL,srch_jtot(ip),PARITY,srch_par(ip),
     x         srch_rterm(ip),4,srch_value(ip)
1501    format(' Var',i4,' for ',2f5.1,2i3,': term #',2i3,' val',f10.5)
            endif ! pralpha

          endif ! kind=4
         enddo !ip

!		Add in any search R-matrix term!
	do ip=1,nvars
	  if(srch_kind(ip)==3) then
	    if(abs(JTOTAL-srch_jtot(ip))<.1
     x         .and.abs(PARITY-srch_par(ip))<.1) then
             damp = 0.0
             term = srch_rterm(ip)
             do ip2=1,nvars
	      if(srch_kind(ip2)==7.and.srch_rterm(ip2)==term)
     x        then
                c1 = EL
                if(srch_Brune(ip2)) then  

	          call BarkerTransform(term,rm_set,mterms,rm_energy,
     x               rm_vec,nch,btype,beta,jtotal,PSIGN(PARITY+2),
     x               LVAL,ETACO,ENFF,coef,RMA,pralpha,final,
     x               E_Brune2,W_Brune2)
                  if(pralpha) write(imp,*) 'ip,term,ip2,brune',
     x              ip,term,ip2,srch_Brune(ip2)," E_B:",E_Brune2(term)

                    E1 = E_Brune2(term)  + ECM(c1,1) - ETOTAL  
                  else
                    E1 = srch_value(ip) + ECM(c1,1) - ETOTAL  
                  endif
                ET(1) = E1        - srch_damp(ip2)  ! pole above threshold
                ET(2) = ECM(c1,1) - srch_damp(ip2)  ! scattering above threshold
                penet(:) = 0.0
                if(pralpha .and. abs(srch_value(ip2))>1e-5) then
                  write(imp,*)
                  write(imp,149) ip, srch_value(ip),
     x                 real(srch_rterm(ip))
                  write(imp,149) ip2,srch_value(ip2),
     x                 real(srch_leff(ip2)),ET,ENLAB
                 endif
                if (srch_damp(ip2) < -1.) then
                  penefac = 1.0
                else if (ET(2) < 0.01 .or. ET(1)<0.01) then
                  penefac = 0.0
                else  ! calculate penefac = P(ET2)/P(ET1) in elastic channel for L = leff
                   ic = part(c1,1)
                   do I=1,2
                   kp = -ET(I)/coef(C1)
                   wn = sqrt(abs(kp))
                   S = ((n-1)*HP(C1)) * wn
                   eta1 = ETACNS * MASS(2+1,IC)*MASS(2+2,IC) *
     x                  SQRT(rmass(ic)/abs(ET(I)))
                   L = int(abs(srch_leff(ip2)+0.5)-0.5)
                   AN = abs(srch_leff(ip2)+0.5)-0.5 - L
                   T = srch_leff(ip2)
                  !write(imp,*) C1,ET(I),kp,wn,(n-1)*HP(C1),S,eta1,an,L,T
                 CALL COULFG(S,eta1,AN,T,fc_s,gc_s,fcp_s,gcp_s,1,0,J,M1)
                    T = fc_s(L+1)**2 + gc_s(L+1)**2
                   penet(I) = S / T
                   enddo
   		  penefac = penet(2)/penet(1)
		endif
                damp = damp + srch_value(ip2) * penefac
                if(pralpha.and.abs(srch_value(ip2))>1e-5) 
     x          write(imp,149) jtotal,parity,ip2,srch_value(ip2),
     x                         penefac,penet,ET
149	        format(f5.1,i2,' Damping from param',i3,' is',f10.5,
     x         ' *',f12.4,', from P=',2f12.3,', from E=',2f12.3)
               endif
	     enddo
              rm_energy(srch_rterm(ip)) = 
     x           cmplx(srch_value(ip),-(damp+srch_damp(ip))*0.5-gameps)

	      if(pralpha) then
                ki = srch_kind(ip)
                ! write(imp,151) ki,srch_value(ip)
        	write(imp,1501) ip,JTOTAL,srch_jtot(ip),PARITY,srch_par(ip),
     x           srch_rterm(ip),3,srch_value(ip)
               write(imp,*) 'Term ',srch_rterm(ip),' set:',JTOTAL,PARITY,ip
               endif
	      endif
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
         if(firstE.and.rm_Brune(ki)) then  ! output to E_Brune(ki),W_Brune(ki)
	  call BarkerTransform(ki,rm_set,mterms,rm_energy,rm_vec,nch,
     x       btype,beta,jtotal,PSIGN(PARITY+2),LVAL,ETACO,ENFF,coef,RMA,
     x       pralpha,final,E_Brune,W_Brune)
	    if(pralpha) write(6,*) ki,'Brune: E:',E_Brune
          endif
	 do I=1,NCH
          SYMR(I) = sqrt(-COEF(I))
         enddo
          EobsShift = 0.0
          if(eobs) then
   	  do i=1,nch
	  AC = (n-1)*HP(i) 
	  EobsShift = EobsShift + rm_vec(i,ki)*rm_vec(i,ki)
     x           * (Pshift(i,ki) - AC*beta(i))
          write(160,1517) i,rm_vec(i,ki),Pshift(i,ki),
     x      (Pshift(i,ki) - AC*beta(i)),
     x           rm_vec(i,ki)**2* (Pshift(i,ki) - AC*beta(i))
          enddo
          write(160,1517) ki,EobsShift
          endif
1517	format(i4,4f12.5)

	do 152 i=1,nch
	do 152 j=1,nch
	 if(.not.chweak(j)) then
!	  T = sqrt(RMASS(PART(i,1))/RMASS(PART(j,1)))
	  T = sqrt(HP(j)*COEF(j)/(HP(i)*COEF(i)))

	  Rmat(i,j) = Rmat(i,j) +  T* (n-1)*HP(i)*
     x          rm_vec(i,ki)*rm_vec(j,ki)/
     x         (rm_energy(ki) + EobsShift -ETOTAL)	

	if(pralpha) ! .or.rm_vec(i,ki)*rm_vec(j,ki)>0.)
     x  write(imp,151) ki,rm_energy(ki),i,j, rm_vec(i,ki),rm_vec(j,ki),
     x   ECM(INITL(1),1),Rmat(i,j),T,
     x       Rmat(i,j)   *SYMR(i)/SYMR(j)
 151	format(' New R term',i4,' at',2f8.4,:,':',2i3,'  with',2f10.5,
     x ' at ',f9.4,' MeV => R =',2f10.4,' C/C =',f8.3,2f10.4)
!       if(i>j.and.abs(rmat(i,j))>1e-20) rtrace = Rmat(i,j)
	endif

 152	continue
 	endif
 	enddo


C****************************************************** Find open channels

	call nopen(NCH,ECM(1,1),iop,nop)
	if(pralpha.and.nop<nch) then
  	   write(KO,*) 'Only ',nop,' open channels'
	   write(imp,*) 'Only ',nop,' open channels'
	endif

C****************************************************** Print final R-matrix
	if(pralpha.and.r_added) then
	   write(imp,*) ' r*R matrix with additional search terms:'
          do 153 i=1,nch
          SYMR(I) = sqrt(-COEF(I))
  153      write(imp,45) i,(Rmat(i,j) , j=1,nch)
	   write(imp,*) ' r*R matrix symmetrized'
          do i=1,nch
           write(imp,45) i,(SYMR(i)*Rmat(i,j)/SYMR(j),j=1,nch)
          enddo
! 	write(155,*) ECM(INITL(1),1),abs(rtrace)
	  endif ! prrm
	
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
!  235	   format(1x,i3,18f8.4,:,/(4x,18f8.4))
!	  if(prrm) then
!	  T = (atan(real(Kmat(1,1)))-HARDSPH(1))*180/pi
!	  if(T<0) T = T + 180.
!	  T4 = (atan(real(Rmat(1,1)/rmrad1)*pmax(1))-HARDSPH(1))*180/pi
!	  if(T4<0) T4 = T4 + 180.
!	  an1 = 1 - Rmat(1,1)/rmrad1 * SHIFT(1)
!	 ac=(atan(real(Rmat(1,1)/rmrad1)*pmax(1)/an1)-HARDSPH(1))*180/pi
!	  if(ac<0) ac = ac + 180.
!        an2=(atan(ggam*0.5/(alpha(1,1)+esh-ECM(1,1)))-HARDSPH(1))*180/pi
!	  if(an2<0) an2 = an2 + 180
!	  T4 = (atan(ggam*0.5/(0.44-ECM(1,1)))-HARDSPH(1))*180/pi
!	  if(T4<0) T4 = T4 + 180.
!	  rho2 = (atan(0.54*0.5/(0.44-ECM(1,1)))-HARDSPH(1))*180/pi
!	  if(rho2<0) rho2 = rho2 + 180.
!	  write(177,'(f8.3,6f10.4)')  ECM(1,1),T,ac,an2,T4,rho2
!	  rmst = alpha(1,1)+esh
!	  write(179,'(f8.3,6f10.4)')  ECM(1,1),ggam,esh
!	  endif
	 endif
C******************************************* S Scattering Matrix
C
C    S = [1 + i K] * [1 - i K]^{-1}
C
	call GETS(Kmat,Smatr,prtm,iop,nop,nch,lval,imp)

        if(prtm) then
	DO IOI=1,NOP
        T = 0.0
         do IOJ=1,NOP
           T = T + abs(Smatr(IOP(IOI),IOP(IOJ)))**2
         enddo
           unitar(IOI) = T
	enddo
        write(IMP,'('' Unitarity ='',4x,9f16.5)') unitar(1:NOP)

	 do IOI=1,NOP
          I = IOP(IOI)
          SYMS(I) = (2. * ecm(I,1)/RMASS(PART(I,1)))**0.25d0
         enddo
	write(IMP,*) ' S matrix symmetrized, but no i^L factors'
	DO IOI=1,NOP
          I = IOP(IOI)
        T = 0.0
         do IOJ=1,NOP
           Ssym(IOI,IOJ) = SYMS(I)*Smatr(I,IOP(IOJ))/SYMS(IOP(IOJ))
     x              / (0.,1.)**(LVAL(IOJ)-LVAL(IOI))
           T = T + abs(Ssym(IOI,IOJ))**2
        enddo
          write(IMP,451) IOI,(Ssym(IOI,IOJ),IOJ=1,NOP)
           unitar(IOI) = T
   57   format(1x,I3,18f8.4,:,/(4x,18f8.4))
        ENDDO
        write(IMP,'('' Unitarity ='',4x,9f16.5)') unitar(1:NOP)
        endif
C******************************************** S Columns for Cross sections
C
      	FUSL(:) = 0.0
	if(WOUT) allocate(psi(N,NCH))
	if(.not.WOUT) allocate(psi(1,1))

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
	  TC  = rhs(ki,j) * WFP(j) * (-coef(j))
!	  if(PRALPHA) write(imp,252) c,j,ib,TC
252	  format('TC for ',3i4,'=',2f10.5)
	do i=1,n
     	   psi(i,c) = psi(i,c) + ybas(i,c,ib) * TC
	enddo
	  rmmat(c,j) = rmmat(c,j) -ybas(n,c,ib)*rhs(ki,j)*coef(j)
        enddo
  260	continue
        if(PRALPHA) write(imp,*) ' R-matrix again:'
          do 261 i=1,nch
  261      if(PRALPHA)write(imp,45) i,(rmmat(i,j) , j=1,nch)
	WFN(:) = 0.0
	do c=1,nch
	WFN(:) = WFN(:) + rmat(:,c)*WFP(c)
!        if(PRALPHA) write(imp,2615) WFP(c),WFN(:)
!2615	format(' Contribution from ',2f9.5,': ',20f10.5)
	enddo
!        if(PRALPHA) write(imp,*)
        if(PRALPHA) write(imp,262) (WFP(c),WFN(c),psi(n,c),c=1,NCH)
262     format(' WFP =',2f10.5,' WF =',2f10.5,' psi='2f10.5)
       endif

         DO C=1,NCH
         SMAT(C) = Smatr(C,EL)  ! With GETS (0.,1.)**(LVAL(EL)-LVAL(C))
         if(WOUT) psi(:,C) = psi(:,C) * (0.,1.)**(LVAL(EL)-LVAL(C))
         ENDDO


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
!1430 FORMAT(5(1X,F11.5,' +i*',F9.5,',') )
 1430 FORMAT(5(1X,F11.8,' +i*',F11.8,',') )
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

	if(allocated(aa)) deallocate(aa,rhs,ipiv)
	if(allocated(kk)) deallocate(kk)
	if(allocated(ybas)) deallocate(ybas,alpha)
	return
	END SUBROUTINE RMATRIX

      SUBROUTINE INTEG(LT,LL1,NMAX,H,KP,VPOT,Y,EM,LOGD,imp)
      implicit real*8(a-h,o-z)
      REAL*8 KP,LL1
      COMPLEX*16 LOGD
      DIMENSION Y(NMAX),VPOT(NMAX),EM(NMAX)
      HP=H*H
      Y(1)=0.
      Y(2)=max(H**LT,100*tiny(1d0))
      G2=0D0
      IF(LT.EQ.1)G2=-2*Y(2)
	W1 = VPOT(2)
      G3=-LL1-(W1-KP*EM(2))*HP+12D0

      DO 30 N=3,NMAX
      G1=G2
      G2=Y(N-1)*G3
	WN = VPOT(N)
      G3=-LL1/real(N-1)**2-(WN-KP*EM(N))*HP
      G3=G3+12D0
      Y(N)=((16D0*Y(N-1)-G2)*9D0-G2-G1)/G3
   30 CONTINUE

      LOGD=147.0*Y(NMAX)-360.0*Y(NMAX-1)+450.0*Y(NMAX-2)
     1 -400.0*Y(NMAX-3)+225.0*Y(NMAX-4)-72.*Y(NMAX-5)+10.*Y(NMAX-6)
      LOGD=LOGD/(60.0*Y(NMAX)*H)

      RETURN
      END
      SUBROUTINE INTEGC(LT,LL1,NMAX,H,KP,VPOT,Y,EM,LOGD,imp)
      implicit real*8(a-h,o-z)
      REAL*8 KP,LL1,EM(NMAX)
      COMPLEX*16 LOGD,Y(NMAX),VPOT(NMAX),W1,WN,G3,G2,G1
      HP=H*H
      Y(1)=0.
      Y(2)=max(H**LT,100*tiny(1d0))
      G2=0D0
      IF(LT.EQ.1)G2=-2*Y(2)
	W1 = VPOT(2)
      G3=-LL1-(W1-KP*EM(2))*HP+12D0

      DO 30 N=3,NMAX
      G1=G2
      G2=Y(N-1)*G3
	WN = VPOT(N)
      G3=-LL1/real(N-1)**2-(WN-KP*EM(N))*HP
      G3=G3+12D0
      Y(N)=((16D0*Y(N-1)-G2)*9D0-G2-G1)/G3
   30 CONTINUE

      LOGD=147.0*Y(NMAX)-360.0*Y(NMAX-1)+450.0*Y(NMAX-2)
     1 -400.0*Y(NMAX-3)+225.0*Y(NMAX-4)-72.*Y(NMAX-5)+10.*Y(NMAX-6)
      LOGD=LOGD/(60.0*Y(NMAX)*H)
!	write(imp,10) Y(1:NMAX)
!10	format(10(2g10.2,','))

      RETURN
      END
      SUBROUTINE SEARCH11(NODES,LT,NMAX,H,CONV,EIG,BETA,
     X     W,Y,EM,PCON,IFAULT,N0,IMP)
      implicit real*8(a-h,o-z)
      PARAMETER(ERROR=1d-8)
      REAL*8 INTT,INTP,KP,LOGD,WFNODES(2,1000+NODES),
     X    G3NODES(1000+NODES),INODES(1000+NODES)
      LOGICAL INNERTP,ISOUT,INVERTED
      INTEGER BC,COUNT,CC,CCP,PCON
      DIMENSION Y(NMAX),W(NMAX),EM(NMAX)
!      NMAXP=NMAX+1
      NMAXP=NMAX
	if(mod(NMAX,2)/=0) write(6,*) ' SEARCH11: ERROR - UNEVEN ',NMAX
!     write(65,*)  'SEARCH11',NODES,LT,NMAX,H,CONV,EIG,BETA,
!    X     PCON,IFAULT,N0,IMP
!     write(65,*)  'W',W
      IFAULT = 0
      small = tiny(H)
      ssmall = small**0.3
      HP=H*H
      LM=(LT+1)*LT
      R0 = N0*H
      Y(:) = 0d0
      Y(N0)=min(ssmall,R0**LT)
      NT=-1
      COUNT=1
      RATIO=1.
      NMISS=0
      INNERTP = .true.
      ISOUT = .false.
      INVERTED = .false.; Q = 1.
        G2 = 0.
      CCP=1
      BC=0
      NLAST = 0
      LNODE = 0
      NUSUCC = 0
      PPI = 1.0
      PPP=EIG
      P=PPP
      PP=PPP
      IF(PCON.GE.3)WRITE(IMP,207) NODES,EIG,LT
      IF(PCON.GE.5)WRITE(63,2071) NODES,EIG,LT
  207 FORMAT(/'      PARAMETERS      MISMATCH    NODES  =>',i4,f10.3,i4)
 2071 FORMAT('# Look  for ',i4,' nodes, starting at ',f10.3,' L=',i4)
  102 KP=-EIG*CONV
C                 Use logarithmic derivative BETA at Rmax   G3=E-V
	Y(NMAXP)= 1d-50
	if(EIG<0.) Y(NMAXP)= ssmall
!	Y(NMAXP)= 1.
	Y(NMAXP-1)=Y(NMAXP) *(1d0 - BETA * H)
      G3=-real(LM)/(NMAX*NMAX*HP)-(W(NMAX)+KP*EM(NMAX))
C			use local energy approximation
	IF(G3.lt.-0.001) then
	    X = sqrt(-G3)
	    T = (X-BETA)/(2d0*X)
	    TP= 1d0 - T	    
	    Y(NMAXP-1)=Y(NMAXP) *(T*exp(X*H) + TP*exp(-X*H)
     X                             -H**3/6d0*real(LM)*2d0/(NMAX*H)**2)
	    T=(Y(NMAXP)-Y(NMAXP-1))/(H*Y(NMAXP))
	else IF(G3.gt.0.001) then
	    X = sqrt(G3)
	    T = BETA/X
	    TP= 1d0 	    
	    Y(NMAXP-1)=Y(NMAXP) *(T*sin(-X*H) + TP*cos(-X*H))
c     X                             -H**3/6d0*real(LM)*2d0/(NMAX*H)**2)
	    T=(Y(NMAXP)-Y(NMAXP-1))/(H*Y(NMAXP))
	endif
	if(COUNT<1) then
	do N=1,NMAXP
	write(150,'(f8.3,4f12.5)') N*H,(W(N)+KP*EM(N)-KP)/CONV!,EM(N)
	enddo
	write(150,'(a,f12.5)') '&',EIG
	endif
      DEL=0.
      NS=NMAXP-1
      INTT=0.
      INTP=0.
      Y(N0)=Y(N0)*RATIO
      RATIO=1.
      NCO=1
C		Find matching point: first allowed (E>V) as integrating in
C		NOT(?) first attractive turning point as integrate in.
C		If no such point found, use NMAX*2/3 (for positive EIG),
C		 and NMAX/10 for negative EIG.

	if(INVERTED.and.NUSUCC>0) then
	   NU = NUSUCC
	   go to 302 ! use previous successful NU !!
	   endif
      NU=0
	G1 = 1.
      T = 0.0
       I = NS-5
      DO 30 N=I,N0+5,-1
      G3=-real(LM)/(N*N*HP)-(W(N)+KP*EM(N))  ! G3 = kinetic energy
       NU = N
        IF(INNERTP) then
	 if(G3>T.and.G1<T) go to 302  ! first allowed r
	else
         if(G3<T.and.G1>T) go to 302  ! first forbidden r
         if(N<I.and.G3>G1) go to 302  ! over peak of Coulomb barrier
	endif
   30   G1 = G3
! 301   continue
!       if(EIG.ge.0) NU = NMAX*2/3
       if(EIG.ge.0) NU = NMAX-10
       if(EIG.lt.0) NU = NMAX/10
       NU = max(NU,N0+5)
  302 T=0.
      TP=0.
!        write(62,*) '# At E =',real(EIG),' NU =',NU,real(H)*NU
      G3=-real(LM)/(NMAXP*NMAXP)-(W(NMAXP)+KP*EM(NMAXP))*HP+12.
        ISOUT = G3>0.
      N=NS
   31 G1=G2
      G2=Y(N-NT)*G3
      G3=-real(LM)/real(N)**2-(W(N)+KP*EM(N))*HP
c     IF(NU.eq.0) THEN
c       IF ((G3.GE.0.OR. N.LT.NUMIN).and.N.LT.NMAX/2)NU=N
c      ENDIF
      G3=G3+12.
      IF(N.EQ.1.AND.LT.EQ.1)G2=-2*Y(1)
      I=N-NT
      IF(N.NE.NS) Y(N)=((16.*Y(I)-G2)*9.0-G2-G1)/G3
      X=Y(N)**2
       DV=-CONV
      T=T+X*DV
      TP=X+TP
	if(PCON>4) write(63,315) N*H,Y(N),W(N),(G3-12.)/HP/CONV,EM(N)
	if(PCON>4) call flush(63)
315 	format(f8.2,1p,e12.3,0p,3f10.4)
	ETHR=0.5
	ETHR=-1e6
      IF(Y(N)*Y(I)<.0.AND.N>N0.and.G3-12.>ETHR*HP*CONV)then
	if(PCON.ge.3) then
	INODES(NCO)  = N*H
	G3NODES(NCO)  = (G3-12.)/HP/CONV
	WFNODES(1,NCO) = Y(N)
	WFNODES(2,NCO) = Y(I)
	endif
	NCO=NCO+1
	endif
      IF(N.EQ.NU)GO TO 41
      N=N+NT
      GO TO 31
   41 TP=TP-0.5*X
	if(PCON>4) write(63,*) '&'
      ICC=NU-NT
      N=ICC-NT
      NS=N-NT
      LTT=NS-NT
      I=LTT-NT
      IFA=I-NT
	if(IFA<1) then
	  write(6,*) 'Error in SEARCH11 as NU= ',NU,' too small'
	  write(6,*) 'LT,LM,H,KP,EIG,CONV =',LT,LM,H,KP,EIG,CONV
       	   DO N=1,NMAX
      	G3=-real(LM)/(N*N*HP)-(W(N)+KP*EM(N))  ! G3 = kinetic energy
	    write(6,42) N*H,W(N),EM(N),G3
42		format(f8.3,3f12.5)
	   enddo
	  stop
	endif
      LOGD=147.0*Y(NU)-360.0*Y(ICC)+450.0*Y(N)
     1 -400.0*Y(NS)+225.0*Y(LTT)-72.*Y(I)+10.*Y(IFA)
      LOGD=NT*LOGD/(60.0*Y(NU)*H)
      NT=-NT
      NS=1
      DEL=-LOGD-DEL
      RATIO=1.0/(RATIO*Y(NU))
      INTT=T/X+INTT
      INTP=TP/X+INTP
      IF(NT.LE.0) GO TO 3
      X=Y(N0)**2
       DV=-CONV
      T=X*DV
      TP=X
      N=N0+1
      G2=0
      IF(LT.EQ.1.and.N0==1)G2=-2*Y(1)
      G3=-real(LM)/real(N0*N0)-(W(N0)+KP*EM(N0))*HP+12.
      GO TO 31
    3 COUNT=COUNT+1
      IF(ABS(INTT).gt.1d300) go to 200
      INTP=INTP*X*RATIO*RATIO

      if(NCO==LNODE.and.NCO/=NODES) NMISS = NMISS+1
      if(NCO==NODES) NUSUCC = NU

      INVERTED = (NCO-LNODE)*Q < 0 
     x    .and. NODES==LNODE .and. COUNT>2.and.BC==1

      IFAULT=0
      IF(COUNT.GT.80)GO TO 200
      IF(PCON. GE.3) then
        WRITE(IMP,2170) COUNT,P,DEL,NMISS,NCO,NU*H
     X			,INNERTP,INVERTED !! ,DEL/(INTT*H)
        WRITE(63,2171)P,DEL,NCO,NU*H
     X			,(INODES(I),i=1,min(6,NCO-1))
     	call flush(imp)
     	call flush(63)
     	endif
C  217 FORMAT(E14.6,E15.3,I6)
C  217 FORMAT(F14.6,G15.6,I6,I4,10i4)
! 217 FORMAT(F14.6,G15.6,I6,f6.1,1x,6f6.1)
 2170 FORMAT(1x,i2,F13.6,G15.6,2i3,f6.1,2l2,f11.6)
 2171 FORMAT('#',F13.6,G15.6,I6,f6.1,1x,6f6.1)
! 218 FORMAT(14x,15x,6x,4x,1p,4e10.2)
! 219 FORMAT(14x,15x,6x,4x,4f10.3)
  
      	if(NMISS>3.and.INNERTP.and.ISOUT) INNERTP=.false.
      IF(NCO.NE.NODES) GO TO 6
      Q=DEL/(INTT*H)
!      IF(COUNT.GT.50.and.ABS(Q).lt.1e-3 )GO TO 7
      IF(ABS(DEL).LT.ERROR)GO TO 7
      IF(ABS(P).GT.ERROR) then
	if(ABS(Q/P).lt.ERROR)GO TO 7
	endif
      QMAX = 10.
      if(P>100.) QMAX = P*0.25
      if(abs(Q).gt.QMAX) Q = SIGN(QMAX,Q)
   50 P = P+Q
      BC=1
      Q=2*Q
   51 NLAST = NCO
   52 EIG=P
   	LNODE = NCO
       GO TO 102
    6 IF(BC.NE.0)GO TO 61
      IF(NLAST.ne.0.and.(NCO-NODES)*(NLAST-NODES).lt.0)GO TO 61
      CC=1
      IF(NCO.ge.NODES)CC=-1
 
c     CC=-CC
c     IF(CC.LT.0)CCP=-IABS(CCP)
c     IF(CCP.GT.0)PP=P
c     IF(CCP.LT.0)PP=0.5*PP
      PPI=PPI*1.5
      Q=PPI*CC
      P = P+Q
      Q=2*Q
      NMISS = max(NMISS-1,0)		! do not count starting trials
      GO TO 51
   61 Q=0.5*Q
      P=P-0.5*Q
      GO TO 52
    7 INTP=INTP*H
      INTP=SQRT(1d0/INTP)
 
      Q = 0d0
      G = 4d0/3d0
      DO 8 N=1,NMAX
      X=N*H
      T=1.d0
      IF(N.LE.NU)T=RATIO
      Y(N)=Y(N)*INTP*T
c	write(4,*)  N,Y(N),T,INTP
C   ------------------- Scale by 1/r if required:
C      Y(N)=Y(N)/X
	if(N.eq.NMAX) G=1d0/3d0
C   ------------------- Scale by sqrt(m*)
      Y(N)=Y(N)*sqrt(EM(N))
	Q = Q + H*G*Y(N)**2
	G=2d0-G
    8 CONTINUE
	Q = sqrt(1d0/Q)
	do 9 N=1,NMAX
    9   Y(N)=Y(N)*Q 

      IF(PCON.gt.0)WRITE(IMP,991)P,NODES,COUNT,NMISS,DEL,Y(NMAX)
      call flush(imp)
C 991 FORMAT(///' EIG ',12H VARIED FROM,F9.5,4H TO ,F9.5)
  991 FORMAT(' Eigenstate found at',G12.5,' (',i3,' nodes after'
     X  ,2I3,' ) DEL,Y(NMAX) =',1p,2e10.2)
      IF(PCON.lt.0.or.MOD(PCON,2).EQ.0)RETURN
      WRITE(IMP,992)
  992 FORMAT(//'   N     R      Y(R)')
      WRITE(IMP,998)(N,N*H,Y(N),N=10,NMAX,10)
  998 FORMAT(4(I4,F8.3,2X,E11.4,3X))
      RETURN
  200 IFAULT = 1
      WRITE(IMP,986) IFAULT
  986 FORMAT(' Eigenstate FAILURE- ROUTINE SEARCH ',I5)
      WRITE(IMP,992)
      do 987 N=1,NMAX*0
  987 WRITE(IMP,*) N*real(H),Y(N)
      RETURN
      END


	Subroutine HEMATRIX(Nstate,alpha,nsturm,SS,VVsturm,numr,
     X		AA,nd,nch,hcm,FORMF,CLIST,NCLIST,NFLIST,NF,LCHL,coef,
     X 		EL,ECM,nbas0,nbmax,symm,prba,MR,CHNO,NL,NLO,NLN,
     X		NLC,MLT,MLM,CUTOFF,ICUTC,ISNONO,kk,kkdim,vsearch,
     X            emass,vemass,PTYPE,CPSO)
	use parameters
	implicit real*8 (a-h,o-z)
	real*8 VVsturm(numr,Nstate),ECM(Nstate),kk(kkdim,kkdim)
	real*8 SS(numr,Nstate,0:nbmax),alpha(Nstate,nbmax)
	real*8 NN,hcm(Nstate),wrad(numr),coef(Nstate),EE,
     x        emass(numr,Nstate),vemass(numr,Nstate),
     X	  FNL(NLN,NLO),VNL(1:numr,1:MLM),EXF(NLN),VR,wwsrch(numr)
	integer LCHL(Nstate),EL,nsturm(Nstate),nstmax,nbas0(Nstate),
     X    NFLIST(MAXCH,MAXCH,MCLIST),NCLIST(MAXCH,MAXCH),
     X    MR,NL,CHNO(MFNL,6),NLO,NLN,NLC,MLT,MLM,ICUTC,C,D,
     X    CUTOFF(Nstate),vsearch,PTYPE(12,MLOC)
	complex*16 FORMF(MAXM,MLOC),CLIST(MAXCH,MAXCH,MCLIST),VV,VV1
	complex*16 AA(nd,nd),ww(numr),SSC(nbmax),wwd(numr),
     X	  FNC(NLN,NLO),VNC(1:numr,1:MLM),ECF(NLN)
	complex*16,allocatable:: AAS(:,:)
	logical pertcent,coupled,blas1,blas2,symm,FFR,SH,NREV,ISNONO,
     X		prba,incl,derd,CPSO
	parameter (pertcent=.false.,blas1=.true.,blas2=.false.,
     X  	   SH=.false.)
!			Adjust blas1 and blas2 for best times in your system

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
	  write(70,*) (i-1)*real(hcm(1)),real(SS(i,1,l1))   !,wrad(i)
	  enddo
	  write(70,*) '&'
	enddo
c --------------------------------------------
c       matrix elements for MM=<SS|SS> delta(lji)
c --------------------------------------------
 	do kv1=1,Nstate
	do l1=1,nsturm(kv1)
	  k1 = nbas0(kv1)+l1
 	  do kv2=1,Nstate
	  do l2=1,nsturm(kv2)
	    k2 = nbas0(kv2)+l2
	  if(k2<=k1) then
 	    AA(k1,k2)=0d0
            if (kv1.eq.kv2) then
	    do i=1,numr
	      AA(k1,k2)=AA(k1,k2)+wrad(i)*hcm(kv1)*
     x		 	SS(i,kv1,l1)*SS(i,kv2,l2)
	    enddo
	    endif
            AA(k2,k1)=AA(k1,k2)
	    endif
          enddo
          enddo
	enddo
	enddo

 	   write(60,*) 'MM normalisation matrix:'
 	   call wrcmat(AA,nch,nch,nd,4,60)
  	 endif	
c --------------------------------------------
c       matrix elements for AA=<SS|V2|SS>*YY2
c --------------------------------------------
	
	AA(:,:) = 0d0
	if(ISNONO.or.vsearch>0) kk(:,:) = 0d0
 	do kv1=1,Nstate
	  ryev = -coef(kv1)
	  kv2l=Nstate
	  if(symm) kv2l=kv1
	do kv2=1,kv2l
	  ww(:) = 0d0; wwd(:) = 0d0
	  if(kv1==kv2) ww(:) = vemass(:,kv1)*ryev
	  if(vsearch>0) wwsrch(:) = 0.
	  coupled=.false.; derd=.false.
	  do NC=1,NCLIST(kv1,kv2)
	   ifm = NFLIST(kv1,kv2,NC)
!	    incl = vsearch==0 .or. vsearch>0.and.ifm.ne.vsearch
	    incl = ifm.ne.vsearch
	   if(abs(CLIST(kv1,kv2,NC))>1e-20) then
	     coupled=.true.
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
	      k20 = nbas0(kv2)
	  do i=1,numr
	    rrad = (i-1)*hcm(kv1)
	       VV = ww(i) * emass(i,kv1)
	      if(pertcent.and.kv1.eq.kv2) 
     & 		 VV = VV + ryev*(LCHL(kv1)*(LCHL(kv1)+1)
     & 		                -LCHL(EL)*(LCHL(EL)+1))/rrad**2
	      VV = VV * wrad(i) * hcm(kv1) 
     & 			* (0.,1.)**(LCHL(kv1)-LCHL(kv2))
	    do l1=1,nsturm(kv1)
	      k1=nbas0(kv1)+l1
	      VV1 = VV * SS(i,kv1,l1)
	      l2max = nsturm(kv2)
	      if(symm) l2max = min(l2max,k1 - k20)
c					l2max so k2.le.k1 if symm.
	    if(l2max>0) then
	    if(blas2) then
!				SS is real not complex!
	        SSC(1:l2max) = SS(i,kv2,1:l2max)
	        call zaxpy(l2max,VV1,SSC,1,AA(k1,k20+1),nd)
	      else
!	    	do l2=1,l2max
!	     	AA(k1,k20+l2) = AA(k1,k20+l2) + SS(i,kv2,l2) * VV1
	     	AA(k1,k20+1:k20+l2max) = AA(k1,k20+1:k20+l2max) + 
     X			SS(i,kv2,1:l2max) * VV1
!             	enddo
	      endif ! blas2
	      endif ! l2max>0
             enddo ! l1
	  enddo ! i
!				acting on WF derivatives to right (less optimised code!)
	 if(derd) then
          do i=2,numr-1
            rrad = (i-1)*hcm(kv1)
               VV = wwd(i) * emass(i,kv1)
              if(pertcent) stop 'percent & derivatives not implemented'
              VV = VV * wrad(i) * hcm(kv1) / hcm(kv2) / 2.0
     &                  * (0.,1.)**(LCHL(kv1)-LCHL(kv2))
            do l1=1,nsturm(kv1)
              k1=nbas0(kv1)+l1
              VV1 = VV * SS(i,kv1,l1)
              l2max = nsturm(kv2)
              if(symm) l2max = min(l2max,k1 - k20)
c                                       l2max so k2.le.k1 if symm.
                AA(k1,k20+1:k20+l2max) = AA(k1,k20+1:k20+l2max) +
     X            (SS(i+1,kv2,1:l2max)-SS(i-1,kv2,1:l2max)) * VV1
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
	       VV = - wwsrch(i)  * emass(i,kv1) ! negative because V is eigenvalue
	      VV = VV * wrad(i) * hcm(kv1) 
     & 			* (0.,1.)**(LCHL(kv1)-LCHL(kv2))
	    do l1=1,nsturm(kv1)
	      k1=nbas0(kv1)+l1
	      VV1 = VV * SS(i,kv1,l1)
	      l2max = nsturm(kv2)
	      if(symm) l2max = min(l2max,k1 - k20)
c					l2max so k2.le.k1 if symm.
	    if(l2max>0) then
	     	kk(k1,k20+1:k20+l2max) = kk(k1,k20+1:k20+l2max) + 
     X			SS(i,kv2,1:l2max) * VV1
	      endif ! l2max>0
             enddo ! l1
	    enddo ! i
	   endif ! vsearch
	  endif ! coupled
	enddo  ! kv2
	enddo  ! kv1
	if(symm) then
	do 35 k1=1,nch
	do 35 k2=1,k1-1
35	AA(k2,k1) = AA(k1,k2) 
	 if(vsearch>0) then
	 do 37 k1=1,nch
	 do 36 k2=1,k1-1
36	  kk(k2,k1) = kk(k1,k2) 
37	  if(abs(kk(k1,k1))<1e-4) kk(k1,k1) = 1e-4*sig
	 endif
	endif
!	   imp=160
!	   write(imp,*) 'V matrix:'
!	   call WRCMAT(AA,nd,nd,nd,4,imp)
c -------------------------------------------------
c       matrix elements for AA=<SS|V(nonlocal)|SS>
c -------------------------------------------------

      HI = 1.0 / DBLE(MR)
      N = numr
       if(NL>0) REWIND 12
      DO 50 INL=1,NL
         FFR = CHNO(INL,3).LT.10
       NLL = CHNO(INL,4)
 !      write(6,*) 'HEMATRIX: read NLL,NLO =',NLL,NLO,'; NLN=',NLN
          IF(FFR) then
	    READ(12) ((FNL(I,J),I=1,NLL),J=1,NLO)
	    VNL(:,:) = 0d0
	   else
            READ(12) ((FNC(I,J),I=1,NLL),J=1,NLO)
	    VNC(:,:) = 0d0
	   endif
         D = ABS(CHNO(INL,1))
         C = CHNO(INL,2)
         NONO = MOD(CHNO(INL,3),10)
         NREV = CHNO(INL,1).LT.0 .OR. C.EQ.D
         IF(SH) WRITE(6,22) INL,D,C,NLL,.NOT.NREV,FFR,NONO
22       FORMAT(' NL coupling #',I4,' to',I4,' from Ch.',2I4,
     &         ', Reverse=',L2,', Real =',L2,', NONO =',I2)
            SQH = SQRT(hcm(C)/hcm(D))
        IF(SH.AND.FFR) CALL DISPLY(FNL,NLN,NLO,NLN,SCALE)
        IF(SH.AND..not.FFR) CALL DISPLR(FNC,NLN,NLO,NLN,SCALE)
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
            VNL(I,J) = VR * emass(i,D)
29           I = I + MR
          ELSE
            IF(II.EQ.1) CALL ECPAND(FNC,NLN,NLL,NLO,ECF,J,MLT)
             I = (KMIN-1)*MR + II
            DO 30 K=KMIN,KMAX
             VV = (-P2*ECF(K-1)+Q*ECF(K+2))*X + (P1*ECF(K)-P*ECF(K+1))*Y
            VNC(I,J) = VV * emass(i,D)
30           I = I + MR
         ENDIF
31      CONTINUE
!       IF(SH.AND.FFR) CALL DISPLY(VNL,numr,MLM,numr,SCALE)
!       IF(SH.AND..not.FFR) CALL DISPLR(VNC,numr,MLM,numr,SCALE)
!          SCALE = SCALE * 1E-6
!        IF(SH.AND..not.FFR) CALL DISPLI(VNC,numr,MLM,numr,SCALE)

!	A(kv1,l1;kv2,l2) = * + <SS(kv1,l1,:) | V | SS(kv2,l2,:)>
!                        = * + int(i,j) <SS(kv1,l1,i) | V(i,j) | SS(kv2,l2,j)>
	 kv1 = D
	 kv2 = C
	      k10 = nbas0(kv1)
	      k20 = nbas0(kv2)
	    do l1=1,nsturm(kv1)
	      k1=k10 + l1
	      l2max = nsturm(kv2)
!	      if(symm) l2max = min(l2max,k1 - k20)
c					l2max so k2.le.k1 if symm.
	    do l2=1,l2max
		VV1 = 0d0
	 	VR = 0d0
		do J=1,MLM
             	JJ = (J - NLC*MLT  - 1)
             	IMIN = 1 + MAX(-JJ ,0)  +  MAX(CUTOFF(D),ICUTC)
             	IMAX = N - MAX(JJ ,0)   -  5
         	 IF(FFR) THEN
	         do I=IMIN,IMAX
	           VR = VR + SS(I,kv1,l1)*VNL(I,J)*SS(I+JJ,kv2,l2)
         	   enddo
          	 ELSE
	         do I=IMIN,IMAX
	          VV1 = VV1 + SS(I,kv1,l1)*VNC(I,J)*SS(I+JJ,kv2,l2)
         	   enddo
	         endif
	        enddo
	         if(FFR) VV1 = VR
		VV1 = VV1 * HCM(kv1)
	    IF(NONO>1) then  			 ! Do normal couplings
	     	AA(k1,k20+l2) = AA(k1,k20+l2) + VV1
	     	if(.not.NREV) AA(k20+l2,k1) = AA(k20+l2,k1) + VV1
	     else				 ! Do NONO overlaps
	       if(vsearch==0) then
	     	kk(k1,k20+l2) = kk(k1,k20+l2) + VV1
	     	if(.not.NREV) kk(k20+l2,k1) = kk(k20+l2,k1) + VV1
               endif
		if(ISNONO.and.vsearch==0) then
		 if(NONO==0) then  	! NONO post contribution to AA
	     	   AA(k1,k20+l2) = AA(k1,k20+l2) +
     x				VV1*(alpha(kv1,l1)-ECM(kv1))
	     	   if(.not.NREV) AA(k20+l2,k1) = AA(k20+l2,k1) + 
     x				VV1*(alpha(kv1,l1)-ECM(kv1))
		 else			! NONO prior contribution to AA
	     	   AA(k1,k20+l2) = AA(k1,k20+l2) + 
     x				VV1*(alpha(kv2,l2)-ECM(kv2))
	write(60,48) k1,k20+l2,VV1,alpha(kv2,l2),ECM(kv2),
     x        VV1*(alpha(kv2,l2)-ECM(kv2))
48	format(' AA(',2i4,'): add ',2f8.4,'*(',f8.4,'-',f8.4,') =',2f8.4)
	     	   if(.not.NREV) AA(k20+l2,k1) = AA(k20+l2,k1) + 
     x				VV1*(alpha(kv2,l2)-ECM(kv2))
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

 	do kv1=1,Nstate
	do l1=1,nsturm(kv1)
	      k1=nbas0(kv1)+l1
	    AA(k1,k1)=AA(k1,k1)+alpha(kv1,l1)-ECM(kv1)  ! * emass(i,kv1)  ! i ?
          enddo
	enddo
!	   write(imp,*) 'T+V-E matrix:'
!	   call WRCMAT(AA,nd,nd,nd,4,imp)
 
c --------------------------------------------
c       matrix elements for NN=<SS|Vo|SS> delta(lji)
c --------------------------------------------
 	do kv1=1,Nstate
	do l1=1,nsturm(kv1)
	      k1=nbas0(kv1)+l1
	    ryev = -coef(kv1)
	  do l2=1,nsturm(kv1)
	    k2 = nbas0(kv1)+l2
 	    NN=0d0; EE=0d0
	    do i=1,numr
	      NN=NN+wrad(i)*SS(i,kv1,l1)*SS(i,kv1,l2)*VVsturm(i,kv1)
	      EE=EE+wrad(i)*SS(i,kv1,l1)*SS(i,kv1,l2)*(1-emass(i,kv1))  ! correction to -E term
	    enddo  ! VVsturm (bpot) already has  * emass(i) factor
	     AA(k1,k2)=AA(k1,k2)-NN*ryev*hcm(kv1) + EE*hcm(kv1)*ECM(kv1)
          enddo
	enddo
	enddo
!	   write(imp,*) 'T+V-V0-E matrix:'
!	   call WRCMAT(AA,nd,nd,nd,4,imp)

	if(CPSO.and..true.) then ! Symmetrise the matrix explicitly with folded spin-orbit terms
	 write(6,*) ' SYMMETRIZING THE H-E MATRIX !'
	 allocate(aas(nd,nd))
	 do k1=1,nd
	 do k2=1,nd
	  aas(k1,k2) = (aa(k1,k2) + aa(k2,k1))*0.5d0
	 enddo
	 enddo
	  aa(:,:) = aas(:,:)
	  deallocate(aas)
	endif

	RETURN
	END
c --------------------------------------------
      SUBROUTINE WRRMAT (A,N,M,L,NCOL,IOUT)
      IMPLICIT INTEGER (A-Z)
      REAL*8 A
      DIMENSION A(L,*)
*
*     WRCMAT : PRINTS A REAL N*M MATRIX STORED IN AN L*L ARRAY
*
      NTIM = M / NCOL
      NF = 0
      IF ( NTIM .GT. 0 ) THEN
         DO 20 I = 1, NTIM
            NI = NF + 1
            NF = NF + NCOL
            DO 10 J = 1, N
               WRITE (IOUT,1000) (A(J,K),K = NI,NF)
   10       CONTINUE
            WRITE (IOUT,1010)
   20    CONTINUE
      ENDIF
      NI = NF + 1
      IF ( NI .LE. M ) THEN
         DO 40 J = 1, N
            WRITE (IOUT,1000) (A(J,K),K = NI,M)
   40    CONTINUE
      ENDIF
      RETURN
*
c1000 FORMAT (7G16.7)
c1000 FORMAT (7F16.8)
c1000 FORMAT (8F13.8)
 1000 FORMAT (8G12.4)
 1010 FORMAT (/)
      END

      SUBROUTINE WRCMAT (A,N,M,L,NCOL,IOUT)
      IMPLICIT INTEGER (A-Z)
      COMPLEX*16 A
      DIMENSION A(L,*)
*
*     WRCMAT : PRINTS A COMPLEX N*M MATRIX STORED IN AN L*L ARRAY
*
      NTIM = M / NCOL
      NF = 0
      IF ( NTIM .GT. 0 ) THEN
         DO 20 I = 1, NTIM
            NI = NF + 1
            NF = NF + NCOL
            DO 10 J = 1, N
               WRITE (IOUT,1000) (A(J,K),K = NI,NF)
   10       CONTINUE
            WRITE (IOUT,1010)
   20    CONTINUE
      ENDIF
      NI = NF + 1
      IF ( NI .LE. M ) THEN
         DO 40 J = 1, N
            WRITE (IOUT,1000) (A(J,K),K = NI,M)
   40    CONTINUE
      ENDIF
      RETURN
*
c1000 FORMAT (7G16.7)
c1000 FORMAT (7F16.8)
c1000 FORMAT (8F13.8)
!1000 FORMAT (4(G12.4,G9.2))
 1000 FORMAT (4(F12.6,F9.6))
 1010 FORMAT (/)
      END


      SUBROUTINE GETK(NC,IOP,NOP,CFMAT,CGMAT,FCWFN,CH,MAXCH,RMAT,
     X			BETA,LVAL,KMAT,IOUT,F,FP,PR)
C
      IMPLICIT INTEGER (A-Z)
      REAL*8 F,FP,ZERO,DF,ONE,CFMAT,CGMAT,BETA(NC)
      COMPLEX*16 AA,BB,RMAT,KMAT,CH,CI,H,HP,PHASE,AAA
      LOGICAL FCWFN,PR
      PARAMETER ( ZERO=0.0D0 , ONE=1d0, CI=(0d0,1d0))
      DIMENSION F(NC,NC,2),FP(NC,NC,2),RMAT(NC,NC),KMAT(NC,NC),
     X          AA(NC,NC),BB(NC,NC),IPVT(NC),IOP(NOP),CH(2,NC),
     &          CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),LVAL(NC)
     x   , AAA(NC,NC)
*
*   Put appropriate Asymptotic wave functions in F and FP:
	if(.not.FCWFN) then
	  F(:,:,:) = 0d0; FP(:,:,:) = 0d0
*				Uncoupled asymptotics
	 do 200 I=1,NC
	   F(I,I,1)  = aimag(CH(1,I))
	   F(I,I,2)  = real(CH(1,I))
	   FP(I,I,1) = aimag(CH(2,I))
  200	   FP(I,I,2) = real(CH(2,I))
!	write(imp,*) ' F,G,FP,GP =',(F(1,1,J),J=1,2),(FP(1,1,J),J=1,2)
!	write(KO,*) ' F,G,FP,GP =',(F(1,1,J),J=1,2),(FP(1,1,J),J=1,2)
	else
*				Coupled asymptotics
	 do 210 I=1,NC
	 do 210 J=1,NC
!	    PHASE = CI**(LVAL(J)-LVAL(I))
	    PHASE = 1d0
	    H = dcmplx(CGMAT(I,J,1),CFMAT(I,J,1)) * PHASE
	    HP= dcmplx(CGMAT(I,J,2),CFMAT(I,J,2)) * PHASE
	   F(I,J,1)  = aimag(H)
	   F(I,J,2)  = real(H)
	   FP(I,J,1) = aimag(HP)
  210	   FP(I,J,2) = real(HP)
	endif

	  if(PR) then
          WRITE (IOUT,*) 'GETK F-matrix:'
	  CALL WRTMAT (F(1,1,1),NC,NC,NC,10,IOUT)
          WRITE (IOUT,*) 'GETK G-matrix:'
	  CALL WRTMAT (F(1,1,2),NC,NC,NC,10,IOUT)
!          WRITE (IOUT,*) 'GETK FP-matrix:'
!	  CALL WRTMAT (FP(1,1,1),NC,NC,NC,10,IOUT)
!          WRITE (IOUT,*) 'GETK GP-matrix:'
!	  CALL WRTMAT (FP(1,1,2),NC,NC,NC,10,IOUT)
	  endif
	   !write(60,'(10f10.5)') BETA
!	if(PR) then
!	write(imp,*) ' F,G,FP,GP =',(F(1,1,J),J=1,2),(FP(1,1,J),J=1,2)
!	 do 220 I=1,NC
! 220    write(imp,225) I,(F(I,I,J),J=1,2),(FP(I,I,J),J=1,2)
! 225    format(' Diagonal F,G,FP,GP in ch ',i3,' :',4f12.6)
!	endif

*	   	  Calculate BB =  [F-RF']
      DO 80 IOJ = 1, NOP
	J = IOP(IOJ)
         DO 70 I = 1, NC
            BB(I,IOJ) = - F(I,J,1)
   70    CONTINUE
   80 CONTINUE
      DO 110 IOJ = 1, NOP
	J = IOP(IOJ)
         DO 100 K = 1, NC
            DF = FP(K,J,1) - BETA(K) * F(K,J,1)
            DO 90 I = 1, NC
               BB(I,IOJ) = BB(I,IOJ) + RMAT(I,K) * DF
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
	  if(PR) then
          WRITE (IOUT,*) 'GETK B-matrix:'
          CALL WRCMAT (BB,NC,NOP,NC,4,IOUT)
	  endif
*
      LHS = 1
      DO 20 J = 1, NC
         K = 1
	 DO 5 IO=1,NOP
    5	 IF(J.eq.IOP(IO)) K = 2
         DO 10 I = 1, NC
            AA(I,J)   = F(I,J,K)
   10    CONTINUE
   20 CONTINUE
      DO 50 J = 1, NC
         KJ = 1
	 DO 25 IO=1,NOP
   25	 IF(J.eq.IOP(IO)) KJ = 2
!	write(IOUT,*) 'Channel ',J,' has KJ =',KJ
         DO 40 K = 1, NC
            DF = FP(K,J,KJ) - BETA(K) * F(K,J,KJ)
            DO 30 I = 1, NC
               AA(I,J) = AA(I,J) - RMAT(I,K) * DF
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      DO 60 I = 1, NC
         IF ( AA(I,I) .EQ. ZERO ) LHS = 0
   60 CONTINUE
	  if(PR) then
          WRITE (IOUT,*) 'GETK A-matrix:'
          CALL WRCMAT (AA,NC,NC,NC,4,IOUT)
	    write(IOUT,*) 'GETK: NC,NOP,LHS =',NC,NOP,LHS
	    AAA(:,:) = AA(:,:)
	    do i=1,NC
	    do j=1,NC
	    AAA(i,j) = AA(i,j)/F(j,j,2)
	    enddo
	    enddo
          WRITE (IOUT,*) 'GETK (1-R.L0)-matrix:'
          CALL WRCMAT (AAA,NC,NC,NC,4,IOUT)
	    if(nc==2) then
	  DF = AAA(1,1)*AAA(2,2) - AAA(1,2)*AAA(2,1)
	  write(IOUT,'('' dets ='',1p,4e12.4)') DF,AAA(1,1)
	    endif
	  endif
      IF ( LHS .GT. 0 ) THEN
         IF ( NC .EQ. 1 ) THEN
!	    write(imp,*) 'K numerator =',BB(1,1)
!	    write(imp,*) 'K denominator =',AA(1,1)
            BB(1,1) = BB(1,1) / AA(1,1)
!	    write(imp,*) 'K =',BB(1,1)
         ELSE
!	    write(KO,*) 'GETK: call ZGETRF @',NC,NOP
            CALL ZGETRF (NC,NC,AA,NC,IPVT,IER)
            IF ( IER .NE. 0 ) THEN
               WRITE (IOUT,1000) IER
               STOP
            ENDIF
!	    write(KO,*) 'GETK: call ZGETRS @',NC,NOP
            CALL ZGETRS('N',NC,NOP,AA,NC,IPVT,BB,NC,IER)
            IF ( IER .NE. 0 ) THEN
               WRITE (IOUT,1010) IER
               STOP
            ENDIF
!	    write(KO,*) 'GETK: ZGETRS done'
         ENDIF
         DO 130 IOJ = 1, NOP
            DO 120 IOI = 1, NOP
               KMAT(IOI,IOJ) = BB(IOP(IOI),IOJ)
  120       CONTINUE
  130    CONTINUE
      ENDIF
      RETURN
 1000 FORMAT (' GETK : ERROR RETURN FROM ZGETRF, IERR = ',I6)
 1010 FORMAT (' GETK : ERROR RETURN FROM ZGETRS, IERR = ',I6)
      END
      SUBROUTINE GETS (KMAT,SMAT,PRTM,IOP,NOP,NC,LVAL,IOUT)
      IMPLICIT INTEGER (A-Z)
      COMPLEX*16 KMAT,SMAT,TMAT,A,TWOI
      REAL*8 ONE
      LOGICAL PRTM
      PARAMETER ( ONE=1.D0, TWOI=(0d0,2d0))
      DIMENSION KMAT(NC,NC),TMAT(NOP,NOP),A(NC,NC),SMAT(NC,NC),
     X    IPVT(NC),IOP(NOP),LVAL(NC)
*
!          WRITE (IOUT,*) 'GETS starts with K-matrix:'
!          CALL WRCMAT (KMAT,NOP,NOP,NC,4,IOUT)

*     TMAT : T-matrix computed from K-matrix
*            S = ( 1 + I * K ) / ( 1 - I * K ) ;  S = DEL + T
*            T = 2 * I * K / ( 1 - I * K )
C
C  Solve T = (1 - i K)^-1 * 2i K   
C  by    (1-i K) * (T) = 2i K
*
      IF ( NC == 1 ) THEN
         TMAT(1,1) = TWOI*KMAT(1,1)/(ONE - (0d0,1d0)*KMAT(1,1))
!	write(imp,*) 'K,T,S =',KMAT(1,1),TMAT(1,1),1+TMAT(1,1)
      ELSE
C
         DO 10 I = 1, NOP
         DO 10 J = 1, NOP
	  A(I,J) = (0d0,-1d0) * KMAT(I,J)
	  if(IOP(I).eq.IOP(J)) A(I,J) = A(I,J) + ONE
	  TMAT(I,J) = TWOI * KMAT(I,J)
10	CONTINUE
!     IF (PRTM) THEN
! 	write(IOUT,*) ' 1-iK:'
!          CALL WRCMAT (A,NOP,NOP,NC,4,IOUT)
! 	write(IOUT,*) ' RHS:'
!          CALL WRCMAT (TMAT,NOP,NOP,NOP,4,IOUT)
!	ENDIF

c   A =P*L*U factorisation:
         CALL ZGETRF (NOP,NOP,A,NC,IPVT,IER)
         IF ( IER .lt. 0 ) THEN
            WRITE (IOUT,1000) IER
            STOP
         ELSE IF ( IER .gt. 0 ) THEN
            WRITE (IOUT,1001) IER
	    A(IER,IER) = ONE
         ENDIF

C   Solve A*T = B
!	    write(KO,*) 'GETS: call ZGETRS @',NC,NOP
         CALL ZGETRS('N',NOP,NOP,A,NC,IPVT,TMAT,NOP,IER)
         IF ( IER .NE. 0 ) THEN
            WRITE (IOUT,1010) IER
            STOP
         ENDIF
      ENDIF

!      IF (PRTM) THEN
!          WRITE (IOUT,1020)
!          CALL WRCMAT (TMAT,NOP,NOP,NOP,4,IOUT)
!      ENDIF
*
*  Convert T to S matrix 
*  and expand open-channels to all-channels matrix
*  and put back the i**L factors
	SMAT(:,:) = 0d0
         DO 30 IOI = 1, NOP
         DO 30 IOJ = 1, NOP
   30      SMAT(IOP(IOI),IOP(IOJ)) = TMAT(IOI,IOJ) 

	 DO 40 I=1,NC
         SMAT(I,I) = SMAT(I,I) + ONE
   40    CONTINUE
C
      IF (PRTM) THEN
	write(IOUT,*) ' S matrix (without i**L factors):'
	DO 55 I=1,NC
   55   write(IOUT,56) I,(SMAT(I,J),J=1,NC)
!  56   format(1x,I3,18f8.5,:,/(4x,18f8.5))
   56   format(1x,I3,1p,10e10.2,:,/(4x,10e10.2))
   57   format(1x,I3,18f8.4,:,/(4x,18f8.4))
      ENDIF
	do 50 i=1,NC
	do 50 j=1,NC
50	SMAT(I,J) = SMAT(I,J) * (0.,1.)**(LVAL(J)-LVAL(I))
      IF (PRTM) THEN
	write(IOUT,*) ' S matrix (WITH i**L factors):'
	DO 65 I=1,NC
   65   write(IOUT,451) I,(SMAT(I,J),J=1,NC)
  451   format(1x,I3,18f10.5,:,/(4x,18f10.5))
      ENDIF

      RETURN
*
 1000 FORMAT (' GETS : Error return from GETS ZGETRF, IER =',I6)
 1001 FORMAT (' GETS : ZGETRF: 1-iK matrix singular after ',I6)
 1010 FORMAT (' GETS : Error return from GETS ZGETRS, IER =',I6)
!1020 FORMAT ('  T - MATRIX :')
!1030 FORMAT ('  S - MATRIX :')
      END
        subroutine nopen(Nstate,ECM,iop,nop)
        implicit real*8 (a-h,o-z)
        integer iop(Nstate)
        real*8  ECM(Nstate)

        nop = 0
        do 50 i =1,Nstate
        if(ECM(i)>0.) then
          nop = nop + 1
          iop(nop) = i
          endif
50      continue
        return
        end

	subroutine BarkerTransform(term,set,mterms,En,vec,nch,
     x       btype,beta,J,p,Lval,Etaco,ENff,coef,RMA,debug,final,
     x       E_Brune,W_Brune)
	! use searchpar, only: E_Brune, W_Brune
        implicit real*8(A-H,O-Z)
	real*8 vec(nch,mterms),new(nch),beta(nch),J,Penetr(nch),
     x          Etaco(nch),ENff(nch),coef(nch),RMA(nch),
     x          E_Brune(0:mterms),W_Brune(0:mterms)
        character p,btype
        complex*16 En(mterms)
 	integer term,terms(mterms),c,lev,lval(nch)
	logical set(mterms),debug,final
        complex*16 CF2,RHO,ETA,PM,LOGD
	real*8,allocatable:: aa(:,:),evec(:,:)
	real*8,allocatable:: ens(:),ensw(:)
	integer,allocatable:: ior(:),iw(:)
	 
 	maxiter = 1000
        eps = 1e-5
	! write(6,*) 'Barker for term=',term
	! write(6,*) 'En:',En

!	 write(6,*) "C,L,B,ETACO,ENff,coef,RMA:"
	do c=1,nch
!	write(6,'(2i4,6f9.4)') c,Lval(c),beta(c)*RMA(c),Etaco(c),
!     x            ENff(c),coef(c),RMA(c)
	enddo

	levels = 0
	do i=1,mterms
         if(set(i)) then
	   levels = levels+1
	   terms(levels) = i
           if(i==term) lev = levels
           !write(6,*) i,'gam:',vec(:,i)
	 endif
 	enddo
	E = En(term)
	!write(6,*) 'lev,terms,levels,E:',lev,terms(1:levels),levels,E
	allocate(aa(levels,levels),evec(levels,levels))
	allocate(ens(levels),ensw(levels),ior(levels),iw(levels))

        new(:) = vec(:,term) ! for btype=A or S, no change to widths

	do iter=1,maxiter

	  aa(:,:) = 0d0
	  do l=1,levels
           lt = terms(l)
           aa(l,l) = En(lt)
           do c=1,nch
            RR = RMA(c)
            LL = Lval(c)
	    Ecm = ENff(c) + E
	    ETA = Etaco(c)/sqrt(abs(Ecm))
            wn = sqrt(abs(-Ecm/coef(c)))
	    RHO = RR * wn
             if(Ecm<0.) then
	        RHO = dcmplx(0d0,+real(RHO))
	        ETA = dcmplx(0d0,-real(ETA))
	  	endif
            PM = (0d0,1d0)
! 	    write(6,*) iter,J,p,E,RHO,ETA
            LOGD=CF2(RHO,ETA,dcmplx(LL),PM,1d-7,AERR,20000,1d-15)*RHO
 	    Shift = real(LOGD)
            Penetr(c)= aimag(LOGD)
            !write(6,*) 'Ch',c
 	
	    do k=1,levels
             kt = terms(k)
 	     aa(k,l) = aa(k,l)-vec(c,lt)*vec(c,kt) * (Shift-beta(c)*RR)
             !write(6,*) lt,kt,'g,g,S-B:',vec(c,lt),vec(c,kt)
	    enddo
             !write(6,*) 'aa(1,2) = ',aa(1,2)
	   enddo
	  enddo
           if(btype=='A' .or. btype=='S') go to 30
		
	  !write(6,*) 'aa:',aa
          call HDIAG(aa,levels,levels,0,evec,nc)
	   do l=1,levels
            ens(l)= aa(l,l)
	    enddo
	    ! write(6,*) 'Evals1:',ens,'from',(aa(l,l),l=1,levels)
            call sortreals(ens,ensw,ior,iw,levels,levels,ifail)
	    ! write(6,*) 'Evals2:',ior,(ens(ior(l)),l=1,levels)
          Enew = aa(ior(lev),ior(lev))
          err = abs(Enew-E)
          if(iter<maxiter/4) then
              E = Enew               ! to start with
          else if( iter<maxiter/3) then
              E = (E + Enew)*0.5     ! getting desparate
          else if( iter<maxiter/2) then
              E = (3*E + Enew)*0.25  ! getting more desparate
          else
              E = (9*E + Enew)*0.1   ! getting most desparate
          endif

           if(debug) write(60,5) lev,iter,E,Enew,err
5	  format(' Barker #',i3,' @',i4,': E,Enew='2f12.6,1p,e10.3)
          if(err < eps.and.iter>2) exit

	enddo

	 new(:) = 0.
	 do k=1,levels
	   do c=1,nch
	     new(c) = new(c) + evec(lev,k) * vec(c,terms(ior(k)))
	   enddo
	 enddo

30	E_Brune(term) = E
	
	t = 0.
	 do c=1,nch
	 t = t + 2*new(c)**2 * Penetr(c)
         enddo
	W_Brune(term) = t
        if(debug.or.final) 
     x     write(6,10) J,p,term,real(En(term)),E,iter,t
10	format(' Barker convergence in J,pi=',f5.1,a1,' #',i2,' from ',
     x   f10.6,' to ',f10.6,' took ',i5,' iterations. GF=',f8.3,' MeV')

	deallocate(aa,evec,ens,ensw,ior,iw)
	end
      subroutine sortreals(a,w,ind,indw,n,nw,ifail)
c sorts a vector of real numbers into ascending order and provides an in
c indicating the position of the sorted numbers in the original array.
      double precision  a(n), w(nw)
      dimension ind(n),indw(nw),ist(20),kount(20)
      k = 0
      if (n-1) 20, 680, 40
   20 k = 1
   40 if (nw) 640, 640, 60
   60 if (k) 660, 80, 660
   80 nn = min(n,nw)
      do 100 k = 1, n
         ind(k) = k
  100 continue
      l3 = 1
      iko = 0
  120 iko = iko + 1
  140 nm = 1
      l1 = l3
  160 l2 = min(l1+nm,n+1)
      l3 = min(l2+nm,n+1)
      kount(iko) = l3 - l1
      ist(iko) = l1
  180 iw = 1
      il1 = l1
      il2 = l2
  200 ndo = min(il2-l1,l3-l2,nn-iw+1)
      if (ndo) 320, 320, 220
  220 do 300 i = 1, ndo
         if (a(l1)-a(l2)) 240, 240, 260
  240    w(iw) = a(l1)
         indw(iw) = ind(l1)
         l1 = l1 + 1
         go to 280
  260    w(iw) = a(l2)
         indw(iw) = ind(l2)
         l2 = l2 + 1
  280    iw = iw + 1
  300 continue
      go to 200
  320 if (il2-l1) 340, 340, 460
  340 iw1 = iw - 1
  360 do 380 i = 1, iw1
         a(il1) = w(i)
         ind(il1) = indw(i)
         il1 = il1 + 1
  380 continue
  400 if (iko-1) 420, 420, 580
  420 if (l3-n) 440, 440, 700
  440 iko = 2
      go to 140
  460 if (l3-l2) 480, 480, 520
  480 j = l3 - 1
      k = il2 - l1
      j1 = il2 - 1
      do 500 i = 1, k
         a(j) = a(j1)
         ind(j) = ind(j1)
         j = j - 1
         j1 = j1 - 1
  500 continue
      iw1 = iw - 1
      if (iw1) 400, 400, 360
  520 j = il2 - 1
      k = l2 - 1
      n21 = il2 - l1
      do 540 i = 1, n21
         a(k) = a(j)
         ind(k) = ind(j)
         k = k - 1
         j = j - 1
  540 continue
      do 560 i = 1, nn
         a(il1) = w(i)
         ind(il1) = indw(i)
         il1 = il1 + 1
  560 continue
      l1 = k + 1
      go to 180
  580 if (l3-n) 600, 600, 620
  600 if (kount(iko-1)-nm-nm) 620, 620, 120
  620 iko = iko - 1
      nm = kount(iko)
      l1 = ist(iko)
      go to 160
  640 k = k + 2
  660 ifail = 1
      return
  680 ind(1) = 1
  700 ifail = 0
      return
      end
