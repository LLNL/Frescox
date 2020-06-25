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
	SUBROUTINE RMATRIXP(JTOTAL,NCH,MINTL,INITL,ECM,COEF,HP,ITC,CHL,
     X	  PART,EXCIT,LVAL,JVAL,JPROJ,JTARG,N,FORMF,NF,ITCM,
     X    CUTVAL,ISOCEN,bndx,weak,nbasi,nrbmin,nd1,KO,CDETR,WOUT,KFUS,
     X    pralpha,PCON,CH,LL1,NCHAN,KS,PARITY,PSIGN,K,RMASS,BLOCKD,
     X    CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,symm,CLIST,NCLIST,NFLIST,
     X    SIGJ,JCOEF,XS,FUSL,SMATS,SMATL,CHANS,DONE,JTMAX,JTMIN,SCALE,
     X    JEX,ABSEND,DONES,NTHRESH,RESM,THRJ,IF1,IF2,TOTJ,FUSJ,meigsi,
     X    MR,CHNO,NL,NLO,NLN1,NLC,MLT,MLM,ICUTC,ISNONO,FLIP,EXCH,GAP,
     X    ETA,FORML,AFRAC,QNF,MASS,DROPPED,NSTEPD,iams,TIME0,TIME0J,
     X    XCOEF,PTYPE,CFUSJ,SMALLCOUP,SMALLCHAN,SMALLS,CHPRES,RTURN,NEX,
     X	  VEFF,FNC,WNM,NSA,WAVES,WDISK,CFG,WRITTEN,ENLAB,phasin,linel,
     x    IEXCH,FUSLL,CPSO,say,ETOTAL,BAND)
	use parameters
	use searchpar
	use drier
	use fresco1, only: vsearch,echan,enodes,jleast,btype,KRM
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
     x    DROPPED,ki1,ki2,ib2,nbasi,mag,mode1,kfn,ifail,kkdim,ich,
     x    IOI,IOJ,ip2
	logical pralpha,prrm,prkm,prtm,prbut,changed,symm,pauli,prham,
     X 		FCWFN,FJSWTCH,ISNONO,tr,BLOCKD(NCH),FLIP,prba,prmats,
     X		WOUT,FAIL,SMALLJ,WRITTEN(299),r_added,
     X          rm_set(mterms),drop(NCH),allorder_exch,nonel(MAXCH),red,
     x          adjusted,em,CPSO,say,Brune
	real*8 jtotal,JCOEF,ECM(MAXCH,3),JVAL(NCH),JPROJ(NCH),
     X		COEF(NCH),HP(NCH),beta(NCH),CRCRAT(MAXCH),ABSEND,
     &           CFMAT(MAXCH,MAXCH,2),CGMAT(MAXCH,MAXCH,2),Eshift,bndx,
     X		 JEX(6,MXP,MXX),TOTJ,OUTJ,FUSJ,EPS,E1,AN,E,ELAST,WN,
     X		 EXCH(MAXCH,MAXCH),PNORMS(MXP),FORML(MAXNLN,MSP,2),
     X		 AFRAC(MXPEX,MXPEX,2,MSP),RINTP,tnorm,rms(2),rmst,
     X		 MASS(4,MXP+1),ac,an1,an2,rho2,rcore,am12,am123,
     X		 TM,TIME0J,SECOND,TIME0,g,XCOEF,PI,T4,SMALLCOUP,
     X		 SMALLCHAN,CHSIZES(MXPEX),RTURN,GAP,ETA(MXP,MXX),
     X           ENLAB,FCOUL(NCH,NCH,2),FCOULP(NCH,NCH,2),fracshift,
     X           rm_vec(nch,mterms),JTARG(NCH),unitar(NCH),rmr(NCH),
     x           rtrace,xsc,PMAX(MAXCH),vmin(NCH),weak,rmrad1,damp,
     x           HARDSPH(MAXCH),SHIFT(MAXCH),ggam,esh,LL1(NCH),EOF,
     x           FUSLL(LMAX1,2),rmdist(NCH),SV,JV,racah,SYMR(NCH),
     x           ETACO(NCH),ENFF(NCH),RMA(nch)
	REAL*8 S,T,kp,EIG,conv,K(MXP,MXX),RMASS(MXP),SIGJ(0:MXPEX),r,
     X 	 XS(NCH),JTMIN,JTMAX,THRJ(NTHRESH),RESM(NTHRESH),FUSL(1+NFUS1),
     X	 CFUSJ(NFUS),AMDSQS,CFG(MAXMUL,4),CF,CG,Z,phasin(MINTL),
     x   fc_s(lmax1),gc_s(lmax1),fcp_s(lmax1),gcp_s(lmax1),lgd,eta1,
     x   sph,cph,ETOTAL,rwamp,pene,gameps,rm_shift(nch,mterms),
     x   SYMS(NCH)

	complex*16 FORMF(MAXM,MLOC),CLIST(MAXCH,MAXCH,MCLIST)
   	COMPLEX*16 CHL(LMAX1,MXPEX,2),CH(2,MAXCH),SOMA(NCH,1),
     x     HMINp,HPLp,WFP(nch),WF,HMIN,HPL,phaze,WFN(nch),
     x     FNC(NLN1,NSA),WNM(NLN1,NSA,2),CI,C6,C7,VPOT(N,3),SRC(N,1),
     x     rm_energy(mterms)
	character PSIGN(3),SCALE(0:MXPEX),COMPARE
	character*5 eigen,wid
**** Local arrays:
!	real*8,allocatable:: ybas(:,:,:),alpha(:,:),kk(:,:),gswf(:,:),
!     x	   pvec(:),qvec(:),vec(:),val(:),phys(:),evec(:,:),aar(:,:),
!     x	   kkt(:,:),wfnxy(:,:),ovl(:),schro(:)
!	complex*16,allocatable:: aa(:,:),rhs(:,:),rvec(:),phic(:),yc(:),
!     x	   rhsb(:)
!	integer,allocatable:: ipiv(:)
!	logical,allocatable:: donell(:,:)
!	real*8 bpot(n,nch),phi3(nch),EB,rhsprob(nbasi*2)
!	real*8 emass(n,nch),rm(n),drm(n),d2rm(n),vemass(n,nch)
	complex*16 Rmat(nch,nch),Kmat(nch,nch),Smatr(nch,nch),TC,
     X		logd(nch),Rmati,Rex0,Rmat0,SMAT(nch),Rmmat(nch,nch)
        complex*16 RmatLS(nch,nch),Rmatsym(nch,nch),SmatLS(nch,nch)
        real*8 svalLS(nch),LSconv(nch,nch)
        integer lvalLS(nch),ichLS,nchLS,partLS(nch),exLS(nch),LV,LPAR
    	integer nop,iop(nch),mfails,M1,BAND(2,MXP,MXX),ismin2,ismax2,is2
!	integer nodmin(nch),info,nop,iop(nch),failed,mfails,nodest(nch),
!     X          nbasis(nch),nbas0(nch),nb0,nbmax,NS,IS,nb,ii,ifout
	parameter(imp=60)
	parameter (fracshift=0.25, mfails=5, tr=.true., prba=.false.,
     x             allorder_exch=.true., Z=0d0)
	parameter(mag=1)   ! increase of basis for elastic
	 TM(I) = SECOND() - TIME0

C
!	open(imp,form='formatted')
	call openif(imp)

          if(pralpha) 
     x  write(imp,'(/''JT ='',f5.1,a1/)') jtotal,PSIGN(PARITY+2)
           ichLS = 0
           LSconv(:,:) = 0.0
          do ic=1,MXP
          do ia=1,NEX(ic)
           ismin2 = nint(2*abs(jex(1,ic,ia)-jex(2,ic,ia)))
           ismax2 = nint(2*abs(jex(1,ic,ia)+jex(2,ic,ia)))
           LPAR = PARITY *  SIGN(1, BAND(1,IC,IA)*BAND(2,IC,IA) )

           do is2 = ismin2,ismax2,2
            SV = is2*0.5
            do LV =nint(abs(JTOTAL-SV)),nint(JTOTAL+SV)
             if((-1)**LV  == LPAR) then
               ichLS=ichLS+1
               lvalLS(ichLS) = LV
               svalLS(ichLS) = SV
               partlS(ichLS) = ic
               exlS(ichLS) = ia
               if(pralpha) write(imp,*) '#',ichLS,'ic,ia:',ic,ia,'LS=',LV,SV
               do c=1,nch
                if(PART(c,1) ==ic .and.
     x             EXCIT(c,1)==ia .and.
     x             lval(c)   ==LV) then
                 
                    rmdist(c) = sqrt((2*SV+1.) * (2*jval(c)+1.)) *
     x            racah(lval(c)+Z,JPROJ(c),JTOTAL,JTARG(c),jval(c),SV)
                 I = 1
                  if(abs(KRM)==3) I = (-1)**nint(JTOTAL-SV-lval(c))
!                 if(abs(KRM)==4) I = (-1)**nint(SV-JPROJ(c)-JTARG(c))
                  if(abs(KRM)==4) I = (-1)**nint(JTOTAL-lval(c)
     x                                        -JPROJ(c)-JTARG(c))
                  LSconv(ichLS,c) = rmdist(c)*I
           if(pralpha.and.PCON>0.and.abs(KRM)>=1)
     x            write(imp,14997) ichLS,c,lval(c)+Z,JPROJ(c),JTOTAL,
     x            JTARG(c),jval(c),SV,LSconv(ichLS,c),I
                  endif
                enddo
             endif
            enddo
            enddo
          enddo
          enddo
        nchLS = ichLS
        if(pralpha.and.PCON>0) then
           write(imp,*) '   S-J conversion matrix:'
          do i=1,nchLS
           write(imp,451) i,(LSconv(i,j),j=1,nch)
          enddo
        endif

	if(PCON>=3) call openif(62)
	if(prba) call openif(70)
!	write(KO,*) ' PRALPHA = ',pralpha
	symm = .false.
	if(say.and.final)write(48,*) 'RMATRIXP, REVC: set symm=',symm

	prrm = pralpha.and.PCON>0 ;prkm = prrm ; prtm = pralpha
	prmats = prrm
	prbut = pralpha
	prham = pralpha .and. PCON>4
        gameps = 1d-10 ! MeV spreading of R-matrix poles

!	write(KO,*) ' PRALPHA,prrm= ',pralpha,prrm
	PI = 4d0*atan(1d0)
	CI = (0d0,1d0)
	EL = INITL(1)
	NICH = NCH
	IEX = NCH
	M = N
	NLN = min(NLN1,(N-1)/MR+1)
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

        Brune = btype == 'A'
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
	if(pralpha) write(imp,*) ' Ch. ',C,': beta =',beta(C),' from '
     x                             ,btype,real(bndx),real(rmrad1)
	 enddo


	Rmat(:,:) = 0.
	if(.not.FJSWTCH) then  ! match asymptotic functions to zero here

	if(prrm) then
	   write(imp,*) ' r*R matrix:'
          do 43 i=1,nch
   43      write(imp,45) i,(Rmat(i,j) , j=1,nch)
   45	   format(1x,i3,18f8.3,:,/(4x,18f8.3))
  451	   format(1x,i3,18f10.5,:,/(4x,18f10.5))
	   write(imp,*) ' r*R matrix diag. (before Buttle correction):'
   	   write(imp,45) nch,(Rmat(j,j) , j=1,nch)
	  endif ! prrm
C
	endif   ! nb0>0   (else Rmat=0 + buttle)

!		Add in any search R-matrix term!
	do ip=1,nvars
	 ki = srch_rterm(ip)
	  if(srch_kind(ip)==3) then
	    if(abs(JTOTAL-srch_jtot(ip))<.1
     x         .and.abs(PARITY-srch_par(ip))<.1) then
             damp = 0.0
             do ip2=1,nvars
              if(srch_kind(ip2)==7.and.srch_rterm(ip2)==srch_rterm(ip))
     x          damp = damp + srch_value(ip2)
             enddo
              rm_energy(srch_rterm(ip)) =
     x            cmplx(srch_value(ip),-(damp+srch_damp(ip))*0.5-gameps)
              rm_set(srch_rterm(ip)) = .true.
           if(pralpha)
     x	      write(60,*) 'J/pi=',JTOTAL,PARITY,' so set',srch_rterm(ip)
	      if(pralpha) then
                write(imp,151) ki,rm_energy(ki)
        	write(imp,1501) ip,JTOTAL,srch_jtot(ip),PARITY,srch_par(ip),
     x           srch_rterm(ip),3,srch_value(ip)
               endif
	      endif
          endif
          enddo

	do ip=1,nvars
	  ich = srch_r_ch(ip)
	  if(srch_kind(ip)==4.and.rm_set(max(1,srch_rterm(ip))).and.
     x            ich<=nch) then
	     rmdist(:) = 0.0

	     if(ich>=1) then ! have channel 'ich'
                rmdist(ich) = 1.0
                C1 = ich
	     else ! search for channel with ic,ia,l,j values
             JV = srch_r_jch(ip)
             SV = srch_r_sch(ip)
                 do i=1,nchLS
                 if(partLS(i)==srch_r_ic(ip).and.
     x              exLS(i)==srch_r_ia(ip).and.
     x              lvalLS(i)==srch_r_lch(ip).and.
     x              abs(svalLS(i)-SV)<0.1) then
                  ichLS = i
                  endif
                  enddo
               
	       do c=1,nch
	        if(PART(c,1) ==srch_r_ic(ip) .and. 
     x             EXCIT(c,1)==srch_r_ia(ip) .and.
     x             lval(c)   ==srch_r_lch(ip)) then
	         if     (JV>-0.1) then   ! J specified
                   if(abs(jval(c)-JV)<0.1) then
                       rmdist(c)=1.0
                       C1 = C
                    endif
	         else if(SV>-0.1) then   ! S specified
                  rmdist(c) = LSconv(ichLS,c) 
      if(pralpha)
     x    write(imp,14997) ichLS,c,lval(c)+Z,JPROJ(c),JTOTAL,JTARG(c),
     x            jval(c),SV,rmdist(c)
14997	 format(2i3,' HH*RACAH(',4f5.1,';',2f5.1,') =',f8.4,i4)
                 else
                    write(6,14998) ip,srch_r_jch(ip),srch_r_sch(ip)
14998               format(' Neither JCH nor SCH given for variable',
     x                       i3,2f5.1)
	            if(abs(srch_value(ip))>1d-30) stop
	         endif  ! JV or SV given
                rm_shift(C,srch_rterm(ip)) = srch_r_shift(ip)
                C1 = C
                endif   ! part,excit, L agree
	       enddo    ! c
                rm_shift(C1,srch_rterm(ip)) = srch_r_shift(ip)  ! be sure
	       if(abs(sum(rmdist(:)**2)-1.) > 1e-5) then
	 write(6,15001) ip,srch_rterm(ip),srch_name(ip),srch_r_ic(ip),
     x     srch_r_ia(ip),srch_r_lch(ip),srch_r_jch(ip),srch_r_sch(ip)
15001   format(' Rmatrix channel not found #',i3,' term ',i3,': ',
     x           A10,' ic,ia,l,J,S =',3i3,2f5.1)
	        if(abs(srch_value(ip))>1d-30) stop
	       endif
	     endif
             

             rwamp = srch_value(ip)
             if(.not.srch_rwa(ip)) then
                E1 = rm_energy(srch_rterm(ip)) + ECM(c1,1) - ETOTAL  ! real part
                kp = -E1/coef(C1)  !  C1 is some channel with correct L and E, at least
                ic = part(c1,1)
                wn = sqrt(abs(kp))
                S = ((n-1)*HP(c1)) * wn
                eta1 = ETACNS * MASS(2+1,IC) * MASS(2+2,IC) *
     x                  SQRT(rmass(ic)/abs(E1))
                L = srch_r_lch(ip)
                T = L
                CALL COULFG(S,eta1,0D0,T,fc_s,gc_s,fcp_s,gcp_s,2,0,I,M1)
                pene = S / (fc_s(L+1)**2 + gc_s(L+1)**2)
                rwamp = sqrt(abs(srch_value(ip)/(2*pene)))
                if(srch_value(ip)<0.0) rwamp = -rwamp
                endif

	     do i=1,nch
              rm_vec(i,srch_rterm(ip)) = rm_vec(i,srch_rterm(ip)) 
     x              + rwamp*rmdist(i)
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
               endif
 	write(imp,1501) ip,JTOTAL,srch_jtot(ip),PARITY,srch_par(ip),
     x         srch_rterm(ip),4,srch_value(ip)
            endif ! pralpha
	  endif
	enddo
1501	format(' Var',i4,' for ',2f5.1,2i3,': term #',2i3,' val',f10.5/)

C****************************************************** Find open channels

           do c1=1,nch
            ic = part(c1,1)
           RMA(c1)   = (n-1)*HP(C1)
           ETACO(c1) = ETACNS * MASS(2+1,IC)*MASS(2+2,IC) *
     x                  SQRT(rmass(ic)) ! for Barker
           ENFF(c1)  = ECM(c1,1) - ETOTAL ! for Barker
           enddo

        do ki=1,mterms
        if(rm_set(ki)) then
         if(firstE.and.rm_Brune(ki)) then  ! output to E_Brune(ki),W_Brune(ki)
          call BarkerTransform(ki,rm_set,mterms,rm_energy,rm_vec,nch,
     x       btype,beta,jtotal,PSIGN(PARITY+2),LVAL,ETACO,ENFF,coef,RMA,
     x       pralpha,final,E_Brune,W_Brune)
          endif
        endif
        enddo
C****************************************************** Find open channels

	call nopen(NCH,ECM(1,1),iop,nop)
	if(pralpha.and.nop<nch) then
  	   write(KO,*) 'Only ',nop,' open channels'
	   write(imp,*) 'Only ',nop,' open channels'
	endif

         do I=1,NCH
          SYMR(I) = sqrt(-COEF(I))
          SYMS(I) = (2. * ecm(I,1)/RMASS(PART(I,1)))**0.25d0
         enddo
C**************************************************************
      if (KRM>0.and..not.Brune)  then  !  Usual R > K > S method using boundary conditions beta

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
!	  T = sqrt(RMASS(PART(i,1))/RMASS(PART(j,1)))
!	  T = sqrt(COEF(j)/COEF(i))
      T = sqrt(HP(j)*COEF(j)/(HP(i)*COEF(i)))
	  Rmat(i,j) = Rmat(i,j) +  T* (n-1)*HP(i)*
     x    rm_vec(i,ki)*rm_vec(j,ki)/(rm_energy(ki)-ETOTAL)

	if(pralpha) ! .or.rm_vec(i,ki)*rm_vec(j,ki)>0.)
     x  write(imp,151) ki,rm_energy(ki),i,j, rm_vec(i,ki),rm_vec(j,ki),
     x   ECM(INITL(1),1),Rmat(i,j)
 151	format(' New R term',i4,' at',2f8.4,:,': ',2i3,' with',2f10.5,
     x ' at ',f9.4,' MeV',:,' => R =',2f10.4)

 152	continue
 	endif
 	enddo

        if(pralpha.and.r_added) then
           write(imp,*) ' r*R matrix with additional search terms:'
          do 153 i=1,nch
          rmr(I) = sqrt((n-1)*HP(i))
  153      write(imp,45) i,(Rmat(i,j) , j=1,nch)
           write(imp,*) ' r*R matrix symmetrized'
          do i=1,nch
           write(imp,45) i,(SYMR(i)*Rmat(i,j)/SYMR(j),j=1,nch)
          enddo
           write(imp,*) '   R matrix symmetrized'
          do i=1,nch
           do j=1,nch
           Rmatsym(i,j) =SYMR(i)*Rmat(i,j)/(rmr(i)*rmr(j)*SYMR(j))
           enddo
           write(imp,451) i,(Rmatsym(i,j),j=1,nch)
!    x       (SYMR(i)*Rmat(i,j)/(rmr(i)*rmr(j)*SYMR(j)),j=1,nch)
          enddo
       if(PCON>0) then
           write(imp,*) '   S-J conversion matrix:'
          do i=1,nchLS
           write(imp,451) i,(LSconv(i,j),j=1,nch)
          enddo
           write(imp,*) '   R matrix symmetrized in S basis'
           RmatLS(:,:) = 0.
          do i=1,nchLS
            do j=1,nchLS
            TC = 0.0
            do c=1,nch
            do c2=1,nch
             TC = TC + LSconv(i,c)*Rmatsym(c,c2)*LSconv(j,c2)
            enddo
            enddo
             RmatLS(i,j)=TC
            enddo
           write(imp,451) i,(RmatLS(i,j),j=1,nch)
          enddo
!       write(155,*) ECM(INITL(1),1),abs(rtrace)
        endif ! PCON>0
          endif ! prrm

C****************************************************** K Scattering Matrix
C
C    K =  - [G - RG~]^-1 [F - RF~] where RF~ = r R F' - r beta F, etc

	call GETK(nch,iop,nop,CFMAT,CGMAT,FCWFN,CH,MAXCH,
     X			Rmat,beta,LVAL,Kmat,imp,FCOUL,FCOULP,prmats)
C
	if(prkm) then
	   write(imp,*) ' K matrix:'
           do 230 i=1,nop
  230      write(imp,231) i,(Kmat(i,j) , j=1,nop)
  231	   format(1x,i3,8g14.4,:,/(4x,8g14.4))
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
           T = T + abs(Smatr(IOI,IOP(IOJ)))**2
        enddo
           unitar(IOI) = T
        ENDDO
        ! write(IMP,'('' Unitarity ='',4x,9f16.5)') unitar(1:NOP)


        write(imp,*) ' S matrix symmetrized, without i^L factors'
        DO IOI=1,NOP
          I = IOP(IOI)
        write(IMP,451) I,(SYMS(I)*Smatr(I,IOP(IOJ))/SYMS(IOP(IOJ))
     x                          * (0d0,1d0)**(LVAL(I)-LVAL(IOP(IOJ)))  ! print after removing i^L factors after GETS
     x                  ,IOJ=1,NOP)
        T = 0.0
         do IOJ=1,NOP
           T = T + abs(Smatr(I,IOP(IOJ)))**2
        enddo
           unitar(IOI) = T
   57   format(1x,I3,18f8.4,:,/(4x,18f8.4))
        ENDDO
        ! write(IMP,'('' Unitarity ='',4x,9f16.5)') unitar(1:NOP)
        if(PCON>0) then
        write(imp,*) ' S matrix symmetrized in S basis, no i^L factors'
          do i=1,nchLS
            do j=1,nchLS
            TC = 0.0
            do IOI=1,NOP
            c = IOP(IOI)
            do IOJ=1,NOP
            c2 = IOP(IOJ)
             C6 = SYMS(c)*Smatr(c,c2)/SYMS(c2)
     x                          * (0d0,1d0)**(LVAL(c)-LVAL(c2))  ! print after removing i^L factors after GETS
             TC = TC + LSconv(i,c)*C6           *LSconv(j,c2)
            enddo
            enddo
             SmatLS(i,j)=TC
            enddo
           write(imp,451) i,(SmatLS(i,j),j=1,nchLS)
          enddo
        endif
        endif

C******************************************** Level-Matrix method (KRM<0 or Brune
    	else !  get S from R using level-matrix method
        do ki=1,mterms
         if(pralpha.and.rm_set(ki)) then
        do 58 i=1,nch
        do 58 j=1,nch
      T = sqrt(HP(j)*COEF(j)/(HP(i)*COEF(i)))
	  Rmat(i,j) = Rmat(i,j) +  T* (n-1)*HP(i)*
     x    rm_vec(i,ki)*rm_vec(j,ki)/(rm_energy(ki)-ETOTAL)
        write(imp,151) ki,rm_energy(ki),i,j, rm_vec(i,ki),rm_vec(j,ki),
     x   ECM(INITL(1),1),Rmat(i,j)
 58    continue
        endif
        enddo



	if(pralpha) write(60,*) "Level Matrix Method:"
       	call levelmatrix(N,ETOTAL,NCH,NOP,LVAL,IOP,HP,CH,beta,ECM,COEF,
     x      mterms,rm_set,rm_vec,rm_energy,rm_shift,Brune,
     x      SmatLS,pralpha)

        if(pralpha) 
     x  write(60,*) ' S matrix nonsymmetric, i^L factors: Level Matrix'
        Smatr(:,:) = 0.0
        DO IOI=1,NOP
          I = IOP(IOI)
          DO IOJ=1,NOP
          J = IOP(IOJ)
          Smatr(I,J) = SYMS(J)*SmatLS(I,J)/SYMS(I)  ! back from symmetrized
     x                          * (0d0,1d0)**(LVAL(J)-LVAL(I))  ! and with i^L factors (like GETS)
          enddo
        if(pralpha) write(60,451) I,(Smatr(I,J),J=1,NCH)
          enddo

	endif  ! methods for R > S
C******************************************** S Columns for Cross sections
C
      	FUSL(:) = 0.0
	do 600 jin=1,MINTL
	 EL = INITL(JIN)
         PEL = PART(EL,1)
         EXL = EXCIT(EL,1)

         SMAT(:) = Smatr(:,EL)  ! With GETS or after levelmatrix (0.,1.)**(LVAL(EL)-LVAL(:)) 

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
         I=0
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


       include 'usescatr.f'

      if(CDETR.GT.1) CDETR = CDETR - 1
C                next JIN    :
600   enddo
	return
	END SUBROUTINE RMATRIXP

	subroutine levelmatrix(N,ETOTAL,NCH,NOP,LVAL,IOP,HP,CH,beta,ECM,
     x      COEF,mterms,rm_set,rm_vec,rm_energy,rm_shift,Brune,Smat,pr)
	implicit real*8(A-H,O-Z)
	integer LVAL(NCH),IOP(NOP)
	logical rm_set(mterms),pr,Brune
	Real*8 HP(NCH),rm_vec(NCH,mterms),beta(NCH),
     x   ECM(NCH),COEF(NCH),VEC(mterms,mterms),rm_shift(NCH,mterms)
	COMPLEX*16 CH(2,NCH),rm_energy(mterms),Smat(NCH,NCH)

	real*8 Pen(NCH),rad(NCH),rho(NCH),unitar(NCH),Pen2(NCH)
	complex*16 SiPmB(NCH),HPOH,S,HSP(NCH),
     x             A(mterms,mterms),BB(mterms,NCH)
	integer c,c2,L,M,LT(mterms),IPVT(mterms)

!	write(60,*) N,ETOTAL,NCH,NOP,mterms,'IOP:',IOP
!	write(60,*) 'A',NOP,LVAL,HP
!	write(60,*) 'B',CH
!	write(60,*) 'C',ECM,COEF,'BETA:',BETA

	do c=1,NCH
	rad(c) = (N-1)*HP(c)
	HPOH = CH(2,c)/CH(1,c)
        XK = sqrt(abs(-ECM(c)/COEF(c)))
        rho(c) = XK*rad(c)
	SiPmB(c) = rad(c)*HPOH 

	Pen(c) = rad(c)*XK/abs(CH(1,c))**2
	Pen2(c)= aimag(rad(c)*HPOH)
 	 if(ECM(c)<0d0) Pen(c) = 0d0
! 	HSP(c) = sqrt(conjg(CH(1,c))/CH(1,c))
	HSP(c) = abs(CH(1,c))/CH(1,c)
	if(pr) write(60,10) c,rad(c),rho(c),Pen(c),Pen2(c),SiPmB(c)
10	format(' c,R,rho,P,phi; S:',i4,4f10.5,';',2f10.5)	
	if(.not.Brune) SiPmB(c) = SiPmB(c) - rad(c)*beta(c)
	enddo

	if(pr) write(60,*) 'Level set:',rm_set
	nterms = 0
	do ki=1,mterms
	 if(rm_set(ki)) then
	   nterms = nterms + 1
	   LT(nterms) = ki
	 endif
	enddo
	if(pr) then
	write(60,20) (LT(L),L=1,nterms)
20	format('Level list:',60I3)
        do L=1,nterms
        write(80,25) LT(L),(rm_shift(c,LT(L)),c=1,NCH)
25      format(' Term:',i3,'  S(:)=',20f9.5)
        enddo
	endif

	do L=1,nterms
	do M=1,nterms
	 S = 0d0
	  do c=1,NCH
	    S = S - rm_vec(c,LT(L)) * SiPmB(c) * rm_vec(c,LT(M))
            if(Brune) then
             if(L /= M) then
               S = S + rm_vec(c,LT(L)) * rm_vec(c,LT(M)) *
     x          (rm_shift(c,LT(L))* (ETOTAL - rm_energy(LT(M))) - 
     x           rm_shift(c,LT(M))* (ETOTAL - rm_energy(LT(L))) ) /
     x            (rm_energy(LT(L)) - rm_energy(LT(M)))
              else  ! L==M
               S = S + rm_vec(c,LT(L))**2 * rm_shift(c,LT(L))
              endif
            endif
	  enddo
	  if(L==M) S=S + rm_energy(LT(L)) - ETOTAL
	A(L,M) = S
	enddo !M
	do c=1,NCH
	 VEC(L,c) =  rm_vec(c,LT(L))
	enddo !c
	enddo !L
	if(pr) write(60,*) 'A matrix'
        do 230 i=1,nterms
  230      if(pr) write(60,231) i,(A(i,j) , j=1,nterms)
  231      format(1x,i3,8g14.4,:,/(4x,8g14.4))

	if(nterms>0) then
	CALL ZGETRF (nterms, nterms,A,mterms,IPVT,IER)
	if(IER/=0) write(60,*) 'ZGETRF:',IER
        if(pr) write(60,*) 'IPVT = ',IPVT(1:nterms)
	
	do io=1,NOP
           c = iop(io)
	BB(1:nterms,io) = VEC(1:nterms,c) *sqrt(2d0*Pen(c))
	enddo
 	CALL ZGETRS('N', nterms,NOP,A,mterms,IPVT,BB,mterms,IER)
	if(IER/=0) write(60,*) 'ZGETRS:',IER
 	endif

	Smat(1:NCH,1:NCH) = 0.0

	do io=1,NOP
	  c = IOP(io)
	  Smat(c,c) = 1d0
	do io2=1,NOP
	  c2 = IOP(io2)

	  S = 0d0
	  do L=1,nterms
	   S = S + VEC(L,c2)*BB(L,io)
	  enddo
           S = (0d0,1d0) * sqrt(2d0*Pen(c2)) * S
	  Smat(c2,c)= HSP(c2) * (Smat(c2,c)+  S) * HSP(c)
	  
	enddo ! io2
	enddo ! io
	  
        if(pr) then
        write(60,*) ' S matrix symmetric, NO i^L factors: Level Matrix'
        DO c2=1,NCH
        write(60,451) c2,(Smat(c2,c),c=1,NCH)
        T = 0.0
         do c=1,NCH
           T = T + abs(Smat(c,c2))**2
         enddo
           unitar(c2) = T
   57   format(1x,I3,18f8.4,:,/(4x,18f8.4))
  451	   format(1x,i3,18f10.5,:,/(4x,18f10.5))
        ENDDO		
        endif
	
!       write(60,'('' Unitarity ='',/,4x,9f16.5)') unitar
	return
	end

!       Smat(c2,:) = Smat(c2,:) * (0d0,1d0)**(LVAL(:)-LVAL(c2))

