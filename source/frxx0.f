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
************************************************************************
	subroutine freadf(ki,ko,koe,TMP,llmax,jjbord,jjump,readvar)
	use parameters
	use factorials
	use drier
	use kcom
	use trace
	use fresco1, dllmax=>lmax, djjump=>jump,djjbord=>jbord
	use parallel
      	use searchpar, only: nvars,boldplot
	use gails, only: cfalloc
        use io, only: grace,pre_is,post_is
!$      use omp_lib
#ifdef MPI
        use mpi
#ifdef ALTIXMPI
	use mpi_struct
#endif /* ALTIXMPI */
#endif /* MPI */
        implicit real*8(a-h,o-z)
	parameter (mp=200,ma=2000,mmxp=20,maxex=5000)
        real jmax,masst,massp,zp,zt,jex(2,maxex,mmxp),jbord(7)
	real et,ep,tp,tt,qval,be,ampl,er
        real def(0:8),str,j,triton(8),jsp(ma),jn,jnp,jp,kkp,jt,kkt
        real jcom,jcomp,jcor,jcorp,j1,sn,kcor,kcorp
        integer qc,qcm,qcmin,qcmax,la,lam,lap,lambda
        integer parity(2,maxex,mmxp),par
	real coef(mp),masses(2,mmxp),p(0:12),mnep(0:12),
     x		defp(0:12),deft(0:12),mnet(0:12),autowid
	real ap,at,rc,ac,p0,p1,p2,p3,p4,p5,p6,p7,v,r0,a,av,w,wr0,aw,wa,
     X	     r0w,rw,rv,vr0,vso,rso,rso0,aso,awd,wda,vsoi,rsoi,asoi,
     X	     wd,wdr0,wdr,vd,vdr,vda
        real x,y,z
        real w3j,w6j,w9j,ep1
        real pcoup0,pcoup1,pcoup2,pcoup3,pcoup4,pcoup5,pcoup6
	complex qscale(0:11)
!	equivalence (def,p),(p0,p(0)),(p1,p(1)),(p2,p(2)),(p3,p(3)),
!    X	 (p4,p(4)),(p5,p(5)),(p6,p(6)),(p7,p(7)),
!    X   (at,p1),(ap,p2),(rc,p3),(ac,p4),(v,p1,vso),(r0,vr0,rv,rso,p2),
!    X   (a,av,aso,p3),(w,p4),(wr0,r0w,rw,p5),(aw,wa,p6)
	equivalence (def,p),(p0,p(0)),(p1,p(1)),(p2,p(2)),(p3,p(3)),
     X	 (p4,p(4)),(p5,p(5)),(p6,p(6)),(p7,p(7)),
     X   (r0,vr0,rv),(a,av),(wr0,r0w,rw),(aw,wa),(rso,rso0),(awd,wda),
     X   (wdr,wdr0)
	equivalence (bandp,ptyp),(bandt,ptyt)
	!real, allocatable:: afrac(:,:,:,:,:,:)
	logical*1, allocatable:: afrac(:,:,:,:,:,:),usedcc(:,:,:)
        logical*1 alogic
        character*120 line,input_file
	character*70 TMP,MASFIL
        character*20 CCSPLIT
        character*14 version
	character*10 hscale(0:2),planes(0:3),buttles(0:4),methods(0:6)
	data hscale/'none      ','projectile','target    '/
	data planes/'none      ','exit      ','entrance  ','entr.&exit'/
	data methods/'default   ','ERWIN     ','ERWINCC   ',
     x  'ERWINRKC  ','ERWINRK   ','ERWINRK+NT','NL-SEARCH '/
	character*8 namet,namep,name,CHME
	character*4 word
	character*30 comp
	character*200 datafile
	data buttles/'complex+sh','real+shift',
     x               'complex-sh','real-shift','none    '/
	character*2 ex(2)
	character cpwf,ch1,citt,cset,jl,ccrealc,gscouplc
	integer cp,ptyp,bandp,bandt,ptyt,copyp,copyt,tnt(4,mp),kp,
     x	   typei,shapei,type,shape,qq,lsp(ma),knb(ma),cpot,
     X     cpots(maxex,mmxp),nexx(mmxp),kindsp(ma),kna(ma),
     X	   usedlow(ma),usedhigh(ma),knused(ma),nchmax(mmxp),jump(6),
     X     infam,outfam,jjump(7,3),usedcore(2,ma),knib(ma),knia(ma),
     X     maxcoup(3),kii,readstates,pcomps(0:maxex),infile,inia(ma)
      integer SYCTXT,NPROW,NPCOL,MYROW,MYCOL
      integer, allocatable:: usermap(:,:)
        real*8 jjbord(8),expand(11),expandi(11)
	logical bigj,used(2,ma),fail3,frac,uu,nml,keep,itt,pwf,nosub,
     X	fexch,ignore,trnl,adjq,fail,given,parall,haso(ma),foldso,cdc,
     x  complexbins,pr,gamma,readvar
	namelist/fresco/hcm,rmatch,rintp,hnl,rnl,centre,hnn,
     X   rnn,rmin,rsp, rasym,accrcy,switch,ajswtch, sinjmax,plane,
!@     X   gailitis,gailacc,
     X   jtmin,jtmax,absend,dry,rela,nearfa,jump,jbord,pset,jset,jleast,
     x   kqmax,pp,thmin,thmax,thinc,koords,cutl,cutr,cutc,complexbins,
     x   ips,it0,iter,fatal,iblock,pade,psiren,iso,nnu,maxl,minl,mtmin,
     X   epc,erange,dk, nosol,nrbases,nrbmin,pralpha,pcon,rmatr,bndx,
     X   meigs,buttle,weak,llmax,expand,hktarg,mpihelp,eigens,nlagcc,
     X   chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,ccreal,
     x   lampl,veff,kfus,wdisk,bpm,melfil,cdcc,nfus, nparameters,
     x	 smallchan,smallcoup,sumform,maxcoup,ompform,initwf,pluto,
     X	 inh,TMP,MASFIL,unitmass,finec,pel,exl,lab,lin,lex,elab,nlab,
     x   vsearch,echan,enodes,gscouplonly,ccbins,sock,btype,boldplot,
     x   KRM, dist0,tcfile,dgam,eobs,grace,elpmax,relref,
     x   tcfilename,numnode, sumform  ! not used: only for reading old inputs

	namelist/partition/namep,massp,zp,nex,pwf,namet,masst,zt,qval,
     x	                   readstates,prmax,lpmax,mixpot
	namelist/states/ jp,  copyp,ptyp,bandp,ep,kkp,tp, cpot,
     X		jt,copyt,ptyt,bandt,et,kkt,tt, fexch,ignore,infam,outfam
	namelist/pot/ kp,type,itt,nosub,shape,def,mnep,mnet,ap,at,rc,ac,
     X		 p,p0,p1,p2,p3,p4,p5,p6,p7,v,r0,rv,vr0,a,av,
     X		 w,wr0,rw,aw,wa,r0w,vso,rso,rso0,aso,vsoi,rsoi,asoi,
     X           wd,wdr,wda,wdr0,awd,defp,deft,vd,vdr,vda,
     x           lshape,jl,xlvary,alvary,citt,datafile
	namelist/potl/ kpi,typei,itt,nosub,shapei,p,lshapei,jl,
     x 			xlvary,alvary,datafile
	namelist/step/ib,ia,k,str
	namelist/overlap/ kn1,kn2,ic1,ic2,in,kind,ch1,nn,l,lmax,sn,
     &         ia,j,ib,kbpot,krpot,be,isc,ipc,nfl,nam,ampl,keep,
     & 	       dm,nk,er,e,ppower,rsmin,rsmax,nlag,phase,autowid
	namelist/dalitz/ triton
	namelist/twont/ tnt,coef
	namelist/coupling/icto,icfrom,kind,ip1,ip2,ip3,ip4,ip5,
     X			  p1,p2,jmax,rmax,kfrag,kcore,nforms,infile
	namelist/inel/ ib,ia,k,no,kp,a
	namelist/cfp/ in,ib,ia,kn,a,keep
	namelist/scale/ qscale

      frac(x) = abs(x-nint(x)).gt.1e-5
      fail3(x,y,z) = frac(x+y+z) .or. x.gt.y+z .or. x.lt.abs(y-z)

#ifdef ALTIXMPI
	if(MPIC.and..not.STDINALL)then
          include 'nlstruct-mpi2.f'
        endif
#endif /* ALTIXMPI */
        
        eps = 1d-10
	z = 0d0
	adjq = .false.
	fail = .false.; cset = ' '
	given = .false.; 
	symm = .true.;   nn2wf = .false.
        ccbins = .false.; complexbins = .false.
        sumccbins = .false.; sumkql = .false.; eigens = 0
	maxcoup(:)=0; sock(:,:)=1.0
	nparameters = 0
         expand(1) = 1.2; expand(2) = 1.00; expand(3) = 0; 
	 expand(4) = 4.0; expand(5) = 0.4; expand(6) = 2.   
	 expand(7) = 1.0; expand(8) = 1.0; expand(9) = 1.0
	 expand(10)= 1.0; expand(11)= 1.0
 	 expandi(:) = expand(:)
        maxit = 0  ; nchbinmax=1
	knused(:) = 0
	maxqp = 0 ; maxqp1 = 0 ; maxqp2 = 0 
	ko3 = 301 + 2*iame
	call compiler(comp)
!	write(6,*) ' I am ',iame
    	pr = iame==0
        n_svn = 0
        call version_number(version)

    	if(pr) write(koe,1002) version,comp
 1002 	format(' FRESCOX - version ',a,
     X      ': Coupled Reaction Channels             on ',a30/)
        call GET_COMMAND_ARGUMENT(1,input_file,I)
!        write(6,*) 'Command argument <'//trim(input_file)//'>',I
	uu = .false.
 	inquire(file=input_file,exist=uu)
	if(uu) then
	      if(pr)
     X        write(koe,*) ' Using as input filE <'//
     X          trim(input_file)//'>  NOT stdin'
  	      open(ki,file=input_file,status='old')
	      rewind ki
              else
              if(MPIC) write(6,*) 
     x    'For MPI, recommend use command line argument for input file'
	      endif
 1003 	inquire (UNIT=ki,NAME=input_file) 
!       write(6,*) 'Reading from <'//trim(input_file)//'>'
        if(iame==0.or.STDINALL)read(ki,1005) headng
!	 write(6,'(a)') 'Heading: <'//headng//'>'
      	if(iame==0.or.STDINALL)read(ki,1005) line
!	 write(6,'(a)') 'Line   : <'//line//'>'
 1005 	format(a120)
#ifdef ALTIXMPI
        if(MPIC.and..not.STDINALL)
     !            call MPI_BCAST(MPI_BOTTOM,1,headng_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
        nml = line(1:1)=='N'.or.line(1:1)=='n'
        cdc = line(1:1)=='C'.or.line(1:1)=='c'
	trnl = line(1:1)=='N' .or. line(1:1)=='C'
#ifdef ALTIXMPI
        if(MPIC.and..not.STDINALL.and..not.nml)
     >   write(koe,*)'MPI WITH ONLY ROOT STDIN REQUIRES NAMELIST INPUT'
#endif /* ALTIXMPI */
! 	write(koe,*) 'nml,cdc =',nml,cdc,' as read ',trim(line)
	if(cdc) then
	   if(pr) write(koe,*) 'Assuming CDC input.'
              I = ichar('0')
       CHME = char(I+mod(iame/100,10))//char(I+mod(iame/10,10))//
     X       char(I+mod(iame,10))//'.301'
        if(iame==0) CHME = '301'
           open(ko3,recl=125,form='formatted',
     x          delim='apostrophe',file='fort.'//trim(CHME))

	   call cdcin(ki,ko3,headng,trnl,pr)
	   rewind ko3
	   if(uu) close(ki)
	   uu = .false.
	   ki = ko3
	   go to 1003
	endif
	if(nml.and.pr)
     >  write(koe,*) ' Using NAMELIST input'
!@	call defaults(TMP,MASFIL,gailitis)
	call defaults(TMP,MASFIL)
	jbord(:) = 0.; jump(:) = 0
	if(nml) then
	 ios=0
         if(trnl) then
	 if(iame==0.or.STDINALL)
     >      read(ki,nml=fresco,IOSTAT=ios,end=2,err=2)
	 else
	 if(iame==0.or.STDINALL)read(ki,nml=fresco)
	 endif
   2     if ( ios .ne. 0 ) then
            write (koe,*) ' Input namelist FRESCO read error: ',ios
            write (koe,fresco)
            stop
         endif
#ifdef ALTIXMPI
         if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nlfresco_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
         endif
      if(pr)write(koe,1009) headng
 1009 format(/1x,a120/) 

	fcwfn=.false.
	if(.not.nml) then
		hcm=0;rmatch=0;rintp=0;hnl=0;
		rnl=0;centre=0;hnn=0;rnn=0;rmin=0;rsp=0
        read(line,1010,iostat=ios,err=10101) hcm,rmatch,rintp,hnl,
     x                            rnl,centre,hnn,rnn,rmin,rsp
	go to 10109
 1010 	format(10f8.3)
10101   write(koe,*) ' **** Error ',ios,' reading line: hcm,rmatch,..' 
	write(koe,1005) line
	write(koe,*) ' **** Data read so far are:'
        write(koe,1010) hcm,rmatch,rintp,hnl,rnl,centre,hnn,rnn,rmin,rsp
        read(line,1010) hcm,rmatch,rintp,hnl,
     x                            rnl,centre,hnn,rnn,rmin,rsp
	stop
	endif
10109 	if(pr)
     >  write(koe,1010) hcm,rmatch,rintp,hnl,rnl,centre,hnn,rnn,rmin,rsp
      	if (rmatch<0.and..not.nml) then
           rasym=0;accrcy=0;switch=0;ajswtch=0;sinjmax=0
      	   read(ki,1005) line
           read(line,1011) rasym,accrcy,switch,ajswtch,sinjmax
!@     x                   ,gailitis,gailacc
        endif
      	if (abs(rasym)>eps) then
         if(pr)
     >   write(koe,1011) rasym,accrcy,switch,ajswtch,sinjmax
!@     x                   ,gailitis,gailacc
 1011 	format(f8.2,f8.5,f8.2,3f8.2,e8.1)
        if(accrcy.LE.0.0) accrcy = 0.01
        if(switch.LE.1.0) switch = 1000.
        if(ajswtch.LE.1.0) ajswtch = 10.
!@	ngail = nint(gailitis)
	ngail = 0
!@	if(ngail>0.and.gailacc<1e-20) gailacc = accrcy**2
	endif
	   m = nint(abs(rmatch)/hcm)
	   m = (m/2)*2 + 1    ! m is odd, so rmatr is even multiple of hcm
	   mrm = m
	      rmatch = (m-1)*hcm
	     n = m + 1
	      md = m - 6
	   mr= nint(rintp/hcm)
	      if(mr.le.0) mr = 4
	      rintp = mr * hcm
	      RIN = 1.0/RINTP
            if(rsp.le..01) rsp = abs(rmatch)
	    mint = nint(rsp/hcm) + 1
	      rsp = (mint-1)/mr*rintp ! + rintp  ! extra step
	    mint = nint(rsp/hcm) + 1
	    nln = (n-1)/mr + 1
	      if(hnl.lt.1e-5) hnl = hcm
	   mlt = nint( hnl/hcm )
	   nlt = nint( hcm/hnl )
			   hnl = hcm
	      if(nlt.gt.1) hnl = hcm/nlt
	      if(mlt.gt.1) hnl = hcm*mlt
	      nlt = max(nlt,1)
	      mlt = max(mlt,1)
	   hsp = hcm * mlt
	   nlo = nint(rnl/hsp) + 1
	      if(nlo.lt.4) nlo = 4
           RNL = (NLO - 1)*HSP
          NLM = (NLO-1)*NLT + 1
          MLM = (NLO-1)*MLT + 1
          NLCN=  CENTRE/HSP
           CENTRE = NLCN*HSP
          NLC = NLO/2 - NLCN
	   if(hnn.lt.eps) hnn = rintp
	   nnn = (rnn-rmin) / hnn
     	   rnn = nnn * hnn + rmin
      if(pr)
     >WRITE(koe,1020) N,HCM,RMATCH,RINTP,NLM,HNL,RNL,CENTRE,NNN,HNN,RNN,
     X               RMIN,RSP
 1020 FORMAT(/' Centre-of-mass Range is',I7,' *',F7.4,' fm., Maximum at'
     X, F7.2,' fm., Interpolating NL forms every',F6.2,' fm.'/
     X ' Non - locality width is',I7,' *',F7.4,' fm., Maximum of',
     X  F7.2,' fm., Centred at',F6.2,' fm.'/
     X ' 2-Nucleon Separation of',I7,' *',F7.4,' fm., Maximum of',
     X  F7.2,' fm., Minimum at',F6.2,' fm.'/
     X ' Maximum single particle bins of',F8.4,' fm.')
      IF(fcwfn) THEN
         IF(RASYM.gt.0.and.pr)
     *   WRITE(koe,1012) RASYM,ACCRCY,SWITCH,AJSWTCH,sinjmax
!@     * 			,ngail,gailacc
         IF(RASYM.lt.0.and.pr)
     *    WRITE(koe,1013) RASYM,ACCRCY,SWITCH,AJSWTCH,sinjmax
!@     * 			,ngail,gailacc
 1012 FORMAT(/' CRCWFN  Rasymptotic = ',F8.2,' fm,  Accur = ',F8.6,
     *  ',   Switch = ',F8.2,' fm.,   L switch = ',F8.2,
     *  ',   SinJmax =',f8.2) !@',  Gailitis =',i3,' gailacc =',1p,e9.1)
 1013 FORMAT(/' CRCWFN  -Theta-min = ',F8.2,' deg,  Accur = ',F8.6,
     *  ',   Switch = ',F8.2,' fm.,   L switch = ',F8.2,
     *  ',   SinJmax =',f8.2) !@',  Gailitis =',i3,' gailacc =',1p,e9.1)
       ENDIF
      	mintm2 = mint -2
	if(pr)write(koe,*) '              M,Mint = ',m,mint

	maxn = n
	maxm = max(n,mint)
	maxnln = nln
	maxnlo = nlo
	maxnnn = nnn

      	if(.not.nml) then
           read(ki,1030) jtmin,jtmax,absend,dry,cset,rela,
     x         nearfa,(jump(i),jbord(i),i=1,6),jleast
	   endif
	JJUMP(2:7,1) = jump(1:6)
	JJBORD(2:7) = jbord(1:6)
        JJUMP(1,1) = 1
        JJBORD(1)= MAX(0D0,JTMIN)
        JJBORD(8)= 0
	bigj = abs(jtmax).gt.999d0 .or. maxval(jjbord).gt.999d0
!	if(.not.bigj) then
!      	write(koe,1030) jtmin,jtmax,absend,dry,cset,rela,nearfa,
!     x         (jump(i,1),jbord(i),i=2,7)
!	else
!      	write(koe,1031) jtmin,nint(jtmax),absend,dry,cset,rela,nearfa,
!     x         (jump(i,1),nint(jbord(i)),i=2,7)
!	endif
 1030 	format(2f4.0,  f8.4,l2,1x,a1,a2,i2,6(i4,f4.0),f4.0)
!1031 	format(f4.0,i4,f8.4,l2,1x,a1,a2,i2,10i4)
      I = MAX(1, MIN(ABS(NEARFA),13))
      NEARFA = ISIGN(I,NEARFA)
	if(.not.nml) then
         jset = 0
         pset = 0
         IF(CSET.EQ.'T') THEN
           JSET = 1
         ELSE IF (LGE(CSET,'0') .AND. LLE(CSET,'9')) THEN
            JSET = ICHAR(CSET) - ICHAR('0')
         ELSE IF(CSET.EQ.'P') THEN  ! positive parity
           PSET = 1
         ELSE IF(CSET.EQ.'M') THEN  ! negative parity
           PSET = -1
         ELSE IF(CSET.eq.'N') THEN  ! natural parity
           PSET = 2
         ENDIF
	endif
      relref = max(relref,1)


      CCSPLIT = ' '
      if (abs(NEARFA)>=10) CCSPLIT = '(splitting Coulomb)'
      if(pr) WRITE(koe,1040) JTMIN,JTMAX,jleast,ABSEND,
     x                       DRY,PSET,JSET,RELA,NEARFA,CCSPLIT
 1040 FORMAT(//' Range of total J is',F8.1,' <= J <=',F8.1,
     X ' (at least',f8.1,') and Absorbtion => ',F10.4,' mb.',/,
     X  ' Dry Run =', L2,', CC set limits =',2I3,', ',
     X  'Relativistic kinematics = ',A2,
     x ', Both/Far/Near Analyses =',I3,1x,a)
!      IF(index(' abAB',RELA)==0) WRITE(*,*) 
!     x  'RELATIVISTIC KINEMATICS '''//rela//''' NOT IMPLEMENTED!!'
      IF(index(RELA,'h')/=0) WRITE(*,*) 
     x     ' Relativistic reference channel:',relref
      NJ = 0
      DO K=7,1,-1
       IF(JJUMP(K,1)/=0) then
	NJ = K
	GO TO 9
	endif
	ENDDO
9     JJBORD(NJ+1) = JTMAX
       DO 10 K=1,NJ,+1
      IF(JJUMP(K,1).EQ.0) GO TO 10
       IF(JJBORD(K+1).GT. JTMAX) JJBORD(K+1) = JTMAX
       IF(JJBORD(K)  .LT. JTMIN) JJBORD(K)   =  JTMIN
       IF((JJBORD(K+1)-JJBORD(K))/JJUMP(K,1) .LT. 3) JJUMP(K,1) = 1
       IF(JJUMP(K,1).gt.0) THEN
         I = MOD(NINT(JJBORD(K+1)-JJBORD(K)),JJUMP(K,1))
         IF(I.ne.0) THEN
	   I = JJUMP(K,1) - I
           if(pr)
     >     WRITE(koe,1048) JJBORD(K),JJBORD(K+1),JJUMP(K,1),I
 1048      FORMAT(/' *** From J =',F6.1,' to ',F6.1,' is NOT a ',
     X     ' multiple of ',I4,'!!. Upper limit increased by',I3)
           JJBORD(K+1) = JJBORD(K+1) + I
         ENDIF
        ENDIF
       IF(JJBORD(K+1)-JJBORD(K).LT.1 .AND. NJ.GT.1) JJUMP(K,1) = 0
   10 CONTINUE
       IF(NJ.GT.1.and.pr)
     > WRITE(koe,1050) (JJUMP(K,1),JJBORD(K),K=1,NJ)
 1050 FORMAT(/' Jump J-total steps',4('  of',I5,' after',F8.1:','),/
     X             19x            ,3('  of',I5,' after',F8.1:','))
      JTMAX = JJBORD(NJ+1)

      if(.not.nml) then
        read(ki,1055) kqmax,pp,thmin,thmax,thinc,koords,
     X				  cutl,cutr,cutc,dist0,dgam
        if (dist0==0.0) dist0 = -1.0
       endif
 1055 format(2i1,f6.2,f8.3,f6.3,i2,4f8.3,i2)
!      write(koe,1055) kqmax,pp,thmin,thmax,thinc,koords,cutl,cutr,cutc
      !PP = mod(PP,4)
      IF(abs(CUTL).LT.EPS) CUTL = -1.6
      ICUTC = MAX(INT(CUTC/HCM)+1,1)
      CUTC = (ICUTC-1)*HCM
      WORD = '    '
      if(CUTR.lt.0.) WORD = 'Turn'
      if(pr)
     >WRITE(koe,1060) KQMAX,PP,PPKIND(PP),THMIN,THMAX,THINC,dgam,grace,
     X KOORDS,'(',COORDS(1),(',',COORDS(I+1),I=1,MIN(KOORDS,3)),')'
      if(pr.and.dist0>-.1) write(koe,1064) dist0
      if(pr)
     >WRITE(koe,1065) CUTL,WORD,CUTR,CUTC
 1060 FORMAT(/' Cross Sections (and up to T',I1,
     x ' for ',I1,'=',A10,') for Theta from',F6.1,' to',F6.1,
     x ' in steps of',F5.1,' degrees, DGAM=',i1,', grace=',l1,
     X ',  Coordinates =',I2,' ',5(A1,A4))
 1064 FORMAT(' Make distributions, with Coulomb from ',f7.1,' deg')
 1065 FORMAT(/' Lower Radial Cutoff = maximum of',F6.2,'*L*h & ',A4,
     X F5.1,' fm.,  Lower Cutoff for Couplings =',F6.1,'  fm.'/)
       if(.not.grace) then
            pre_is = 'legend string '
            post_is= ''
        else
            pre_is = '             s'
            post_is= ' legend '
        endif

      if(DGAM>0.and. (PP/=2 .and. PP/=3.or.kqmax<2)) then
 	write(koe,*) 'For decay gammas need PP=2 or 3, and KQMAX>>0'
    	stop
 	endif

        if(readvar) call readvars(nparameters,koe,rterms)

	dk = 0d0
      if(.not.nml) read(ki,1070)
     x     ips,it0,iter,iblock,pade,iso,nnu,maxl,minl,mtmin,
     X     epc,erange,dk,inh,plane,smallchan,smallcoup,initwf,mpihelp,
     x     ccrealc  !,gscouplc
 1070 format(f6.4,i2,i4,2i2,a1,i3,2i4,i2,f6.2,2f8.4,2i2,2f8.6,i4,i2,2a1,
     x      f8.4)
!     write(koe,1070) ips,it0,iter,
!    x     iblock,pade,iso,nnu,maxl,minl,mtmin,epc,erange,dk,inh,plane
!    X     ,smallchan,smallcoup,initwf,mpihelp,ccreal !,gscouplonly
      if(.not.nml) fatal = ITER.GT.0
      ITER = ABS(ITER)
 	line =  ' Convergence failures not fatal'
	if(fatal) line = 'Convergence failures are fatal'
	if(ITER<=3) line = ' '
      if(.not.nml) nosol = ITER < IT0
!      IT0 = MAX(IT0,1)
      ITMIN = MIN(IT0,ITER)
      ITER = MAX(IT0,ITER)
!      IF(IPS.NE.0 .AND. ITER.GE.4) ITMIN = MAX(ITMIN,1)
      IF(ABS(PADE).GE.1) THEN
       maxit = iter
       if(.not.nml) psiren = PADE.GT.0
       PADE = ABS(PADE)
      ENDIF
      NSOL = ITER + 1
      IF(NNU.LT.18)  NNU = 18
         NNU = ((NNU+5)/6)*6
	  maxnnu = nnu
      ISOCEN = 0
      IF(ISO.eq.'A'.or.ISO.eq.'J') ISOCEN = 1
      IF(ISO.eq.'B'.or.ISO.eq.'L') ISOCEN = 2
      IF(MTMIN.EQ.0) MTMIN = 6
      IF(EPC.lt.1e-5) EPC = (30./NNU)**2
      IF(EPC.LT..001) EPC = 0.001
      if(EPC.GT.5.) EPC = 5.0
	maxpw = jtmax+20
	if(llmax>=0) maxpw = llmax
	lmax1 = max(nint(jtmax)+20,MAXL+1)
!	  if(maxpw<lmax1) write(koe,*) '  Maximum L value =',maxpw
      IF(MAXL.EQ.0) MAXL = JTMAX + Max(6,mtmin) + 4
      IF(ABS(JTMIN).GT.6..AND.MINL.LT.0) MINL = ABS(JTMIN) - 6.
      IF(MINL.LT.0) MINL = 0
C     IF(JTMAX.LE.-JTMIN+.1) MINL = MAXL + 1
C
	if(.not.allocated(fact))then
	 mfact = 2*lmax1+200
!	   write(370,*) ' lmax1 =',lmax1,' so  mfact = ',mfact
	 allocate (fact(mfact))
	 fact(1) = 0d0
         CALL LOGFAC(MFACT)
	endif

	if(NOSOL) then
       if(pr)
     > WRITE(koe,*) ' Coupled equations only set up, NOT solved'
	else
       if(pr)
     > WRITE(koe,1074) ITMIN,ITER,IPS,trim(line),
     X IBLOCK,PADE,ISO,ISOCEN,NOSOL,ccreal,initwf !,gscouplonly 
 1074 FORMAT(
     X/' Iterate Couplings between',I4,' and',I4,' times, to',F7.3,
     X' % if sooner.  ',A40
     X /,'  Block solved exactly =',I4,' chs., with Pade =',I2,
     X   ' & Isocen = ',A1,'=',I1,', NOSOL = ',L1,', CCREAL =',L2,
     X   ', initwf =',i3)
!    X   ', GScouplOnly =',L2)
       if(abs(smallchan+smallcoup)>0.and.pr)
     > WRITE(koe,1075) smallchan,smallcoup
 1075 FORMAT('    Small channels are ',1p,e9.2,' and small couplings ',
     X       'are ',e9.2,' of unitarity')
        if(initwf/=0) write(koe,10751) initwf
10751  format('   READ FIXED CHANNEL WFS from file #',i4)
	endif
	if(.not.nml) ccreal = ccrealc == 'T'
	if(.not.nml) gscouplonly = gscouplc  == 'T'
	gscouplonly = .false.  ! not accurate to use gscouplonly=T !

       if(iblock<0.and.iblock>-9.and..not.nml) then
	 read(ki,1076) eigens,nrbases,nrbmin,buttle,pralpha,pcon,
     x			meigs,rmatr,bndx,weak,btype,eobs
 1076    format(i1,i3,i4,i1,l1,I2,i2,f6.2,3f8.4,1x,A1,L1)
	 write(48,1076) eigens,nrbases,nrbmin,buttle,pralpha,pcon,
     x			meigs,rmatr,bndx,weak,btype,eobs
	 if(eigens>0) then
  	  read(ki,10761) vsearch,echan,enodes
10761     format(3i4,l2,f6.0)
  	  write(48,10761) vsearch,echan,enodes
	  call flush(48)
	  endif
	endif
	   if(rmatr<eps) rmatr = rmatch
	   rmatr = min(rmatch,rmatr)
	   mrm = nint(abs(rmatr)/hcm)
	   mrm = (mrm/2)*2 + 1    ! mrm is odd, so rmatr is even multiple of hcm
!			        for Simpsons rule in Rmatrix
	   rmatr = (mrm-1)*hcm
!why?	   m = mrm
!why?	   md = mrm - 6
	   buttle = max(0,min(buttle,4))
	   rterms = rterms .or. nrbases/=0
	   if(nrbases<0) then
	 	nlagcc = -nrbases
		nrbases= 0
	      endif
       if(btype=='A') KRM=min(KRM,-1)  ! allow KRM=1 later, when Barker transform coded
       if(rterms) then
	   if(nrbmin==0) nrbmin = nrbases
	 if(pr.and.nrbases>0)
     >   write(koe,1077) nrbases,nrbmin,buttle,buttles(buttle),
     x	pralpha,pcon,meigs,rmatr,bndx,weak,eigens,btype,eobs,KRM
 1077    format(/'  R-MATRIX SOLUTIONS with <=',i3,' basis',
     x  ' states/channel , and at least',i3,' states.',/
     X  '   buttle=',i1,'=',a10,', pralpha=',L1,', pcon=',i2,
     X     ', meigs =',i4,', Rmatr =',f8.4,', bndx =',2f8.4,
     x     ', weak =',f8.5,', eigens =',i2,', btype =',A1,1x,L1,
     x     ', KRM =',i2)
	 if(pr.and.nlagcc>0) write(koe,10772)nlagcc,pralpha,pcon,
     x                meigs,rmatr,bndx,weak,eigens,btype
10772    format(/'  LAGRANGE MESH R-matrix SOLUTIONS with ',i3,
     x           ' basis states/channel '/,
     X  '   pralpha=',L1,', pcon=',i2,', meigs =',i4,', Rmatr =',f8.4,
     x  ', bndx =',2f8.4,', weak =',f8.5,', eigens =',i2,', btype =',A1)
	 if(eigens>0) then
  	  write(koe,10775) vsearch,echan,enodes
10775    format('   Bound states: search on potential form at ',i4,
     x    ' so ch. ',i3,' has ',i3,' nodes.')
	 endif
	endif
      	if (abs(rasym)>eps) then
         if(0.0.lt.rasym.and.rasym.le.min(rmatr,abs(rmatch))+5*hcm) then
           if(pr)
     >     write(koe,*) 'RASYM < RMATCH,RMATR: SO NO CRCWFN'
         else
          fcwfn = .true.
         IF(RASYM.gt.0.and.pr) 
     *   WRITE(koe,1012) RASYM,ACCRCY,SWITCH,AJSWTCH,sinjmax
!@     * 			,ngail,gailacc
         IF(RASYM.lt.0.and.pr)
     *    WRITE(koe,1013) RASYM,ACCRCY,SWITCH,AJSWTCH,sinjmax
!@     * 			,ngail,gailacc
         endif
         endif
      if(pr)
     >WRITE(koe,1080) NNU,MAXL,MINL,MTMIN,EPC
 1080 FORMAT(/'  NL quadrature with ',I4,' Gaussian points, ', 
     X  'Calculate multipoles up to ',I5,' from',I4,/,
     X  '  M-transfers for lp+lt greater than or equal to ',I2,
     X  ', Angular Integration Cutoff below',F8.4,' %')
      WORD = '    '
      IF(erange.LT.0) WORD = ' MeV'
      IF(ABS(erange).GT.EPS.and.pr)
     >WRITE(koe,1090) erange,WORD,DK
 1090 FORMAT(/' Range of successive Continuum BINS =',F9.4, A4,
     X ', formed by integration with DK <',F7.4,' / fm.')
      if(inh>2) inh = 0
      if(INH>0.and.pr) write(koe,1095) inh,hscale(inh)
 1095 format('  Scale partition radii according to inh =',i2,': ',a10)
      if(plane>0.and.pr)
     >            write(koe,1096) plane,planes(mod(plane,4))
 1096 format(/' Plane wave treatment:',i2,'=',a10)
 
 !! PARALLELISM
      maxthread = 0; numthread = 0
!$      numthread = OMP_GET_MAX_THREADS()
!$      maxthread = OMP_GET_NUM_PROCS()
	mpihelp = max(1,mpihelp)
!# serial tests 	mpihelp = min(mpinodes,mpihelp)
 	mpisets = max(mpinodes/mpihelp,1)
	iams = iame/mpihelp ! set number, 0 to mpisets-1
        mpinodes = min(mpinodes,mpisets*mpihelp)
      if(pr.and.mpihelp+numthread>1) 
     x       write(koe,1098) numthread,maxthread,mpisets,mpihelp
 1098 format(/' Requested number of OpenMP threads:',i4,
     x  ' (out of node limit of ',i4,' processors) in each of ',
     x  i4,' MPI sets with ',i3,' helpers per set')
      if(numthread>1 .and. ompform<0)then
!       ompform = 2  ! use from /fresco/ namelist only
      elseif(numthread<2)then
!       ompform = 0
      endif
!$    if(pr)write(koe,1099)ompform
 1099 format(/' Calculating formf using threaded version ',i1)
	if(numthread==0) numthread=mpihelp; ! parallel MPI
      isparallel = parall()  
      IF(lnbl(TMP).gt.1)then
       ilen=len(TRIM(TMP))
       if(TMP(ilen:ilen)/="/") TMP(ilen+1:ilen+1)="/"
      endif
      IF(lnbl(TMP).gt.1.and.MACH>=3) then
          if(pr) write(koe,1004) mpinodes,trim(TMP)//'*'
 1004  format(' Calculation with',I3,' MPI nodes and',
     X        ' temporary files in ',A)
      endif
!################################################################## 
#ifdef MPI
! Supernumeries can exit
      if(iame+1 > mpisets*mpihelp) then
      	 write(6,*) ' MPI node ',iame,' unneeded'
	 call MPI_finalize(ierr)
      	 stop
      	 endif
! Helpers go away here
	ihelp = mod(iame,mpihelp)
	master = iame/mpihelp
	call MPI_COMM_SPLIT(MPI_COMM_WORLD,master,ihelp,commgroup,ierr)
	!write(6,*) ' iame,ihelp =',iame,ihelp, ', err =',ierr
	
#ifdef SCP
! set up Scalapack grid
!	CALL BLACS_PINFO( IAM, NPROCS )
!	write(500+iame,*) ' BLACS_PINFO gives ',IAM,NPROCS	
            CALL BLACS_GET( -1, 0, SYCTXT )
	write(500+iame,*) ' BLACS_PINFO got SYCTXT ',SYCTXT
	allocate(usermap(1,mpihelp),ICNTXT(mpisets))
	
	do i=1,mpisets  ! different J/pi sets, each with scalapack context
	usermap(:,:) = 0
	 do j=1,mpihelp
	 usermap(1,j) = (i-1)*mpihelp + j-1  
	 enddo
	 write(500+iame,1104) i,usermap(1,1:mpihelp)
1104	 format(' Jpi set ',i3,' has usermap =',20i4)

	 NPCOL = 1; NPROW = mpihelp  
	 NPROW = 1; NPCOL = mpihelp
            ICNTXT(i) = SYCTXT
	 call BLACS_GRIDMAP(ICNTXT(i),usermap,1,NPROW,NPCOL)
	 
	write(500+iame,*) ' BLACS_GRIDINIT done, giving',ICNTXT(i)
            CALL BLACS_GRIDINFO( ICNTXT(i), NPROW,NPCOL,MYROW,MYCOL)
	write(500+iame,*) ' BLACS_GRIDINFO gives R,C =',MYROW, MYCOL

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)	
	enddo
	deallocate(usermap)
#endif /* SCP */		

	if(ihelp > 0) call erwinhelper
#endif /* MPI */
!################################################################## 

      if(.not.nml) read(ki,1100) chans,listcc,treneg,cdetr,smats,
     x     xstabl,nlpl,waves,lampl,veff,kfus,wdisk,bpm,melfil,cdcc,nfus,
     x     tcfile
 1100 format(40i2)
!      write(koe,1100) chans,listcc,treneg,cdetr,smats,xstabl,nlpl,waves,
!     x        lampl,veff,kfus,wdisk,bpm,melfil,cdcc,nfus,tcfile
      if(pr)
     >WRITE(koe,1110) CHANS,LISTCC,TRENEG,CDETR,SMATS,XSTABL,NLPL,WAVES,
     X        LAMPL,VEFF,KFUS,WDISK,BPM,melfil,cdcc,nfus,tcfile
 1110 FORMAT(//' Trace switches are',
     X' :  CHANS =',I3,', LISTCC =',I3,', TRENEG =',I3,', CDETR =',I3
     X,', SMATS =',I3,', XSTABL =',I3,', NLPL =',I3//23x,'WAVES =',I3,
     X', LAMPL =',I3,', VEFF =',I3,', KFUS =',I3,', WDISK =',I3,
     X', BPM =',I3,', MELFIL =',I3//23x,'CDCC =',I3,', NFUS =',I3,
     x', TCFILE =',I3)
      IF(PADE.GE.1.AND.(VEFF.NE.0 .OR. WAVES.NE.0)) then
	if(pr)WRITE(koe,1120)
 1120 FORMAT(/' CAUTION: WAVE FUNCTIONS ARE NOT IMPROVED BY PADE ',
     X 'ACCELERATION - ONLY THE S-MATRIX ELEMENTS.'/)
	if(psiren.and.pr) write(koe,11201) 
11201 format('          WAVE FUNCTIONS SIMPLY ',
     X'RESCALED BY S-MATRIX CORRECTIONS')
	endif
      IF(VEFF.ne.0.and.mpinodes.gt.1) THEN
       write(koe,*) '*** ERROR: Parallel VEFF calculation NOT supported'
       VEFF = 0
       ENDIF
      say = chans+listcc+treneg+cdetr+smats>0 .or. rterms
C     IF(PADE.GE.1) VEFF=0
	if(KFUS==0) NFUS=0
	if(KFUS/=0) NFUS=NFUS+1		! NFUS now includes elastic channel
	 nfus1 = max(1,NFUS)
	cfalloc = melfil/=0 .or. fcwfn

 	HBC =   197.32705d0		! hcross * c
!	finec = 137.03599d0		! 1/alpha (fine-structure constant)
	AMU   = 931.49432d0		! 1 amu/c^2 in MeV
	COULCN = HBC/finec		! e^2
!		use unitmass from user input (1 or 1.008 or whatever) in amu
	FMSCAL = 2d0 * unitmass * AMU / HBC**2
	ETACNS = COULCN * sqrt(FMSCAL) * 0.5d0
	pmass = 1.0078  ! amu
        muN   = sqrt(COULCN)*HBC/(2.0D0*pmass*amu)  
	   ! nuclear magneton=e h/(2 mp c) = e hc/(2mp c^2)

*				OLD-STYLE CONSTANTS
!     ETACNS = 0.157454
!     FMSCAL = 1d0/20.736d0
!     FMSCAL = 1d0/20.748d0
*     FMSCAL = 0.0478326  ! == unitmass=0.999740 and finec=137.0455
!     FMSCAL = 0.048196   ! == unitmass=1.007335 and finec=137.5648
!     COULCN = 2. * ETACNS / SQRT(FMSCAL)

      if(pr)
     >WRITE(koe,1121) unitmass,finec,finec/137.03599d0,hbc
 1121 FORMAT(/' Using unit mass =',F9.6,' amu ',
     X   'and 1/fine-structure constant =',f10.5,
     X   ' (',f9.6,' * true ) with hc =',f10.5,' MeV.fm,')
      if(pr)
     >WRITE(koe,1123) FMSCAL,1d0/FMSCAL,ETACNS,COULCN,muN
 1123 FORMAT( '  thus 2*amu/hbar^2 =',F10.7,' = 1/',f7.4,
     x        	' and Coulomb constant =',F10.7,' so e^2 =',f12.8,
     x   ', and nuclear magneton=',f10.7/)
	
      if(INH/=0.and.pr) then
	write(koe,*) ' Rescaling partition radii '
      if(INH<=2) write(koe,*)          '   in order to include',
     X  '  longitudinal recoils in zero-range transfer couplings'
      if(INH==1) write(koe,*) '   according to the projectile masses.'
      if(INH==2) write(koe,*) '   according to the target masses.'
      ! if(INH==3) write(koe,*) '   according to specified PRMAX values'
	write(koe,*) 
	endif

	if(maxval(abs(sock-1.0))>1e-3) write(koe,1124) sock
1124	format('   Factors for spin-orbit folding=',5f6.2,',',5f6.2)
	if(maxval(maxcoup)>0) write(koe,1125) maxcoup
1125	format('   Adjust limits ''maxcoup'' =',3i10)
	if(maxval(expand-expandi)>1e-3) write(koe,1126) expand
1126	format('   Expansion factors for arrays=',11f8.3)
!	                            write(koe,1124) sock

      if(pr)
     >write(koe,'('' Now pre-scan input and save to file 3''/)')
	   call flush(koe)
      ic = 1
	it = 1
	mxx = 0
   30 if(nml) then
	 namep='        ';massp=0; zp=0; nex=0; pwf=.true.; prmax=0.
	 namet='        ';masst=0; zt=0; qval=0; readstates = 0
         mixpot = 0
	 ios=0;  lpmax=-1
         if(trnl) then
         if(iame==0.or.STDINALL)
     >            read(ki,nml=partition,IOSTAT=ios,end=32,err=32)
	 else
         if(iame==0.or.STDINALL)
     >            read(ki,nml=partition) 
	 endif
   32     if ( ios .ne. 0 ) then
            write (koe,*) ' Input namelist PARTITION read error: ',ios
            write (koe,partition)
            stop
         endif
#ifdef ALTIXMPI
         if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nlpartition_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
	else
	 read(ki,1130) namep,massp,zp,nex,cpwf,namet,masst,zt,qval,
     x                 readstates,prmax,lpmax1,mixpot
         lpmax = lpmax1-1
	 pwf = cpwf.eq.'T'.or.cpwf.eq.' '
	endif
!      if(massp+masst.lt.1e-5) go to 62
      if(namep=='        '.or.namet=='        ') go to 62
	name = namep; call ucase(name)
        gamma = name(1:3)=='GAM'.or.name(1:3)=='PHO'
	if(abs(massp)<eps.and..not.gamma) then
	 if(.not.given) then
	   write(koe,*) ' Using mass table ',trim(MASFIL),'for',namep
	   given = .true.
	   endif
        call GETMASS(RM,AM,EM,IA,IZ,IFLAG,namep,MASFIL)
	if(iflag==0) then
	  zp = iz
	  massp = rm
 	  write(koe,33) namep,rm,iz
33 	  format('  Found ',A8,' in mass tables: M,Z =',f8.4,i4)
	  adjq = .true.
	  else
	  write(koe,*) ' Could not get mass of ',namep,': error =',iflag
	  fail = .true.
	  endif
	endif
	if(abs(masst)<eps) then
	 if(.not.given) then
	   write(koe,*) ' Using mass table ',trim(MASFIL),'for',namet
	   given = .true.
	   endif
        call GETMASS(RM,AM,EM,IA,IZ,IFLAG,namet,MASFIL)
	if(iflag==0) then
	  zt = iz
	  masst = rm
 	  write(koe,33) namet,rm,iz
	  adjq = .true.
	  else
	  write(koe,*) ' Could not get mass of ',namet,': error =',iflag
	  fail = .true.	  
	  endif
	endif
	if(ic==1) then
	  totmass = massp + masst + qval*amu
	endif
	if(adjq.and.ic>1) then
	  qval  = (totmass - massp - masst)*amu 
	  write(koe,1129) ic,qval
 1129     format('   Q-value in partition ',i2,' set to',f15.4,' MeV')
	endif

	kii = ki
	if(readstates>0) then
	   write(ko,1131) ic,readstates
 1131      format('  Read states for partition',i2,' from file ',i4)
	   read(readstates,*)  nex
	   write(ko,'(a,1x,i5,1x,a/)') '   to read ',nex,' states'
	   kii = readstates
	   endif

!      write(ko,1130) namep,massp,zp,nex,pwf,namet,masst,zt,qval,
!    x                readstates,prmax,lpmax1,mixpot
      write(ko,nml=partition) 
 1130 format(a8,2f8.4,i4,a1,1x,a8,2f8.4,f8.4,i4,f8.2,2i4)
	call flush(ko)
	mxx = max(mxx,abs(nex))
	if(ic>mmxp) then
	 write(koe,*) ' Need to increase frxx0 parameter mmxp =',mmxp
	 stop
	 endif
        if(mxx>maxex) then
         write(koe,*) ' Need to increase frxx0 parameter maxex =',maxex
         stop
         endif
	masses(:,ic) = (/ massp , masst /)
	nchmax(ic) = 0
	keep = .false.
      do 60 ia=1,abs(nex)
      if(nml.or.readstates>0) then
        if(.not.keep) then
	jp=0;   copyp=0; bandp =0; ep=0;   kkp=0; tp=0;  cpot=0; 
     	jt=0;   copyt=0; bandt =0; et=0;   kkt=0; tt=0;
	fexch=.false.; ignore=.false.; infam=0; outfam=0
	endif

	 ios=0
       if(trnl) then
        if(iame==0.or.STDINALL)
     >            read(kii,nml=states,IOSTAT=ios,end=54,err=54)
       else
        if(iame==0.or.STDINALL)
     >            read(kii,nml=states) 
       endif
   54  if ( ios .ne. 0 ) then
            write (koe,*) ' Input namelist STATES read error: ',ios
            write (koe,states)
            stop
         endif
#ifdef ALTIXMPI
       if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nlstates_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
        else
       read(kii,1150)  jp,  copyp,bandp ,ep,  kkp,tp, cpot,
     X		jt,  copyt,bandt ,et,  kkt,tt,  ex,infam,outfam
	fexch = ex(1)(1:1).eq.'T' .or. ex(1)(2:2).eq.'T'
	ignore = ex(2)(1:1).eq.'T' .or. ex(2)(2:2).eq.'T'
	endif
	if(cpot==0) cpot = ic
        if(ia==1.and.bandp==0) bandp = 1
        if(ia==1.and.bandt==0) bandt = 1
!       write(ko,1150) jp,  copyp,bandp ,ep,  kkp,tp, cpot,
!     X		jt,  copyt,bandt ,et,  kkt,tt,  ex, infam,outfam
        write(ko,nml=states) 
	call flush(ko)
 1150 format(f4.1,2i2,f8.4,2f4.1,i4,2x,f4.1,2i2,f8.4,2f4.1,2a2,2i4)
	it = it+1
	jex(1,ia,ic) = jp
	parity(1,ia,ic) = ptyp
	if(copyp>0) jex(1,ia,ic) = jex(1,copyp,ic)
	if(copyp>0) parity(1,ia,ic) = parity(1,copyp,ic)
	jex(2,ia,ic) = jt
        parity(2,ia,ic) = ptyt
	if(copyt>0) jex(2,ia,ic) = jex(2,copyt,ic)
	if(copyt>0) parity(2,ia,ic) = parity(2,copyt,ic)
        ichmax = nint(2*jex(1,ia,ic)*jex(2,ia,ic) + 
     x		        jex(1,ia,ic)+jex(2,ia,ic)+1)
	nchmax(ic) = nchmax(ic)+ichmax
	cpots(ia,ic) = cpot
	nexx(ic) = nex
   60	continue
	ic = ic+1
	go to 30
!   62	write(ko,nml=partition)
   62	write(ko,'('' &partition /'')')
	call flush(ko)
	mxp = ic-1
	mxpex = it-1
	 NFUS = min(NFUS,mxpex)         ! the total number of channels!
	if(iblock==-9) then
	   iblock=it
	   if(iame==0) write(koe,*) ' ALL channels in CC set'
	   endif
	if(fail) then
	   write(koe,*) ' STOPPING, as some masses unknown'
	   stop
	  endif
	if(mxp>1) then
	iaf =  real(mxp*mxx)**2*2*ma*sizeof(alogic)/1048576.
	 if(pr.and.iaf>200)write(6,*)' Allocate afrac in ',iaf,' MB'
        if(.not.allocated(afrac))allocate(afrac(mxp,mxx,mxp,mxx,2,ma))
	afrac = .false.    !afrac = 0.0 
	endif
         lsp = 0; jsp = 0.
	 nkmax = 0; haso(:) = .false.
	 knia(:) = 0; knib(:) = 0; inia(:)=0
	if(abs(DK)>1e-9) then
	 nkmax = erange/DK
	 if(erange<0) nkmax = sqrt(fmscal*abs(erange))/DK
	endif
	do kn=1,ma
	  kna(kn)=0
	  knb(kn)=0
	  enddo

	call flush(ko)
	kpmax = 0
	nf = 0; mpwcoup = 0
	ityp = -1
	nix = 1
	pcomps(:) = 0
	kdepfac = 1
70       if(nml) then
	  loop=-1
	  kp=0; type=0; itt=.false.;shape=0; nosub=.false.
	  def(:)=0; defp(:)=0; deft(:)=0; p(:)=0; mnep(:)=0; mnet(:)=0
	  ap=0;at=0;rc=0;ac=0;v=0;rv=0;av=0;
	  w=0;rw=0;aw=0;wd=0;wdr=0;wda=0;vd=0;vdr=0;vda=0
          vso=0;rso=0;aso=0;vsoi=0;rsoi=0;asoi=0
	  jl=' ';xlvary=0;alvary=0;lshape=0;lshapei=0
	  datafile=' '
	 ios=0
         if(trnl) then
          if(iame==0.or.STDINALL)
     >            read(ki,nml=pot,IOSTAT=ios,end=705,err=705)
	 else
          if(iame==0.or.STDINALL)
     >            read(ki,nml=pot) 
	 endif
  705     if ( ios .ne. 0 ) then
            write (koe,*) ' Input namelist POT read error: ',ios
            write (koe,pot)
            stop
	  endif
#ifdef ALTIXMPI
         if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nlpot_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
	  if(abs(ap)+abs(at)<eps) at=1d0
	  if(maxval(abs(p))+maxval(abs(def))>0) loop=-2  ! Ordinary input only
	else
	  loop=-2
	  datafile=' '
	  p(0)=0.;p(8)=0.
	  read(ki,72) kp,type,citt,shape,(p(k),k=1,7)
72        format(i3,i2,a1,i2,f8.3,6f8.4)
	    jl=' ';xlvary=0;alvary=0;lshape=0;lshapei=0
	  if(shape>=30.and.shape<40) then
	    read(ki,73) jl,lshape,xlvary,alvary
73	    format(4x,a1,i3,2f8.4)
	  endif
	  itt = citt=='1'.or.citt=='3'
	  nosub=citt=='2'.or.citt=='3'
	endif
	kpi = kp
	typei = type
	shapei = shape
	lshapei = lshape
	kpmax = max(kp,kpmax)
c			Find partition using this coupling:
	    ick = 0
	    do 1295 ic=1,mxp
	    do 1295 ia=1,abs(nexx(ic))
1295	    if(cpots(ia,ic)==kp) ick = ic
        if(kpi.eq.0) then
!	  write(ko,nml=potl)
   	  write(ko,'('' &potl /'')')
	  call flush(ko)
	  go to 80
	endif
71	loop=loop+1
        if(loop==0) then
          typei=0; p1=at; p2=ap; p3=rc; p4=ac; p(5:7)=0
                  if(abs(rc)<eps) go to 71
        else if(loop==1) then
          typei=10; p = mnep
                  if(sum(abs(p))<eps) go to 71
        else if(loop==2) then
          typei=11; p = mnet
                  if(sum(abs(p))<eps) go to 71
        else if(loop==3) then
          typei=1; p1=v; p2=rv; p3=av; p4=w; p5=rw; p6=aw; p7=0
                  if(abs(p1)+abs(p4)<eps) go to 71
        else if(loop==4) then
          typei=2; p1=vd; p2=vdr; p3=vda; p4=wd; p5=wdr; p6=wda; p7=0
                  if(abs(p1)+abs(p4)<eps) go to 71
        else if(loop==5) then
          typei=10; p = defp
                  if(sum(abs(p))<eps) go to 71
        else if(loop==6) then
          typei=11; p = deft
                  if(sum(abs(p))<eps) go to 71
        else if(loop==7) then
          typei=3; p1=vso; p2=rso; p3=aso; p4=vsoi;
                   p5=rsoi; p6=asoi;p7=0
                  if(abs(p1)+abs(p4)<eps) go to 71
        else if(loop>7) then
          go to 70
        endif
        type = abs(typei)
        shape = abs(shapei)
	kp = abs(kpi)
	if(shape==45) kdepfac = max(kdepfac,2)
!	 if(say) write(48,*) ' kdepfac = ',kdepfac
	haso(kp) = haso(kp).or.(type==3.and.abs(p1)+abs(p4)>1e-9)
 	  if(say) write(48,711) kp,loop,type,haso(kp),p
711	  format(' kp,loop =',2i3,' is type,so,p(:) =',i4,l2,13g12.4)
! 721     format(i3,i2,a1,i2,2f8.4,2f8.2,2f8.1,f8.4)

      if(type.ge.10.and.type.le.16) then
               if(shape.eq.0.or.shape.gt.13) shape = 10
c                        remove next line when coul quadrature written.
               if(ityp.eq.0) shape = min(shape,10)
      		write(ko,nml=potl) 
	        call flush(ko)
	    nf0 = nf
            do 1290 k=0,7
            if(k.ne.0 .and. abs(p(k)).lt.1e-9) go to 1290
            if(k.eq.0 .and. shape.ne.12.and.shape.ne.13) go to 1290
               nf= nf + 1
      		if(type.ge.14.and.type.le.16) nf=nf+1
1290        continue
	    knused(kp) = knused(kp) + (nf-nf0+1)
	    pcomps(ick) = pcomps(ick) + (nf-nf0+1)
	    
       if(type.lt.12.or.type.ge.18) go to 1350
1301    if(nml) then
	  ib=0; ia=0; k=0; str=0
	  if(iame==0.or.STDINALL)
     >            read(ki,nml=step)
#ifdef ALTIXMPI
          if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nlstep_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
	  else
!	  read(ki,*) ib,ia,k,str	! must be able to read blank line
	  read(ki,1302) ib,ia,k,str
	  endif
1302     format(4x,3i4,f8.4)
!        write(ko,*) ib,ia,k,str
        write(ko,nml=step) 
	   call flush(ko)
         symm = .false.
         if(ib.eq.0) go to 1350
            nix = nix + 1
            if(ib.gt.0) go to 1301
1350        continue
	else  if(type.eq.17) then
	 if(ick>0) then
	  nf = nf + nchmax(ick)*(nchmax(ick)+1)/2
	 else
	     write(koe,*) 'TYPE=17: cannot find partition using ',
     x		'potential KP =',kp
	     write(koe,*) (nexx(ic),(cpots(ia,ic),ia=1,abs(nexx(ic))),
     x				ic=1,mxp)
	     stop
	 endif
      	  write(ko,nml=potl) 
	   call flush(ko)
1361    if(nml) then
          ib=0; ia=0; k=0; str=0
          if(iame==0.or.STDINALL)read(ki,nml=step)
#ifdef ALTIXMPI
          if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nlstep_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
          else
          read(ki,*) ib,ia,k,str
          endif
        write(ko,nml=step)
         if(ib.eq.0) go to 1380
            nix = nix + 1
         if(ib.gt.0) go to 1361
1380     continue
	else  if(type.eq.20 .or. type.eq.21) then
!      	write(ko,72) kpi,typei,citt,shapei,(p(k),k=1,7)
      	write(ko,nml=potl) 
         shape = max(1,min(12,shape))
         nf = nf + shape
	 knused(kp) = knused(kp) + 4
      	else
!      	write(ko,72) kpi,typei,citt,shapei,(p(k),k=1,7)
      	write(ko,nml=potl) 
         if(typei.ge.0) nf = nf + 1	        
	 knused(kp) = knused(kp) + 1
	 pcomps(ick) = pcomps(ick) + 1
	endif

      if(type.le.9.or.type.ge.18) ityp = type
	if(loop>=0.and.loop+1<=9) go to 71
	if(kpi>0) go to 70
80	mpair = nix-1
	
!	write(KO,*) 'HASO(1:10) =',HASO(1:10)

	call flush(ko)
	msp = 0	
	tnt=0; coef=0
        lamax = 0
        qcmax = 0
	keep = .false.

!	write(KO,*) ' Pluto(:) =',pluto(:)
	npluto=0
	do i=1,10
	if(pluto(i)>0) npluto=i
	enddo
	if(npluto>0) 
     x	  write(koe,*) ' Pluto will use potentials #',pluto(1:npluto)

85	if(nml) then
	  if(.not.keep) then
          kn1=0;kn2=0;ic1=0;ic2=0;in=0;kind=0;ch1=' ';nn=0;l=0;lmax=0;
          sn=0;ia=0;j=0;ib=0;dm=0;nk=0;er=erange;e=0
          kbpot=0;krpot=0;be=0;isc=0;ipc=0;nfl=0;nam=0;ampl=0
	  tnt=0; coef=0; ppower=0.; rsmin=0.; rsmax=rsp; nlag=0; phase=0
	  autowid=0.
	  endif

	 ios=0
         if(trnl) then
         if(iame==0.or.STDINALL)
     >            read(ki,nml=overlap,IOSTAT=ios,end=851,err=851)
	 else
         if(iame==0.or.STDINALL)
     >            read(ki,nml=overlap) 
	 endif
  851     if ( ios .ne. 0 ) then
            write (koe,*) ' Input namelist OVERLAP read error: ',ios
            write (koe,overlap)
            stop
	  endif
#ifdef ALTIXMPI
          if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nloverlap_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
	 if(abs(e)>eps) be = -e
 	else
	dm=0;er=0;nk=0;e=0; datafile=' '
        read(ki,852) kn1,kn2,ic1,ic2,in,kind,ch1,nn,l,lmax,sn,ia,j,ib,
     &         kbpot,krpot,be,isc,ipc,nfl,nam,ampl,nk,er,rsmin,rsmax,
     &         ppower,nlag,phase,autowid
        endif
	ini = in
	im = abs(in)
852   format(2i3,4i2,1x,a1,3i2,f4.1,i2,f4.1,i2,2i3,f8.4,4i3,f8.4,
     X		i3,f8.3,3f8.4,i4,2f8.3)
	if(max(kn1,kn2)>ma) then
	 write(koe,*) ' Need to increase frxx0 parameter ma =',ma,
     x	   ' to at least ',max(kn1,kn2)
	  stop
	 endif
      if(kn1*ic1*ic2*ini.eq.0) then
!	kn1 = 0
!	write(ko,nml=overlap)
   	write(ko,'('' &overlap /'')')
	go to 105
	endif
!       write(ko,852) kn1,kn2,ic1,ic2,ini,kind,ch1,nn,l,lmax,sn,ia,j,ib,
!     &         kbpot,krpot,be,isc,ipc,nfl,nam,ampl,nk,er,ppower
          dmm = masses(im,ic1) - masses(im,ic2)
      	  icr = ic2
      	  icp = ic1
      	  if(dmm<0) then
	    icr = ic1
      	    icp = ic2	 
	  endif
!********************************************************
        kn2 = max(kn2,kn1)
        lsp(kn1:kn2) = l
        ssp = sn
        if(ib==0)then
!	   if(ccbins.and.kind.ne.3.and.im==2) ib=nexx(icp)+
           if(.not.ccbins.or.(ccbins.and.kind.ne.3)) go to 820
	   ib=kn1
	   endif
        jcom = jex(im,ib,icp)
        il1=0
        lmin=0
        if(lmax==0)lmin=l
        if(lmax==0)lmax=l
        nch = kn2-kn1+1
        do iia=1,abs(nexx(icr))
         jcor=jex(im,iia,icr)
         par=parity(im,ib,icp)*parity(im,iia,icr)
         qcmax = max(qcmax,nint(2*jcor))
         do l1=lmin,lmax
          if((-1)**l1.ne.par)cycle
          njn=nint((lmax+sn)*2.)
          do ij=0,njn
           j1=ij*0.5
           if(fail3(l1+z,j1,sn))cycle
           if(fail3(j1,jcor,jcom))cycle
           il1=il1+1
           if(il1<=nch)then
            lsp(kn1-1+il1)= l1
            jsp(kn1-1+il1)= j1
            knia(kn1-1+il1) = iia
            inia(kn1-1+il1) = im
            lamax = max(lamax,nint(2*j1))
           endif
          enddo
         enddo
        enddo
!********************************************************
!         if(nam>0.and.ia>0.and.ib>0) then
          if(nam>0) then
           do kn=kn1,kn2
             if(abs(ampl)*sqrt(real(nam))>eps)then
              afrac(icp,max(ib,1),icr,max(ia,1),im,kn)=.true.
             endif
!             afrac(icp,max(ib,1),icr,max(ia,1),im,kn)
!     X          = ampl*sqrt(real(nam))
!            write(6,*) 'afrac(',icp,max(ib,1),icr,max(ia,1),im,kn,
!    X          ') =',ampl*sqrt(real(nam))
	   enddo
	  endif
820	continue
        if(ib>0)knib(kn1:kn2) = ib
        knb(kn1:kn2) = kn2
        kindsp(kn1:kn2) = kind
        if(kind/=3)jsp(kn1:kn2) = j
        kna(kn1:kn2) = kn1

!	write(koe,*) 'Overlap ',kn1,': DMM,DM =',real(dmm),real(dm)
!	write(koe,*) 'masses:',in,ic1,ic2,((masses(i,k),i=1,2),k=1,mxp)
	if(abs(dm)<eps) then
	dm = abs(dmm)
      	IF(INI.LT.0) then
	   DM = DM + BE/AMU
*	   BE correction should depend on gs energy differences, 
*	    and not the energy of the current (maybe excited) state!:
C	   Yes: but we do not know gs separation energy!
c     	   DM = DM + (QVAL(ICP)-QVAL(ICR))/AMU
	 endif
	endif
       if(nam.eq.-1) then
         dm = ampl
       else if(nam.le.-2) then
          nk = abs(nam)
          if(nk>2) nk = nk * 5
          er = ampl
          if(abs(er)<eps) er = erange
          if(abs(er)<eps) er = 1.2
       endif
	if(er < 0.and.nk==0) nk = 10
	  nkmax = max(nkmax,nk,10)
	  if(rsmax<hcm) rsmax=rsp
	  rsmax = min(rsmax,rsp)

#ifdef corex
!      if(be<0.and.kind>=1.and.kind<=3) ccbins = .true.  ! coupled channels bins
      if(be<0.and.(kind>=1.and.kind<=3.or.
     &            (kind==0.and.isc<0))) then
         ccbins = .true.  ! coupled channels bins
!	 write(48,*) ' CCBINS, because',kn1,be,kind,isc
	 if(kind==3) nchbinmax = max(nchbinmax, kn2-kn1+1)
	  call flush(48)
	 endif
	 write(48,*) ' CCBINS=',ccbins,', as',kn1,be,kind,isc
#endif /* corex */
      if(kind>=6) nn2wf = .true.   			! 2-nucleon bound states
      write(ko,nml=overlap) 
	msp = max(msp,kn2)
      if(kind.eq.4.and.nn.ge.2) then
	if(nml) then
	 triton=0
	 if(iame==0.or.STDINALL)
     >            read(ki,nml=dalitz)
#ifdef ALTIXMPI
         if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nldalitz_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
	else
	 if(iame==0.or.STDINALL)
     >            read(ki,*) triton
	endif
	write(ko,nml=dalitz)
	endif
      if(kind.gt.5) then
      	nk = nn
        mpair = max(mpair,nk)
      	if(nk.gt.0) then
	  if(nk.le.mp) then
	   if(nml) then
            if(iame==0.or.STDINALL)read(ki,nml=twont)
#ifdef ALTIXMPI
            if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nltwont_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
	   else
            read(ki,8626) ((tnt(i,jj),i=1,4),coef(jj),jj=1,nk)
8626  		format(3(4i3,f8.4))
           endif
!            write(ko,*) ((tnt(i,jj),i=1,4),coef(jj),jj=1,nk)
            write(ko,nml=twont)
	    else
	        write(koe,*) ' Need to increase frxx0 parameter mp to ',
     x			' >=',nk
	  	stop
	    endif
	endif
	endif
	go to 85
	
  105	cp = 1
	cxwf = complexbins .or. ccbins
        sumccbins = ccbins
        sumkql = ccbins
	call flush(ko)
	maxqrn=0 ; maxqrn1=0 ; maxqrn2=0
	mnpair=0; npair = 0
  110 if(nml) then
	icto=0;icfrom=0;kind=0;ip1=0;ip2=0;ip3=0;ip4=-1;ip5=-1;infile=4
	p1=0;p2=0;jmax=0;rmax=0;kfrag=0;kcore=0;qcm=0;lam=0;nforms=0

	 ios=0; nlcouplingend=0
         if(iame==0.or.STDINALL)read(ki,nml=coupling,IOSTAT=ios)
         if ( ios > 0 .and. trnl) then
            write (koe,*) ' Input namelist COUPLING read error: ',ios
            write (koe,coupling)
            stop
	  endif
  118    if ( ios < 0 ) then
            nlcouplingend=1
          endif 
#ifdef ALTIXMPI
        if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(nlcouplingend,1,MPI_INTEGER,
     >                           0,MPI_COMM_WORLD,ierr)
        if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nlcoupling_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
        if(nlcouplingend==1)goto 119
      else
	kfrag=0;kcore=0
	read(ki,1220) icto,icfrom,kind,ip1,ip2,ip3,p1,p2,
     x                   jmax,rmax,ip4,ip5,nforms,infile
      endif
      if(icto.eq.0.or.icfrom.eq.0) go to 119
!      if(kind>8) kind=kind-8
	qq = ip1
	irem = ip2
	lam = lamax
        qcm = qcmax
        ikmax = lamax+qcmax
        ikm = ikmax
        if(ip4>=0.or.ip5>=0) ikm = abs(qq)
	if(ip4>=0) qcm = ip4
        if(ip5>=0) lam = ip5
        if(qcm==0) ikm = abs(qq)
        if(ip4<0.and.ip5<0) lam = abs(qq)
        if(kfrag>0) p1 = kfrag
	if(kcore>0) p2 = kcore
*	kfrag = nint(p1)
*	kcore = nint(p2)
	foldso = .false.
	if(p1>0.and.p2>0)
     x  	foldso = haso(nint(p1)) .or. haso(nint(p2))
	foldso = .false.  !!! NOT PROPERLY IMPLEMENTED YET!
	nty = 1
	if(foldso) then
	   kopt = cpots(1,abs(icto))
	   if(haso(kopt))     nty=nty+1 ! vso(opt)
	   if(haso(nint(p1))) nty=nty+4 ! vso(frag)
	   if(haso(nint(p2))) nty=nty+4 ! vso(core)
	   nty1=nty
	   if(ccbins) nty = nty+nint(expand(4)) ! margin!
		write(ko,*) ' FOLDSO: using kopt =',kopt,' nty=',nty
	   if(sumccbins)write(ko,*)'  No summing of cc bin formfactors'
           sumccbins= .false.
          endif
        ! if(kind>=5.and.kind<=8) sumccbins= .false.
        if(sumform>1 .and. .not.sumccbins)sumform=1
        if(.not.ccbins)sumform=0
        if(sumform<2)sumccbins=.false.
        if(sumform<1)sumkql=.false.
        if(.not.allocated(usedcc))allocate(usedcc(mxx,mxx,0:4*iabs(qq)))
	usedcc = .false.
!      write(ko,1220) icto,icfrom,kind,ip1,ip2,ip3,p1,p2,jmax,rmax,ip4,ip5
      write(ko,nml=coupling) 
 1220 format(3i4,3i2,2f8.2,2f4.0,4i4)
      symm = symm .and. icto>0
	  used = .false.
	  usedcore(:,:) = 0
       	  if(kind.eq.5.or.kind.eq.2) used(1,1) = .true.
       	  if((kind.eq.5.or.kind.eq.6).and.INH.ne.2) then 
	    write(koe,*) ' WARNING: target longitudinal recoil NOT',
     X	     ' included in the ZR coupling # ',cp
	    endif

  	 if(kind.le.8.and.kind>1) then
	   keep = .false.
8735  	   if(nml) then
	     if(.not.keep) then
	     in=0;ib=0;ia=0;kn=0;a=0
	     endif
	     if(iame==0.or.STDINALL)read(ki,nml=cfp,end=8738) 
8738         continue
#ifdef ALTIXMPI
             if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nlcfp_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
	    else
	     read(ki,8739) in,ib,ia,kn,a
	    endif
           if(in.eq.0) then
	      write(ko,'('' &cfp /'')')
	      go to 876
	      endif
!   	   write(ko,8739) in,ib,ia,kn,a
       	   write(ko,nml=cfp) 
8739  	   format(4x,4i4,f8.4)
!	   call flush(ko)
	    ini = in
	    in = abs(ini)
		  ic1 = abs(icto)
		  ic2 = icfrom
	      	  dmm = masses(in,ic1) - masses(in,ic2)
 	     	  icr = ic2
 	     	  icp = ic1
 	     	  if(dmm<0) then
		    icr = ic1
 	     	    icp = ic2	 
		  endif
		do kk=kn,knb(kn)
	       if(icp>mxp.or.ib>mxx.or.icr>mxp.or.ia>mxx.or.kk>ma) then
	       write(0,*) ' One of input afrac parameters too large'
             write(6,'(5(2i8,1h,))') icp,mxp,ib,mxx,icr,mxp,ia,mxx,kk,ma
	       stop
	       endif
		  !afrac(icp,ib,icr,ia,abs(in),kk) = a
		  if(abs(a)>eps)afrac(icp,ib,icr,ia,abs(in),kk)=.true.
        	   knia(kk) = ia   ! may not have been specified yet
        	   knib(kk) = ib   ! ditto
        	   inia(kk) = in
                enddo
	   if(ini.gt.0) go to 8735
876	   continue
	if(kind.le.4) in1 = mod(abs(kind)-1,2)+1
	if(kind.eq.2) in1 = 2
	usedlow(:) = mxx+1
	usedhigh(:) = 0
	  do in=1,2
		  ic1 = abs(icto)
		  ic2 = icfrom
	      	  dmm = masses(in,ic1) - masses(in,ic2)
 	     	  icr = ic2
 	     	  icp = ic1
 	     	  if(dmm<0) then
		    icr = ic1
 	     	    icp = ic2	 
		  endif
	       if(abs(kind)>4) in1 = in
	  do kk=1,msp
	    do ia=1,mxx
	    do ib=1,mxx
	    ! a     = afrac(icp,ib,icr,ia,in,kk)
	     alogic = afrac(icp,ib,icr,ia,in,kk)
             !uu = abs(a).gt.1e-10 .or. used(in,kk)
             uu = alogic .or. used(in,kk)

             if((kind.eq.5.or.kind.eq.2).and.in==1) then
                 uu = kk.eq.1
                 kna(1) = 1
                 knb(1) = 1
                 endif

	     used(in,kk) = uu
       	     if(in==in1.and.alogic) then
		      usedlow(kk) = min(ib,usedlow(kk))
		      usedhigh(kk) = max(ib,usedhigh(kk))
	        usedcore(in,kk) = 
     x                max(1,nint(2*jex(in,ia,icr)+expand(5)))
      	      endif
            enddo
            enddo
	  enddo
	  enddo
	 endif
	 if((kind.eq.3.or.kind.eq.4).and.ip3.ge.10) then
	 if(nml) then
	 qscale(:)=(0.,0.)
         if(iame==0.or.STDINALL)read(ki,nml=scale) 
#ifdef ALTIXMPI
         if(MPIC.and..not.STDINALL)
     >            call MPI_BCAST(MPI_BOTTOM,1,nlscale_struct,
     >                           0,MPI_COMM_WORLD,ierr)
#endif /* ALTIXMPI */
	 else
         read(ki,232) (qscale(i),i=max(0,-qq),abs(qq))
232   		format(6e12.4)
	 endif
!         write(ko,232) (qscale(i),i=max(0,-qq),abs(qq))
         write(ko,nml=scale) 
	 endif
	if(kind==1.and.qq==0) then
	   npair = 30   ! we do not know how much data in file 4!
	   npair = max(npair,nint(expand(6)*mxx**2))   ! we do not know how much data in file 4!
	   if(nforms==0) nf = nf+npair
 	   if(nforms/=0) nf = nf+nforms
	   endif
	if(kind==1.and.qq==1) npair = 1
	if(kind==9) then
	   npair = 30   ! we do not know how much data in file 4!
	   npair = max(npair,nint(expand(6)*mxx**2))   ! we do not know how much data in file 4!
	   if(nforms==0) mpwcoup = mpwcoup+npair*abs(jtmax)
 	   if(nforms/=0) mpwcoup = mpwcoup+nforms
	   endif
	mnpair = mnpair + npair
	 write(48,*) ' mnpair increased by ',npair,' to ',mnpair
	nq = 0 ; nq1=0 ; nq2 = 0
         in1 = 1
         in2 = 2 
         icp=1 ; icr=2
         if(abs(kind).le.4) in1 = mod(abs(kind)-1,2)+1
         if(abs(kind).le.4) in2 = in1
	do 879 knt1=1,msp
	do 879 knp1=1,msp
	 if(kna(knp1)<knp1) go to 879
	 if(kna(knt1)<knt1) go to 879
	 do knt=knt1,knb(knt1)
	 do knp=knp1,knb(knp1)
           nkq = 0 ; nkq1 = 0
         if(knt==knt1.and.knp==knp1)then
           nkq2 = 0
         endif
         ln = lsp(knt)
         lnp= lsp(knp)
         sn = ssp
         jn = jsp(knt) 
         jnp= jsp(knp)
         ib = knib(knt)
         ibp= knib(knp)
         ia = knia(knt)
         iap= knia(knp)
         if(in1*icp*icr*ia*iap*ib*ibp/=0)then
          jcom =jex(in1,ib ,icp)
          jcomp=jex(in1,ibp,icp)
          jcor =jex(in1,ia ,icr)
          jcorp=jex(in1,iap,icr)
          if(mod(nint(jcor*2),2)==0)kcor =0.
          if(mod(nint(jcor*2),2)==1)kcor =0.5
         else
          jcom =0.;jcomp=0.;jcor =0.;jcorp=0.;kcor=0.
         endif
         qcmax=nint(jcor+jcorp)
         qcmin=nint(abs(jcor-jcorp))
!@@
	 if(kind.eq.2) then
	       nq = nq + 1 + nint((jsp(knp)+1))  + nint(expand(3))
	       nf = nf + abs(qq) 
	       if(irem>=4) nf = nf + abs(qq) + pcomps(icfrom) ! 2 more
	 endif

      	  if(used(in1,knp).and.used(in2,knt)) then
	    if(kind.ge.5.and.kind.le.8) then
	       if(.not.(kindsp(knt1)>=6.and.knp-knp1/=knt-knt1)) then
	       nq = nq + nint((jsp(knp)+1)*(jsp(knt)+1)*expand(7))
	       endif
	    else
c        kind = 3 & 4 :    s.p. inelastic form factors of multipole 'ik'
             if(.not.sumkql)then
               do 8782 ik=0, max(abs(qq),min(lamax,lam))
                ! if(mod(ik+lnp+ln,2).ne.0) go to 8782
                 if(fail3(ik+z,jn+z,jnp+z).and.
     x             (.not.foldso.or.(fail3(ik+1.0,jn,jnp).and.
     x                              fail3(ik-1.0,jn,jnp)) )
     X             .and.kindsp(knt)==0.and.kindsp(knp)==0) go to 8782
                 if(qq.lt.0 .and. ik.ne.abs(qq)) go to 8782
!                                         Impose order to use symmetry
!                                         when the wfs are real:
                 if(.not.cxwf.and.usedlow(knt) > usedhigh(knp)
     X                 .and..not.foldso) go to 8782
                 do 8781 qc=qcmin,min(qcmax,qcm)
                 do 8781 lambda=0,qc
                   if(mod(ln+lnp+ik+qc-lambda,2)/=0) go to 8781
                   if(mod(qc,2)/=0) go to 8781
                 nq = nq + nty
	if(listcc>10)write(6,*) 'knp,knt,ik=',knp,knt,ik,
!     x                          ' qc,lambda=',qc,lambda,
     x                          ' nq =',nq
                 if(foldso) nq = nq + nty
                 if(ccbins) nkq = nkq + nty1**2
 		 if(.not.(ccbins.or.foldso))then
		   nkq = nkq + nty**2
     x               * max(usedcore(in1,knp),usedcore(in2,knt))
		  endif
8781             continue
8782           continue
             else !if(sumkql)then
        ep1=1e-5
               do 920 la=nint(abs(jcom-jcomp)),min(nint(jcom+jcomp),lam)
                   if(fail3(la+z,jcom,jcomp)) goto 920
                 ikql1=0
                 do 910 ik=0, min(ikmax,ikm)
                 do 910 qc=qcmin,min(qcmax,qcm)
                   if (fail3(qc+z,jcor,jcorp)) goto 910
      np0=1
      if(abs(real(ROTORNC(qc,dble(jcor),dble(jcorp),
     >dble(kcor),dble(kcor))))<ep1)np0=0
                 do 905 lambda=0,qc
                   if(mod(ln+lnp+ik+qc-lambda,2)/=0) goto 905
                   if (mod(ik+lambda+la,2)/=0) goto 905
                   if (fail3(la+z,ik+z,lambda+z)) goto 905
        np1=1
         if(W3J(ik+z,lambda+z,la+z,z,z,z)<ep1)np1=0
                 ikql=0
                 do 904 lap=abs(ln-lnp),lnp+ln
                   if(mod(ln+lnp+lap,2)/=0) goto 904
                   if(fail3(lap+z,jn,jnp)) goto 904
                   if(mod(ik+qc-lambda+lap,2)/=0) goto 904
                   if(fail3(lap+z,ik+z,qc-lambda+z)) goto 904
                   if(fail3(qc+z,la+z,lap+z)) goto 904
        np2=1
         if(W3J(ik+z,qc-lambda+z,lap+z,z,z,z)<ep1)np2=0
        np3=1
         if(W3J(lap+z,ln+z,lnp+z,z,z,z)<ep1)np3=0
        np4=1
         if(W6J(jn,jnp,lap+z,lnp+z,ln+z,sn)<ep1)np4=0
        np5=1
         if(W6J(lap+z,la+z,qc+z,lambda+z,qc-lambda+z,ik+z)<ep1)np5=0
        np6=1 
         if(W9J(jcom,jcomp,la+z,jn,jnp,lap+z,jcor,jcorp,qc+z)<ep1)np6=0
!       write(401,'(2(i5,i2,4f4.1),f4.1,5i2,7i2)')
!     >knt,ln,jn,jcor,kcor,jcom,knp,lnp,jnp,jcorp,kcor,jcomp,sn,
!     >la,ik,qc,lambda,lap,
!     >np0,np1,np2,np3,np4,np5,np6
       if(np0*np1*np2*np3*np4*np5*np6>0)then
         ikql=ikql+1
         ikql1=ikql1+1
       endif
!        write(2,'(7i4,7i5,2x,2f4.1,3i4,f4.1,2x,i4,2f4.1)')
!     >knt,knp,la,ik,qc,lambda,lap,np0,np1,np2,np3,np4,np5,np6
!     >,jn,jnp,lap,lnp,ln,sn,qc,jcor,jcorp
    ! x                * max(usedcore(in1,knp),usedcore(in2,knt))
904              continue
               if(ikql>0)then
                 nq = nq + nty
         	 nkq = nkq + nty**2
                 if(.not.usedcc(ib,ibp,la))then
                  nq2 = nq2 + nty
                  nkq2 = nkq2 + nty**2
                  usedcc(ib,ibp,la)=.true.
                 endif
               endif
905              continue
910              continue
               if(ikql1>0)then
                 nq1 = nq1 + nty
                 nkq1 = nkq1 + nty**2
!		 write(6,'(a,3i4,l2,i3)') 'la,ib,ibp,usedcc,nty',
!     x                             la,ib,ibp,usedcc(ib,ibp,la),nty
               endif
920            continue
             endif !sumkql
		nkq  = nkq *expand(1)			! fudge factor!
                nkq1 = nkq1*expand(8)
                nkq2 = nkq2*expand(10)	
	        maxqp  = max(maxqp ,abs(qq)/2+1,nkq )
	        maxqp1 = max(maxqp1,abs(qq)/2+1,nkq1)
	        maxqp2 = max(maxqp2,abs(qq)/2+1,nkq2)
	    endif
	  endif
	  enddo !knt,knf
	  enddo !knt,knf
!  	write(6,'(a,5i6)') 'knt,knp,nq,nq1,nq2=',knt,knp,nq,nq1,nq2
879 	 continue !knt,knf
        nq  = nq *expand(2)      ! fudge factor!
	nq1 = nq1*expand(9)      ! fudge factor!
	nq2 = nq2*expand(11)     ! fudge factor!
	maxqrn  = max(nq ,maxqrn )
        maxqrn1 = max(nq1,maxqrn1)
        maxqrn2 = max(nq2,maxqrn2)
	if(kind>=3.and.kind<=6)then
         if(maxcoup(1)==0)then
          if(sumform==0)then
           nf = nf + nq
          elseif(sumform==1)then
           nf = nf + nq1
          elseif(sumform==2)then
           nf = nf + nq2
          endif
         endif
        endif
	if(nq==0.and.kind/=1.and.kind/=9) then
	  write(koe,*) ' Missing fractional overlaps! for coupling ',cp
	  if(kind/=4.and.kind/=2.and.kind/=5.and.kind/=9) 
     x              write(koe,*) ' Projectile uses:',used(1,1:msp)
	  if(kind/=3) write(koe,*) ' Target     uses:',used(2,1:msp)
	endif
      cp = cp + 1
      if(kind==7.and.abs(irem)==2) cp=cp+1
	if(allocated(usedcc)) deallocate(usedcc)
	go to 110
  119  write(ko,'('' &coupling /'')')
	write(ko,*) 'End of couplings'
       maxcpl = cp-1
	call flush(ko)
	if(.not.nml) PEL = ABS(ICFROM)
         IF(PEL.LE.0) PEL = 1
        if(.not.nml) EXL = KIND
         IF(EXL.LE.0) EXL = 1
        if(.not.nml) LAB = ip1
         IF(LAB.EQ.0) LAB = PEL
        if(.not.nml) LIN = ip2
         IF(LIN.EQ.0) LIN = 1
        if(.not.nml) LEX = ip3
         IF(LEX.LE.0) LEX = 1
	mloc = nf

      if(.not.nml) read(ki,1255) (elab(i),nlab(i),i=1,3),elab(4)
!      write(koe,1255) (elab(i),nlab(i),i=1,3),elab(4)
 1255 format(3(f8.4,i8),f8.4)

	maxcpl = max(maxcpl,1)
	msp = max(msp,1)
	mpair = max(mpair,mnpair,1)
	mloc = max(mloc,1)
	maxqrn  = max(maxqrn ,1)
	maxqrn1 = max(maxqrn1,1)
	maxqrn2 = max(maxqrn2,1)
	maxnnu = max(maxnnu,1)
	mclist  = maxqp +maxval(knused(:))*kdepfac+npair + 2
	if(say) write(48,*) 'mclist =',mclist,' from ', maxqp ,
     x            maxval(knused(:)),npair,',',kdepfac
	mclist1 = maxqp1+maxval(knused(:))+npair + 2
	mclist2 = maxqp2+maxval(knused(:))+npair + 2
      write(ko,1480) mxx,mxp,mxpex,maxit,maxcpl,maxnnu
 1480 format(//' parameters:    mxx    mxp  mxpex  maxit ',
     x'maxcpl maxnnu ' / '   required:',15i7/)
      write(ko,1490)
     x maxn,mint,maxnln,maxnlo,msp,lmax1,mpair
 1490 format(/' parameters:   maxn  mint  maxnln maxnlo ',
     x'   msp  lmax1  mpair' / '   required:',16i7/)
      write(ko,1495) mloc,mclist,maxqrn,mpwcoup
!      write(6,1495) mloc,mclist,maxqrn,mpwcoup
 1495 format(/' parameters:     mloc   mclist   maxqrn  mpwcoup',
     X	/ '   required:',4i9/)
#ifndef corex

#else /* corex */
      if(ccbins)write(ko,1496) maxqrn1,maxqrn2,mclist1,mclist2
!      if(ccbins)write(6,1496) maxqrn1,maxqrn2,mclist1,mclist2
 1496 format(/' parameters:  maxqrn1  maxqrn2  mclist1  mclist2' 
     X	    / '   required:',4i9/)
        if(sumccbins)then
         maxqrn1 = maxqrn2
         maxqrn2 = maxqrn
         maxqrn  = maxqrn1
         mclist  = mclist2
        elseif(sumkql)then
         maxqrn2 = maxqrn
         maxqrn  = maxqrn1
         mclist  = mclist1
        endif
#endif /* corex */
        if(maxcoup(1)>0)then
          maxqrn =max(maxqrn ,maxcoup(1))
          mloc = nf + maxqrn 
        elseif(maxcoup(1)<0)then
          maxqrn =min(maxqrn ,abs(maxcoup(1)))
          mloc = nf + maxqrn 
        endif
        if(maxcoup(2)>0)then
          maxqrn2=max(maxqrn2,maxcoup(2))
        elseif(maxcoup(2)<0)then
          maxqrn2=min(maxqrn2,abs(maxcoup(2)))
        endif
        if(maxcoup(3)>0)then
          mclist =max(mclist ,maxcoup(3))
        elseif(maxcoup(3)<0)then
          mclist =min(mclist ,abs(maxcoup(3)))
        endif
      if(ccbins)write(ko,1497) mloc,maxqrn,maxqrn2,mclist
      write(ko,1497) mloc,maxqrn,maxqrn2,mclist,nchbinmax
 1497 format(/' parameters:     mloc   maxqrn  maxqrn2   mclist' 
     x       ,' nchbinmax',
     X	    / '   required:',5i9/)
	write(ko,1500) headng
 1500   format(/a120,/,'NAMELIST')
        line = ''
        write (ko,fresco, iomsg=line)
!       write(6,*) line
	call flush(ko)
	pi = 4d0 * atan(1d0);r4pi = sqrt(0.25/pi)
	acc8 = epsilon(acc8);fpmax = huge(acc8)**0.8d0
	maxmul = max(lmax1 + 6,MAXL+1)
      	llmax = lmax1 - 1
	write(KO,*) 'cxwf=',cxwf,' ccbins=',ccbins,' sumccbins =',
     x     sumccbins,' sumkql =',sumkql,' ccbins=',ccbins,
     x               ' sumform =',sumform
	if(cxwf) write(ko,*)  '  COMPLEX bins used.   TWONN =',nn2wf
#ifdef corex
        if(ccbins)write(ko,*)'  Formfactor reduction level: SUMFORM='
     &                        ,sumform
        if(sumccbins)then
         write(ko,*)'  summing coupled channels formfactors into',
     &             ' composite projectile formfactor '
        elseif(sumkql)then
         write(ko,*)'  summing deformed core multipole formfactors'
        elseif(ccbins)then
        write(ko,*)'  coupled channel bins with no formfactor reduction'
        endif
#endif /* corex */

	if(allocated(afrac)) then
	 if(pr.and.iaf>10) write(ko,*)' Deallocate afrac'
	 deallocate (afrac)
	endif
	if(say) write(48,*) 'FRXX0: MINL,MAXL =',MINL,MAXL
	if(cdc.and.iame>0) close(ko3,status='delete')
	return
	end


	subroutine defaults(TMP,MASFIL)
!@	subroutine defaults(TMP,MASFIL,gailitis)
	use parameters
	use drier
	use kcom
	use trace
	use fresco1
	use parallel
      	use searchpar, only: boldplot,maxleg
        use io, only: grace


      	implicit real*8(a-h,o-z)
    	character*70 TMP,MASFIL
	 hcm=0.1; rmatch=0; rintp=0.5; hnl=0; rnl=0; centre=0; hnn=0; 
         rnn=0; rmin=0;rsp=0
	 rasym=0; accrcy=0.001; switch=0; ajswtch=0
	 sinjmax=0 ; plane=0; lmax=-1; ccreal=.false.; dist0=-1.
!@	 gailitis=0; gailacc=0; ngail=0
         jtmin=0; jtmax=0; absend=0.001; pset=0; jset=0; rela='  '
	 jump(:,:)=0; jbord(:)=0; jleast=0.; nearfa=0; 
         kqmax=0; pp=0; thmin=0; thmax=180; thinc=1; koords=0; dgam=0
	 cutl=0; cutr=0; cutc=0; smallchan=0; smallcoup=1d-12
         ips=0; it0=0; iter=0; fatal=.true.; iblock=0; pade=0; 
	 psiren=.false.; iso=' '; nnu=0; maxl=0; minl=0; mtmin=0; 
         epc=0; erange=0; dk=0;  nosol=.false.; mpihelp = 1
	 nrbases=0; nrbmin=0; pcon=0; pralpha=.false.; rmatr=0; nlagcc=0
         buttle=0; meigs=0; bndx(:)=0.; weak=1e-3; llmax=-1
  	 vsearch=0; echan=0; enodes=0; 
         chans=0; listcc=0; treneg=0; cdetr=0; smats=0; xstabl=0; 
	 nlpl=0; waves=0; lampl=0; veff=0; kfus=0; wdisk=0
	 cdcc=0; nfus=0; bpm=0; melfil=0; btype = ' '; KRM=0
	 inh = 0; sumform=2; maxqrn=0; mclist=0; initwf=0
         ompform = 0; gscouplonly=.false.; pluto(:)=0; eobs=.false.
	 dry = .false.; boldplot=.false.; tcfile=0; maxleg=-1
         tcfilename = 'fort.5420'
         elpmax = 0.01; relref=1
	 TMP='/tmp/'; grace=.true.
	 MASFIL = '/netopt/PHS1IT/lib/m88lrb'
	 hktarg = 0.2

!				NEW-STYLE CONSTANTS:
	 unitmass = 1d0
	 finec = 137.03599d0	! 1/alpha (fine-structure constant)
!				OLD-STYLE CONSTANTS:
!	 unitmass=1.007335d0; FINEC=137.5648d0
     	 pel=0;exl=0;lab=0;lin=0;lex=0
	 elab=0; nlab=0
	 return
	 end subroutine defaults


      FUNCTION W3J(A,B,C,D,E,F)
      IMPLICIT NONE
      REAL W3J,A,B,C,D,E,F
      REAL*8 WIG3J
      W3J=real(WIG3J(DBLE(A),DBLE(B),DBLE(C),DBLE(D),DBLE(E),DBLE(F)))
      W3J=abs(W3J)
      RETURN
      END
      FUNCTION W6J(A,B,C,D,E,F)
      IMPLICIT NONE
      REAL W6J,A,B,C,D,E,F
      REAL*8 WIG6J
      W6J=real(WIG6J(DBLE(A),DBLE(B),DBLE(C),DBLE(D),DBLE(E),DBLE(F)))
      W6J=abs(W6J)
      RETURN
      END
      FUNCTION W9J(A,B,C,D,E,F,G,H,Z)
      IMPLICIT NONE
      REAL W9J,A,B,C,D,E,F,G,H,Z
      REAL*8 WIG9J
      W9J=real(WIG9J(DBLE(A),DBLE(B),DBLE(C),DBLE(D),DBLE(E),DBLE(F)
     >,DBLE(G),DBLE(H),DBLE(Z)))
      W9J=abs(W9J)
      RETURN
      END
      
        subroutine readvars(nparameters,koe,rterms)
	use parameters
	use factorials
	use io
	use searchpar
	implicit real*8(a-h,o-z)
	integer kp,pline,col,nafrac,channel,q,par,nparameters,
     x          term,type,maxrank,points,pline2,col2,dataset(2)
	real*8 jtot,energy,afrac,potential,width,ratio2,jch
	real*8 value,step,valmin,valmax,srch_error(mvars), stepi,error
	character*10 name
	character*5 wid,red
        character psign(3)
	logical rterms,nopot,rwa
        data psign / '-','?','+' /

	namelist /variable/ name,kind,step,valmin,valmax, nul,
     x                      kp,pline,col,potential, dataset,datanorm,
     x     	            nafrac,afrac, energy,jtot,par,channel,width,
     x			    term,nopot,pline2,col2,ratio2,rwa,ivar,B,
     x                      icch,iach,lch,jch,sch,damp,dataEshift
	
        mterms = 0
       	if(nparameters<=0) go to 99
	nvars = nparameters

	nul = -124578   ! 'undefined'
        srch_B(:) = nul
!
!####### Read in specification of search parameters
       	write(koe,'(/''  Define SEARCH VARIABLES:''/)')
        i0 = ichar('0')
	rterms=.false.
	stepi = 0.01
        do ip=1,nparameters
        name='Var'//char(mod(ip/10,10)+i0)//char(mod(ip,10)+i0)
        kind=0;valmin=0;valmax=0; kp=0;pline=0;col=0;
	pline2=0; col2=0; ratio2=0.0; damp=0.; rwa=.true.
        potential=nul; step=stepi; width=0; rwa=.true.; B=nul
        nafrac=0;afrac=nul; energy=nul;jtot=0;par=0;channel=0;term=1
        dataset(1)=1;dataset(2)=0; datanorm=1.0;    nopot=.false.
        icch=1; iach=1; lch=-1; jch=-1.0; sch=-1.0; dataEshift=0.
        read(ki,nml=variable)

        srch_kind(ip) = kind; srch_name(ip) = name

         if(kind==1) then !   potential parameter
 	  write(koe,1010) ip,name,kp,pline,col
1010	  format('   Variable',i3,'=',a10,' is potential KP=',i2,
     X       ', line=',i2,' col=',i2)
          srch_value(ip) = potential; 
          srch_kp(ip) = kp; srch_pline(ip)=pline; srch_col(ip)=col
          srch_pline2(ip)=pline2; srch_col2(ip)=col2
          srch_ratio2(ip) = ratio2;
          if(pline2*col2>0) then
            write(6,1011) ip,name,kp,pline2,col2,ratio2
1011        format('   Variable',i3,'=',a10,' in potential KP=',i2,
     X       ' is linked to variable2 at line=',i2,' col=',i2,
     x       ' with ratio',f8.4)
          endif


         else if(kind==2) then ! spectroscopic amplitude
          write(koe,1012) ip,name,nafrac
1012	  format('   Variable',i3,'=',a10,' is Afrac #',i3)
          srch_value(ip) = afrac; srch_nafrac(ip) = nafrac; 

         else if(kind==3) then ! R-matrix energy
 	  write(koe,1013) ip,name,term,jtot,psign(par+2),nopot
1013	  format('   Variable',i4,'=',a10,' is energy of R-matrix term'
     X       ,i4,' at J/pi =',f5.1,a1,' [NoPot=',L1,']')
          srch_value(ip) = energy; srch_rterm(ip) = term
	  srch_jtot(ip) = jtot; srch_par(ip) = par
	  srch_nopot(ip) = nopot
	  rterms=.true.

         else if(kind==4) then ! R-matrix width
          if(channel>0)  then
           write(6,1014) ip,name,term,channel,1000.0*width**2
1014       format('   Variable',i4,'=',a10,' is width of R-matrix term',
     X       i4,' in channel',i3,' (rw =',f10.2,' keV)')
           srch_r_ch(ip) = channel
          else
!           if(jch<0.0) jch = lch
           if(rwa) then
              wid = 'r.w.a';  red='r'
              widk = 1000*width**2
            else
	      wid = 'width';  red=' '
              widk = 1000*width
            endif
           write(6,1015) ip,name,wid,term,icch,iach,lch,jch,sch,red,widk
1015       format('   Variable',i4,'=',a10,' is ',a5,' of R-matrix term',
     X       i4,' in channel ic=',i1,' ia=',i2,' l=',i2,' j=',f5.1,
     x       ' s=',f5.1,' (',a1,'w =',f10.2,' keV)')
           srch_r_ic(ip) = icch
           srch_r_ia(ip) = iach
           srch_r_lch(ip) = lch
           srch_r_jch(ip) = jch
           srch_r_sch(ip) = sch
          endif

          srch_value(ip)=width; srch_rterm(ip)=term; srch_rwa(ip)=rwa

         else if(kind==5) then ! dataset normalisation
          write(koe,1018) ip,name,dataset
1018	  format('   Variable',i3,'=',a10,
     X       ' is normalisation for dataset ',2i3,' >> IGNORED')
          srch_value(ip) = datanorm; srch_datanorm(:,ip) = dataset(:)

         else if(kind==6) then ! dataset energy shift up
          write(koe,10182) ip,name,dataset
10182	  format('   Variable',i3,'=',a10,
     X       ' is energy shift up for dataset ',2i3,' >> IGNORED')
          srch_value(ip) = dataEshift; srch_datanorm(:,ip) = dataset(:)
          srch_damp(ip) = damp;        srch_B(ip) = B
        endif

         if(abs(srch_value(ip)-nul)>.001) then
          write(koe,1019) srch_value(ip)
1019	  format('     value ',f10.4)
	 else
          write(koe,1020)  
1020	  format('     value from Fresco input')
	 endif
	mterms = max(mterms,term)
       enddo

99     allocate(rm_Brune(0:mterms),E_Brune(0:mterms),W_Brune(0:mterms))
       rm_Brune(:) = .false.

       end
