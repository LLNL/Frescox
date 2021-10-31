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
****frxx1.f*********************************************************
	MODULE arrays

C    FORM FACTORS AND THEIR PARAMETERS
C    ---------------------------------
      COMPLEX*16,allocatable:: FORMF(:,:),ELCOEF(:),FORMC(:,:,:)
      COMPLEX*16,allocatable:: BSMAT(:,:,:,:)
      REAL*8,allocatable:: FORML(:,:,:),FORMDEF(:,:),RMAS(:)
      REAL*8,allocatable:: AFRAC(:,:,:,:),OO(:,:),EVEC(:,:)
      REAL*8,allocatable:: CCFRAC(:,:,:)
C
C    PARTIAL WAVES AND THEIR PARAMETERS
C    ----------------------------------
      COMPLEX*16,allocatable:: CLIST(:,:,:),SRC(:,:)
      COMPLEX*16,allocatable:: PSI(:,:),WNM(:,:,:),FNC(:,:)
      integer,allocatable:: NFLIST(:,:,:),NCLIST(:,:)
      integer,allocatable:: levelp(:),levelt(:)
      CHARACTER*1,allocatable:: INCOME(:)
C   new format:
      COMPLEX*16,allocatable:: CLIS(:)
      INTEGER,allocatable:: NFLIS(:),PFLIS(:,:)
      INTEGER MFLI,NFLI

C
C    COULOMB FUNCTIONS
C    -----------------
      COMPLEX*16,allocatable:: CH(:,:)
      REAL*8,allocatable:: CFMAT(:,:,:),CGMAT(:,:,:),CRCRAT(:)
C
C    SOLVING THE COUPLED EQUATIONS
C    -----------------------------
      COMPLEX*16,allocatable:: SMAT(:),SOMA(:,:),SPAD(:,:,:),SRATIO(:),
     X		FIM(:,:),FIMD(:,:)
      complex*16,allocatable:: ferw(:)
      INTEGER,allocatable:: LVAL(:),PART(:,:),EXCIT(:,:),CUTVAL(:),
     X        INITL(:),CHNO(:,:),PTYPE(:,:),JPWCOUP(:,:)
      LOGICAL,allocatable:: finishedcc(:),READWF(:)
      LOGICAL,allocatable:: BLOCKD(:),SAME(:),EMPTY(:),
     X			    SKIP(:,:),SIMPLE(:)
      REAL*8,allocatable:: BPHASE(:,:,:),EXCH(:,:),alphas(:,:,:),LAM(:),
     x                     GDER(:),LAMERR(:),PWCOUP(:,:),
     X			   FIMR(:,:,:),FIMDR(:,:,:),TCOEFS(:,:)

      REAL*8,allocatable:: XS(:),SIGFUS(:,:),TCOEF(:,:)
      COMPLEX*16,allocatable:: FUSUM(:,:)
      REAL*8,allocatable:: JVAL(:),JPROJ(:),JTARG(:),ECM(:,:),LL1(:)
	End MODULE arrays
	

****fr*********************************************************
      subroutine fr
	use parameters
	use factorials
	use drier
	use kcom
	use trace
	use io
	use fresco1
	use gails, ipmax=>numax, cff=>cf
	use parallel
	use searchdata
	use searchpar
	use arrays
!$      use omp_lib
#ifdef MPI
        use mpi
#endif /* MPI */
!      implicit real*8(a-h,o-z)
      implicit none
C
C    PARTIAL WAVES AND THEIR PARAMETERS
C    ----------------------------------
      COMPLEX*16 PHI(MAXN,2),VOPT(MAXNLN),VPOT(MAXN,3)
      integer NCHPART(MXP),CHBASE(MXP),nphases,linel(100)
      real*8 phaze(40),phasin(100),phaslst(40),phaseadd(40)
      COMPLEX*16 CHL(LMAX1,MXPEX,2),SMATI
      LOGICAL DERIV,INITWFE
      logical, EXTERNAL::refer
      REAL*8 CF,CG,R,RH, JTOTALI,JVALI,JVALE,HPI,ENLABI,ENCOM
      integer NAI,PARITYI,LVALI,LVALE,ITI,NRF
C
C    SOLVING THE COUPLED EQUATIONS
C    -----------------------------
      COMPLEX*16 C6,C7,S,CI,CRho,CEta,CF2
      INTEGER NEX(MXP+1),BAND(2,MXP,MXX),COPY(2,MXP,MXX,2),
     X        ITC(MXP,MXX),CPOT(MXP,MXX),M1,NBL,LPMAX(MXP+1)
      INTEGER QNF(19,MSP),BHEAD(2,2,0:9), NSTEPD
      INTEGER CP,ICTO(MAXCPL+1),ICFROM(MAXCPL+1),KIND(MAXCPL+1),
     X        LOCF(MAXCPL),FPT(7,MAXQRN,MAXCPL),NKP(2,MAXCPL),
     X        IREM(MAXCPL+1),KPCORE(MAXCPL+1),ICOR(MAXCPL,2),
     X        GPT(2,MAXQRN,MAXCPL),FILE(MAXCPL),NFI(3),
     X        NLL(MAXCPL),KLT(MAXCPL),NIB(MAXCPL),
     X        QQ(MAXCPL+1),ICOM(MAXCPL,2),PARITY,C,C1,C2,
     X        IP4(MAXCPL+1),APFLG(5),NLN1,ip,id
	real*8 AERR,GMAX,DEGENY,APEPS,WERRMAX,WERR,DIFF,FRACT
	real*8 WID(MAXCPL+1),CENTR(MAXCPL+1),GVAL,Shift,Pen
	logical prwf
        INTEGER EL,CUTOFF,NCH,NQ,NF,IERFLG,NPAR,NCH2,C2I,C2F,MAXLAM0
	LOGICAL SHSMAT,SMATEL,FJSWTCH
C
C    CONTROL VARIABLES
C    -----------------
      INTEGER DONE,DONES,IOFAM(2,MXP,MXX),KOO,MDONE
      LOGICAL GIVEXS(MXP),EXTRA(2,MXP,MXX),FLIP
      LOGICAL LOCAL(MAXCPL),REV(MAXCPL),WREQ,COUPLE(MAXCPL),FALLOC,
     X        NPRIOR,NPOST,OPN(2),ffreal(MAXCPL),REW14,OPEN12,WOUT,WOUT0
      LOGICAL LTRANS(MAXQRN,MAXCPL),MLCALL(2,MAXCPL),CPSO(MAXCPL),
     X        SCHONP,REPEAT,OPEN18,STOP19
      LOGICAL FAIL,DISC8,DISCIT,PWFLAG(MXP),SMALLJ,HERE,SHFL,OP
	integer NF0,LAMAX,NC,NCLREQ,F51,MIXPOT(MXP+1)
	integer NCP,IEXCH,NSA,NJA,NSB,NJB
        real*8 TCC,TCC0,TIME0W,TIME1W,JUMPER,RELS
        real*8 R0,XCOEF,TIME0J,RASYMAX,RTURN1,EE,AL,BEST,RERR,R1
	real*8 T,TH1,TH2,ALOSSM,TSTD,RENORM,RTURN,RS,GAP,AMSA,XLMAX,X
	integer JPSET,LUSED,JCCSET,IEX,MINTL,MEL,NTHRESH,MCCSET
	integer NIX,IEN,ILEN,NCHAN,MAXF,DROPPED,LEN,NLEN
	integer IC,IA,NA,MAM,IF1,IF2,ITCM,IMA,JBLOCK,L,IPARIT,KS
	integer JF,JIN,IAM,J,NLJ,IN,IPUT,IGET,LAP,IF,I,IT
	integer KN,NSP,L1,IE,NPROJ,NTARG,KINTL,nparit,JUMP2
	integer NICH,MAXXCH,IB,JFT,JFTT,NFL,NR,IR0,IR,IC1
	integer IC2,IN1,ICH,KP,KPB,IX,KK,II,NBEST,NSOLI,ITNL,IPARI
	integer NJDONE,LL,JJ,I0,nbas,CHANSI,JF0,NFD,NEXK,KPA,KPA2
	integer LAST,J2LAST,J2,KNLAST,KPB2
	integer I1,I2,IX1,IX2,LAM1,LAM2
        integer NJTOTAL,IJTOTAL
C
C    DEFINING THE MASS PARTITIONS AND THEIR EXCITED STATES
C    -----------------------------------------------------
      PARAMETER(NTHRESH=5)
      REAL*8 MASS(4,MXP+1),RMASS(MXP),HP(MXP),JEX(6,MXP,MXX)
      REAL*8 ENEX(2,MXP,MXX),QVAL(MXP+1),T4,RESM(NTHRESH),THRJ(NTHRESH)
      REAL*8 BE(MSP,5),D0(MSP,2),BPROJ(MSP,MXP),HPOT(MLOC),BEPROJ,BETG
      REAL*8 ENEX0(2,MXP,MXX),BETARG
C
C    BPM and VEFFPOT
C    ---------------
      REAL*8 BARE,TLE,SIGL,HOM,EKL,EKL1,DER,RSO,POS,BAR,RTURN2
      REAL*8 FUSIO,RSOO
C
C    COUPLINGS
C    ---------
      REAL*8 BETAR(MAXCPL+1),BETAI(MAXCPL+1),JMAX(MAXCPL+1),SCALEL,
     X	     RMAX(MAXCPL+1),MEK(MPAIR),SPINTR(2,MPAIR),VARYL(2,MLOC),
     X       XA(MAXCPL),XB(MAXCPL),XP(MAXCPL),XQ(MAXCPL),STREN(MLOC),
     X       EMID(MSP),KMINX(2,MSP),SPINV(2),MASSV(2),CHRGV(2),
     x       TENSOR,TENSOR2,ETAS,ETATG,ETAPR,ETATG2,WPA,WPB,WPOT
      INTEGER MATRIX(6,MPAIR),LAMBDA(MLOC),LSHAPE(MLOC),LDEP(MLOC),
     X   NKBIN(2,MSP),IBIN(MXX),NBINS,NNJMAX,PARV(2),NORBIN(MSP),
     X   CHSIGN(MXX),BSIGN(0:MSP),NFUSCH,mag,POTCAP(MLOC,2,MAXCPL),
     x   INFILE(MAXCPL),NPWCOUP
	LOGICAL TWOW,ISNONO,J40P,J14T16
C
C    INCOMING ENERGIES
C    -----------------
      REAL*8 K(MXP,MXX),ETA(MXP,MXX),ETOTAL,ENLAB,EOFF,EPOLE,ETOTALR,
     X     RMK,RMKD,CFG(MAXMUL,4),DE,ELAST,CSIG(LMAX1,MXPEX),
     x     GAM(MXP,MXX),ECMC(MXP,MXX),dataEshift
      REAL*8 JTOTAL,JAP,JSWITCH,JAL,JN,LJMAX,SSWITCH,JNLAST
C
C    ARRAYS FOR NON-LOCAL COUPLING FORM FACTORS
C    ------------------------------------------
      REAL*8 DNL(MAXNLO)
      REAL*8 NLOC(NLM,MAXCPL)
      COMPLEX*16 FFC
      EXTERNAL FFC
C
C    CALCULATION OF CROSS SECTIONS
C    -----------------------------
      REAL*8 JCOEF,TOTFUS(3),CORFUS(3,NFUS1),SIGT(3),FUSL(1+NFUS1),
     x       SIGTOT(3),SIGEL(3),datanorm,FUSLL(LMAX1,2),FUSLLT(2)
      REAL*8 SIGJ(0:MXPEX),SIGR(MXP,MXX),OUTJ,SFAC(MXPEX),rtrace,
     X		CHSIZES(MXPEX),XSC,AMDSQS,FUSJ,TOTJ,CFUSJ(NFUS1),sfa,
     x          strength(0:2),R0PSTR,RSCATTR,SIGR2(MXP,MXX),SIGREAC,
     x          ELAS,SIGELE(MXPEX),DSPINS(MXX),FRACFUS,SIGCH
      integer LCROSS,LFAM,LXSEC,NANGL,SMALLS(MXPEX),NSMALL,
     x        CHPRES(MXPEX),LEG,MAXPLM,DMULTIES(MXX),NMULTIES,LGAM,
     x        DLEVEL(MXX)
C
    	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN  
      COMMON /TIMER/ START
      real*8 START,TIME1,SECOND,memuse
      real*4 TM,TM0,TSYNC
      real*8 Z,HALF,EPS,DDEXP
C
      CHARACTER*80 CHME*3
      CHARACTER*10 COMMAND
      CHARACTER*8 NAME(2,MXP+1),NAMEV(2),DNAME
      CHARACTER*4 W
      CHARACTER*1 SCALE(0:MXPEX),PSIGN(3),COMPARE
C
      integer NFNL,MMXF,MMXCH,MMXQRN,MMXIT,MAL1,MMXICH,MAXCHT
      DATA NFNL,MMXF,MMXCH,MMXQRN,MMXIT,MAL1,MMXICH/7*0/
      DATA PSIGN / '-','?','+' /
      DATA RESM / 0.01, 0.1, 0.5, 0.90, 0.99 /
      parameter(mag=1)
   
	namelist/nlparams/ mxx,mxp,mxpex,maxit,maxcpl,maxnnu,
     X		maxn,maxm,maxnln,maxnlo,msp,lmax1,mpair,nkmax,
     X          mloc,maxch,maxqrn,mfnl,mclist,maxpw,
     x   	inh,maxich,maxcch,melfil,
     x		unitmass,finec,fmscal,coulcn,amu,etacns,
     X  	dry,psiren,rintp,epc,
     X		nl,maxl,mlt,nln,nlo,nlc,
     X          n,mr,minl,jtest,nlt,nlm,nnu,
     X		chans,listcc,treneg,cdetr,smats,smatl,xstabl,nlpl,waves,
     X         	lampl,veff,kfus,wdisk,bpm,cdcc,nfus

	namelist/nlfresco1/ headng,iso,rela,TMP,
     x	hcm,rmatch,hnl,rnl,centre,hnn,rnn,rmin,rasym,accrcy,rsp,
     X   switch,ajswtch,sinjmax,jtmin,jtmax,absend,erange,dk, 
     X   thmin,thmax,thinc,cutl,cutr,cutc,ips,jbord,elab,
     X 		nearfa,jump,koords,kqmax,pp,isocen,nnn,ngail,
     X		it0,iter,iblock,pade,mtmin,nlab,m,md,lmax,meigs,
     X		pel,exl,lin,lab,lex,nrbases,nrbmin,pcon,mrm,plane,
     X		fatal,nosol,fcwfn,buttle,pralpha,symm,locfil,mcalls,
     x  	rterms,rin,bndx,rmatr,smallchan,smallcoup,!@ gailacc,
     x		mlm,nlcn,mint,mintm2,nj,jset,pset,icutc,itmin,nsol,
!    X   	scalx,rscalx,
     X		nforml1,nforml2,weak,allpart,number_calls,
     x          sumccbins,cxwf

      namelist/nlsearch/nvars,srch_kind, srch_name, srch_value,
     x  srch_step,nul,
     x	srch_minvalue,srch_maxvalue,srch_kp,srch_pline,srch_col,
     x  srch_pline2,srch_col2,srch_ratio2,
     x	srch_nafrac,numafrac, srch_afrac_overlap, 
     x	srch_jtot,srch_par,srch_r_ch, srch_datanorm,
     x	datasets,datalen, datangles,datavals,dataerr,
     x	data_lab,data_matched,data_ic,data_ia,data_thfile,ndof,
     x	data_idir,data_rank_k,data_rank_q,data_type,data_energies

C scalars:
	namelist /nlreadin/ WOUT,WOUT0,NF0,NF,ISNONO,FLIP,
     X     NCHAN,ITCM, NIX ,NSP,
     x  NFI,TWOW,OPN,NPRIOR,NPOST,REW14,IEXCH,
     x  NSA, NJA, XCOEF, DISCIT, DISC8, LAMAX
C variable arrays:
!	namelist /nlreadin/ STREN,
!     X     NAME,MASS,NEX,QVAL,RMASS,HP,COPY,EXTRA,IOFAM,GIVEXS,
!     X     JEX,CPOT,ENEX,BAND,ITC,
!     X     PWFLAG,FORMF,FORML,FORMC,PTYPE,FILE,ffreal,
!     X     HPOT,LAMBDA,MATRIX,MEK, QNF,D0,BE,AFRAC,
!     x  ICTO,ICFROM,KIND,QQ,IREM,KPCORE,BETAR,BETAI,JMAX,RMAX,REV,
!     x  LOCAL,NIB,LTRANS,FORMDEF,
!     x  MLCALL,BPROJ,XA,XB,XP,XQ,SPINTR,FPT,NKP,VARYL,LSHAPE,LDEP,
!     x  NLL,LOCF,KLT,COUPLE,ICOM,ICOR,NCP,GPT	

      DDEXP(T) = EXP(MIN(100D0,MAX(-100D0,T)))
      TM(I) = real(SECOND() - TIME0)
C
      START = 0.0
      X=0d0
      TIME0 = SECOND()
!$    TIME1W = OMP_GET_WTIME()
	TM0=0.
      I=0; T = dble(TM(I))
!      write(KOI,*) 'Starting  @ ',real(T)
	koo = 93
        TSYNC=0.0
C
      Z = 0.0
      HALF = 0.5D0
      CI = (0D0,1D0)
      EPS = 1E-10
      R0PSTR=1.35 ! fm
!	scalx = 1d-100; rscalx = 1d0/scalx

      NFNL=0 ; MMXF=0 ; MMXCH=0 ; MMXQRN=0 ; MMXIT=0 ; MAL1=0
      NCLREQ=0; LUSED = 0; MMXICH=0 ; MAXXCH=0
      NIX=1; NPWCOUP = 0
	CHSIGN(:) = 0
	lanecoup = .false.
	numafrac=0
	rewind 3
	call rewop(4)
	if(KO/=6) rewind KO  ! keep only final Fresco run!
	if(rterms.and.final) call rewop(60)
	if(rterms.and.pralpha.and.final) call rewop(61)
	NANGL = (abs(THMAX) - THMIN)/abs(THINC) + 1 + 0.5	
      data_chisq(:) = 0.0
      STOP19 = .false.
	allocate (PTYPE(12,MLOC))
C
      if(final.and.eigens==0) then
	   inquire(7,opened=op); if(op) close(7)
      OPEN(7,form='formatted',file='fort.7',recl=110)
!      OPEN(80,form='formatted',file='fort.80',recl=90)
	   inquire(16,opened=op); if(op) close(16)
      OPEN(16,form='formatted',file='fort.16')
      call rewop(16)
      call rewop(35)
      call rewop(7)
      endif
C
**MPI*******************************************************************
#ifdef MPI
c     MPI inititialized in fresco.f
      if(MPIC)then
!	write(500+iame,*) ' FRXX1 call barrier at start'
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!       write(6,*)' This is ',iame+1,' of ',mpisets,' MPI threads'
	allocate (finishedcc(mpinodes))
      endif
#endif /* MPI */
************************************************************************
       I = ichar('0')
      CHME = char(I+mod(iams/100,10))//char(I+mod(iams/10,10))//
     X       char(I+mod(iams,10))
      TMPN = trim(TMP)//'fort.'//CHME//'.'
!      TMPN = 'fort.'//CHME//'.'   ! for debugging: files to current directory
!      TMPN = trim(TMP)//'/node.'//CHME//'/fort.'
**MPI*******************************************************************
#ifdef MPI
      IF(MPIC)THEN
        IF(KO.ne.KOI) THEN
          OPEN(KO,FILE=trim(TMPN)//'51',form='formatted')
          T = TM(I)
          write(KO,*) ' Node ',iame,' started, o/p file ',KO,T
        else
           OPEN(51,FILE=trim(TMPN)//'51',form='formatted')
        ENDIF
	F51 = 51
	F51 = 500+iame ! temp debugging
        OPEN(50,FILE=trim(TMPN)//'50',form='formatted')
        OPEN(F51,FILE=trim(TMPN)//'F51',form='formatted')
C   The following distributed files are opened elsewhere
*    OPEN(1 OPEN(8 OPEN(9 OPEN(11 OPEN(12 OPEN(14 OPEN(18
        SMATL = SMATS
        if(.not.final.or..not.say) SMATL=0
        if(SMATL>=2)then
         if(iams==0)then
          OPEN(38,form='formatted')
         else
          OPEN(38,FILE=trim(TMPN)//'38',form='formatted')
         endif
        endif
        if(iams==0)then
         OPEN(45,form='formatted')
        else
         OPEN(45,FILE=trim(TMPN)//'45',form='formatted')
        endif
        OPEN(48,FILE=trim(TMPN)//'log',form='formatted')
	written(48) = .true.
        OPEN(49,FILE=trim(TMPN)//'49',form='unformatted')
	written(49) = .true.
        rewind 49
         T = TM(I)
        write(48,*) ' Node ',iame,' started, o/p file ',KO,real(T)
      ENDIF ! if MPI
#endif /* MPI */
!	rewind 48
      if(TRENEG>0) then
	open(34,form='formatted',access='sequential')
	rewind 34
	endif
      OPEN12 = .false.
      OPEN18 = .false.
	INFILE(:) = 0
      XA=0.;XB=0.;XP=0.;XQ=0.;MATRIX=0;MEK=0.
        ENEX(:,:,:) = 0.; BE(:,:) = 0.; GPT= 0.; FILE(:)=0;FPT=0
	CPSO(:) = .false.

      if(.not.allocated(AFRAC)) then

      memuse = (MAXM/1024.)*(MLOC/1024.)*(16/1024.)
      if(memuse>.01) write(KO,'(a,i5,i8,a,f7.3,a)')'Allocating FORMF(',
     &MAXM,MLOC,') in ',memuse,' GB' 
      call flush(KO)

      allocate(FORMF(MAXM,MLOC),ELCOEF(MLOC),FORMDEF(12,MLOC))
      allocate(RMAS(MSP))
      if(.not.cxwf) then
        allocate(FORML(MAXNLN,MSP,2),FORMC(1,MSP,2))
        maxnlr = MAXNLN; maxnlc = 1
	FORML(:,:,:) = 0.; FORMC(:,:,:) = 0.
        endif
      if(     cxwf) then
          allocate(FORMC(MAXNLN,MSP,2))
          maxnlc = MAXNLN
          if(nn2wf) then
                allocate(FORML(MAXNLN,MSP,2))
                maxnlr = MAXNLN
              else
                allocate(FORML(1,MSP,2))
                maxnlr = 1
          endif 
	FORML(:,:,:) = 0.; FORMC(:,:,:) = 0.
        endif

      nforml1 = maxnlr; nforml2 = MSP

      if(.not.sumccbins)then
      
!$    memuse = (MXPEX/1024.)*(MXPEX*2/1024.)*(MSP*8/1024.)
!$    write(KO,'(a,2i5,i2,i5,a,f7.3,a)')
!$   &'Allocating AFRAC(',MXPEX,MXPEX,2,MSP,') in ',memuse,' GB' 

      allocate(AFRAC(MXPEX,MXPEX,2,MSP))
      allocate(CCFRAC(1,1,1))
      
      else

!$    memuse = (MXPEX/1024.)*MXPEX*8/dble(1024**2)
!$    write(KO,'(a,i7,i2,a,f7.3,a)')'Allocating CCFRAC(',
!$   &MXPEX,2,') in ',memuse,' GB' 
      
      allocate(AFRAC(1,1,1,1))
      allocate(CCFRAC(MXPEX,2,MXPEX))
      endif
      if(MPWCOUP>0) allocate(PWCOUP(3,MPWCOUP),JPWCOUP(9,MPWCOUP))
      if(MPWCOUP==0)allocate(PWCOUP(3,1),JPWCOUP(9,1))

      if(NKMAX>=0) then
	allocate(BPHASE(2,max(NKMAX,1),MSP))
        BPHASE(:,:,:) = 0.
        if(CDCC>1) then
 	  allocate(BSMAT(NCHBINMAX,NCHBINMAX,max(NKMAX,1),MSP))
         else
          allocate(BSMAT(NCHBINMAX,NCHBINMAX,1,1))
	 endif
	endif

      endif
 
C
C   READ IN MAIN INPUT VARIABLES
C   ----------------------------
C
	CALL READIN(WOUT,WOUT0,NF0,NF,STREN,ISNONO,LPMAX,MIXPOT,
     X     NAME,MASS,NEX,QVAL,RMASS,HP,COPY,EXTRA,IOFAM,GIVEXS,
     X     JEX,CPOT,ENEX,BAND,ITC,FLIP,NCHAN,ITCM,PSIGN,BHEAD,
     X  PWFLAG,FORMF,FORML,FORMC,PTYPE,CPSO,EPS,MMXQRN,FILE,ffreal,
     X  HPOT,LAMBDA,MATRIX,MEK,NIX ,NSP, QNF,D0,BE,AFRAC,CCFRAC,BPHASE,
     x  ICTO,ICFROM,KIND,QQ,IREM,KPCORE,BETAR,BETAI,IP4,JMAX,RMAX,REV,
     x  NFI,TWOW,OPN,NPRIOR,NPOST,REW14,LOCAL,NIB,LTRANS,FORMDEF,
     x  MLCALL,BPROJ,XA,XB,XP,XQ,SPINTR,FPT,NKP,VARYL,LSHAPE,LDEP,
     x  NLL,LOCF,KLT,COUPLE,ICOM,ICOR,NCP,IEXCH,GPT,POTCAP,INFILE,
     x  NSA, NJA, XCOEF, DISCIT, DISC8, LAMAX,NBINS,EMID,KMINX,NKBIN,
     x  NORBIN,STOP19,WID,CENTR,NPWCOUP,PWCOUP,JPWCOUP,BSMAT,RMAS)
!     if(.not.cxwf) write(KO,*) 'WF:',FORML(1:5,2,1)
!     if(.not.cxwf) write(KO,*) 'PT:',FORML(1:5,2,2)
!     if(     cxwf) write(KO,*) 'WF:',FORMC(1:5,2,1)
!     if(     cxwf) write(KO,*) 'PT:',FORMC(1:5,2,2)
	J40P = .false.; J14T16=.false.; MAXLAM0 = 0
	 DO JF=1,NF0
	  if(PTYPE(5,JF)>=40) J40P=.true.
	  if(14<=PTYPE(2,JF).and.PTYPE(2,JF)<=16) J14T16=.true.
	  MAXLAM0 = max(MAXLAM0,PTYPE(3,JF))
	 enddo
	if(say.and.final)
     x    write(48,*) ' MAXLAM0,J14T16,J40P =',MAXLAM0,J14T16,J40P
	if(NOSOL) WOUT=.false.
	FALLOC = WOUT.or.ITER>0
	LOCFIL = .not.rterms 
	LOCFIL = LOCFIL .and. .not.FCWFN   ! eventually fix PWISER!
	LOCFIL = LOCFIL .and. NF>0		! set some efficiency minimum
	LOCFIL = LOCFIL .and. BPM==0.and.VEFF==0.and.NFUS==0
     x		        .and. melfil==0 .and. .not.STOP19
	LOCFIL = .false.			! ALWAYS SEEMS TO SLOW, NOT SPEED UP!
	if(LOCFIL) write(KO,*)  ' Putting FORMF in file 19'
	 NLOC(:,:) = 0.0
       DO CP=1,NCP
       IF( (KIND(CP)==7.OR.KIND(CP)==8.or.
     X     (KIND(CP)==1.and.KPCORE(CP)<0).or.
     X     (KIND(CP)==2.and.IP4(CP)>0).or.
     x     (KIND(CP)==9.and.QQ(CP)==1))      .AND..NOT.OPEN12) THEN
	   inquire(12,opened=op)
	   if(op) then
	     WRITE(KO,*) ' FILE ',12,' ALREADY OPEN!!?'
	     endif
           IF(MACH.eq.3.or.MACH.eq.4) THEN
!	   write(F51,*) 'Open file 12 as ',trim(TMPN)//'12'
!	   call flush(F51)
           OPEN(12,ACCESS='SEQUENTIAL',FILE=trim(TMPN)//'12',
     X           FORM='UNFORMATTED',status='unknown')
             ELSE
           OPEN(12,ACCESS='SEQUENTIAL',STATUS='SCRATCH',
     X           FORM='UNFORMATTED')
             ENDIF
           OPEN12 = .TRUE.; written(12) = .true.
           ENDIF
	ENDDO
        T = TM(I)
        if(final.and.iams==0)
     >   write(KOI,*) 'Finished all Couplings @ ',real(T)
!        write(KOI,*)'NF0,NF,OPEN12,NCHAN,EPS,MMXQRN,NSP,NNN,
!     X               NFI,TWOW,OPN,REW14,NLL,IEXCH,NCP,NSA,NJA,
!     X		     XCOEF,DISCIT,DISC8,LAMAX='
!     X              ,NF0,NF,OPEN12,NCHAN,EPS,MMXQRN,NSP,NNN,
!     X               NFI,TWOW,OPN,REW14,NLL,IEXCH,NCP,NSA,NJA,
!     X		     XCOEF,DISCIT,DISC8,LAMAX

	if(final) then
	if(symm) write(ko,'(/'' Symmetric Hamiltonian'')')
	if(.not.symm) write(ko,'(/'' Non-symmetric Hamiltonian'')')
	if(PLANE>0) write(ko,'(/'' Plane wave treatment ='',i4)') PLANE
	endif
	ryev = 1.
      	if(RMASS(PEL)>1e-5)ryev = 1./(FMSCAL*RMASS(PEL))
	ipmax = LAMAX+1
 	if(say.and.final)write(48,*) ' ipmax =',ipmax
	
	if(melfil/=0) then
	   IC = PEL
           allocate (levelp(NEX(IC)),levelt(NEX(IC)))	
           ECMC(PEL,EXL) = ELAB(1) * 
     x                MASS(3-LIN,LAB) / (MASS(2,LAB)+MASS(1,LAB))
	   call WRITESPEC(LAMAX+1,HEADNG,RMASS,PEL,NLN,NAME,symm,
     X	    RINTP,NEX,COPY,MASS,JEX,BAND,ENEX,ntarg,nproj,levelp,levelt,
     X      ECMC(PEL,EXL),LMAX+1)
	  endif

**MPI*******************************************************************
        if(MPIC)then ! if MPI
! no sending of messages necessary, each thread initialized already
      IF(LISTCC.GT.0.and.LISTCC.le.90) LISTCC = LISTCC - mpisets
        endif ! MPI
************************************************************************

	if(LOCFIL) then
	  DO I=1,MAXM
	  open(19,form='unformatted',status='scratch')
	  write(19) (FORMF(I,IF),IF=1,NF)
	  enddo
	  deallocate(FORMF)
	  allocate (FORMF(MLOC,1))
	endif
!	if(.not.rterms.and..not.mcalls) then
!	   write(48,*)  'FORML(1,1,1) =',FORML(1,1,1)
!	   write(48,*) 'FORML:',allocated(FORML); call flush(48)
!	   deallocate(FORML)
!	   allocate(FORML(1,1,2))
!           nforml1 = 1; nforml2 = 1
!	   write(48,*)  'FORML(1,1,1) =',FORML(1,1,1); call flush(48)
!	   endif
	
         ELAST = 0.
	CHANSI = CHANS
       if(say.and.final.and.(abs(NLAB(1))>1.or.num_energies>0)) then
	  call rewop(71)
	  write(71,'(''@legend ON'')')
	  call rewop(171)
	  write(171,'(''@legend ON'')')
	  endif
	do IC=200,210  ! see IFOUT and IFO in CRISS for the reason for these numbers
	  if(written(IC)) rewind IC
	  enddo
	LFAM=37
!!        if(written(LFAM)) rewind LFAM
!	write(6,*) ' DISCIT, DISC8, WREQ =',DISCIT,DISC8,WREQ
	phaseadd(:) = 0.0
	phaslst(:) = 0.0
	 MDONE = NSA*NJA * 3
	 ENLAB = 0.0
      if(say.and.final)write(48,*) 'MINL,MAXL,MDONE =',MINL,MAXL,MDONE
	LEN = 0; firstE=.true.
      DO 950 IEN = 1,4
      IF(num_energies==0.and.abs(ELAB(IEN)).lt.1e-12) GO TO 960
         IF(NLAB(IEN).EQ.0) NLAB(IEN) = 1
         DE = (ELAB(IEN+1) - ELAB(IEN)) / abs(NLAB(IEN))
         NLEN=abs(NLAB(IEN))
         if(num_energies>0) NLEN=num_energies
      DO 950 ILEN=0,NLEN
       if(num_energies==0) then
           IF(ILEN.GT.0 .AND. ELAB(IEN+1).lt.1e-20) GO TO 950
	     ENLAB = ELAB(IEN) + ILEN*DE   ! linear
           if(NLAB(IEN)<0) then            ! log
	     CF = (log(ELAB(IEN+1))-log(ELAB(IEN)))
     x               /(ELAB(IEN+1)-ELAB(IEN))
             CG =  log(ELAB(IEN)) - CF * ELAB(IEN)
	     ENLAB = exp(CF*ENLAB + CG)
	     endif
           IF(ABS(ENLAB-ELAST).LT.1E-12) GO TO 950
	 else
	   if(IEN>1.or.ILEN==0) go to 950
	   ENLAB = energy_list(ILEN)
           PEL   = pel_list(ILEN)
           LAB = PEL
           NSA = 2*JEX(1,PEL,EXL) + 1.1
           NJA = 2*JEX(2,PEL,EXL) + 1.1
           XCOEF = 10.0/(NSA * NJA)
	 endif
         ELAST = ENLAB
      if(final) WRITE(KO,1260) NAME(LIN,PEL),NAME(LIN,LAB),ENLAB
 1260 FORMAT('1',131('*'),/132('*')//8X,'INCOMING ',A8,' ;',
     X        4X,'LABORATORY ',A8,  ' ENERGY =',
     X G13.5,' MeV.' //       1X,131('*'),/1X,131('*') /)
!     if(.not.final) WRITE(KO,1261) NAME(LIN,PEL),NAME(LIN,LAB),ENLAB
!1261 FORMAT(/'  *** In ',A8,' ; Lab ',A8,' at',G13.5,' MeV.' )
      if(.not.final) WRITE(KO,1262) ENLAB
 1262 FORMAT(/'*** ',G15.7,' MeV.' )
      if(LAMPL>0) WRITE(LFAM-1,1263) ENLAB
 1263 format(' 0Scattering Amplitude: PL coefficients. Elab =',g12.5)
	LEN = LEN+1
	CHANS = CHANSI
!	if(.not.final) CHANS = 0
	call flush(KO)
!	if(KO/=6) write(6,*)  'Scattering at lab energy ',real(ENLAB)
	
	if(allocated(finishedcc)) finishedcc(:) = .false.
C
      EOFF= - QVAL(LAB) + ENEX(1,LAB,LEX) + ENEX(2,LAB,LEX)
      ecmrat = MASS(3-LIN,LAB) / (MASS(2,LAB)+MASS(1,LAB))
      ENCOM = ENLAB * ecmrat
	if(index(rela,'b')>0) then
	 TLE = MASS(2,LAB)+MASS(1,LAB)
	 X = ENLAB*2.0*MASS(3-LIN,LAB) / TLE**2 / AMU
	 ENCOM = TLE * (sqrt(1.0+x)-1.0) * AMU
	 if(abs(x)<0.1) ENCOM = TLE*(X/2.-X**2/8. + X**3/16.) * AMU
	 if(say.and.final)write(48,*) 'ECMNR,TLE,X,ENCOM =',
     x              ENLAB * ecmrat,TLE,X,ENCOM
        endif
      ETOTAL = ENCOM + EOFF
      etarat = ETACNS * MASS(2+1,pel)*MASS(2+2,pel) * SQRT(RMASS(pel))
C
C    CHANNEL ENERGIES
C    ----------------
      MAL1 = MIN(LMAX,INT(JTMAX+20.)) + 1
            XLMAX = MAL1 - 1
	hktarg = min(hktarg,0.20d0)
        T = MASS(1,relref) + MASS(2,relref)   ! use reference partition = 'relref'
        ! RELS = T**2 * AMU + 2.*MASS(3-LIN,relref)*ENLAB ! for rela='h'
        RELS = T**2 * AMU + 2.*MASS(3-LIN,PEL)*ENLAB ! for rela='h'
     x            + 2.*T * EOFF
	FAIL = .false.
	BEST = HCM
        ETOTALR = ETOTAL
      DO 145 IC=1,NCHAN
      NA = NEX(IC)
      DO 145 IA=1,NA
      IT = ITC(IC,IA)
         ECMC(IC,IA) = ETOTAL + QVAL(IC) - ENEX(1,IC,IA)-ENEX(2,IC,IA)
	 GAM(IC,IA) = 1.0
!@@
         IF(RMASS(IC).lt.1e-5) then
            K(IC,IA) = ECMC(IC,IA)/HBC
         ELSE
            EE = FMSCAL*RMASS(IC) * ABS(ECMC(IC,IA))
            K(IC,IA) = SQRT(EE)
	  if(index(rela,'a')>0.or.index(rela,'b')>0) then
	   if(index(rela,'a')>0)    ! Ingemarsson eq(16)
     x	    X = (1. + ECMC(IC,IA)/(2.*AMU*RMASS(IC))) 
     x        / (1. + 2*ECMC(IC,IA)/((MASS(1,IC)+MASS(2,IC))*AMU))
 	   if(index(rela,'b')>0)    ! Ingemarsson eq(17)
!     x	    X = (1. + ECMC(IC,IA)/(2.*AMU*MASS(1,IC)))      ! Antonio queries this when LIN=2
     x	    X = (1. + ECMC(IC,IA)/(2.*AMU*MASS(LIN,IC))) 
     x        / (1. + ECMC(IC,IA)/((MASS(1,IC)+MASS(2,IC))*AMU))
	    K(IC,IA) = sqrt(EE*X)
	    ECMC(IC,IA) = ECMC(IC,IA)*X
	    if(index(rela,'g')>0) then ! Ingemarrson eq(21)   ! needs X from a or b
		GAM(IC,IA)  = X* 
     x	        (1. + ECMC(IC,IA)/(AMU*MASS(LIN,IC))) 
     x        / (1. + ECMC(IC,IA)/(AMU*MASS(LIN,IC)*2.))
	    endif  ! g
	  endif  ! a or b
	  if(index(rela,'f')>0) then   ! independent of a or b
		RS = (MASS(1,IC)+MASS(2,IC))*AMU
		EE = sqrt(RS**2 + 2*RS*ECMC(IC,IA))
		T = (MASS(2,IC)**2 - MASS(1,IC)**2)/(2*EE)*AMU**2
		TH1 = EE*0.5 - T
		TH2 = EE*0.5 + T
		GAM(IC,IA)  = TH1*TH2/EE / (RMASS(IC)*AMU)
!		write(190,*) IC,IA,EE,TH1,TH2,
!     x                       TH1*TH2/EE/AMU,RMASS(IC)
	   endif ! f or g
	  if(index(rela,'r')>0) then   ! getting reduced mass from ratios of total energies
	     th1 = K(IC,IA)*HBC/(MASS(1,IC)*AMU)
	     th2 = K(IC,IA)*HBC/(MASS(2,IC)*AMU)
	     GAM(IC,IA) = (MASS(1,IC)+MASS(2,IC))/
     x        (MASS(1,IC)/sqrt(1.+TH2**2) + MASS(2,IC)/sqrt(1.+TH1**2))
	  endif
          if(index(rela,'h')>0) then   ! get K and ETA from Hale method
             TH1 = MASS(1,IC) + ENEX(1,IC,IA)/AMU ! projectile in this partition
             TH2 = MASS(2,IC) + ENEX(2,IC,IA)/AMU ! target
             ECMC(IC,IA) = (RELS - (TH1+TH2)**2*AMU)/(2.*(TH1+TH2))
             if (IC==relref.and.IA==1) then
              T = ECMC(IC,IA) - ETOTAL
              ETOTALR = ECMC(IC,IA)
              endif
             XSC = (RELS-(TH1+TH2)**2*AMU)*(RELS-(TH1-TH2)**2*AMU)
     x               /(4.*RELS)
             K(IC,IA) = sqrt(abs(XSC)*AMU)/HBC
             EE = TH1*AMU
             GAM(IC,IA) = (1. + ENLAB/EE) / sqrt(1. + ENLAB/(2*EE))
          endif

         ENDIF     ! not gamma
         ETOTAL = ETOTALR
!@@
         ETA(IC,IA) = ETACNS * MASS(2+1,IC) * MASS(2+2,IC)
     X                       * SQRT(RMASS(IC)/ ABS(ECMC(IC,IA)))
     x                       * GAM(IC,IA)
	if(IC.eq.PEL.and.IA.eq.EXL) then
	 if(mod(PLANE,2)==1) ETA(IC,IA) = 0.  ! elastic channel
	else
	 if(    PLANE/2 ==1) ETA(IC,IA) = 0.  ! nonelastic channels
	endif
         SIGR(IC,IA) = 0.0
         NSB = 2*JEX(1,IC,IA) + 1.1
         NJB = 2*JEX(2,IC,IA) + 1.1
            MAM = NSA * NJA * NSB * NJB
            MMXF = MAX(MMXF,MAM)
!		write(161,*) ic,ia,mam,mmxf
!		call flush(161)
 	t = hcm*K(IC,IA)
	NAMEV(1) = 'OK'
	if(t > hktarg) then
	  NAMEV(1) = 'NOT OK'
	  write(ko,1445) ENLAB,IT,t,NAMEV(1),hktarg
 1445 format(/' Accuracy analysis at',f10.3,' MeV in ch.',i3,':',
     x        '  Elastic h*k =',f6.3,' so ',a8,' compared with',f6.3/)
          endif
	FAIL = FAIL .or. t > hktarg*1.2
	  BEST = min(BEST,hktarg/K(IC,IA))
  145   CONTINUE
	if(FAIL) then
	  write(KO,146) HCM,BEST
  146 	format(' Step size HCM',f8.5,' is too large: REDUCE TO <',
     x              f8.5,', so STOP!')
	  stop
	 endif
	MAXF = MMXF
      if(.not.allocated(SIGFUS)) allocate(SIGFUS(MAXF,3))
C
!      IF(FCWFN.and.final.and.JTMAX>100.) THEN
      IF(          final.and.JTMAX>50.) THEN
C              give warning of limits for Coulomb excitations:
       T = 2. * ETA(PEL,EXL) * 180./PI
        TH1 = T / (JTMAX+0.5)
            JTOTAL = -1.
	    X = -1
           IF(RASYM.lt.-0.01) THEN
		X = abs(RASYM)
                JTOTAL = T/(-RASYM)
                RASYM = max(JTOTAL/K(PEL,EXL),RMATCH+10*HCM)
		RS = RASYM
           ELSE IF(RASYM.gt.0.01) THEN
		RS = RASYM
	   ELSE
		RS = RMATCH
           ENDIF
       TH2 = T / (K(PEL,EXL) * RS)
!@       if(NGAIL>0) TH2 = 0.
C       IF(TH1.gt.THMIN .or. TH2.gt.THMIN)
            WRITE(KO,162) TH1,TH2
  162  FORMAT(/'  NOTE: Coulomb excitation cut off ',
     X 'below',F8.3,' deg by JTMAX',:,', and below'
     X ,F8.3,' deg by Rmax'/)
!@        IF(X.gt.0..and.NGAIL==0) WRITE(KO,163) TH2,RS,JTOTAL
        IF(X.gt.0.) WRITE(KO,163) TH2,RS,JTOTAL
  163  FORMAT(' To get scattering to',F5.1,' deg, integrating to',F7.0,
     X   ' fm., and JTMAX should be at least',F6.0/)
!@        IF(X.gt.0..and.NGAIL>=0) WRITE(KO,164) X,JTOTAL
!@  164  FORMAT(' To get scattering to',F5.1,' deg',
!@     X   ', JTMAX should be at least',F8.0/)
 
	RASYMAX = RASYM
       ENDIF
C
C    COULOMB FUNCTIONS
C    -----------------
	DERIV = rterms

!	IF(melfil.eq.0) then
****	CALL COULOMB
	include 'coulomb.f'
!	endif
C
C
C    CDCC OUTPUT (1)
C    -----------------
 	BEPROJ=0.
	KN=1; 
	if(MSP>0) BEPROJ = BE(KN,1)
	KPB = 1 !	partition assumed for deuteron
	ENEX0(:,:,:) = ENEX(:,:,:)
      if(CDCC/=0) then
	I0 = 0
	do IC=1,NCHAN
	  if(MASS(1,IC)<MASS(1,KPB)) then
		I0 = IC
	        go to 152
	   endif
	enddo
	write(KO,*) 'CDCC: Could not find core partition! Stop'
	stop
152	continue
	do 153 KN=1,MSP
	 if(sumccbins)then
          if(abs(CCFRAC(ITC(KPB,1),1,MAX(QNF(5,KN),1)))<1e-9) go to 153
         else
          if(abs(AFRAC(ITC(KPB,1),ITC(I0,1),1,KN))<1e-9) go to 153
	 endif
          BEPROJ = BE(KN,1)
	  SPINV(1) = QNF(10,KN)*0.5
	  PARV(1) = 1
	  go to 154
153	continue
	write(KO,*) 'CDCC: Could not find projectile g.s! Stop'
	stop
154	continue
	BETARG = 0.
        do 1531 KN=1,MSP
         if(sumccbins)then
          if(abs(CCFRAC(ITC(I0,1),2,MAX(QNF(5,KN),1)))<1e-9) go to 1531
         else
          if(abs(AFRAC(ITC(I0,1),ITC(KPB,1),2,KN))<1e-9) go to 1531
         endif
          BETARG = BE(KN,1)
          go to 1541
1531     continue
1541     continue

!				Try to make bin phases continuous:
        	LAST=-1; J2LAST=-1;KNLAST=-1
	do 156 ib=NBINS,1,-1
	 KN = NKBIN(1,IB)
	 L = QNF(9,KN)
	 J2 = QNF(11,KN)
	 if(ib==NBINS.or.L/=LAST.or.J2/=J2LAST) then
	   BSIGN(IB) = 1
	  else
	   I = nint((BPHASE(1,NKBIN(2,IB),KN)-BPHASE(1,1,KNLAST))/PI)
	     BSIGN(IB) = (-1)**I 
	     BPHASE(1,1:NKBIN(2,IB),KN)=BPHASE(1,1:NKBIN(2,IB),KN)-PI*I
	     if(I/=0) 
     X	     write(KO,155) ib,-I,BSIGN(IB)
155          format('  Increase phase of bin ',i4,' by pi*',i3,
     X			', so  overall sign change = ',i2)
	  endif
	 LAST = L
	 J2LAST = J2
	 KNLAST = KN
156	continue    

	MASSV(1) = MASS(1,KPB)-MASS(1,I0)
	CHRGV(1) = MASS(1+2,KPB)-MASS(1+2,I0)
	NAMEV(1) = 'Valence'
	if(abs(MASSV(1)-1.)<.01.and.abs(CHRGV(1)-0.)<.01) 
     x         NAMEV(1)='Neutron '
	if(abs(MASSV(1)-1.)<.01.and.abs(CHRGV(1)-1.)<.01)
     x          NAMEV(1)='Proton  '
	if(abs(MASSV(1)-2.)<.01.and.abs(CHRGV(1)-1.)<.01)
     x          NAMEV(1)='Deuteron'
	if(abs(MASSV(1)-3.)<.01.and.abs(CHRGV(1)-1.)<.01)
     x          NAMEV(1)='Triton  '
	if(abs(MASSV(1)-4.)<.01.and.abs(CHRGV(1)-2.)<.01)
     x          NAMEV(1)='Alpha   '
	MASSV(2) = MASS(2,I0)
	CHRGV(2) = MASS(2+2,I0)
	NAMEV(2) = NAME(2,I0)
	SPINV(2) = JEX(2,I0,1)
	PARV(2) = sign(1,BAND(2,I0,1))
	IBIN(1) = 0
	BSIGN(0) = 1
	NNJMAX = 0
	do 159 IA=2,NEX(KPB)
	 NNJMAX = max(NNJMAX, nint(2.*JEX(1,KPB,IA)+1.))
	 IB = 0
	 do 157 IB=1,NBINS
	 KN = NKBIN(1,IB)
	 if(sumccbins)then
          if(abs(CCFRAC(ITC(KPB,IA),1,MAX(QNF(5,KN),1)))>1e-9) go to 158
         else
          if(abs(AFRAC(ITC(KPB,IA),ITC(I0,1),1,KN))>1e-9) go to 158
         endif
157	 continue
	 write(KO,*) 'CDCC: Could not bin for projectile state ',IA
!	 stop
	 IB = 0
	 KN = 1
158	 IBIN(IA) = IB
	 CHSIGN(IA) = BSIGN(IB)
!	 write(KO,*) 'Excited state ',ia,' has bin ',ib,
!     X		' with phase ',BSIGN(IB)
         ENEX0(1,KPB,IA) = -BE(KN,1)+BEPROJ
!	 write(48,*) ' State ',KPB,IA,' has energy ',ENEX0(1,KPB,IA)
159	continue
         
	IC = 0
	if(PEL/=KPB.and.PEL/=I0) IC=1
	if(CDCC==1) then
        write(57,'(i2)') CDCC
        write(57,'(a120)') HEADNG
        write(57,'(F12.4,4F8.4,i4,f8.4)') ENLAB,BEPROJ,1./FMSCAL,COULCN,
     x               BETARG,IC,(BEPROJ+QVAL(KPB)-QVAL(PEL),IB=1,IC)
	write(57,'(7f8.4)') (MASS(I,KPB),I=1,2),MASS(1,I0),MASSV	! Masses
     x                     ,((MASS(I,PEL),I=1,2),IB=1,IC)
	write(57,'(7f8.1)') (MASS(I+2,KPB),I=1,2),MASS(1+2,I0),CHRGV  ! Charges
     x                    ,((MASS(I+2,PEL),I=1,2),IB=1,IC)
	write(57,'(7A8)') (NAME(I,KPB),I=1,2),NAME(1,I0),NAMEV  	! Names
     x                  ,((NAME(I,PEL),I=1,2),IB=1,IC)
	write(57,'(7f8.1)') (JEX(I,KPB,1),I=1,2),JEX(1,I0,1),SPINV  	! Spins
     x                    ,((JEX(I,PEL,1),I=1,2),IB=1,IC)
	write(57,'(7i8)') (sign(1,BAND(I,KPB,1)),I=1,2),
     X		           sign(1,BAND(1,I0,1)),PARV  	! Parities
     x                  ,((sign(1,BAND(I,PEL,1)),I=1,2),IB=1,IC)
	write(57,'(4i4)') NBINS,NKMAX,NEX(KPB)-1,NNJMAX
	write(57,'(i4,2f8.4)') max(0,NANGL),THMIN,THINC
	do IA=1,NBINS
	  KN = NKBIN(1,IA)
          write(57,'(i2,f4.1,3f8.4,3i4)') QNF(9,KN),QNF(11,KN)*0.5,
     X		 EMID(IA),(KMINX(I,IA),I=1,2),NKBIN(2,IA),KN,NORBIN(IA)
	  write(57,'(10f8.4)') (BPHASE(1,I,KN),I=1,NKBIN(2,IA))
	  enddo
	else ! CDCC=2
        write(57,'(i2)') CDCC
        write(57,'(a120)') HEADNG
        write(57,'(F12.4,4F8.4,i4,f8.4)') ENLAB,BEPROJ,1./FMSCAL,COULCN,
     x               BETARG,IC,(BEPROJ+QVAL(KPB)-QVAL(PEL),IB=1,IC)
        write(57,'(7f8.4)') (MASS(I,KPB),I=1,2),MASS(1,I0),MASSV        ! Masses
     x                     ,((MASS(I,PEL),I=1,2),IB=1,IC)
        write(57,'(7f8.1)') (MASS(I+2,KPB),I=1,2),MASS(1+2,I0),CHRGV  ! Charges
     x                    ,((MASS(I+2,PEL),I=1,2),IB=1,IC)
        write(57,'(7A8)') (NAME(I,KPB),I=1,2),NAME(1,I0),NAMEV          ! Names
     x                  ,((NAME(I,PEL),I=1,2),IB=1,IC)
        write(57,'(7f8.1)') (JEX(I,KPB,1),I=1,2),JEX(1,I0,1),SPINV      ! Spins
     x                    ,((JEX(I,PEL,1),I=1,2),IB=1,IC)
        write(57,'(7i8)') (sign(1,BAND(I,KPB,1)),I=1,2),
     X                     sign(1,BAND(1,I0,1)),PARV    ! Parities
     x                  ,((sign(1,BAND(I,PEL,1)),I=1,2),IB=1,IC)
        write(57,'(5i4)') NBINS,NKMAX,NEX(KPB)-1,NNJMAX,NCHBINMAX
        write(57,'(i4,2f8.4)') max(0,NANGL),THMIN,THINC
	write(57,'(i4)') NEX(I0)-1
	 do IA=2,NEX(I0)
	 write(57,'(I4,2f8.4)') BAND(1,I0,IA),JEX(1,I0,IA),ENEX(1,I0,IA)
	 enddo
       do IA=1,NBINS
          KN = NKBIN(1,IA)
!	   IC=QNF(3,KN)         ! core state in partition I0
           IB=max(1,QNF(5,KN))  ! composite state in partition KPB
          write(57,'(f4.1,3I4,3f8.4,3i4)') 
     x           JEX(1,KPB,IB),BAND(1,KPB,IB),QNF(17,KN),QNF(18,KN),     ! Jex,Pex,nch,IL
     X           EMID(IA),(KMINX(I,IA),I=1,2),NKBIN(2,IA),KN,NORBIN(IA)
	  do I=1,QNF(17,KN)
          write(57,'(i4,f4.1,i4)') QNF(9,KN+I-1),
     x                             QNF(11,KN+I-1)*0.5,QNF(3,KN+I-1)      ! l,j,IC
	  enddo

!	 write(6,*) 'KN,QNF(17,KN)=',KN,QNF(17,KN),' max =',NCHBINMAX
	 if(QNF(17,KN)>NCHBINMAX) then
	 write(6,*) ' *** NOT ENOUGH ROOM FOR BIN OF',QNF(17,KN),' chs!'
 	 write(6,*) ' NCHBINMAX =',NCHBINMAX
	 write(6,*) ' For XCDCC, need to compile with -Dcorex option'
	 stop
	 endif
	  do I=1,NKBIN(2,IA)
            write(57,'(10f10.6)') BPHASE(1:2,I,KN)
 	   if(maxval(abs(BSMAT(1:QNF(17,KN),1:QNF(17,KN),I,KN)))<10.) 
     x 	   then
            write(57,'(10f10.6)') BSMAT(1:QNF(17,KN),1:QNF(17,KN),I,KN)
	   else
            write(57,'(1p,10e10.3)') 
     x                            BSMAT(1:QNF(17,KN),1:QNF(17,KN),I,KN)
	   endif
            ! write(157,*) IA,I
            ! write(157,*) BSMAT(1:QNF(17,KN),1:QNF(17,KN),I,KN)
	  enddo
          enddo
	endif   ! CDCC types

	call flush(57)
	endif   ! CDCC > 0

C
      DO 165 I=1,NSA*NJA
  165 SIGFUS(I,1) = 0.0
      DO 166 I=1,3
      TOTFUS(I) = 0.0
      if(NFUS>0) CORFUS(I,:) = 0.0
      SIGEL(I) = 0.
      SIGTOT(I) = 0.
  166 SIGT(I) = 0.0
      DO 167 I=1,NTHRESH
  167 THRJ(I) = 0.0
	strength(:) = 0.0; RSCATTR=0.0
        inquire(10,opened=op)
!       if(op) close(10,status='delete')
        if(.not.op) open(10,form='unformatted',status='scratch')
       do CP=1,NCP
       if(INFILE(CP)>0) then
        rewind INFILE(CP)
       if(KIND(CP)==1.and.QQ(CP)==1.and.KPCORE(CP)<=-4) 
     x            read(INFILE(CP),*)  ! skip header for FNLREADs
	endif
       enddo

      JSWITCH = 0.0
      SSWITCH = 0.0
      ALOSSM = 0.0
      JTEST = 0
      TSTD = 1.0
      RENORM = 1.0
      IF1 = 1
      IF(ITCM.EQ.1) IF1 = 2
      IF2 = 1+NFUS
	if(final.and.say) write(38,*) ITCM,IF1,IF2
      IF(MACH.le.2.or.MACH==5) THEN
        	inquire(10,opened=op)
!	 write(KO,*) ' File 10: ILEN=',ILEN,' opened =',op
!       if(ILEN==0) 
	if(.not.op) 
     x open(10,access='sequential',status='scratch',form='unformatted')
       rewind 10
      ELSE IF(MACH.eq.6.or.MACH.eq.7) THEN
       if(ILEN==0) open(10,access='sequential',form='unformatted')
       rewind 10
      ELSE IF(MACH==3.or.MACH==4) THEN ! IF(MPIC)
       IF(iams.eq.0.and.ILEN==0) THEN
        T = TM(I)
        if(say.and.final)write(48,*) 'Open file 10 ',real(T)
        open(10,file=TRIM(TMP)//'fort.10',form='unformatted',
     X                status='unknown')
        rewind 10
       ELSE

       ENDIF
      ELSE
	write(6,*) ' Here is old IPSC code. Stop'; STOP 'IPSC'
	ENDIF
	allocate (FNC(NLN,NSA),WNM(NLN,NSA,2))
      IF(VEFF.ne.0) then
        IF(mpisets.gt.1) THEN
         write(KO,*)'*** ERROR: Parallel VEFF calculation NOT supported'
	 stop
        ENDIF
      DO 170 I=1,NLN
      FNC(I,1) = 0.0
  170 WNM(I,1,2) = 0.0
      if(MPIC)then
       open(15,file=TRIM(TMP)//'fort.15',
     X access='sequential',form='unformatted')
      else
       open(15,access='sequential',status='scratch',form='unformatted')
      endif
      REWIND 15
      DO 180 IMA=1,NSA
  180 WRITE(15) (FNC(I,1),I=1,NLN)
      DO 190 IMA=1,NSA
  190 WRITE(15) (WNM(I,1,2),I=1,NLN)
      written(15) = .true.
      endif
	if(INITWF/=0) rewind abs(INITWF)
 	INITWFE = .false.
C
C    FIND EACH J-TOTAL/PARITY COMBINATION
C    ------------------------------------
!$      TCC0 = OMP_GET_WTIME()
      JPSET = 0
      DONE = 0
      LJMAX = 0.0
      SMATL = SMATS
      if(.not.final) SMATL=0
      SHFL = SMATS.ge.6
      JCCSET = 0; MCCSET=0
	MAXB = 0
	MAXCHT = 0
      nphases = 0
      SMALLS(:) = 0
      SMALLJ = .false.
      FUSLL(:,:) = 0.0

      CALL FLUSH(KO)
	if(TRENEG>0) write(89,211) 0,0,0,0
	if(eigens==0) then
	if(SMATL>0) write(56,195) ENLAB,ETOTAL
  195 format('# Fusion cross sections for ELAB=',f12.6,', ECM=',f12.6)
	written(56) = .true.
	endif

	if(say.and.final)write(48,*) 'JBORD:',JBORD(1:NJ+1)
	if(say.and.final)write(48,*) 'JUMP:',JUMP(1:NJ,1)
	if(say.and.final)write(48,*) '740 M,N = ',M,N
      DO 740 JBLOCK=1,NJ
C     IF(CHANS+LISTCC+SMATL.GE.3) WRITE(KO,*) JBLOCK,
C    #       JBORD(JBLOCK),JBORD(JBLOCK+1),JUMP(JBLOCK,1)
         JUMP(JBLOCK,2) = 0
         JUMP(JBLOCK,3) = 0
      IF(JUMP(JBLOCK,1) .LT. 1) GO TO 740
         JAP=JBORD(JBLOCK)
         IF(JBLOCK.GT.1.AND.JUMP(MAX(JBLOCK-1,1),1).GT.1)
     X                                      JAP=JAP+JUMP(JBLOCK,1)
         JAL=JBORD(JBLOCK+1) + 0.1
         IF(JBLOCK.LT.NJ.AND.JUMP(JBLOCK,1).EQ.1) JAL = JAL - 1
C
!      DO 730 JTOTAL=JAP,JAL,JUMP(JBLOCK,1)+Z
      NJTOTAL=NINT((JAL-JAP)/(JUMP(JBLOCK,1)+Z))
      DO 730 IJTOTAL=0,NJTOTAL
      JTOTAL = JAP + IJTOTAL*JUMP(JBLOCK,1)+Z
	if(say.and.final)write(48,*) 'Try JTOTAL =',JTOTAL
         SCHONP = .FALSE.
      JCOEF = (2*JTOTAL+1)* PI  * XCOEF / ABS(K(PEL,EXL))**2
        JUMP2 = JUMP(JBLOCK,1)
        if(IJTOTAL==NJTOTAL.and.JBLOCK<NJ) JUMP2=JUMP(JBLOCK+1,1)
        JUMPER = (JUMP(JBLOCK,1)+JUMP2)*0.5
        T = ETA(PEL,EXL)
        T4= K(PEL,EXL)
        L = JTOTAL
        RTURN =(T+SQRT(T**2 + L*(L+1d0)))/T4
        DONES=0
	FUSJ = 0.0
	TOTJ = 0.0
	SIGJ(:) = 0.0; FUSL(:)=0.0
	if(NFUS>0) CFUSJ(:) = 0.0
	JNLAST = JTOTAL
	SSWITCH = SWITCH
	if(JTOTAL>SINJMAX.and.abs(SINJMAX)>.1) SSWITCH = 1e6
	nparit = 0
      DO 720 IPARIT=2,3
         PARITY = (-1)**IPARIT
        if(abs(PSET).eq.1.and.PARITY.ne.PSET) go to 720
        if(abs(PSET).eq.2.and.
     x     PARITY.ne.(-1)**nint(JTOTAL+.1)*isign(1,PSET)) go to 720
      TIME0J = TM(I)
!$    TIME0W = OMP_GET_WTIME()
C
C     IN CONCURRENT MACHINES, REDIRECT KO and KS outputs
C     ----------------------------------------------------
       KS = 10
!        WRITE(48,*) 'A0',JTOTAL,PSIGN(PARITY+2)
**MPI*******************************************************************
      IF(MPIC) THEN
       CALL FLUSH(KO)
        if(iams==0)then
          KS = 10
          KO = KOI
        else
          KS = 49
          KO = 50
          REWIND KS
          REWIND KO
          if(JCCSET==iams)REWIND 51
        endif
      T = TM(I)
      if(say.and.final)write(48,*) 'Node ',iame,' to consider ',JTOTAL,
     X                    PSIGN(PARITY+2),real(T)
      ENDIF
!        WRITE(48,*) 'A',JTOTAL,PSIGN(PARITY+2)
************************************************************************

      call flush(ko)
C
C   MAKE SET OF COUPLED EQUATIONS
C   -----------------------------
C
***      CALL MAKESET
	include 'makeset.f'


	IF(mpisets<=1) THEN
  	  HERE = .TRUE.
	 ELSE
      	  HERE = MOD(JCCSET-1,mpisets)==iams
	 ENDIF
      if(NCH==0.or.NOSOL.or..not.HERE) go to 710
      if(MINTL==0.and.melfil==0) go to 710
	nparit = nparit + 1
	NCLREQ = maxval(NCLIST(1:NCH,1:NCH))
	if(rela /= '  ') then
	  DO 85 C1=1,NCH
	  DO 85 C2=1,NCH
	  DO 85 NC=1,NCLIST(C1,C2)
85	  CLIST(C1,C2,NC) = GAM(PART(C1,1),EXCIT(C1,1))*CLIST(C1,C2,NC)
	endif
C
	if(.not.rterms) then

C   SOLVE SET OF COUPLED EQUATIONS BY CDE
C   -------------------------------------
C
***	CALL SOLVESET
	include 'solveset.f'
	if(numthread>1) 
     x   write(500+iame,*) ' Node ',iame,'  done solveset'
C
	else
		
C   SOLVE SET OF COUPLED EQUATIONS BY R-MATRICES
C   --------------------------------------------
C
	IF(JSET.GT.0 .AND. JPSET.GT.JSET) GO TO 710

      IF(FLIP) then
	CALL XCH(PSI,N,NCH,NICH,SMAT,JTOTAL,LVAL,
     X      JVAL,JPROJ,JTARG,PART,EXCIT,COPY,MXP,MXX,MAXN,SMATL,SAME,
     X      .true.,EXCH,MAXCH)
       else
 	 EXCH(:,:) = 0.0    ! oopennly allocated (1,1) now!
       endif

	if(nrbases>0) then
	nbas = nrbases

	if(FJSWTCH) nbas = 0
        if(KRM==0) then
	CALL RMATRIX(JTOTAL,NCH,MINTL,INITL,ECM,ECM(1,2),ECM(1,3),ITC,
     X	  CHL,PART,EXCIT,LVAL,JVAL,JPROJ,JTARG,MRM,FORMF,NF,ITCM,EOFF,
     X    CUTVAL,ISOCEN,bndx(IPARIT-1),weak,nbas,nrbmin,nbas*NCH,
     X    KO,CDETR,WOUT,KFUS,
     X    pralpha,PCON,CH,LL1, NCHAN,KS,PARITY,PSIGN,K,RMASS,BLOCKD,
     X    CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,symm,CLIST,NCLIST,NFLIST,
     X    SIGJ,JCOEF,XS,FUSL,SMATS,SMATL,CHANS,DONE,JTMAX,JTMIN,SCALE,
     X    JEX,ABSEND,DONES,NTHRESH,RESM,THRJ,IF1,IF2,TOTJ,FUSJ,meigs,
     X    MR,CHNO,NL,NLO,NLN,NLC,MLT,MLM,ICUTC,ISNONO,FLIP,EXCH,GAP,ETA,
     X    FORML,AFRAC,QNF,MASS,DROPPED,NSTEPD,iams,TIME0,TIME0J,
     X    XCOEF,PTYPE,CFUSJ,SMALLCOUP,SMALLCHAN,SMALLS,CHPRES,RTURN,NEX,
     X	  VEFF,FNC,WNM,NSA,WAVES,WDISK,CFG,WRITTEN,ENLAB,phasin,linel,
     x    IEXCH,FUSLL,CPSO,say,ETOTAL)
         else
        CALL RMATRIXP(JTOTAL,NCH,MINTL,INITL,ECM,ECM(1,2),ECM(1,3),ITC,
     X    CHL,PART,EXCIT,LVAL,JVAL,JPROJ,JTARG,MRM,FORMF,NF,ITCM,
     X    CUTVAL,ISOCEN,bndx(IPARIT-1),weak,nbas,nrbmin,nbas*NCH,
     X    KO,CDETR,WOUT,KFUS,
     X    pralpha,PCON,CH,LL1, NCHAN,KS,PARITY,PSIGN,K,RMASS,BLOCKD,
     X    CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,symm,CLIST,NCLIST,NFLIST,
     X    SIGJ,JCOEF,XS,FUSL,SMATS,SMATL,CHANS,DONE,JTMAX,JTMIN,SCALE,
     X    JEX,ABSEND,DONES,NTHRESH,RESM,THRJ,IF1,IF2,TOTJ,FUSJ,meigs,
     X    MR,CHNO,NL,NLO,NLN,NLC,MLT,MLM,ICUTC,ISNONO,FLIP,EXCH,GAP,ETA,
     X    FORML,AFRAC,QNF,MASS,DROPPED,NSTEPD,iams,TIME0,TIME0J,
     X    XCOEF,PTYPE,CFUSJ,SMALLCOUP,SMALLCHAN,SMALLS,CHPRES,RTURN,NEX,
     X    VEFF,FNC,WNM,NSA,WAVES,WDISK,CFG,WRITTEN,ENLAB,phasin,linel,
     x    IEXCH,FUSLL,CPSO,say,ETOTAL,BAND)
         endif
	else ! Lagrange mesh
	nbas = nlagcc
	if(FJSWTCH) nbas = 0

	CALL LAGRANGE(JTOTAL,NCH,MINTL,INITL,ECM,ECM(1,2),ECM(1,3),ITC,
     X	  CHL,PART,EXCIT,LVAL,JVAL,JPROJ,JTARG,MRM,FORMF,NF,ITCM,EOFF,
     X    CUTVAL,ISOCEN,bndx(IPARIT-1),weak,nbas,nbas*NCH, 
     x    KO,CDETR,WOUT,KFUS,
     X    pralpha,PCON,CH, NCHAN,KS,PARITY,PSIGN,K,RMASS,BLOCKD,
     X    CFMAT,CGMAT,FCWFN,FJSWTCH,CRCRAT,symm,CLIST,NCLIST,NFLIST,
     X    SIGJ,JCOEF,XS,FUSL,SMATS,SMATL,CHANS,DONE,JTMAX,JTMIN,SCALE,
     X    JEX,ABSEND,DONES,NTHRESH,RESM,THRJ,IF1,IF2,TOTJ,FUSJ,meigs,
     X    MR,CHNO,NL,NLO,NLN,NLC,MLT,MLM,ICUTC,ISNONO,FLIP,EXCH,GAP,ETA,
     X    FORML,AFRAC,QNF,MASS,DROPPED,NSTEPD,iams,TIME0,TIME0J,
     X    XCOEF,PTYPE,CFUSJ,SMALLCOUP,SMALLCHAN,SMALLS,CHPRES,RTURN,NEX,
     X	  VEFF,FNC,WNM,NSA,WAVES,WDISK,CFG,WRITTEN,ENLAB,phasin,linel,
     x    IEXCH,FUSLL,IT0,ITER,IPS,say,ETOTAL)

	endif
        JPSET = JPSET + MINTL
	endif
	     if(eigens>0) then
		JCCSET = 0
 		go to 710    
		endif
	     if(abs(NLAB(1))>1.or.num_energies>0) then

	     do JIN=1,MINTL
	     I = nphases+JIN
	     if(I<=40) then
		phaze(I) = phasin(JIN)
		if(LEN==1.and.final.and.say) then
		
                do II=71,171,100
		write(II,88) pre_is,I-1,post_is,
     x  	 JTOTAL,PSIGN(PARITY+2),mod(LVAL(linel(JIN)),100),JIN
88	        format('@',a14,i0,a,' "',f0.1,a1,'/',i2,i2,'"')
                enddo
		
		else
		 if(phaze(I)<phaslst(I)-90) phaseadd(I)=phaseadd(I)+180
		 if(phaze(I)>phaslst(I)+90) phaseadd(I)=phaseadd(I)-180
		endif
		phaslst(I)=phaze(I)
		phaze(I) = phaze(I) + phaseadd(I)
	     endif ! I<=40
	     enddo

		do id=1,datasets
		if(data_type(id)==4
     x	          .and.abs(JTOTAL-data_jtot(id))<.1 
     x            .and.PARITY==data_par(id)) then
                dataEshift = 0.0
          ip = data_shiftvar(id)
          if(ip>0) dataEshift=dataEshift + srch_value(ip)
!                do ip=1,nvars
!                  if(refer(ip,id).and.srch_kind(ip)==6) 
!     x              dataEshift=dataEshift + srch_value(ip)
!                enddo
		 do ip=1,datalen(id)
 		  if(abs(enlab-data_energies(ip,id)-dataEshift)<1e-5) then
		   JIN = data_ch(id)
		   JIN = min(JIN,MINTL); JIN=max(JIN,1)
		   T = phaze(nphases+JIN)
	 	   EE = (T-datavals(ip,id))/dataerr(ip,id)
                   data_chisq(id) = data_chisq(id) + EE**2
                   theoryvals(ip,id) = T
!                  theoryplot(ip,ILEN,id) = T
		   endif
	          enddo
		  endif
	         enddo

	     nphases = nphases+MINTL
	     endif

	if(numthread>1)
     x  write(500+iame,*) ' Node ',iame,'  done cc set 3: ',mpisets
**MPI*******************************************************************
#ifdef MPI
        IF(MPIC.and.mpisets>1) THEN
*--master---------------------------------------------------------------
!  get S from each node to write to 10
      if(iams==0)then
           write(F51,*) ' DONE  from ',DONE,' to',DONE+DONES
	   call flush(F51)
      DONE = DONE + DONES
      T = TM(I)
      if(say.and.final)write(48,*) 'DONE incremented a by ',DONES,
     x     ' to ',DONE,real(T)

         I = (JTOTAL-JAP)/JUMP(JBLOCK,1) + 0.5
         IF(DONE.GE.MAX(MDONE,MINTL+1).and.I.gt.4) DONE=DONE+1000
	      if(say.and.final)write(48,*) 'DONE now a ',DONE
        
        IPUT = MOD(iams+1,mpisets)*mpihelp
	write(F51,*) ' Node ',iams,' to send DONE=',DONE,' to ',IPUT
	call flush(F51)

      call MPI_send(DONE,1,MPI_INTEGER,
     >               IPUT,10,MPI_COMM_WORLD,ierr) 

        T = TM(I)
        write(F51,*) ' DONE=',DONE,' to node',IPUT,' at J=',JCCSET
          call flush(F51)

        T = TM(I)
!        write(48,*) ' Node ',iame,' into sync region:',real(T),
!     X				real(T)-TSYNC

!WRONG    DO 714 i=1,mpinodes-1,mpihelp
        DO 714 i=mpihelp,mpinodes-1,mpihelp
!			Receive SMATs from fresco nodes i to me at 0
      if(.not.finishedcc(i)) then
        write(F51,*) ' Ask for S  from node ',i
	call flush(F51)
        call MPI_recv(JN,1,MPI_DOUBLE_PRECISION,
     >                i,11,MPI_COMM_WORLD,status,ierr) 
      else
        write(F51,*) ' Node ',i,' done'
        goto 714
      endif
	write(F51,*) ' Received JN=',JN,' from Node ',i
	call flush(F51)

	if(JN<0.) then
	  finishedcc(i) = .true.
	!  DONE = DONE+1000 ! dont skip out of loop otherwise
                            ! jccset not counted to final JTOTAL
                            ! and wrong number of JF read back in
        	write(F51,*) ' End-marker received from node ',i
		call flush(F51)
	  go to 714		! nothing more from there!
	  endif

      call MPI_recv(NCH,1,MPI_INTEGER,
     >              i,12,MPI_COMM_WORLD,status,ierr) 
       if(NCH>MAXCH) then
	deallocate (LVAL,PART,EXCIT,JVAL,SMAT)
        deallocate (JPROJ,JTARG)
	allocate(LVAL(NCH),PART(NCH,3),EXCIT(NCH,3),SMAT(NCH),JVAL(NCH))
        allocate(JPROJ(NCH),JTARG(NCH))
	MAXCHT = NCH
	endif
	if(say.and.final)write(48,*) 'JN,JNLAST,SMATL =',JN,JNLAST,SMATL
      if(abs(JN-JNLAST)>0.1.and.SMATL>0) then  ! print out sums for previous J & reset
	write(156,1446) JNLAST,FUSJ,TOTJ,-1,JN,real(iams)
        if(NFUS==0) then
	WRITE(56,1446) JNLAST, FUSJ,TOTJ,TOTJ-FUSJ !,real(JUMP(JBLOCK,1))
	else if(NFUS==1) then
	WRITE(56,1446) JNLAST, FUSJ,TOTJ,TOTJ-FUSJ,CFUSJ(1)
!    x                 ,real(JUMP(JBLOCK,1))
	else
	WRITE(56,1446) JNLAST, FUSJ,TOTJ,TOTJ-FUSJ,
     x   CFUSJ(1:min(7,NFUS)),sum(CFUSJ(1:NFUS)) !,real(JUMP(JBLOCK,1))
	endif
	TOTJ=0.; FUSJ=0.;
	if(NFUS>0) CFUSJ(:) = 0.0
	JNLAST = JN
	endif
      call MPI_recv(MINTL,1,MPI_INTEGER,
     >              i,13,MPI_COMM_WORLD,status,ierr) 
      call MPI_recv(IPARI,1,MPI_INTEGER,
     >              i,14,MPI_COMM_WORLD,status,ierr) 
      call MPI_recv(JUMPER,1,MPI_DOUBLE_PRECISION,
     >              i,15,MPI_COMM_WORLD,status,ierr) 
      call MPI_recv(LVAL,NCH,MPI_INTEGER,
     >              i,16,MPI_COMM_WORLD,status,ierr) 
      call MPI_recv(PART(1,1),NCH,MPI_INTEGER,
     >              i,17,MPI_COMM_WORLD,status,ierr) 
      call MPI_recv(EXCIT(1,1),NCH,MPI_INTEGER,
     >              i,18,MPI_COMM_WORLD,status,ierr) 
      call MPI_recv(JVAL,NCH,MPI_DOUBLE_PRECISION,
     >              i,19,MPI_COMM_WORLD,status,ierr) 
         WRITE(10) JN,NCH,MINTL,IPARI,JUMPER
         WRITE(10) (LVAL(C),JVAL(C),PART(C,1),EXCIT(C,1),C=1,NCH)
         
      DO JIN=1,MINTL 	
      call MPI_recv(EL,1,MPI_INTEGER,
     >              i,20,MPI_COMM_WORLD,status,ierr) 
      call MPI_recv(FUSL(1:1+NFUS1),NFUS+1,MPI_DOUBLE_PRECISION,
     >              i,21,MPI_COMM_WORLD,status,ierr) 
      call MPI_recv(OUTJ,1,MPI_DOUBLE_PRECISION,
     >              i,22,MPI_COMM_WORLD,status,ierr) 
      call MPI_recv(SMAT,NCH,MPI_DOUBLE_COMPLEX,
     >              i,23,MPI_COMM_WORLD,status,ierr) 
         WRITE(10) EL,(SMAT(C),C=1,NCH),FUSL(2:1+NFUS1)
        FUSJ = FUSJ + FUSL(1)
        TOTJ = TOTJ + OUTJ
	if(NFUS>0) CFUSJ(1:NFUS) = CFUSJ(1:NFUS) + FUSL(2:1+NFUS)
	write(156,1446) JN,FUSJ,TOTJ,real(iams)
      ENDDO

       call flush(10)
	write(F51,*) '  ALL SMATS received ',JN,' from ',i
	call flush(F51)     	
	call FILERECV(i,30,KOI)
	call FILERECV(i,31,38)
	call FILERECV(i,32,45)
	write(F51,*) '  ALL stdout received from ',i,';',real(JN)
	call flush(F51)     	
      LUSED = max(LUSED,MAXVAL(LVAL(1:NCH)))
     
714	continue  ! do i=1,mpisets-1
     
C
CCCCC   Included for iPSC/860 & other concurrent
C       Receive baton to allow writing to files KS->10 & KO -> KOI
C
        IGET = MOD(iams-1+mpisets,mpisets)*mpihelp
	if(.not.finishedcc(IGET)) then
           write(F51,*) 'Ask for done from ',IGET
	   call flush(F51)
C       Receive baton from last Node
      call MPI_recv(DONE,1,MPI_INTEGER,
     >              IGET,10,MPI_COMM_WORLD,status,ierr)
           write(F51,*) 'Last DONE = ',DONE,' from ',IGET
           call flush(F51)
	endif

*--slave----------------------------------------------------------------
      elseif(iams/=0)then
! send S from each thread to master
        T = TM(I)
        write(F51,*) ' Node ',iame,' into SMAT region:',real(T),
     X				real(T)-TSYNC
	call flush(F51)
        TSYNC=TM(I)
        i = MOD(iams-1,mpisets)*mpihelp
        write(F51,*) ' Node ',iame,' asks for done from node ',i
	call flush(F51)

C       Receive baton to ask for smats
      call MPI_recv(DONE,1,MPI_INTEGER,
     >              i,10,MPI_COMM_WORLD,status,ierr) 
           write(F51,*) 'Received DONE=',DONE,' from node',i
           call flush(F51)
      IF(DONE.GE.1000) GO TO 750

        IF(MINTL*NCH.gt.0) THEN
         REWIND KS
         READ(KS) JN,NCH,MINTL,IPARI,JUMPER
         READ(KS) (LVAL(C),JVAL(C),PART(C,1),EXCIT(C,1),C=1,NCH)
          call MPI_send(JN,1,MPI_DOUBLE_PRECISION,
     >                   0,11,MPI_COMM_WORLD,ierr)
          call MPI_send(NCH,1,MPI_INTEGER,
     >                   0,12,MPI_COMM_WORLD,ierr)
          call MPI_send(MINTL,1,MPI_INTEGER,
     >                   0,13,MPI_COMM_WORLD,ierr)
          call MPI_send(IPARI,1,MPI_INTEGER,
     >                   0,14,MPI_COMM_WORLD,ierr)
          call MPI_send(JUMPER,1,MPI_DOUBLE_PRECISION,
     >                   0,15,MPI_COMM_WORLD,ierr)
          call MPI_send(LVAL,NCH,MPI_INTEGER,
     >                   0,16,MPI_COMM_WORLD,ierr)
          call MPI_send(PART(1,1),NCH,MPI_INTEGER,
     >                   0,17,MPI_COMM_WORLD,ierr)
          call MPI_send(EXCIT(1,1),NCH,MPI_INTEGER,
     >                   0,18,MPI_COMM_WORLD,ierr)
          call MPI_send(JVAL,NCH,MPI_DOUBLE_PRECISION,
     >                   0,19,MPI_COMM_WORLD,ierr)
           T = TM(I)
           if(say.and.final)
     x    write(48,*) 'to write ',nint(JN),' to node 0 ',real(T)-TSYNC
           call flush(48)
           TM0 = TM(I)

         DO JIN=1,MINTL
          READ(KS) EL,(SMAT(C),C=1,NCH),FUSL(2:1+NFUS1),OUTJ
           call MPI_send(EL,1,MPI_INTEGER,
     >                    0,20,MPI_COMM_WORLD,ierr)
           call MPI_send(FUSL(1:1+NFUS1),NFUS+1,MPI_DOUBLE_PRECISION,
     >                    0,21,MPI_COMM_WORLD,ierr)
           call MPI_send(OUTJ,1,MPI_DOUBLE_PRECISION,
     >                    0,22,MPI_COMM_WORLD,ierr)
           call MPI_send(SMAT,NCH,MPI_DOUBLE_COMPLEX,
     >                    0,23,MPI_COMM_WORLD,ierr)
         ENDDO

         CALL FLUSH(10)
           T = TM(I)
           write(F51,*) ' Sent ',nint(JN),' to file Node 0:',real(T)-TM0
           call flush(F51)
        ENDIF !IF(MINTL*NCH.gt.0)
C --------------------- COPY OVER 50 -> stdout
	CALL FILESEND(KO,0,30)
        CALL FILESEND(38,0,31)
        CALL FILESEND(45,0,32)
C
           write(F51,*) ' DONE  from ',DONE,' to',DONE+DONES
           call flush(F51)
      DONE = DONE + DONES
      T = TM(I)
      if(say.and.final)write(48,*) 'DONE incremented b by ',DONES,
     x                    ' to ',DONE,real(T)
      call flush(48)
         I = (JTOTAL-JAP)/JUMP(JBLOCK,1) + 0.5
         IF(DONE.GE.MAX(MDONE,MINTL+1).and.I.gt.4) then
	   DONE=DONE+1000
	      if(say.and.final)write(48,*) 'DONE now b ',DONE
           write(F51,*) ' DONE  now ',DONE
           call flush(F51)
	   endif
C
CCCCC   Included for iPSC/860 & other concurrent
C       Pass on baton to allow writing to files 6 & 10 from 50 & 49 resp.
C
        IPUT = MOD(iams+1,mpisets)*mpihelp
      call MPI_send(DONE,1,MPI_INTEGER,
     >               IPUT,10,MPI_COMM_WORLD,ierr) 
        T = TM(I)
        if(say.and.final)write(48,*) ' Baton sent to ',IPUT,' at',JCCSET,real(T)-TSYNC
          call flush(48)
C
      endif !if(slave)
      ENDIF !IF(MPIC)
#endif /* MPI */
	if(numthread>1)
     x 	 write(500+iame,*) ' Node ',iame,'  done cc set 4'
************************************************************************
      DONE = DONE + DONES
      T = TM(I)
      if(say.and.final)write(48,*) 'DONE incremented c by ',DONES,' to '
     x      , DONE,real(T)
         I = (JTOTAL-JAP)/JUMP(JBLOCK,1) + 0.5
         IF(DONE.GE.MAX(MDONE,MINTL+1).and.I.gt.4) DONE=DONE+1000
	      if(say.and.final)write(48,*) 'DONE now c ',DONE

      NJDONE = NJ 
      IF(DONE.GE.1000) then
         if(say.and.final)write(48,*) 'DONE ',DONE,' so exit  to 750'
	 GO TO 750
	 endif
C                next PARITY :
  710 if(allocated(CLIS)) deallocate(CLIS,NFLIS,PFLIS)
  720 CONTINUE   ! parity
      if(iams==0.and.HERE.and.SMATL>0) then
	write(156,1446) JNLAST,FUSJ,TOTJ,real(iams)
      if(NFUS==0) then
	WRITE(56,1446) JNLAST, FUSJ,TOTJ,TOTJ-FUSJ !,real(JUMP(JBLOCK,1))
	else if(NFUS==1) then
	WRITE(56,1446) JNLAST, FUSJ,TOTJ,TOTJ-FUSJ,CFUSJ(1) 
!    x        ,real(JUMP(JBLOCK,1))
	else
	WRITE(56,1446) JNLAST, FUSJ,TOTJ,TOTJ-FUSJ,
     x		CFUSJ(1:min(7,NFUS)),sum(CFUSJ(1:NFUS)) !,real(JUMP(JBLOCK,1))
	endif
	call flush(56)
       endif
 1446 	FORMAT(f8.1,11G12.4)
	if(SMATS>5) write(600+nint(JNLAST),*) ETOTAL,
     x                   max(1d-30,FUSJ)/JCOEF/nparit,JNLAST
      IF(MPIC)THEN
       IF(LISTCC.GT.0.and.LISTCC.le.90) LISTCC = LISTCC - mpisets
      ELSE
       IF(LISTCC.GT.0.and.LISTCC.le.90) LISTCC = LISTCC - 1
      ENDIF
      IF(DRY) LISTCC = 0
      IF(DONE.GE.1000)GO TO 750
C                next JTOTAL :
  730 CONTINUE
C                next JBLOCK :
  740 CONTINUE
	     if(ECMC(PEL,EXL).lt.0) go to 948
          T = TM(I)
      if(MPIC.and.final) write(48,*) 'To finish all CC sets @',real(T)
  750 continue
      if(allocated(CLIS)) deallocate(CLIS,NFLIS,PFLIS)
**MPI*******************************************************************
#ifdef MPI
      if(MPIC) then
      if(iams==0)then

            T = TM(I)
            write(48,*) 'Finished all CC sets @ ',real(T)
            call flush(48)

      else !if(iams/=)then
	write(F51,*) 'Send end-marker to Node 0'
	JN = -1.
	call MPI_send(JN,1,MPI_DOUBLE_PRECISION,
     >                 0,11,MPI_COMM_WORLD,ierr)

        IPUT = MOD(iams+1,mpisets)*mpihelp
      call MPI_send(DONE,1,MPI_INTEGER,
     >               IPUT,10,MPI_COMM_WORLD,ierr) 

            T = TM(I)
	    write(48,*) 'Finished all CC sets @ ',real(T)
            call flush(48)

        if(OPEN18) then
		close(18,status='delete')
		OPEN18 = .false.
	  endif
        if(OPEN12) then
		close(12,status='delete')
		OPEN12 = .false.
	  endif
	close(KS,status='delete')

      endif !if(iams/0)
      endif !if(MPIC)
!      if(MPIC)call MPI_barrier(MPI_COMM_WORLD,ierr)
#endif /* MPI */
************************************************************************
      REWIND 10
        if(OPEN18) then
		close(18,status='delete')
		OPEN18 = .false.
		written(18) = .false.
	  endif
!       if(OPEN12) then
!	close(12,status='delete')
!	OPEN12 = .false.
!  	endif
      IF(iams/=0) GOTO 772
        T = TM(I)
      if(final) write(KOI,*) 'Finished all CC sets @ ',real(T)
      call flush(KOI)
C  Serial operation for this section:
         KS=10
         REWIND KS
         MAXCH=MAX(MAXCH,MAXCHT)
         allocate (FUSUM(MAXCH,NSA*NJA))
         if(tcfile>0) then
      	    allocate(TCOEF(LMAX1,NSA),TCOEFS(LMAX1,NSA))
      	    TCOEF(:,:) = 0.0
            TCOEFS(:,:) = 0
           endif
         JN=-1.0
        DO 761 JF=1,JCCSET
!        write(48,*)' reading ',JF,'/',JCCSET,' from file 10',real(T)-TM0
         READ(KS,END=761) JTOTAL,NCH,MINTL,IPARIT,JUMPER
          if(NCH>MAXCH) then
           deallocate (LVAL,PART,EXCIT,JVAL,SMAT)
           deallocate (JPROJ,JTARG,FUSUM)
        allocate(LVAL(NCH),PART(NCH,3),EXCIT(NCH,3),SMAT(NCH),JVAL(NCH))
           allocate(JPROJ(NCH),JTARG(NCH),FUSUM(NCH,NSA*NJA))
           MAXCHT = NCH
           MAXCH = NCH
           endif

         READ(KS) (LVAL(C),JVAL(C),PART(C,1),EXCIT(C,1),C=1,NCH)
	   DO C=1,NCH
      		JPROJ(C)= JEX(1,PART(C,1),EXCIT(C,1))
      		JTARG(C)= JEX(2,PART(C,1),EXCIT(C,1))
	   ENDDO
          RS = 0.0
           IF(abs(JTOTAL-JN).gt.0.1) THEN
              RS=1.0
           ENDIF
           DO 751 C=1,NCH
           DO 751 I=1,NSA*NJA
  751      FUSUM(C,I) = 0.0
         DO 753 JIN=1,MINTL
          READ(KS,end=753) EL,(SMAT(C),C=1,NCH),FUSL(2:1+NFUS1)
C
C     FIND CUMULATIVE SUMS OF REACTION CROSS SECTIONS
C     AND GIVE SUMMARY OF S-MATRIX ELEMENTS ON FILE 7
C     -----------------------------------------------
         DO 752 IT=0,ITCM
  752    SIGJ(IT) = 0.0
      JCOEF = (2*JTOTAL+1)* PI  * XCOEF / ABS(K(PEL,EXL))**2
      CALL SUMX(SMAT,NCH,JIN,EL,JCOEF,XS,SIGT,SIGJ,SIGR,PART,EXCIT,
     X           JUMPER,FUSL,CSIG,JEX,ITC,K,RMASS,PEL,EXL,LVAL,
     X         NSA,NJA,JVAL,JPROJ,JTARG,JTOTAL,TOTFUS,FUSUM,CORFUS,
     X         SIGEL,SIGTOT,SMATS,ITCM,XSTABL,IF1,IF2,FRACFUS)

      if(tcfile>0) then
        RH = (2*JTOTAL+1.)/((2*JVAL(EL)+1.)*(2*JEX(2,PEL,EXL)+1.))
        I = nint(JVAL(EL) - (LVAL(EL) - JEX(1,PEL,EXL)) ) + 1
        L1 = LVAL(EL)+1

        TCOEF(L1,I) = TCOEF(L1,I) +   RH * max(FRACFUS,0d0)
        TCOEFS(L1,I) = TCOEFS(L1,I) + RH

 	endif

!					STRENGTH FUNCTIONS
 	if(LVAL(EL)<=2) then		! S_L
	  strength(LVAL(EL)) = strength(LVAL(EL)) + 
     x        (JTOTAL+0.5)*(1.-Abs(SMAT(EL))**2)/(2*LVAL(EL)+1)
 	endif
 	if(LVAL(EL)==0) 
     x      RSCATTR = imag((1-SMAT(EL))/(1+SMAT(EL)))/K(PEL,EXL)
  753    CONTINUE
!      write(48,*) ' read ',JF,' = ',JTOTAL,' from file 10',real(T)-TM0
      DO 755 I=1,NSA
      SIGFUS(I,2) = 0.0
  755 SIGFUS(I,3) = 0.0
      T = JCOEF*NSA*NJA*JUMPER
      IAM = 0
      DO 758  I=1,NSA
      T4 = 0.0
         DO 756 J=1,NJA
         IAM = IAM + 1
          DO 756 C=1,NCH
  756    T4 = T4 + ABS(FUSUM(C,IAM))**2
  758 SIGFUS(I,IPARIT) = T4 / NJA
      DO 759 I=1,NSA
        AMSA = I-1 -JEX(1,PEL,EXL)
        T4 = 0.0D0
        IF(JTOTAL+0.2.GT.ABS(AMSA)) 
     X   T4 = T * (RS - SIGFUS(I,2) - SIGFUS(I,3))
  759   SIGFUS(I,1) = SIGFUS(I,1) + T4
C
       JN=JTOTAL
  761  CONTINUE
	if(mod(tcfile,2)==1) then
        open(5420,file=tcfilename,access='sequential',form='formatted')
!	 write(5420,*) 0.,0,0,0.,-1, 'END full LSJ,JT',NSA,LUSED
           T = 0.
           do I=1,NSA
              LAP=-1
             do L1=1,LUSED+1
              if(TCOEFS(L1,I)>1e-3) LAP=L1-1
              enddo
	 write(5420,'(a1,2i4,a,f6.1,10x,f10.4,'' MeV cm'')') '#',LAP,NSA, 
     x            '  TC: J-L =',I-1-JEX(1,PEL,EXL),ECMC(PEL,EXL)
  	   RS = JEX(1,PEL,EXL)
           do L=0,LAP
            L1 = L+1
            JVALI = L - JEX(1,PEL,EXL) + I-1
           if(JVALI>=0.0) then
            RH = 0.0
            if(TCOEFS(L1,I)>0.) RH = TCOEF(L1,I)/TCOEFS(L1,I)
            T4 = RH* (2*JVALI+1.)/(2*RS+1.) * 10*PI/K(PEL,EXL)**2
            T = T+T4
           if(TCOEFS(L1,I)>1e-3)
     x      write(5420,'(I4,1p,e14.7,0p,f6.1,2f14.5,f8.4,a5,f10.4)') 
     x     L,TCOEF(L1,I)/TCOEFS(L1,I),JVALI,T4,T,TCOEFS(L1,I)
     x     , '  Ecm',ECMC(PEL,EXL)
 	   endif
	
            enddo ! L
	   write(5420,*) '&'
            enddo ! I


	if(tcfile<2) then
	   write(5420,762) 0
762	  format('#',I4,' states cross sections')
          else
	   write(5420,762) ITCM
          endif

         endif
        if(tcfile>0) deallocate(TCOEF,TCOEFS)
        call flush(5420)
       

       deallocate (FUSUM)
       if(MAXCH>0.and.number_calls<0.and..false.) then ! dealloc if not needed more
       if(final) write(KO,*) 'Deallocate channel arrays'
         deallocate (CLIST,NCLIST,NFLIST,CH,ECM,XS,INCOME,READWF,
     X    CFMAT,CGMAT,CRCRAT,CUTVAL,JVAL,JPROJ,JTARG,CHNO,CFF,
     x   LVAL,PART,EXCIT,INITL,BLOCKD,LL1,LAM,GDER,SMAT,EXCH,PSI,SAME)

         if(.not.rterms) then
		deallocate(SRC,SRATIO,SOMA,SPAD,FIM,FIMD,SIMPLE,EMPTY,SKIP)
         	deallocate(ferw)
		endif
       
       MAXCH=0
       MAXICH=0
       endif
C
C    GIVE ACCUMULATED AND DIFFERENTIAL CROSS SECTIONS
C    ------------------------------------------------
!	MCCSET = JCCSET
	MCCSET = JPSET

	if(num_energies>3) then
!	if(LEN==1) allocate (alphas(MCCSET,MAXXCH,nrbmax*mag))

	endif
      DO 764 CP=1,NCP
  764 IF(MLCALL(2,CP).or.KIND(CP)==7.and.QQ(CP)==2) 
     x       CALL NLSTAT(CP,NLOC(1,CP),NLM,EPC*.01,HNL,CENTRE)
C
      DO 765 I=2,3
      SIGEL(I) = SIGEL(I)/(SIGEL(1) + 1E-20)
      SIGTOT(I) = SIGTOT(I)/(SIGTOT(1) + 1E-20)
      SIGT(I) = SIGT(I)/(SIGT(1) + 1E-20)
      TOTFUS(I) = TOTFUS(I)/(TOTFUS(1) + 1E-20)
	DO 765 IA=1,NFUS
  765 CORFUS(I,IA) = CORFUS(I,IA)/(CORFUS(1,IA) + 1E-20)
      if(.not.final) then
      if(SIGT(1)>1e-5) WRITE(KO,1450) SIGT(1),SIGTOT(1),SIGEL(1)
 1450 FORMAT('  Sig-R,T,E =',3F12.4)
      DO 770 IC=1,NCHAN
      T = maxval(abs(SIGR(IC,1:NEX(IC))))
      if(T>1e-20) then
        WRITE(KO,14611) IC,(SIGR(IC,IA),IA=1,NEX(IC))
	endif
  770 continue
      else
      WRITE(KO,1451) 'REACTION',SIGT
 1451 FORMAT('0CUMULATIVE ',a8,' cross section               =',
     X  F17.5,:,'  <L> =',F9.2,'  <L**2> =',F9.1)
      if(ABS(ETA(PEL,EXL))<1e-9) then
         WRITE(KO,1451) ' TOTAL  ',SIGTOT
         WRITE(KO,1451) ' ELASTIC',SIGEL
       else
         WRITE(KO,1451) ' ELASTIC',SIGELE(ITC(PEL,EXL))
       endif
      DO 771 IC=1,NCHAN
      T = maxval(abs(SIGR(IC,1:NEX(IC))))
      if((T>1e-3.or.T.eq.0d0).and.T<1e4) then 
        WRITE(KO,1460) IC,(SIGR(IC,IA),IA=1,NEX(IC))
	else
        WRITE(KO,1461) IC,(SIGR(IC,IA),IA=1,NEX(IC))
	endif
  771 continue
      endif
	call flush(KO)
 1460 FORMAT('0CUMULATIVE outgoing cross sections in partition',I2,
     X  ' :',7F11.5,/,(52X,7F11.5))
 1461 FORMAT('0CUMULATIVE outgoing cross sections in partition',I2,
     X  ' :',1p,7e11.3,/,(52X,7e11.3))
14611 FORMAT('  Out in',I2,' :',1p,7e11.3,/,(14X,7e11.3))
      if(final.and.say) then
	BETG = BEPROJ + QVAL(MXP)-QVAL(1)
      write(13,*) 'Integrated cross sections for all states'
      write(13,'(3i3,1p,e12.4,0p,2f8.4)') 
     x           PEL,EXL,NCHAN,ENLAB,BEPROJ,BETG
      do IC=1,NCHAN
      write(13,'(2i4)') IC,NEX(IC)
        do IA=1,NEX(IC)
         IT = ITC(IC,IA)
	if(IC==PEL.and.IA==EXL) then
          write(13,'(2(f5.1,i3,f8.4),1p,3e12.4)') 
     X    (JEX(j,IC,IA),BAND(j,IC,IA),ENEX(j,IC,IA),j=1,2),SIGT(1),
     X    SIGTOT(1),SIGEL(1)
        else
	  if(abs(ENEX(1,IC,IA))+abs(ENEX(2,IC,IA))<.01.and.CDCC/=0) then
          write(13,'(2(f5.1,i3,f8.4),1p,e12.4)')   ! special adiabatic case
     X    (JEX(j,IC,IA),BAND(j,IC,IA),ENEX0(j,IC,IA),j=1,2),SIGR(IC,IA)
	  else
          write(13,'(2(f5.1,i3,f8.4),1p,e12.4)') 
     X    (JEX(j,IC,IA),BAND(j,IC,IA),ENEX(j,IC,IA),j=1,2),SIGR(IC,IA)
          endif
	endif
	call flush(13)
	enddo
      enddo
      written(13) = .true.
      WRITE(40,14651) ENLAB,TOTFUS(1),SIGT(1),
     X               ((SIGR(IC,IA),IA=1,NEX(IC)),IC=1,NCHAN),
     x               SIGTOT(1),SIGEL(1)
      written(40) = .true.
      call flush(40)
      endif
      if(final) WRITE(KO,1465) TOTFUS,NAME(1,PEL),(SIGFUS(I,1),I=1,NSA)
 1465 FORMAT('0Cumulative ABSORBTION by Imaginary Potentials  ',
     X '  =',F16.5,'  <L> =',F9.2,'  <L**2> =',F9.1,
     X/,'   Fusion for specific ',A8,' M-states :',4F15.6,/,1X,
     * 8F15.6,/,6F15.6)
        WRITE(KO,1451) 'OUTGOING',SIGT(1)-TOTFUS(1)
	call flush(KO)
!14651 FORMAT(F8.4,16G12.4/(8X,16G12.4))
!14651 FORMAT(G12.5,30G12.5/(12X,30G12.5))
!14651 FORMAT(G12.5,30G13.6/(12X,30G13.6))  
!14651 FORMAT(G12.5,1P,E13.6,0P,38G13.6/(12X,40G13.6))  
!14651 FORMAT(G11.5,1P,E14.6,0P,38G13.6/(12X,40G13.6))  
14651 FORMAT(G12.7,1P,E14.6,0P,38G13.5/(12X,40G13.5))  
!  Put total fusion & reaction cross section in file 39:
!!      WRITE(39,14651) ETOTAL,TOTFUS(1),SIGT(1)
!  Fusion & reaction cross section in file 39:
!lab      IF(NFUS.eq.0) WRITE(39,14651) ENLAB,TOTFUS(1),SIGT(1)
!lab      IF(NFUS.NE.0) WRITE(39,14651) ENLAB,TOTFUS(1),SIGT(1),CORFUS(1,:)
	! write(480,*) say,final,eigens,say.and.final.and.eigens==0
       ELAS = SIGT(1)-TOTFUS(1) 
      if(say.and.final.and.eigens==0) then
      IF(NFUS.eq.0) WRITE(39,14651) ENLAB,TOTFUS(1),SIGT(1),SIGTOT(1),
     x         ELAS,(sum(SIGR(IC,1:NEX(IC))),IC=1,NCHAN)
      IF(NFUS.NE.0) WRITE(39,14651) ENLAB,TOTFUS(1),SIGT(1),SIGTOT(1),
     x         ELAS,(sum(SIGR(IC,1:NEX(IC))),IC=1,NCHAN),CORFUS(1,:)
      endif
      SIGR2(:,:) = SIGR(:,:)
      SIGREAC = SIGT(1) - SIGR(PEL,EXL)
      ! ELAS = SIGT(1)-TOTFUS(1) + SIGR(PEL,EXL)
      ELAS = SIGTOT(1)-SIGREAC
        if(abs(ETA(PEL,EXL))>1e-9) ELAS = SIGELE(ITC(PEL,EXL))   ! Coulomb elastic: integral of deviation from Rutherford
      SIGR2(PEL,EXL) = 0.0  !  All elastic-reorientation is ELASTIC, not REACTION !!
      if(say.and.final.and.eigens==0) then
      WRITE(239,14651) ENLAB,TOTFUS(1),SIGREAC,SIGTOT(1),
     x ELAS,   (SIGR2(IC,1:NEX(IC)),IC=1,NCHAN)
      written(39) = .true.
      call flush(39)
      endif
	T = 0
      DO IC=1,MXP
      if(CPOT(IC,1)==KFUS) then
	 IN =1  !  assume: projectile states 
	if(NFUS>1) write(92,14653) min(NFUS,NEX(IC)),KFUS,IC,IN
14653 	format('#',i3,' states: Partial fusion cross sections ',
     x        'from potential',i3,' (as for partition',i2,' pt#',i2,')')
	if(NFUS>1) write(94,14653) LUSED+1,KFUS,IC,IN
      do IA=1,NEX(IC)
	IT = ITC(IC,IA)
	if(IT<=NFUS) then
      WRITE(KO,14655) IT,KFUS,CORFUS(:,IT)
	T = T + CORFUS(1,IT)
14655 FORMAT('0Cumulative ABSORBTION in state ',i3,' by Imaginary ',
     X  'Potential',I3,'   =',F11.5,:,'  <L> =',F9.2,'  <L**2> =',F9.1)
	R1  = CORFUS(1,IT)/(CORFUS(1,IT)+SIGR(IC,IA))  ! spreading ratio
      WRITE(92,14656) IA,CORFUS(1,IT),SIGR(IC,IA),R1,
     x                JEX(IN,IC,IA),ENEX(IN,IC,IA)-QVAL(IC),ECMC(IC,IA)
14656 format(i3,3f10.5,f5.1,2f8.3)
	endif ! IT
      enddo ! IA

	FUSLLT(:) = 0
	do IA=1,LUSED+1
	 T4 = FUSLL(IA,2) / (FUSLL(IA,1)+FUSLL(IA,2)+1e-10)
	 FUSLLT(:) = FUSLLT(:) + FUSLL(IA,:)
	write(94,14658) IA-1,FUSLL(IA,1:2),T4
14658	format(I4,3f12.5)
	enddo
	write(94,'(a,2f12.5)') '#Tot',FUSLLT(:)
	
	endif
      enddo ! IC
      if(T>0.) WRITE(KO,1466) KFUS,T,FUSLLT(:)
1466  FORMAT('0Cumulative ABSORBTION in all states : Imaginary ',
     X  'Potential',I3,'   =',F11.5, ', XO,PF=',2f11.5)

C		Strength functions
	sfa = sqrt(100./ENCOM)/(2*PI)
	strength(:) = strength(:)*sfa
	
	T = ( K(PEL,EXL)*R0PSTR*MASS(3-LIN,LAB)**(1./3.) )**2
	strength(1) = strength(1) * (1. + 1./T)
	strength(2) = strength(2) * (1. + 3./T + 9./T*T)
        if(final) WRITE(KO,14666) strength(:),R0PSTR,RSCATTR
14666   FORMAT(' Strength functions * 10^4 for L=0-2  =',3f10.5,
     x         ' (with r0=',f5.2,' fm).  R'' =',f8.3,' fm')

      sfa = 1.
      if(abs(ETA(PEL,EXL))<50.) then
      sfa = exp(2d0*PI*ETA(PEL,EXL)) * ENCOM
!     x    ENLAB * MASS(3-LIN,LAB) / (MASS(2,LAB)+MASS(1,LAB))
      if(final) WRITE(KO,14661) sfa
14661 FORMAT(' To convert to S-factors (MeV.mb = keV.b), ',
     X  ' multiply by',1P,E12.4/)
      I=0
       DO 14662 IC=1,NCHAN
       DO 14662 IA=1,NEX(IC)
       if(.not.(IC==PEL.and.IA==EXL)) then
	 I=I+1
	 SFAC(I) = SIGR(IC,IA)*sfa
         endif
14662 	continue
      if(I>0.and.final.and.SMATS+CHANSI+XSTABL>0) then
	   call openif(35); call openif(75)
	 	T = sqrt(SIGEL(1)/(40*PI))
         WRITE(35,14651) ENCOM,SFAC(1:I),strength(:),RSCATTR,T
         WRITE(75,14651) ENLAB,SFAC(1:I),strength(:),RSCATTR,T
	 call flush(35); call flush(75)
         written(35) = I>0
         written(75) = I>0
       endif
      endif

		do id=1,datasets
		if(data_type(id)==3) then

!               Adjust any datanorm search parameter!
                datanorm=1.0; dataEshift = 0.0
          ip = data_normvar(id)
          if(ip>0) datanorm = datanorm * srch_value(ip)
          ip = data_shiftvar(id)
          if(ip>0) dataEshift=dataEshift + srch_value(ip)
!           do ip=1,nvars
!           if(refer(ip,id)) then
!           if(srch_kind(ip)==5) datanorm = datanorm * srch_value(ip)
!           if(srch_kind(ip)==6) dataEshift=dataEshift + srch_value(ip)
!            endif
!           enddo
!           WRITE(KO,*) ' dataset ',id,kq1,lq1,real(datanorm),real(dataEshift)

		 do ip=1,datalen(id)
 		  if(abs(enlab-data_energies(ip,id)-dataEshift)<1e-5) then
	           if(data_ic(id)>0) then
		    T = sigr(data_ic(id),data_ia(id))
		   else
		    if(data_ia(id)==0) T = SIGT(1)
		    if(data_ia(id)==1) T = TOTFUS(1)
		    if(data_ia(id)>=2) T = CORFUS(1,data_ia(id)-1)
		    if(data_ia(id)==-1) T = SIGEL(1)   ! neutron elastic
		    if(data_ia(id)==-2) T = SIGT(1)+SIGEL(1)  ! neutron total=elastic+reaction
		    if(data_ia(id)==-3) T = RSCATTR           ! s-wave potential scattering radius R'
		    if(data_ia(id)==-4) T = strength(0)       ! s-wave strength function
		    if(data_ia(id)==-5) T = strength(1)       ! p-wave strength function for r0=1.35 fm
		    if(data_ia(id)==-6) T = strength(2)       ! d-wave strength function for r0=1.35 fm
	 	   endif
		   if(data_idir(id)==-1) then  ! convert to absolute, first time
		       datavals(ip,id) = datavals(ip,id)/sfa  
		       dataerr(ip,id) = dataerr(ip,id)/sfa  
		     endif
	 	   EE = (T/datanorm-datavals(ip,id))/dataerr(ip,id)
                   data_chisq(id) = data_chisq(id) + EE**2
                   theoryvals(ip,id) = T
		   endif
	          enddo
		  endif
	         enddo
      SIGR(PEL,EXL) = SIGR(PEL,EXL) + SIGT(1)
	if(final.and.SMATL>0) write(56,*) '& '

C
  772 IF(abs(THINC).lt.1e-9) GO TO 785
      rewind 10
C
      IF(MPIC.and.iams.gt.0) GOTO 948
C
      NLJ = NINT(2. * LJMAX + 1.)
      LCROSS = 16
      LFAM=37
      LXSEC=0 ! 199  Not used now, but you can set to unused file no.

      JCCSET = 0
      LEG = 0
      DO 780 IC=1,NCHAN
      DO 780 IA=1,NEX(IC)
      IF(BAND(1,IC,IA) .EQ. 0) GO TO 780
      
      if(final) WRITE(KO,1470) (NAME(IN,IC),IN=1,2),IA,(JEX(IN,IC,IA),
     X      PSIGN(SIGN(1,BAND(IN,IC,IA))+2)    ,IN=1,2),CHSIGN(IA)
 1470 FORMAT(/' CROSS SECTIONS FOR OUTGOING ',A8,' & ',A8,' in state #',
     X I4,' with spins & parities',F4.1,1X,A2,' & ',F4.1,1X,A2,';',i3/)
         IF(DRY .AND..NOT.(IC.EQ.PEL .AND. IA.EQ.EXL)) GO TO 780
         IF(.NOT.GIVEXS(IC)) GO TO 780
         IF(eigens>0) GO TO 780
         IF(BAND(1,IC,IA)*BAND(2,IC,IA).EQ.0) GO TO 780
         IP = SIGN(1,BAND(1,IC,IA))
         MAXPLM = max(1,nint(JEX(1,IC,IA)+JEX(2,IC,IA)+
     x                 JEX(1,PEL,EXL)+JEX(2,PEL,EXL)))
       NMULTIES=0
       if(DGAM>0) then
         IN = PP-1
         DNAME = NAME(IN,IC)
         do 1471 IB=1,IA-1
         IF(COPY(IN,IC,IB,1)>0) go to 1471
           LGAM=abs(JEX(IN,IC,IA)-JEX(IN,IC,IB))
           LGAM = max(LGAM,1)   !  no L=0 photons
           if((-1)**LGAM * SIGN(1,BAND(IN,IC,IA))*SIGN(1,BAND(IN,IC,IB))
     x        <0) LGAM=LGAM+1
           if(LGAM > JEX(IN,IC,IA)+JEX(IN,IC,IB)) go to 1471
             NMULTIES=NMULTIES+1
	     DSPINS(NMULTIES) = JEX(IN,IC,IB)
           DMULTIES(NMULTIES) = LGAM
           DLEVEL(NMULTIES) = IB
1471	 continue
        endif
       SIGCH = SIGR(IC,IA)
       if(IC==PEL.and.IA==EXL) SIGCH = ELAS
	
      CALL CRISS(IC,IA,JEX(1,IC,IA),NSA,NJA,PEL,EXL,JEX(1,PEL,EXL),
     X           LUSED,K,ETA,RMASS,CSIG,KOORDS,MMXCH,MAXF,ENLAB,LEN,
     X           THMIN,THMAX,THINC,EXTRA(1,IC,IA),LAMPL,NEARFA,KQMAX+1,
     X           SIGR(IC,IA),XSTABL,LJMAX,NLJ,NJDONE,JBORD,JUMP,
     X           ITC(IC,IA),ITC(PEL,EXL),IEXCH,PP,PPKIND(PP),CDCC,IP,
     X           IBIN(IA),ENEX(1,IC,IA),PI,HEADNG,IOFAM(1,IC,IA),
     X		 IOFAM(2,IC,IA),LCROSS,LFAM,LXSEC,CHSIGN(IA),
     X           ETOTAL-EOFF,MASS(LIN,LAB),MASS(3-LIN,LAB),
     X           ECMC(IC,IA),MASS(LIN,IC),MASS(3-LIN,IC),
     x           dist0,SIGELE(ITC(IC,IA)),LEG,MAXPLM,ILEN,TCFILE,SIGCH,
     x           DSPINS,DMULTIES,NMULTIES,DLEVEL,DNAME)
	call flush(16)
C
  780 CONTINUE
	if(lampl/=0) write(36,*) 0,0,0
**MPI*******************************************************************
!      IF(MPIC)THEN
!            CALL FILEMV(KO,KOI)
!      ENDIF
************************************************************************
          T = TM(I)
          if(final) write(KOI,*) 'Finished all xsecs @ ',real(T)
	if(mpisets.gt.1) then
          T = TM(I)
          write(48,*) 'xsecs  done @ ',real(T)
          call flush(48)
	  endif
  785 IF(iams.gt.0) GOTO 948
C                            Serial operation now.
      KO = KOI
***      IF(VEFF.NE.0) CALL VEFFPOT
***      IF(BPM.NE.0) CALL BPMFUS
      IF(VEFF.NE.0) then
	include 'veffpot.f'
	close(15)
	endif
      IF(BPM.NE.0) then
	include 'bpmfus.f'
	endif

!  948 if(VEFF.ne.0) deallocate (FNC,WNM)
!  948 deallocate (FNC,WNM)
  948 continue
      if(MACH.ne.6)deallocate (FNC,WNM)
C
C ----------------- NEXT INCIDENT ENERGY
       if(say.and.final.and.(abs(NLAB(1))>1.or.num_energies>0)) then
	 write(71,949) ENLAB,phaze(1:min(20,nphases))
	 call flush(71)
	 write(171,9491) ECM(EL,1),phaze(1:min(40,nphases))
	 call flush(171)
	 endif
! 949   format(f12.8,f8.2,19f7.1)
  949   format(f12.8,20f9.3)
 9491   format(1p,e15.8,',',40(e12.4,','))
	firstE = .false.
  950 CONTINUE
C
  960 IF(iams.gt.0) GO TO 970
       call flush(6)
      close(10,status='delete')
      written(10) = .false.
  970 if(MPIC)then
        if(iams>0)then
         close(38,status='delete')
         close(45,status='delete')
        endif
        close(8,status='delete')
        close(48,status='delete')
        close(49,status='delete')
        close(50,status='delete')
        close(51,status='delete')
	written(48)=.false.
        written(8)=.false.; written(49)=.false.
	written(50)=.false.; written(51)=.false.
      endif
	if(say.and.final.and.iams==0) call fkind(written,KO)
	close(40); close(71); close(171)
	if(OPN(1)) close(11)
	if(OPN(2)) close(9)
	if(OPEN12) then
	   close(12,status='delete')
	   written(12) = .false.
	   endif

      do id=1,datasets
           if(data_type(id)==6) then
            t = (srch_value(data_par(id))-datavals(1,id))/dataerr(1,id)
            data_chisq(id) = t**2
            endif
	   if(data_idir(id)==-1) then
	     data_idir(id) = 0
	     endif
	 ! write(6,*) 'id,type:',id, data_term(id)
           if(data_type(id)==7.or.data_type(id)==8) then
            I = data_term(id)
	     if(pralpha) write(6,*) 'Term',I,'E,W:',E_Brune(I),W_Brune(I)
	     if(pralpha) write(6,*) 'E_Brune(;)',E_Brune
	     if(data_type(id)==7) 
     x           theoryvals(1,id) = E_Brune(data_term(id)) ! Brune energy
	     if(data_type(id)==8) 
     x           theoryvals(1,id) = W_Brune(data_term(id)) ! Brune total observed width
            t = (theoryvals(1,id)-datavals(1,id))/dataerr(1,id)
            data_chisq(id) = t**2
            endif
      enddo

      if(say.and.final.and.iams==0) then
       WRITE(KO,1485) MAXQRN,MLOC,LMAX1,MCLIST,MFNL,MPWCOUP,
     X                MMXQRN,NF,LUSED+1,NCLREQ,NFNL,NPWCOUP
 1485 FORMAT(/' PARAMETERS : MAXQRN   MLOC  LMAX1 MCLIST     MFNL',
     X  '  MPWCOUP'/'    ALLOWED :',4I7,2I9/'    REQUIRED:',4I7,2I9/)

	if(ENLAB>1e-5) then
 	t = hcm*K(PEL,EXL)
	NAMEV(1) = 'OK'
	if(t > hktarg) NAMEV(1) = 'NOT OK'
	write(ko,1486) ENLAB,t,NAMEV(1),hktarg
 1486 format(/' ACCURACY ANALYSIS at',f10.3,' MeV :'//
     x        '  Elastic h*k =',f6.3,' so ',a8,' compared with',f6.3/)
      DO 1467 I=1,NTHRESH
1467  IF(THRJ(I).gt..1) WRITE(KO,1468) RESM(I),THRJ(I)
1468  FORMAT('   Real(S-el) >',F5.2,' first at J =',F10.1)
      IF(JSWITCH.gt.0.1) WRITE(KO,1469) RMATR+GAP,JSWITCH
      WRITE(KO,1469) RTURN,JTOTAL
1469  FORMAT('   R-turn = ',F8.2, ' fm    at J =',F10.1)
C              give warning of limits for Coulomb excitations:
       T = 2. * ETA(PEL,EXL) * 180./PI
        TH1 = T / (JTOTAL+0.5)
        TH2 = T / (K(PEL,EXL) * max(abs(rmatch),rasym))
        WRITE(KO,1621) TH1,JTOTAL,TH2
 1621  FORMAT(/' Forward-angle excitation cut off ',
     X 'below',F8.3,' deg by max JT =',f10.1/
     x '                              and below',F8.3,' deg by max R.'/)

	do CP=1,NCP
	if(WID(CP)>0.) then
	NAMEV(1) = 'OK'
	if(WID(CP) > RNL*1.2) NAMEV(1) = 'NOT OK'
	WRITE(KO,220) CP,WID(CP),rnl,NAMEV(1),CENTR(CP),centre
220   FORMAT(/i2,
     &          ': Recommended RNL: non-local width >',F7.2,
     &              ' cf:',f7.2,' fm: ',a8/
     &        '    Recommended CENTRE: centration   ~',F7.2,
     &              ' cf:',f7.2,' fm'/)
	endif
     	enddo

	endif  ! acc
	endif  ! final

	deallocate(PTYPE)
#ifdef MPI
!	write(500+iame,*) ' FRXX1 call barrier at end'
!      if(MPIC)call MPI_barrier(MPI_COMM_WORLD,ierr)
         DO IT=2,mpihelp    ! terminate helpers
          i = IT-1
          write(500+iame,*) ' Node ',iame,' to terminate ',i
          call MPI_send(0,1,MPI_INTEGER,i,1,commgroup,ierr)
 	enddo
#endif /* MPI */
	call flush(KOI)


      IN  = 1
!$    IN  = OMP_GET_MAX_THREADS()
      TIME1 = SECOND() - TIME0
!$    TIME1W  = OMP_GET_WTIME() -TIME1W
      IF(TIME1.GT.0.01) WRITE(KOI,1495) iame,TIME1
!$   &  ,TIME1W,TIME1/(TIME1W*IN)*100.
 1495 FORMAT(' Total CPU',I3,' time = ',F10.2,' seconds'
!$   x      ',  Wall time =',f10.2,' seconds:',f8.3,' %'
     x  )
	if(numthread>1)
     x  write(500+iame,*) ' Node ',iame,'  done cc set end'
      if(CDCC/=0) write(57,*) MCCSET,MAXCH
	inquire(5420,opened=op); if(op) close(5420)
      RETURN
CUNI  DEBUG SUBCHK
      END
        logical function intrac()
        use searchpar, only: interactive
        intrac = interactive
        return
        end

