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
	program sfrescox
	use parameters
	use factorials
	use io
	use searchpar
	use searchdata
	use trace
	use fresco1, only:lmax,jbord,jump,thmin,thmax,thinc,elab,rterms,
     x    peli=>pel,exli=>exl,labi=>lab,lini=>lin,lexi=>lex
	implicit real*8(a-h,o-z)
	integer kp,pline,col,nafrac,channel,q,par,nparameters,ndatasets,
     x      dataset(2),term,type,maxrank,points,pline2,col2,dof,
     x      icch,iach,lch,pel,exl,pelsave(maxen),variables(mvars),pelii
	real*8 jtot,energy,afrac,potential,datanorm,width,ratio2,B,leff
	real*8 value,step,valmin,valmax,srch_error(mvars),esave(maxen),
     x         stepi,error,jch,sch,srch_vinit(mvars),dataEshift
        real*8,allocatable::     emat(:,:)
	real*8 data_delta(mds),data_xmin(mds),data_value(mds),
     x         data_angle(mds),data_error(mds),data_energy(mds)
        integer data_points(mds),data_iscale(mds),data_leg(mds)
	logical data_abserr(mds),data_detinfo(mds)
        logical, EXTERNAL:: refer
	character*80 input_file,output_file,data_file,search_file,
     x       plot_file,tag,param_file,reffile,n_reffile,s_reffile,
     x       in_file
	character*3 cmlab
	character*8 stype(0:4),errortype
	character*15 name,dir
	character*14 info,version
	character*500 cmd,nums
	character*5 wid,type_tag,pl_suffix
        character psign(3),red,shf
	logical xvals,abserr,lab,undef,noerror,ranks(0:4),dplotted,
     x		nopot,nodat,op,loge,detinfo,eplots,logy,eplot,xmgr,rwa
	logical plotted(mdl,mds),valsonly,all_legends,Mflip,stdin,brune
	external fcn
        data psign / '-','?','+' /

	namelist /variable/ name,kind,step,valmin,valmax, reffile,
     x                  kp,pline,col,potential, dataset,datanorm,
     x     	        nafrac,afrac, energy,jtot,par,channel,width,
     x			term,nopot,pline2,col2,ratio2,param_file,B,
     x                  icch,iach,lch,jch,sch,damp,rwa,ivar,dataEshift,
     x                  leff,brune
	namelist /data/ type,data_file,points,xmin,delta,lab,energy,
     X    idir,iscale,abserr,ic,ia,k,q,angle,jtot,par,channel,leg,idat,
     x    pel,exl,labe,lin,lex,value,error,detinfo,info,eplots,logy,ib,
     x    reffile,kind,term,Mflip
	
        data stype/ ' ','in Fm**2','in barns', 'in mb', 'in micbn'/
! DATA:
!    type      = -2       Legendre coefficient for input energy
!              = -1       Legendre coefficient for fixed energy
!              = 0       angular distribution for fixed energy
!              = 1       excitation & angular cross sections 
!              = 2       excitation cross section for fixed angle
!              = 3       excitation total cross section
!			   ic=0: ia=0 is fusion; 
!				 ia=1 is reaction xs,
!				 ia=-1 is angle-integrated elastic
!				 ia=-2 is total cross section (elastic+reaction)
!				 ia=-3 is potential scattering radius R' = - tan(delta(L=0))/k
!				 ia=-4 is s-wave strength function S_0
!				 ia=-5 is p-wave strength function S_1 for r0=1.35 fm
!				 ia=-6 is d-wave strength function S_2 for r0=1.35 fm
!              = 4       excitation phase shift for fixed partial wave
!              = 5 	 search factor for bound state
!              = 6 	 value and error for a search parameter 
!              = 7 	 value and error for pole energy in Brune basis
!              = 7 	 value and error for total formal width of pole in Brune basis 
!
!    idir      = 0       cross-section data are given in absolute units.
!              = 1       cross-section data are ratio to rutherford.
!              = 2       cross sections are given in absolute units but will
!                             converted to ratio to rutherford.
!              =-1       cross sections are given as astrophysical S-factors but will
!                             converted to absolute

!    iscale    =-1       dimensionless (eg ratio to rutherford if idir>0)
!              = 0       absolute cross-section units are fermi-squared/sr.
!              = 1       absolute scale is barn/sr.    (MeV-barn for idir=-1)
!              = 2       absolute scale is mb/sr. (MeV-mb for idir=-1)
!              = 3       absolute scale is micro-b/sr. (MeV-microbarn for idir=-1)
	
    	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN  
!				Change stdout recl on some machines
	call stdout(6,140)
	nul = -124578   ! 'undefined'
	fine = 10000.    ! this will be the penalty for FAILs, eg in EIGCC
	maxleg = -1
        valsonly = .true.
	hasEshift = .false.
        all_legends = .true.   !  legends also for LINEs as well as PLOTs

	written = .false.; written(3) = .true.
        call version_number(version)
1	write(6,1000)  version
 1000   format('SFRESCOX - version ',a,
     x       ': Search Coupled Reaction Channels'/)
	esmin=0 ; esmax = 1e30
!
!####### Read in search input
	write(6,*) 'Please give search file name'
	read(5,'(a)',end=999,err=1) cmd
	read(cmd,*,end=1001,err=1) search_file,esmin,esmax
 1001	write(6,1002) search_file
 1002   format(' SFrescox search file = ',a80/)
	if(esmax<1e29 .or. esmin>0.) write(6,10021) esmin,esmax
10021   format(' Limit data eneriges to [',1p,g12.4,',',g12.4,']')
        open(303,file=search_file,status='old')

	read(303,*) input_file,output_file,nparameters,ndatasets
	nvars=nparameters ; datasets=ndatasets
        if(nvars>mvars) then
          write(6,*) 'ONLY ROOM FOR ',mvars,' SEARCH VARIABLES'
          stop
          endif
        if(datasets>mds) then
          write(6,*) 'ONLY ROOM FOR ',mds,' DATASETS'
          stop
          endif
	write(6,1003) input_file,output_file,nparameters,ndatasets
 1003   format(' Frescox input file = ',a/,:,
     x         '    and output file = ',a/,
     X         ' Search on',i4,' variables for ',i3,' datasets,')

         stdin = .false.
	if(trim(input_file)=='=') then
         stdin = .true.
	 input_file = trim(output_file)//'.in'
	 write(6,'(/a)') 'Rewriting Frescox input file: '//
     x                   trim(input_file)
         open(309,file=input_file,status='unknown')
	 read(303,*) ! skip line
130       read(303,'(a)') cmd
          write(6,'(a)') '>>'//trim(cmd)
          if(cmd(1:3)=='   '.or.cmd(1:1)=='#'.or.cmd(1:1)=='!') goto 130
          if(cmd(1:3)=='EOF') go to 134
          write(309,'(a)') trim(cmd)
          go to 130
134       close(309)
 	 endif
        open(309,file=input_file,status='old')
	write(6,1003) trim(input_file)//' rewritten'
	call MNINIT(5,6,33)

!
!####### Read in specification of search parameters
       	if(nparameters>0) write(6,'(/''  Define SEARCH VARIABLES:''/)')
        i0 = ichar('0')
	undef = .false.
	rterms=.false.
	stepi = 0.01
	nfparam=303
        mterms = 0;  srch_name(0) = 'None'; srch_Brune(:) = .false.
	 ! ndata_offset =  0;  lastset = 0

        do ip=1,nparameters
        name='Var'//char(mod(ip/10,10)+i0)//char(mod(ip,10)+i0)
        kind=0;valmin=0;valmax=0; kp=0;pline=0;col=0;
	pline2=0; col2=0; ratio2=0.0
        potential=nul; step=stepi; width=0; rwa=.true.; B=nul
        nafrac=0;afrac=nul; energy=nul;jtot=0;par=0;channel=0;term=1
        dataset(1)=0;dataset(2)=0; datanorm=1.0; dataEshift= 1.0
        nopot=.false.; param_file=' '; ivar=0; reffile=' ';brune=.false.
        icch=1; iach=1; lch=-1; jch=-1.0; sch=-1.0; damp=0.0; leff=0.0
        read(nfparam,nml=variable)

	 if(lnblnk(param_file)>0) then  ! redirect to new parameter file
          nfparam=305
          open(nfparam,file=param_file,status='old')
	  read(nfparam,nml=variable)  ! re-read from here
	 endif
	

        srch_kind(ip) = kind; srch_name(ip) = name
        srch_minvalue(ip)=valmin; srch_maxvalue(ip)=valmax
	if(abs(step-stepi)<1e-10.and.step>abs(width) ! make default steps small for small widths
     x             .and.abs(width)>1e-20) step=step*abs(width)
        srch_step(ip)=step
	srch_error(ip)=0  ! initially
	srch_damp(ip)=0.  ! initially
	srch_B(ip)=B

         if(kind==1) then !   potential parameter
 	  write(6,1010) ip,name,kind,kp,pline,col
1010	  format('   Variable',i4,'=',a15,' is',i2,' potential KP=',i2,
     X       ', line=',i2,' col=',i2)
          srch_value(ip) = potential; 
          srch_kp(ip) = kp; srch_pline(ip)=pline; srch_col(ip)=col
	  srch_pline2(ip)=pline2; srch_col2(ip)=col2
          srch_ratio2(ip) = ratio2; 
	  if(pline2*col2>0) then
 	   write(6,1011) ip,name,kind,kp,pline2,col2,ratio2
1011	   format('   Variable',i4,'=',a15,i2,' in potential KP=',i2,
     X       ' is linked to variable2 at line=',i2,' col=',i2,
     x       ' with ratio',f8.4)
	  endif

         else if(kind==2) then ! spectroscopic amplitude
          write(6,1012) ip,name,kind,nafrac
1012	  format('   Variable',i4,'=',a15,' is',i2,' Afrac #',i3)
          srch_value(ip) = afrac; srch_nafrac(ip) = nafrac; 

         else if(kind==3) then ! R-matrix energy
 	  write(6,1013) ip,name,kind,term,jtot,psign(par+2),nopot,damp
1013	  format('   Variable',i4,'=',a15,' is ',i2,
     X       ' energy of R-matrix term',i4,' at J/pi =',f5.1,a1,
     x        ' [NoPot=',L1,']',:,' damping=',1p,e10.2)
          srch_value(ip) = energy; srch_rterm(ip) = term
	  srch_jtot(ip) = jtot; srch_par(ip) = par
	  srch_nopot(ip) = nopot; 
          srch_damp(ip) = damp    ! TEMP 
	  rterms=.true.; mterms = max(mterms,term)

         else if(kind==4) then ! R-matrix width
           if(rwa) then
              wid = 'r.w.a';  red='r'
              widk = 1000*width**2
            else
              wid = 'width';  red=' '
              widk = 1000*width
            endif
          if(channel>0)  then
 	   write(6,1014) ip,name,kind,wid,term,channel,red,widk
1014	   format('   Variable',i4,'=',a15,' is',i2,' ',a5,
     x       ' of R-matrix term',
     X       i4,' in channel',i3,' (',a1,'w =',g10.2,' keV)')
	   srch_r_ch(ip) = channel
          else
	   ! if(jch<0.0) jch = lch
 	   write(6,1015) ip,name,kind,wid,term,icch,iach,lch,jch,sch,
     x                   red,widk
1015	   format('   Variable',i4,'=',a15,' is ',i2,' ',a5,
     x       ' of R-matrix term',
     X       i4,' in channel ic=',i1,' ia=',i2,' l=',i2,' j=',f5.1,
     x       ' S=',f5.1,' (',a1,'w =',g10.2,' keV)')
	   srch_r_ic(ip) = icch
	   srch_r_ia(ip) = iach
	   srch_r_lch(ip) = lch
	   srch_r_jch(ip) = jch
	   srch_r_sch(ip) = sch
	  endif
           if(abs(B-nul)>1e-3) write(6,1016) B
1016       format('             Override:  B=',f10.4)

          srch_value(ip) = width; srch_rterm(ip) = term
          srch_rwa(ip) = rwa

         else if(kind==5) then ! dataset normalisation
         if(dataset(2)==0) then
          write(6,1018) ip,name,kind,dataset(1),trim(reffile)
1018	  format('   Variable',i5,'=',a15,' is',i2,' normalisation',
     x      ' for dataset ',i3,' or filenames ',a)
         else
          write(6,10182) ip,name,kind,dataset(1:2),trim(reffile)
10182	  format('   Variable',i5,'=',a15,' is',i2,' normalisation',
     X       ' for datasets',i4,' through',i4,' (inclusive)',
     x       ' or filenames ',a)
         endif
!          lastset = maxval(dataset(1:2))
          srch_value(ip) = datanorm; srch_datanorm(:,ip) = dataset 
          srch_reffile(ip) = reffile

         else if(kind==6) then ! dataset shift
	  hasEshift = .true.
         if(dataset(2)==0) then
          write(6,10184) ip,name,kind,dataset(1),trim(reffile)
10184      format('   Variable',i5,'=',a15,' is',i2,
     X     ' energy shift up for dataset ',i3,' or filenames ',a)
         else
          write(6,10186) ip,name,kind,dataset(1:2),trim(reffile)
10186     format('   Variable',i5,'=',a15,' is',i2,' energy shift up',
     X       ' for datasets',i4,' through',i4,' (inclusive)',
     x       ' or filenames ',a)
         endif
          srch_value(ip) = dataEshift; srch_datanorm(:,ip) = dataset
          srch_reffile(ip) = reffile

         else if(kind==7) then ! R-matrix damping
          write(6,10188) ip,name,kind,term,damp
10188      format('   Variable',i4,'=',a15,' is',i2,' damping width of '
     X       ,'R-matrix term',i4, ': damping=',1p,e10.2,' MeV')
          if(abs(energy-nul)>.1) then
             cmd = ''
             if(brune) cmd = ', in Brune basis'
             write(6,10189) leff,energy,trim(cmd)
10189      format('   Damping rescaled by penetrability for',
     x     ' Leff=',f6.1,' above threshold energy',f15.2,' MeV',a)
            else
             brune = .false.
            endif
          srch_rterm(ip) = term
          srch_value(ip) = damp
          srch_leff(ip) = leff
          srch_damp(ip) = energy
          srch_Brune(ip) = brune

        endif

         if(abs(srch_value(ip)-nul)>.001) then
         if(abs(srch_value(ip))<1e4.and.abs(srch_value(ip))>.01) then
          write(6,1019) srch_value(ip),step,valmin,valmax
1019	  format('     value ',f12.4,8x,' step ',f7.4,
     X           ', min,max ',2f8.4/)
         else
          write(6,1020) srch_value(ip),step,valmin,valmax
1020	  format('     value ',g14.4,4x,' step ',f7.4,
     X           ', min,max ',2f8.4/)

         endif
	 else
          write(6,1021)  step,valmin,valmax
1021	  format('     value from Frescox input',
     X           ', step ',f7.4,', min,max ',2f8.4/)
	  undef = .true.
	 endif
       enddo
       if(nfparam==305) close(305)
!	if(allocated(rm_Brune)) deallocate(rm_Brune)
       allocate(rm_Brune(0:mterms),E_Brune(0:mterms),W_Brune(0:mterms))
       rm_Brune(:) = .false.; E_Brune(:)=0.0; W_Brune(:)=0.0
!
!####### Read in experimental data sets
       	inf=0
	ranks(:) = .false.; maxrank=-1
	num_energies=0; energy_count(:) = 0
        scatmin=1e30 ; scatmax = 0.

        do id=1,ndatasets
	  type=0; angle=0; jtot=-1;par=0;channel=1;ib=0
          data_file="="; xmin=0;delta=-1;idir=0;iscale=-1;lab=.false.
          ic=1;ia=1;k=0;q=0; abserr=.false.;  points=-1; energy=-1.
	  pel=peli; exl=exli; labe=labi; lin=lini; lex=lexi; leg=-3
	  detinfo=.false.; info=' '; eplots=.false.;logy=.false.
          reffile= ' '; kind=-1; Mflip=.false.

	  write(6,*) 
!	  write(6,*) '** Read definition of data set ',id
          read(303,nml=data)

          if(idir.eq.1) iscale=-1
          errortype='relative'; if(abserr) errortype='absolute'
          cmlab='CM '; if(lab) cmlab='LAB'
	  dat_file(id) = data_file
	  dat_info(id) = info
	  dat_eplots(id) = eplots
	  dat_logy(id) = logy
	    data_energy(id) = energy
	    data_delta(id) = delta
	    data_points(id) = points
            data_xmin(id) = xmin
	    data_angle(id) = angle
	    data_error(id) = error
	    data_abserr(id) = abserr
	    data_iscale(id) = iscale
	    data_detinfo(id) = detinfo
	    data_leg(id) = leg
	    data_value(id) = value
          write(6,1025) id,data_file,type,stype(iscale+1),errortype,
     x      	cmlab,ic,ia,idir,detinfo,info,eplots,logy,Mflip
1025      format('** Read DATA SET ',i3,' from file ',a80/
     x     '    of type ',i2,'  ',a8,' with ',a8,' errors,',
     x     1x,a3,' for state/excit',2i3,',  idir=',i2,/
     x     30x,' Detailed info=',l1,';   set info =',a11,
     x         ' eplots=',L1,' with log(Y)=',L1,' Mflip=',L1/)
	  ranks(k) = .true.
	  maxrank = max(maxrank,k)
	  if(type==-3) write(6,10322) k,q,leg
	  if(type==-2) write(6,10312) k,q
	  if(type==-1.and.energy<0) write(6,1027) k,q
	  if(type==-1.and.energy>0) write(6,1028) k,q,energy
	  if(type==0.and.energy<0) write(6,1029) k,q
	  if(type==0.and.energy>0) write(6,1030) k,q,energy
	  if(type==1) write(6,1031) k,q
	  if(type==2) write(6,1032) k,q,angle,cmlab
          if(type<=3.and.ib>0) write(6,1033) 'gamma decay to state',ib
	  if(type==3.and.ic>=1) write(6,1033) 
	  if(type==3.and.ic==0) then
		if(ia==0) write(6,1033) 'fusion'
		if(ia==1) write(6,1033) 'reaction'
		if(ia==-1) write(6,1033) 'elastic'
		if(ia==-2) write(6,1033) 'total'
		if(ia==-3) write(6,1033) 'R'''
		if(ia==-4) write(6,1033) 'S_0'
		if(ia==-5) write(6,1033) 'S_1'
		if(ia==-6) write(6,1033) 'S_2'
		if(ia>=2) write(6,1033) 'outgoing',ia
		endif
	  if(type==4) write(6,1034) jtot,psign(par+2),channel
	  if(type==5) write(6,1035) 
	  if(type==6) write(6,1036) par,value,error,abserr
1027	  format('  Legendre coefficient of cross section T',2i1)
1028	  format('  Legendre coefficient of cross section T',2i1,
     x				' for energy',f8.3,' MeV')
1029	  format('  Angular cross section T',2i1)
1030	  format('  Angular cross section T',2i1,
     x				' for energy',f8.3,' MeV')
1031	  format('  Excitation and Angular cross sections T',2i1)
10312	  format('  Excitation and Legendre cross sections T',2i1)
1032	  format('  Excitation cross section T',2i1,
     x				' for angle',f8.3,' deg ',a3)
10322	  format('  Excitation cross section T',2i1,
     x				' for Legrendre order',i4)
1033	  format('  Excitation cross section',:,' for ',a,:,' in',i3)
1034	  format('  Excitation phase shift in channel',f5.1,a1,' #',i2)
1035	  format('  Target search parameters for bound states')
1036	  format('  Value and error for search parameter',i3,' are'/
     x        2f10.5,' (abserr=',l1,')')
1037	  format('  Value and error for Brune energy of term ',i3,
     x        ' are',/2f10.5,' (abserr=',l1,')', i4)
1038	  format('  Value and error for Brune total formal width of',
     x        'term ',i3,' are',/,2f10.5,' (abserr=',l1,')', i4)
	  if(type<-3.or.type>8) write(6,*) 'Unrecognised data type ',type
	  data_type(id) = type
	if(type==6) then
	      datavals(1,id)=value
	      if(.not.abserr) error = error*value
 	      dataerr(1,id)=error
 		datalen(id) = 1
 	      ip = 2
	else if(type==7 .or. type==8) then
	      datavals(1,id)=value
	      if(.not.abserr) error = error*value
 	      dataerr(1,id)=error
              data_term(id) = term
              rm_Brune(term) = .true.
              par = 0
              do ip=1,nparameters
               if(srch_kind(ip)==3 .and. srch_rterm(ip)==term) par=ip
	      enddo
              data_par(id) = par
              
 		datalen(id) = 1
 	      ip = 2
	  if(type==7) write(6,1037) term,value,error,abserr,par
	  if(type==8) write(6,1038) term,value,error,abserr,par
	else
              factor=1.0
	      if(iscale<0) factor=1. 
	      if(iscale==0) factor=10.
	      if(iscale>0) factor=1000.0/10.0**(3*(iscale-1)) 
	  inf=306
	  if(data_file=="=") inf=303
	  if(data_file=="<") inf=5
          if(inf==306) then
             write(6,1039) id,data_file
1039        format(' To open file for dataset ',i4,': ',a)
             open(inf,file=data_file,status='old')
             endif
          xvals = delta<=0.; x = xmin
          if(points<0) points=99999
	  write(6,*) ' Read data set ',id
      do 10 ip=1,points
	 	info = ''
	    if(type==1.or.type==-2) then
          if(detinfo) then 
		     read(inf,end=111,fmt=*,err=11) x,a,val,err,info
	      else
		     read(inf,end=111,fmt=*,err=11) x,a,val,err
	      endif ! detinfo
        else if(xvals.or.type==5) then
              if(detinfo) then 
               read(inf,end=111,fmt=*,err=11) x,val,err,info
	      	  else
               read(inf,end=111,fmt=*,err=11) x,val,err
	      	  endif
        else 
              read(inf,end=111,fmt=*,err=11) val,err
              x = x + delta            
        endif !type
              if(x<0) go to 11
	      val =factor*val
	      if(abserr) then
	         err =factor*err
	        else
	         err = val * err
	        endif  ! Now all errors are absolute (and mb, except for r/ruth)
	      datavals(ip,id)=val
 	      dataerr(ip,id)=abs(err)
 	      data_info(ip,id) = info

            if(type==0) then        ! angular distribution for fixed energy
                  datangles(ip,id) = x
	      	  data_energies(ip,id) = energy
             else if(type==-1) then  ! Legendre coefficient of cross section for fixed energy
                  datangles(ip,id) = nint(a)  ! integer order
 	          data_energies(ip,id) = energy
             else if(type==1.or.type==-2) then  ! excitation and angular cross sections
 	          data_energies(ip,id) = x
                  datangles(ip,id) = a
             else if(type==2) then  ! excitation cross section for fixed angle
 	      	  data_energies(ip,id) = x
	      	  datangles(ip,id) = angle
             else if(type==-3) then  ! excitation cross section for fixed Legendre order
 	      	  data_energies(ip,id) = x
	    	  datangles(ip,id) = real(leg)
	     else if(type==3) then  ! excitation total cross section 
 	     	  data_energies(ip,id) = x
             else if(type==4) then  ! excitation phase shift
 	          data_energies(ip,id) = x
             else if(type==5) then  ! bound state search
 	          bs_no(ip,id) = nint(x)
	         endif
  10      continue
          ip = points+1
	  go to 111
  11      backspace inf
	endif
 111      datalen(id) = ip-1
          if(inf==306) then
!             write(6,*) ' Close file for dataset ',id,': ',data_file
             close(inf)
             endif
          data_idir(id)=idir; data_idir1(id)=idir; data_lab(id)=lab; 
          data_ic(id) = ic; data_ia(id) = ia; ; data_ib(id) = ib; 
          data_rank_k(id) = k; data_rank_q(id)=q; data_Mflip(id)=Mflip
	  data_ch(id) = channel; data_jtot(id)=jtot; data_par(id)=par
          data_reffile(id) = reffile; data_kind(id) = kind
	  if(data_type(id)/=6) then
          if(energy<0) write(6,112) datalen(id)
112	  format(' ',i4,' data points',:,': for lab energy',f10.3)
          if(energy>0) write(6,112) datalen(id),energy
     	  endif
	  if(datalen(id)==0) then
	    write(6,*) '   NO DATA POINTS !! Stop'
	    stop
	    endif
            if(pel.le.0) pel = 1
            if(exl.le.0) exl = 1
            if(labe.eq.0) labe = pel
            if(lin.eq.0) lin = 1
            if(lex.le.0) lex = 1
	    data_pel(id) = pel
	    data_exl(id) = exl
	 	if(exl>1) then
		 write(6,*) ' EXL=',EXL,'>1 not implemented'
		 stop 'EXL'
		endif
	    data_labe(id) = labe
	    data_lin(id) = lin
	    data_lex(id) = lex
            if(pel+exl+labe+lin+lex>5) then
	      write(6,1250) pel,exl,labe,lin,lex
 1250       format('     Incoming partition',I3,' in state #',I2,
     X             ',   Lab Energy for part.',I3,' Nucleus',I2,
     X             ' in Excitation pair',I2/)
              endif
       if(type==-2) then
          write(6,*)'   Energy  Order  Datum '//cmlab//' Absolute error'
          do ip=1,datalen(id)
	  		info = data_info(ip,id)
          write(6,118) data_energies(ip,id),nint(datangles(ip,id)),
     x                      datavals(ip,id),dataerr(ip,id),info
  118 	   format(1x,f10.5,i8,2g12.4,8x,a)
	   maxleg = max(maxleg,nint(datangles(ip,id)))
          enddo
       else if(type==-1) then
          write(6,*) '   Order  Datum '//cmlab//'  Absolute error'
          do ip=1,datalen(id)
          write(6,119) nint(datangles(ip,id)),
     x                      datavals(ip,id),dataerr(ip,id)
 119       format(1x,i8,2g12.4)
	   maxleg = max(maxleg,nint(datangles(ip,id)))
          enddo

	 else if(type==0) then
          write(6,*) ' Angle '//cmlab//' Datum     Absolute error'
	  do ip=1,datalen(id)
	  write(6,12) datangles(ip,id),datavals(ip,id),dataerr(ip,id)
  12 	   format(1x,f10.5,1x,2g12.4)
	  enddo
	 else if(type==1) then
          write(6,*) '   Energy   Angle '//cmlab//
     x               ' Datum     Absolute error'
	  do ip=1,datalen(id)
	  info = data_info(ip,id)
	  write(6,13) data_energies(ip,id),datangles(ip,id),
     x			datavals(ip,id),dataerr(ip,id),info
  13 	   format(1x,f10.5,f8.3,2g12.4,8x,a)
	  enddo
	 else if(type==5) then
          write(6,*) '   Bound state   Target     Absolute error'
	  do ip=1,datalen(id)
	  write(6,131) bs_no(ip,id),datavals(ip,id),dataerr(ip,id)
  131 	   format(1x,i8,3x,2f12.4)
	  enddo
	 else if(type==6) then

          if(data_par(id)==0) then  ! have to use data_reffile and data_kind to find variable
            reffile = data_reffile(id)
            kind    = data_kind(id)
            do ip=1,nvars
             if(srch_kind(ip)==kind.and.srch_reffile(ip)==reffile) then
               data_par(id) = ip
              endif
            enddo
          endif

          do ip=1,datalen(id)
          write(6,*) '   Search param         Target     Absolute error'
     x      //'  reffile'
	  write(6,132) data_par(id),srch_name(data_par(id)),
     x                 datavals(ip,id),dataerr(ip,id),trim(reffile)
  132 	   format(1x,i3,': ',a15,2f12.4,4x,6x,a)
	   enddo
	 else if(type==7.or.type==8) then
          write(6,*) '   Search term          Target     Absolute error'
	  write(6,132) data_term(id),srch_name(data_par(id)),
     x                 datavals(ip,id),dataerr(ip,id)

	 else
          write(6,*) '  Energy   Datum     Absolute error'
	  do ip=1,datalen(id)
	  write(6,12) data_energies(ip,id),datavals(ip,id),
     x			dataerr(ip,id)
	  enddo
	   if(type==-3) maxleg = max(maxleg,leg)
	 endif
          close(1)
	 if(type<5) then
	  scatmin = min(scatmin,
     x                  minval(data_energies(1:datalen(id),id)))
	  scatmax = max(scatmax,
     x                  maxval(data_energies(1:datalen(id),id)))
	 endif
        enddo      ! id
        
       do ip=1,nvars   ! add in any missing reffile links from norm/shift vars to data
       if(srch_kind(ip)==5.or.srch_kind(ip)==6) then

       dataset=srch_datanorm(:,ip) 
       reffile = srch_reffile(ip)
       if (len(trim(reffile))<1) then
             reffile = dat_file(dataset(1)) 
             if(reffile/='=' .and. reffile/='<') then
               if (dataset(2)>0) then
                i = index(reffile,'@')
                if(i>0) reffile = reffile(1:i+1)//'*'       
               endif
             srch_reffile(ip) = reffile
             endif
         endif
       endif
       enddo

	data_normvar(:) = 0;  data_shiftvar(:) = 0
	do id=1,ndatasets
	do ip=1,nvars
        if(refer(ip,id)) then
        if(srch_kind(ip)==5) then
	  if(data_normvar(id)/=0) then
	     write(6,122) 'norm',id,data_normvar(id),ip
	     stop
	   endif
	  data_normvar(id) = ip
	  endif
        if(srch_kind(ip)==6) then
	  if(data_shiftvar(id)/=0) then
	     write(6,122) 'shift',id,data_shiftvar(id),ip
	     stop
	   endif
	  data_shiftvar(id) = ip
	  endif
         endif
  	enddo !ip
  	enddo !id
122 	format('*** Duplicate ',a,' variables for set',i4,':',2i4,
     x         '!!   STOP')

       
       do id=1,ndatasets ! add in reffile for type=6 data if existing for variable target
       if(data_type(id)==6) then
         ip = data_par(id)
         if(ip==0) then  ! have to use data_reffile and data_kind to find variable
           reffile = data_reffile(id)
           kind    = data_kind(id)
	   do ip=1,nvars
             if(srch_kind(ip)==kind.and.srch_reffile(ip)==reffile) then
	      data_par(id) = ip
             endif
	   enddo
         else
          reffile = srch_reffile(max(1,ip))
          if(len(trim(reffile))>1) then
           if(srch_kind(ip)==5) data_reffile(id) = reffile
           name = 'link:'//trim(srch_name(ip))  ! for information only: names are not unique yet
           if(srch_kind(ip)==3) data_reffile(id) = name
           if(srch_kind(ip)==4) data_reffile(id) = name
           if(srch_kind(ip)==7) data_reffile(id) = name
          endif
         endif
        endif
       enddo
       
	  if(scatmax>=scatmin) write(6,133) scatmin,scatmax
133	  format(' Scattering energies:',f10.5,' to',f10.5,' MeV(lab)')
!	  esmin = max(esmin,scatmin)  ! global
!	  esmax = min(esmax,scatmax)
          call  energylist ! (ndatasets,esmin,esmax)

        close(2)  
	srch_vinit(1:nvars) = srch_value(1:nvars)
        write(6,*) 

!			Pre-read input to find SOME array-parameter limits
!
!####### Read in main fresco input
	ki = 309
	ko = 3
	koe = 307
	call machine(mach)
	open(ko,form='formatted',delim='apostrophe')
	open(koe,file=trim(output_file)//'-init')
	write(6,*) ' FRESCOX output to ',trim(output_file)//'*'
        write(6,*) 

	call freadf(ki,ko,koe,TMP,lmax,jbord,jump,.false.)	
        
	close(ki); close(koe)
        NANGL = (abs(THMAX) - THMIN)/abs(THINC) + 1 + 0.5
	ntheory_pts = max(mdl,NANGL)
	gettheoryplot = .false.
        pelii = peli
	  if(num_energies==0.and.elab(1)>0.) then
	    num_energies = 1
	    energy_list(1) = elab(1)
            pel_list(1) = pelii
	    endif
	  if(num_energies>0) then
	     write(6,15) num_energies,
     x              (energy_list(ien),pel_list(ien),ien=1,num_energies)
   15	  format(' Calculate scattering at ',i5,' energies'/
     x             (1x,10(f11.6,',',i1)))
	     write(6,151) (energy_count(ien),ien=1,num_energies)
 !151	  format('           Number of  data points'/(1x,10i8))
 151	  format(1x,10i8)
             else
	     write(6,*) ' No scattering energies'
	     endif
        write(6,*) 
        allocate (theoryplot(ntheory_pts,num_energies,datasets))
        write(6,*) 'grace:',grace
        if(.not.grace) then
           type_tag = '@TYPE'
           pl_suffix = ''
        else
           type_tag = '@type'
           pl_suffix = '.agr'
        endif

C                   Find dof for chisq fit
        ndof = dof(datalen,datasets,srch_step,nvars)
	
	ki = 3
	ko = 308
	open(ko,file=output_file)
	koi= ko
	written(ko) = .true.
!	rewind ki
	noerror = .true.
	final = .true.
	final = .false.
	MAXCH = 0; MAXICH = 0  	 !	no arrays allocated
       	
C			DO IT!

	call fr

	final = .false.
	if(undef) then
	  do ip=1,nvars
	  do iof=ko,6,6-ko
           if(abs(srch_value(ip)-nul)<1e-3) then
	      write(iof,10482) ' ',ip,srch_name(ip)
	    else
	      write(iof,1048) ' ',ip,srch_name(ip),srch_value(ip)
	    endif
	  enddo
	  enddo
	 endif
	 totchisq = sum(data_chisq(1:datasets))
	if(abs(totchisq/ndof)<1e5) then
	write(6,1040) totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	else
	write(6,10401) totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	endif
	write(6,10402) totchisq,
     x              (data_chisq(id),id=1,datasets)
!	write(6,10401) totchisq,
!    x              (data_chisq(id),id=1,datasets)
1040	format(/'  Total ChiSq/df=',f11.4,:' from ',10f10.3/
     x      (:'#',33x,10f10.2))
10401	format(/'  Total ChiSq/df=',1p,e11.4,:' from ',10e10.3/
     x      (:'#',33x,10e10.3))
!10402	format(/'  Total ChiSq   =',1p,e11.4,:' from ',10e10.3/
!     x      (:'#',33x,10e10.3))
10402	format(/'  Total ChiSq   =',f11.1,:' from ',10f10.1/
     x      (:'#',33x,10f10.1))
1041	format(a1,' Variable ',i5,':',a15,'=',g13.5,
     x         ': ChiSq/df=',f10.3,:' from ',10f10.3,/
     x                             (:'#',26x,6f10.3))
10411	format(a1,' Variable ',i5,':',a15,'=',g13.5,
     x         ': ChiSq/df=',1p,e10.3,:' from ',10e10.3,/
     x                                (:/'#',26x,6e10.3))
1042	format(/a1,' ChiSq/df=',f11.3,:' from ',6f10.3,/
     x                              (:'#',27x,6f10.3))
10421	format(/a1,' ChiSq/df=',1p,e11.3,:' from ',6e10.3,/
     x                                 (:'#',27x,6e10.3))
10422	format(/a1,' ChiSq/df=',f11.3,' part ',:,6f10.3,/
     x                              (:'#',27x,6f10.3))
!      write(6,10425) nparameters,nvarsv,ndatapoints,ndof
10425 format('  ',i4,' parameters (',i4,' variable) and',i7,
     x       ' data points, so dof=',i7)

        do id=1,ndatasets   ! Record first-time conversions of data to absolute units from relative
         idir = data_idir(id)
           if(idir==2) then  !  convert to ratio to Rutherford, the first time
              data_idir(id) = 1
              endif
           if(idir==-1.and.data_type(id)<=2) then  !  convert cross sections to absolute, the first time
              data_idir(id) = 0
              write(191,*) ' set idir ',id,' to zero',IK
              endif
        enddo


	do ip=1,nvars
!	if(srch_kind(ip)==2) then
!	  write(6,1043) ip,srch_name(ip),srch_afrac_overlap(ip)
!1043	  format(' Note: variable ',i3,'=',a15,' is overlap ',a80)
!	endif
!			Give variable names to Minuit
	call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ierflg)
	if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
	enddo
!
!####### INTERACTIVE SECTION

100	continue
	write(6,advance='no',fmt='(''sfrescox> '')')

	read(5,'(a200)',end=999) cmd
	if(cmd(1:2)=='EX'.or.cmd(1:2)=='ex') then
          write(6,'(a)') '>'//trim(cmd)
	  stop

	else if(cmd(1:1)=='Q'.or.cmd(1:1)=='q'
     x	    .or.cmd(1:1)=='V'.or.cmd(1:1)=='v') then
          write(6,'(a)') '>'//trim(cmd)
	  do ip=1,nvars
	  do iof=ko,6,6-ko
	   if(srch_error(ip)>1e-20.or.noerror) then
	    t = abs(srch_value(ip))
            if(abs(srch_step(ip))>1e-10
     x         .or. cmd(1:1)=='Q'.or.cmd(1:1)=='q') then
	     if(t<1e-3.or.t>1e3) then
	      write(iof,10481) ' ',ip,srch_name(ip),srch_value(ip),
     x		srch_step(ip),srch_error(ip),srch_rwa(ip)
	     else
	      write(iof,1048) ' ',ip,srch_name(ip),srch_value(ip),
     x		srch_step(ip),srch_error(ip),srch_rwa(ip)
	     endif
	    endif
     	   else
           if(cmd(1:1)=='Q'.or.cmd(1:1)=='q')
     x	    write(iof,1049) ' ',ip,srch_name(ip),srch_value(ip)
	   endif
1048	  format(a1,'   Var ',i5,'=',a15,' value ',f12.6,:,
     x          ', step ',f8.4,', error ',f8.4,' rwa =',L1)
10481	  format(a1,'   Var ',i5,'=',a15,' value ',1p,e12.4,:,
     x          ', step ',g8.1,', error' ,g9.2,' rwa =',L1)
10482     format(a1,'   Var ',i5,'=',a15,' HAS NO VALUE ASSIGNED!')
1049	  format(a1,'   Var ',i5,'=',a15,' value ',f12.6,' fixed')
	  enddo
	  enddo

        else if(cmd(1:3)=='CHA'.or.cmd(1:3)=='cha') then
          write(6,'(a)') '>'//trim(cmd)
          do iof=ko,6,6-ko
 	  write(iof,1042) '#',totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	  enddo
          do ip=1,nvars
          if(abs(srch_value(ip)-srch_vinit(ip))>1e-20) then
          do iof=ko,6,6-ko
	      sigdev = (srch_value(ip)-srch_vinit(ip))/srch_error(ip)
	   if(srch_kind(ip)/=4)  then
              write(iof,10484) ' ',ip,srch_name(ip),srch_value(ip),
     x          srch_error(ip),srch_vinit(ip),sigdev
           else
	      icch=srch_r_ic(ip); iach=srch_r_ia(ip)
              jch=srch_r_jch(ip); sch=srch_r_sch(ip); 
              lch=srch_r_lch(ip); channel=srch_r_ch(ip)
              write(iof,11484) ' ',ip,srch_name(ip),srch_value(ip),
     x          srch_error(ip),srch_vinit(ip),sigdev
     x          ,icch,iach,sch,lch
            endif
10484	  format(a1,'   Var ',i5,'=',a15,' value ',f15.6,' +/-',f12.6,
     x          ', initially ',f15.6,' (dev =',f8.2,' sig)')
11484	  format(a1,'   Var ',i5,'=',a15,' value ',f15.6,' +/-',f12.6,
     x          ', initially ',f15.6,' (dev =',f8.2,' sig)',
     x          ', IC,IA =',2i2,'; S,L =',f4.1,i3)
          enddo
           endif
          enddo

	else if(cmd(1:3)=='SET'.or.cmd(1:3)=='set') then
          write(6,'(a)') '>'//trim(cmd)
	  read(cmd(4:100),*,err=980,end=981) ip,val
	  if(ip<1.or.ip>nvars) go to 982
	  do iof=ko,6,6-ko
	  write(iof,1050) ip,srch_name(ip),val
1050	  format(/' SET',i3,'=',a15,' to ',g12.4,:,' from',g12.4)
	  enddo
	  srch_value(ip) = val
	  call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ierflg)
	  if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
          if(hasEshift) call energylist  ! remake union list of energies
          ndof = dof(datalen,datasets,srch_step,nvars)

	 call fr
	 totchisq = sum(data_chisq(1:datasets))
         if(datasets>1) then
	  if(abs(totchisq/ndof)<1e6) then
 	   write(6,1041) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	   else
 	   write(6,10411) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	  endif
	 else
	  if(abs(totchisq/ndof)<1e6) then
 	    write(6,1041) ' ',ip,srch_name(ip),val,totchisq/ndof
	    else
 	    write(6,10411) ' ',ip,srch_name(ip),val,totchisq/ndof
	    endif
	 endif
	  if(abs(totchisq/ndof)<1e6) then
 	   write(ko,1041) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	   else
 	   write(ko,10411) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	   endif


	else if(cmd(1:3)=='REA'.or.cmd(1:3)=='rea') then
          write(6,'(a)') '>'//trim(cmd)
	  read(cmd(5:100),*,err=980,end=981) plot_file
	  write(6,*) ' Replace parameters by those in plot file ',
     x        '<'//trim(plot_file)//'>'
	  open(304,file=plot_file)
	  if(index(plot_file,'snap')>0) then
10377	    read(304,10382,end=10380) ic,ip,chi
	    write(6,10379) ic,chi
10379	    format('  Reading snap at #',i6,' with Chisq =',1p,e13.4)
	    read(304,10381,end=10380) (srch_value(ip),ip=1,nvars)
	    go to 10377
10380	    do ip1=1,nvars
!      	    write(6,10485) ip1,name,srch_name(ip1),srch_value(ip1)
	      call MNPARM(ip1,srch_name(ip1),
     x          srch_value(ip1),srch_step(ip1),
     x		srch_minvalue(ip1),srch_maxvalue(ip1),ierflg)
	      if(ierflg>0) write(6,*) 'MNPARM for ',ip1,
     x        			      ': IERFLG =',ierflg
            enddo
10381	 	format(5e12.5)
10382	 	format(2i6,1p,e12.4,1x,1p,10e10.2)
	  else
	  do ip=1,nvars
	  read(304,10483,err=102,end=102)  ip1,name,val
          if(ip1==0) exit
10483	  format(1x,7x,i5,1x,a15,7x,e12.4)
!	  do iof=ko,6,6-ko
	   iof=ko
	  if(abs(srch_value(ip1)-val)>1d-12) 
     x 	    write(iof,10485) ip1,name,srch_name(ip1),val,srch_value(ip1)
10485	  format(' READ',i3,'=',2a15,' to ',g12.4,:,' from',g12.4)
!	  enddo
          !!write(6,*) ' changing var',ip1,'=',srch_value(ip1),'to',val
	  srch_value(ip1) = val
	  call MNPARM(ip1,srch_name(ip1),srch_value(ip1),srch_step(ip1),
     x		srch_minvalue(ip1),srch_maxvalue(ip1),ierflg)
	  if(ierflg>0) write(6,*) 'MNPARM for ',ip1,': IERFLG =',ierflg
	  enddo
	  endif
	  go to 103
102	  write(6,*) ' Read plot file ended'
103	 close(304)
          if(hasEshift) call energylist  ! remake union list of energies
          ndof = dof(datalen,datasets,srch_step,nvars)
	 call fr
	 totchisq = sum(data_chisq(1:datasets))
!	 close(1)
         if(datasets>1) then
	  if(abs(totchisq/ndof)<1e6) then
 	   write(6,1042) ' ',totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	   else
 	   write(6,10421) ' ',totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	  endif
	 else
	  if(abs(totchisq/ndof)<1e6) then
 	    write(6,1042) ' ',totchisq/ndof
	    else
 	    write(6,10421) ' ',totchisq/ndof
	    endif
	 endif
	  if(abs(totchisq/ndof)<1e6) then
 	   write(ko,1042) ' ',totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	   else
 	   write(ko,10421) ' ',totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	   endif

	else if(cmd(1:3)=='FIX'.or.cmd(1:3)=='fix') then
          write(6,'(a)') '>'//trim(cmd)
	  read(cmd(4:100),*,err=980,end=981) ip
	  if(ip<1.or.ip>nvars) go to 982
	  do iof=ko,6,6-ko
	  write(iof,10502) ip,srch_name(ip)
10502	  format(/'  FIX',i3,'=',a15)
	  enddo
	  srch_step(ip) = 0.0
	  call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ierflg)
	  if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg

        else if(cmd(1:5)=='FXWID'.or.cmd(1:5)=='fxwid') then
          write(6,'(a)') '>'//trim(cmd)
          write(6,*) 'Fix all R-matrix widths:'
         do ip=1,nvars
          if(srch_kind(ip)==4) then
          do iof=ko,6,6-ko
          write(iof,10502) ip,srch_name(ip)
          enddo
          srch_step(ip) = 0.0
          call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x          srch_minvalue(ip),srch_maxvalue(ip),ierflg)
          if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
         endif
         enddo

        else if(cmd(1:5)=='FXPBG'.or.cmd(1:5)=='fxpbg') then
          write(6,'(a)') '>'//trim(cmd)
          write(6,*) 'Fix all background poles (with BG in name):'
         do ip=1,nvars
          if(index(srch_name(ip),'BG')>0.and.srch_kind(ip)==3) then
          do iof=ko,6,6-ko
          write(iof,10502) ip,srch_name(ip)
          enddo
          srch_step(ip) = 0.0
          call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x          srch_minvalue(ip),srch_maxvalue(ip),ierflg)
          if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
         endif
         enddo

        else if(cmd(1:4)=='FXBG'.or.cmd(1:4)=='fxbg') then
          write(6,'(a)') '>'//trim(cmd)
          write(6,*) 'Fix all background widths (with BG in name):'
         do ip=1,nvars
          if(index(srch_name(ip),'BG')>0) then
          do iof=ko,6,6-ko
          write(iof,10502) ip,srch_name(ip)
          enddo
          srch_step(ip) = 0.0
          call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x          srch_minvalue(ip),srch_maxvalue(ip),ierflg)
          if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
         endif
         enddo

        else if(cmd(1:5)=='FXNOR'.or.cmd(1:5)=='fxnor') then
          write(6,'(a)') '>'//trim(cmd)
          write(6,*) 'Fix all dataset scaling norms:'
         do ip=1,nvars
          if(srch_kind(ip)==5) then
          do iof=ko,6,6-ko
          write(iof,10502) ip,srch_name(ip)
          enddo
          srch_step(ip) = 0.0
          call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x          srch_minvalue(ip),srch_maxvalue(ip),ierflg)
          if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
         endif
         enddo

        else if(cmd(1:5)=='FXSHF'.or.cmd(1:5)=='fxshf') then
          write(6,'(a)') '>'//trim(cmd)
          write(6,*) 'Fix all dataset scaling norms:'
         do ip=1,nvars
          if(srch_kind(ip)==6) then
          do iof=ko,6,6-ko
          write(iof,10502) ip,srch_name(ip)
          enddo
          srch_step(ip) = 0.0
          call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x          srch_minvalue(ip),srch_maxvalue(ip),ierflg)
          if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
         endif
         enddo

        else if(cmd(1:5)=='FXRPE'.or.cmd(1:5)=='fxrpe') then
          write(6,'(a)') '>'//trim(cmd)
          write(6,*) 'Fix all R-matrix pole energies:'
         do ip=1,nvars
          if(srch_kind(ip)==3.and. srch_step(ip)>0.0) then
          do iof=ko,6,6-ko
          write(iof,10502) ip,srch_name(ip)
          enddo
          srch_step(ip) = 0.0
          call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x          srch_minvalue(ip),srch_maxvalue(ip),ierflg)
          if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
         endif
         enddo

	else if(cmd(1:4)=='STEP'.or.cmd(1:4)=='step') then
          write(6,'(a)') '>'//trim(cmd)
	  step = 0.01
	  read(cmd(5:100),*,err=980,end=981) ip,step
	  if(ip<1.or.ip>nvars) go to 982
	  do iof=ko,6,6-ko
	  write(iof,10503) ip,srch_name(ip),step
10503	  format(/' STEP ',i3,'=',a15,' to ',g12.4)
	  enddo
	  srch_step(ip) = step
	  call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ierflg)
	  if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg

       else if(cmd(1:4)=='ELIM'.or.cmd(1:4)=='elim') then
          write(6,'(a)') '>'//trim(cmd)
	  read(cmd(5:100),*,err=980,end=981) esmin,esmax
	  esmin = max(esmin,scatmin)  ! global
	  esmax = min(esmax,scatmax)
          call  energylist ! (ndatasets,esmin,esmax)
	  write(6,10508) esmin,esmax,num_energies
10508     format(' Restricting data energy range:',f10.4,' to',f10.4,
     x            ' MeV(lab), so ',i5, ' energies')

       else if(cmd(1:5)=='POLES'.or.cmd(1:5)=='poles') then
          write(6,'(a)') '>'//trim(cmd)
          read(cmd(6:100),*,err=980,end=981) pmin,pmax
	  pmincm=pmin*ecmrat; pmaxcm=pmax*ecmrat
          write(6,10509) pmin,pmax,pmincm,pmaxcm
10509     format(' Vary only poles in energy range:',f10.4,' to',f10.4,
     x           ' MeV(lab),'/
     x           '           which is energy range:',f10.4,' to',f10.4,
     x           ' MeV(cm)')
	   id = 0; ifree=0
	  do ip=1,nvars
	   if(srch_kind(ip)==3) then
	     e = srch_value(ip)
	     if(e<pmincm .or. e>pmaxcm) then
 		srch_step(ip)=0
	 	id=id+1
          call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x          srch_minvalue(ip),srch_maxvalue(ip),ierflg)
          if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
	        rterm = srch_rterm(ip)
	        do iq=1,nvars
	         if(srch_kind(iq)==4.and.srch_rterm(iq)==rterm) then
 		   srch_step(iq)=0
          call MNPARM(iq,srch_name(iq),srch_value(iq),srch_step(iq),
     x          srch_minvalue(iq),srch_maxvalue(iq),ierflg)
          if(ierflg>0) write(6,*) 'MNPARM for ',iq,': IERFLG =',ierflg
	            
	           endif  ! need to fix width
	         enddo  ! iq
               else
                ifree = ifree+1
	       endif ! need to fix pole
	   endif   ! pole variable
	  enddo  ! ip
	  write(6,*) id,' poles outside range now fixed, and',
     x               ifree,' still free'

       else if(cmd(1:4)=='FREE'.or.cmd(1:4)=='free') then
          write(6,'(a)') '>'//trim(cmd)
	  write(6,*) ' Free all pole energies and widths'
	  do ip=1,nvars
	   if((srch_kind(ip)==3.or.srch_kind(ip)==4).and.
     x        abs(srch_step(ip))<1e-20) then
             srch_step(ip)=0.01
          call MNPARM(ip,srch_name(ip),srch_value(ip),srch_step(ip),
     x          srch_minvalue(ip),srch_maxvalue(ip),ierflg)
          if(ierflg>0) write(6,*) 'MNPARM for ',ip,': IERFLG =',ierflg
	    endif
	  enddo

       else if(cmd(1:5)=='ESCAN'.or.cmd(1:5)=='escan') then
          write(6,'(a)') '>'//trim(cmd)
	  read(cmd(6:100),*,err=980,end=981) emin,emax,de
	  ne = nint((emax-emin)/abs(de)) + 1
	  loge = de < 0.0
	  if(ne>maxen) write(6,*) ' Only room for ',maxen,' points'
	  ne = min(ne,maxen)
	  ! ne = max(ne,2)
	  if(ne>1) de = (emax-emin)/(ne-1)
	  write(6,1051)  emin,emax,ne,de,pelii
1051	  format(' EXCITATION function from ',f8.4,' to',f8.4,
     X		 ' in ',i5,' steps of',f8.5,' pel=',i2)
          is = num_energies; 
          esave(:) = energy_list(:); pelsave(:) = pel_list(:)
	  num_energies=ne
	  do i=1,ne
          energy_list(i) = emin + (i-1)*de ! linear
          pel_list(i) = pelii
           if(loge) then            ! log
             CF = (log(emax)-log(emin))/(emax-emin)
             CG =  log(emin) - CF * emin
             energy_list(i) = exp(CF*energy_list(i) + CG)
             endif            
	  enddo
	  write(6,15) num_energies,(energy_list(ien),pel_list(ien),
     x                                 ien=1,num_energies)
          ndof = dof(datalen,datasets,srch_step,nvars)
	    final = .true.; say=final
	  call fr
            num_energies=is; 
                 energy_list(:) = esave(:); pel_list(:) = pelsave(:)
	    final = .false.; say=final
	  write(6,*) ' File 71 has all phase shifts, and '
          write(6,*) ' file 40 the ',
     x      'fusion, reaction and non-elastic cross sections'
          write(6,*) ' file 35 the S-factors (CM energies)'
          write(6,*) ' file 75 the S-factors (LAB energies)'

       else if(cmd(1:4)=='SCAN'.or.cmd(1:4)=='scan') then
          write(6,'(a)') '>'//trim(cmd)
	  step=-1
	  read(cmd(5:100),*,err=980,end=981) ip,val1,val2,step
	  if(ip<1.or.ip>nvars) go to 982
	  if(step<=0) step=srch_step(ip)
	  write(6,1052) ip,srch_name(ip),val1,val2,step
1052	  format(/' SCAn',i3,'=',a15,' from ',g13.5,
     X       ' to ',g13.5,' in steps of ',g12.4)
	  oldval = srch_value(ip)
          ns = (val2-val1)/step+1
          do is=1,ns
	  val = val1 + (is-1)*step
	  srch_value(ip) = val
          if(hasEshift) call energylist  ! remake union list of energies
          ndof = dof(datalen,datasets,srch_step,nvars)

  	  call fr

	  totchisq = sum(data_chisq(1:datasets))
          write(2018,*) val,totchisq/ndof
          call flush(2018)
           if(datasets>1) then
	    if(abs(totchisq/ndof)<1e6) then
 	     write(6,1041) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	     else
 	     write(6,10411) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	     endif
	   else
	    if(abs(totchisq/ndof)<1e6) then
 	     write(6,1041) ' ',ip,srch_name(ip),val,totchisq/ndof
	     else
 	     write(6,10411) ' ',ip,srch_name(ip),val,totchisq/ndof
	     endif
	   endif
	    if(abs(totchisq/ndof)<1e6) then
  	     write(ko,1041) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	     else
  	     write(ko,10411) ' ',ip,srch_name(ip),val,totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	     endif
          enddo
	 srch_value(ip) = oldval
         if(hasEshift) call energylist  ! remake union list of energies
         ndof = dof(datalen,datasets,srch_step,nvars)
	 call fr
	 totchisq = sum(data_chisq(1:datasets))

        else  if(cmd(1:3)=='CHI'.or.cmd(1:3)=='chi') then
          write(6,'(a)') '>'//trim(cmd)
	 if(hasEshift) call energylist  ! remake union list of energies
         ndof = dof(datalen,datasets,srch_step,nvars)
          call fr
	 totchisq = sum(data_chisq(1:datasets))
	
         write(6,1054) totchisq,totchisq/ndof,ndof
1054     format(//' Total ChiSq=',1p,e14.4,0p,': ',f12.5,' per dof',
     x     ' from dof =',i6//
     x    ' Data # Type     E/A  V:Norm,Shift  ',
     x    ' ChiSq      ChiSq/N  ChiSq/Tot      N')
         s = 0.
         do id=1,datasets
         s = s +  data_chisq(id)
        if(data_type(id)==0) then
          write(6,1055) id, data_type(id),'E',data_energy(id),
     x  data_normvar(id),data_shiftvar(id),
     x  data_chisq(id),data_chisq(id)/max(1,datalen(id)),
     x  data_chisq(id)/ndof, datalen(id),trim(dat_file(id)),s
1055      format(' ',i5,i4,1x,A1,F9.3,2i5,f14.3,2f10.4,i8,2x,A40,f12.4)

        else if(data_type(id)==1) then
          write(6,1055) id, data_type(id),'m',0.,
     x  data_normvar(id),data_shiftvar(id),
     x  data_chisq(id),data_chisq(id)/max(1,datalen(id)),
     x  data_chisq(id)/ndof,datalen(id),trim(dat_file(id)),s

        else if(data_type(id)==2) then
          write(6,1055) id, data_type(id),'A',data_angle(id),
     x  data_normvar(id),data_shiftvar(id),
     x  data_chisq(id),data_chisq(id)/max(1,datalen(id)),
     x  data_chisq(id)/ndof,datalen(id),trim(dat_file(id)),s

        else if(data_type(id)==3) then
          write(6,1055) id, data_type(id),'A',data_energy(id),
     x  data_normvar(id),data_shiftvar(id),
     x  data_chisq(id),data_chisq(id)/max(1,datalen(id)),
     x  data_chisq(id)/ndof,datalen(id),trim(dat_file(id)),s

        else if(data_type(id)==6)then
          write(6,1056) id, data_type(id),data_value(id),
     x    data_par(id),data_chisq(id),data_chisq(id)/ndof,
     x    'for '//trim(data_reffile(id)) ! srch_name(data_par(id))
1056      format(' ',i5,i4,1x,F10.3,i5,' targ',f14.3,10x,f10.4,10x,A)

        else if(data_type(id)>=7 .and. data_type(id)<=8) then
          write(6,1056) id, data_type(id),data_value(id),
     x    data_term(id),data_chisq(id),data_chisq(id)/ndof,
     x    srch_name(data_par(id))

	endif
	 enddo

        else  if(cmd(1:4)=='SHOW'.or.cmd(1:4)=='show') then
          write(6,'(a)') '>'//trim(cmd)
	 if(hasEshift) call energylist  ! remake union list of energies
         ndof = dof(datalen,datasets,srch_step,nvars)
	 call fr
	 totchisq = sum(data_chisq(1:datasets))
	  do iof=ko,6,6-ko
	  do id=1,datasets
	  if(data_type(id)<6) then 
	  datanorm=1.0; dataEshift = 0.0
!		Adjust any datanorm search parameter!
          ip = data_normvar(id)
          if(ip>0) datanorm = datanorm * srch_value(ip)
          ip = data_shiftvar(id)
          if(ip>0) dataEshift=dataEshift + srch_value(ip)
!	   do ip=1,nvars
!           if(refer(ip,id)) then
!           if(srch_kind(ip)==5) datanorm = datanorm * srch_value(ip)
!           if(srch_kind(ip)==6) dataEshift=dataEshift + srch_value(ip)
!            endif
!	   enddo

          write(iof,*) 
          write(iof,*) ' Dataset ',id,' <'//trim(dat_file(id))//'>'
	  if(data_type(id)==-1) then
          write(iof,*) '   Order   Datum      Abs. error  Theory',
     x               '         Chi'
	  else if(data_type(id)==0) then
          write(iof,*) '   Angle   Datum      Abs. error  Theory',
     x               '         Chi'
	  else if(data_type(id)==1) then
          write(iof,*) '   Energy   Angle     Datum      Abs. error  ',
     x               '  Theory         Chi'
	  else if(data_type(id)==5) then
          write(iof,*) '   Bound state    Datum      Abs. error  ',
     x               'Theory         Chi'
     	  else
          write(iof,*) '  Energy   Datum      Abs. error  Theory',
     x               '         Chi'
          endif
	  tot = 0.0
	  do ip=1,datalen(id)
	   energy = data_energies(ip,id)+dataEshift
 	   if(esmin<=energy .and. energy<=esmax) then
            chi=(theoryvals(ip,id)-datanorm*datavals(ip,id))/
     x			(datanorm*dataerr(ip,id))
	  tot = tot + chi**2
	  if(data_type(id)==-1) then
	  write(iof,1061) nint(datangles(ip,id)),
     x          datanorm*datavals(ip,id),
     X          datanorm*dataerr(ip,id),theoryvals(ip,id),chi**2
	  else if(data_type(id)==0) then
	  write(iof,1062) datangles(ip,id),datanorm*datavals(ip,id),
     X          datanorm*dataerr(ip,id),theoryvals(ip,id),chi**2
     	  else if(data_type(id)==1) then
	  write(iof,1062) data_energies(ip,id)+dataEshift,
     x          datangles(ip,id),datanorm*datavals(ip,id),
     X          datanorm*dataerr(ip,id),theoryvals(ip,id),chi**2
	  else if(data_type(id)==5) then
	  write(iof,10621) bs_no(ip,id),datanorm*datavals(ip,id),
     X          datanorm*dataerr(ip,id),theoryvals(ip,id),chi**2
     	  else
	  write(iof,1062) data_energies(ip,id)+dataEshift,
     x          datanorm*datavals(ip,id),
     X          datanorm*dataerr(ip,id),theoryvals(ip,id),chi**2
	  endif
1061 	   format(1x,i8  ,3g12.5,8f12.4,3x,'"',a,'"  ',f8.2)
1062 	   format(1x,f10.5,3g12.5,8f12.4,3x,'"',a,'"  ',f8.2)
10621 	   format(1x,i8,7x,3g12.5,f10.4)
	  endif ! e in range
	  enddo ! ip
            write(iof,10622) tot,tot/max(1,datalen(id)),penalty
10622       format(' Tot Chisq=',1p,e12.4,0p,' Chisq/N=',f12.4,
     x             ',  Penalty=',f12.4)
          else if( data_type(id)==6 ) then
          if(data_type(id)/=data_type(max(id-1,1)) .or. id==1)  then
          write(6,*)
          write(6,*) '   Search param         Target     Absolute error'
     x           ,'   Parameter      chi^2'
            endif

          value = srch_value(data_par(id))
          chi = (datavals(1,id)-value)/dataerr(1,id)
	  write(6,13211) data_par(id),srch_name(data_par(id)),
     x                 datavals(1,id),dataerr(1,id),value,chi**2

          else if( data_type(id)==7.or.data_type(id)==8 ) then ! type>=7
          write(6,*)
          write(6,*) '   R-matrix term        Target     Absolute error'
     x           ,'   Parameter      chi^2'

          value = theoryvals(1,id)
          chi = (datavals(1,id)-value)/dataerr(1,id)
	  write(6,13211) data_term(id),srch_name(data_par(id)),
     x                 datavals(1,id),dataerr(1,id),value,chi**2
13211 	   format(1x,i3,': ',a15,2f12.4,4x,2f12.4)


	  endif 
	  enddo
	if(abs(totchisq/ndof)<1e6) then
	write(iof,1040) totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	else
	write(iof,10401) totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	endif
	write(iof,10402) totchisq,
     x              (data_chisq(id),id=1,datasets)
	  enddo	

        else  if(cmd(1:4)=='DONE'.or.cmd(1:4)=='done') then
          write(6,'(a)') '>'//trim(cmd)
          tag = ' '
          read(cmd,'(5x,a80)',err=980,end=981) tag
!         write(6,*) ' tag = <'//trim(tag)//'>'
          plot_file = 'done.sfresco'
          plot_file = trim(search_file)//'ed'
          if(tag(1:1).ne.' ') plot_file = trim(tag)
          open(304,file=plot_file,delim='apostrophe')
	  write(304,549) '=',trim(output_file)
549 	  format(' ''',a,''' ''',a,'''')
	  write(304,*)  nparameters,ndatasets
	  write(304,*) 
          open(305,file=input_file,status='old')  ! original FRESCO input file
	  rewind 305
550	  read(305,'(a)',END=551) cmd
	  write(304,'(a)') trim(cmd)
	  go to 550
551	  write(304,'(a/)') 'EOF'
	
	  do ip=1,nvars  ! write parameters for fit
	   name =  srch_name(ip	); kind=srch_kind(ip)
           term = srch_rterm(ip); step = srch_step(ip)
           reffile = srch_reffile(ip)
           valmin = srch_minvalue(ip); valmax = srch_maxvalue(ip)
	   t = srch_value(ip)
	      pline2 = srch_pline2(ip); col2 = srch_col2(ip)

           if(srch_kind(ip)==1 .and. pline2*col2==0)  then
              kp = srch_kp(ip); pline = srch_pline(ip); 
              col = srch_col(ip)
              call write_var1(ip,304,name,kind,kp,pline,col,t,step)

           else if(srch_kind(ip)==1 .and. pline2*col2>0)  then
              ratio2 = srch_ratio2(ip) 
              kp = srch_kp(ip); pline = srch_pline(ip); 
              col = srch_col(ip)
              call write_var12(ip,304,name,kind,kp,pline,col,t,step,
     x                        pline2,col2,ratio2)

           else if(srch_kind(ip)==2)  then
          	nafrac = srch_nafrac(ip)
              call write_var2(ip,304,name,kind,nafrac,t,step)

           else if(srch_kind(ip)==3)  then
              nopot = srch_nopot(ip) ; damp=srch_damp(ip)
	      jtot=srch_jtot(ip); par=srch_par(ip) 
              call write_var3(ip,304,name,kind,term,jtot,par,t,
     x                            nopot,damp,step)

	   else if(srch_kind(ip)==4)  then
	      icch=srch_r_ic(ip); iach=srch_r_ia(ip)
              jch=srch_r_jch(ip); sch=srch_r_sch(ip); 
              lch=srch_r_lch(ip); channel=srch_r_ch(ip)
              rwa=srch_rwa(ip)
              call write_var4(ip,304,name,kind,term,icch,iach,lch,t,
     x                      sch,jch,channel,step,rwa)

	   else if(srch_kind(ip)==5)  then
              dataset=srch_datanorm(:,ip) 
              call write_var5(ip,304,name,kind,dataset,reffile,t,step)

	   else if(srch_kind(ip)==6)  then
              dataset=srch_datanorm(:,ip) 
              call write_var6(ip,304,name,kind,dataset,reffile,t,step)

	   else if(srch_kind(ip)==7)  then
 	      damp=srch_value(ip)
              energy=srch_damp(ip)
              leff = srch_leff(ip)
              brune = srch_Brune(ip)
            if(abs(energy-nul)<.1) then
              call write_var7a(ip,304,name,kind,term,damp,step,valmin,
     x                        valmax)
            else
              call write_var7b(ip,304,name,kind,term,damp,step,valmin,
     x                        valmax,leff,energy,brune)
	    endif

           endif
	  enddo
          nvariables = 0
          do ip=1,nvars
           if(srch_step(ip)/=0.0) then
             nvariables = nvariables+1
             variables(nvariables)= ip
           endif
          enddo
          
        do id=1,ndatasets
        call print_data(304,id,data_type(id),dat_file(id),
     x    data_points(id),data_xmin(id),data_delta(id),data_lab(id),
     X    data_energy(id),data_idir(id),data_iscale(id),data_abserr(id),
     x    data_ic(id),data_ia(id),data_ib(id),data_kind(id),
     x    data_rank_k(id),data_rank_q(id),
     x    data_angle(id),data_jtot(id),data_par(id),data_ch(id),
     x    data_leg(id),data_pel(id),data_exl(id),data_labe(id),
     x    data_lin(id),data_lex(id),data_value(id),
     x    data_error(id),data_detinfo(id),dat_info(id),
     x    dat_eplots(id),dat_logy(id),data_reffile(id),data_term(id))
        enddo

       allocate (emat(nvariables,nvariables))
       emat(:,:) = 0.0
       call mnemat(emat,nvariables)
        if(maxval(abs(emat(:,:)))>0.0) then
          write(304,561) nvariables,nvariables
          write(304,562) variables(1:nvariables)
561       format(/' &Cov nvariables =',i5,'  variables(1:',i5,') = ')
562       format( 20i5 )
          write(304,564)
          write(6,*) nvariables,'variables of Covariance matrix'
          write(6,5625) (variables(ip),srch_step(variables(ip)),
     x                  ip=1,nvariables)
5625      format(1p,10(i4,e11.3,','))
          do ip=1,nvariables
           write(304,563)  ip,nvariables
563        format(' &Cov row=',i5,' emat(1:',i5,') = ')
           write(304,'(4x,1p,5e13.5)') (emat(ip,id1),id1=1,nvariables)
           write(304,564)
564        format('     /')
          enddo
        endif
       deallocate(emat)
        close(303); close(304)
	write(6,*) ' New search file written: ',trim(plot_file)
 
       data_file = trim(plot_file)//'-datafit'
       open(310,file=data_file,delim='apostrophe')
      dplotted = .false.
       do ip=1,nvars
       if(srch_kind(ip)>=5) then
       dataset=srch_datanorm(:,ip) 
       reffile = srch_reffile(ip)
       if(srch_kind(ip)==5) 
     x     call write_varlink5(310,dataset,reffile,srch_value(ip))
       if(srch_kind(ip)==6) 
     x     call write_varlink6(310,dataset,reffile,srch_value(ip))
       dplotted = dplotted .or. srch_kind(ip)==5 .or. srch_kind(ip)==6
        endif
       enddo
       close(310)


       data_file = trim(plot_file)//'-datafit.csv'
       ! open(310,file=data_file,delim='apostrophe')
       open(310,file=data_file,recl=len(cmd))
       write(310,'(a)') 'segment,dataset,filedir,source,EXFOR,points,'//
     x 'low angle bound,high angle bound,aframe,'//
     x 'low energy bound,high energy bound,eframe,'//
     x 'norm,shift,chi^2,chi^2/pt'
       write(310,"('int,str,str,str,int,int,',2('float,'),'str,',
     x          2('float,'),'str,',4('float,'))")
       write(310,565)
 565   format(',,,,,,',2('deg,'),'str,',2('MeV,'),'str,',',MeV,,,')
       do id=1,datasets
        if(data_type(id)<=2) then
        emin = minval(data_energies(1:datalen(id),id))
        emax = maxval(data_energies(1:datalen(id),id))
        amin = minval(datangles(1:datalen(id),id))
        amax = maxval(datangles(1:datalen(id),id))
        cmlab='com'; if(data_lab(id)) cmlab='lab'
        if(data_type(id)==2) then
          amin = data_angle(id)
          amax = data_angle(id)
	endif
       iof = index(dat_file(id),'/')
       if(iof>0) then 
        dir = dat_file(id)(1:iof)
        in_file = dat_file(id)(iof+1:)
       else
        dir = ''
        in_file = dat_file(id)
       endif
        tag = in_file
        I = index(in_file,'@')
        if(I>0) tag = in_file(1:I-1)//'.dat'

        datanorm = 1; dataEshift = 0.
        ipn = data_normvar(id) 
        ips = data_shiftvar(id) 
        n_reffile =''; s_reffile =''
        if(ipn>0) datanorm = srch_value(ipn)
        ! if(ipn>0) n_reffile= srch_reffile(ipn)
        if(ips>0) dataEshift= srch_value(ips)
        ! if(ips>0) s_reffile = srch_reffile(ips)
       
      write(cmd,*) id,',',trim(in_file),',',trim(dir),',',tag,',',
     x             ' ',',',datalen(id),',  '  ! EXFOR unknown here
!     x             ,n_reffile,',',s_reffile,',   '
 
      write(nums,566) amin,amax,cmlab,emin,emax,'lab',
     x    datanorm,dataEshift,data_chisq(id),data_chisq(id)/datalen(id)
 566  format(2(f10.4,',',f10.4,',',a3,','),f10.6,',',f10.6,',',
     x       f10.3,',',f9.4)
       do I=len(nums)-2,2,-1   ! remove trailing 0
         if(nums(I:I+1)=='0,') nums(I:len(nums)-1) = nums(I+1:len(nums))
       enddo
       cmd = trim(cmd)//trim(nums)
       do I=len(cmd)-1,2,-1   ! remove blanks
         if(cmd(I:I)==' ') cmd(I:500-1) = cmd(I+1:500)
       enddo
         if(cmd(1:1)==' ') cmd=cmd(2:500)
        write(310,'(a)') trim(cmd)
        endif
       enddo !id
       close(310)
      write(6,*) ' Data adjustment csv file written: ',trim(data_file)

        else  if(cmd(1:4)=='PLOT'.or.cmd(1:4)=='plot'
     x       .or.cmd(1:4)=='XMGR'.or.cmd(1:4)=='xmgr'
     x       .or.cmd(1:4)=='LINE'.or.cmd(1:4)=='line') then
          write(6,'(a)') '>'//trim(cmd)
     	  nodat = cmd(1:4)=='LINE'.or.cmd(1:4)=='line'
     	  xmgr =  cmd(1:4)=='XMGR'.or.cmd(1:4)=='xmgr'
	  tag = ' '
!	  read(cmd,'(5x,a80)',err=980,end=981) tag
 	  id1 = 1; id2=ndatasets
821	  read(cmd,*,err=822,end=822) wid,tag,id1,id2
!          go to 823
822	  write(6,*) ' file = <'//trim(tag)//'> for sets',id1,id2
	  if(wid(1:1)=='#') go to 100
!823	  continue
	  do k=0,min(2,maxrank)
	  plot_file = 'search.plot' // pl_suffix
          plot_file = trim(search_file)//'-plot' // pl_suffix
	  if(tag(1:1).ne.' ') plot_file = trim(tag)
	  if(k>0) plot_file = trim(plot_file)//char(k+ichar('0'))
	  open(304,file=plot_file)
	  do ip=1,nvars  ! remind parameters for fit
	   if(srch_error(ip)>1e-20.or.noerror) then
	    t = abs(srch_value(ip))
	    if(t<1e-3.or.t>1e3) then
	    write(304,10481) '#',ip,srch_name(ip),srch_value(ip),
     x		srch_step(ip),srch_error(ip),srch_rwa(ip)
	    else
	    write(304,1048) '#',ip,srch_name(ip),srch_value(ip),
     x		srch_step(ip),srch_error(ip),srch_rwa(ip)
	    endif
     	   else
	    write(304,1049) '#',ip,srch_name(ip),srch_value(ip)
	   endif
	  enddo
	 if(hasEshift) call energylist  ! remake union list of energies
		if(k==0) then
		gettheoryplot = .true.; final = .false.; say=final
           ! ndof = dof(datalen,datasets,srch_step,nvars)
	  call fr
		gettheoryplot = .false.; final = .false.; say=final
          totchisq = sum(data_chisq(1:datasets))
		endif
	if(totchisq/ndof+
     x      maxval(data_chisq(1:datasets)/datalen(1:datasets))<1e6) then
 	  write(304,1042) '#',totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	  else
 	  write(304,10421) '#',totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)
	  endif
 	  write(6,1042) ' ',totchisq/ndof,
     x              (data_chisq(id)/max(1,datalen(id)),id=1,datasets)

        if(id1/=1 .or. id2/=ndatasets) then   ! Calculate chisq/pt for this range
	 T = 0; I = 0
          do id=id1,id2
	   T = T + data_chisq(id)
           I = I + datalen(id)
	  enddo
 	  write(304,10422) '#',T/I,
     x              (data_chisq(id)/max(1,datalen(id)),id=id1,id2)
 	  write(6,10422) '#',T/I
	 endif
     
     !!!!!!!!!!!!!!!!!!!!!
	xmin=180; xmax=0; ymin=1e5; ymax=0.
	ig=-1
	maxien = num_energies  +1     ! excitation functions last
	plotted(:,:) = .false.
	eplots=.false.
	do id=1,datasets
	 eplots = eplots .or. dat_eplots(id)
	enddo
!	 write(6,*) ' maxien,eplots =',maxien,eplots

	  do 1068 ILEN=1,maxien
	  eplot = eplots !.and. energy_count(ILEN)>1  ! do not do angular distn for 1 point !
     x             .and. ILEN<=num_energies		      ! keep last ILEN for 1-point/energy data
! 	  write(6,*) ' ILEN =',ILEN,energy_list(ILEN),pel_list(ILEN),
!    x             ' eplot =',eplot,NANGL
	  if(eplot)  then
	   ig = ig+1
	 write(304,81) energy_list(ILEN),pel_list(ILEN)  ,ig,ig ,ig
81	  format(/'#  Graph for energy ',f8.3,' MeV @',i1,/'#'/
     x            '@WITH G',i1,/'@G',I1,' ON'/
     x            '@G',i1,' autoscale type AUTO'/
     x            '@    yaxis  tick major 1',/
     x            '@    yaxis  tick minor 1',/
     x            '@    yaxis  ticklabel prec 0',/
     x            '@    yaxis  ticklabel format power')
!	   write(6,81) energy_list(ILEN),pel_list(ILEN),ig,ig
	   write(304,1060) energy_list(ILEN),pel_list(ILEN)
1060	    format('@subtitle "Energy = ',f8.3,' @',i2,'"',/,
     x             '@subtitle size 0.7'/'@legend OFF')
	  else
	   if(ILEN<=num_energies) go to 1068
	   ig = ig+1
	 write(304,812)  !ig,ig
812	  format(/'#  Graph for single-energy points'/'#')
!    x            '@WITH G',i1,/'@G',I1,' ON')
        if(.not.stdin) then
	 write(304,1063) 'Search file: '//trim(search_file),
     x 	                  'Frescox input: '//trim(input_file)
        else
	 write(304,1063) 'Search file: '//trim(search_file),''
        endif
	  endif
1063	    format('@subtitle "',a,'; ',a,'"',/,
     x             '@subtitle size 0.7'/'@legend ON'/,
     x             '@legend char size 0.7')
!	  write(6,*) ' ILEN =',ILEN,' eplot,ig =',eplot,ig

	  is = -1; idp = -1; logy=.false.
!	  do 899 id=1,datasets
	  do 899 id=id1,id2
	    if((data_rank_k(id)==k.or.(k==0.and.data_type(id)==5))
     x       .and.data_type(id)/=6
     x       .and.(ILEN>num_energies.or.data_type(id)<3)
     x           ) then
	  dplotted=.false.
	  datanorm=1.0; dataEshift = 0.
!		Adjust any datanorm or dataEshift search parameter!
          ip = data_normvar(id)
          if(ip>0) datanorm = datanorm * srch_value(ip)
          ip = data_shiftvar(id)
          if(ip>0) dataEshift=dataEshift + srch_value(ip)
!	   do ip=1,nvars
!           if(refer(ip,id)) then
!           if(srch_kind(ip)==5) datanorm = datanorm * srch_value(ip) 
!           if(srch_kind(ip)==6) dataEshift=dataEshift + srch_value(ip)
!            endif
!	   enddo

	  if(.not.nodat) then
	   write(304,'(a)') type_tag //' xydy'
	   elast = data_energies(1,id)+dataEshift
	   do 83 ip=1,datalen(id)
	    e = datanorm*dataerr(ip,id)
	    d = datanorm*datavals(ip,id)
	    if(e>abs(d)) e = .9*abs(d)
	    x = datangles(ip,id)
	    if(data_type(id)==5)  x = bs_no(ip,id)

 	    if(eplot.and.abs(data_energies(ip,id)+dataEshift-
     x                        energy_list(ILEN))>1e-6) go to 83
	    if(ILEN==num_energies+1.and.plotted(ip,id)) go to 83  ! skip points that are plotted for ILEN >0 !
	
	    plotted(ip,id) = .true.
	      if(data_type(id)==1) then
              if(abs(elast-data_energies(ip,id)-dataEshift)>1e-6) then
!			if(ILEN==1) write(304,*) 'NEW energy'
                elast = data_energies(ip,id) + dataEshift
 	      endif
 	      endif
	 
	    info = data_info(ip,id)
	    if(data_type(id)==-3) then
	 		x = data_energies(ip,id)+dataEshift 	! x-axis = data energy
	    		write(304,1062) x,d,e,info ! x-axis = energy
	  		dplotted=.true.
	    	endif
	    if(data_type(id)==-2) then
	 		x = data_energies(ip,id)+dataEshift 	! x-axis = data energy
	    		write(304,1062) x,d,e,info ! x-axis = energy
	  		dplotted=.true.
	    	endif
	    if(data_type(id)==-1) then
	    		write(304,1061) nint(x),d,e,info ! x-axis = Legendre order
	  		dplotted=.true.
	    	endif
	    if(data_type(id)==0) then
	    		write(304,1062) x,d,e,info
	  		dplotted=.true.
	    	endif
	    if(data_type(id)==1) then
	       if(eplot) then
	       	write(304,1062) x,d,e,info ! x-axis = data angle
	       	dplotted=.true.
	       	endif
	       if(ILEN==num_energies+1)  then
	 		x = data_energies(ip,id)+dataEshift 	! x-axis = data energy
			write(304,1062) x,d,e,info,datangles(ip,id)
	  		dplotted=.true.
	        endif
	       endif
	    if(data_type(id)==5) then
	    		write(304,1062) bs_no(ip,id)+0.,d,e,info
	    		dplotted=.true.
		endif
	    if(data_type(id)>1.and.data_type(id)<5) then
	      energy = data_energies(ip,id)+dataEshift
 	      if(esmin<=energy .and. energy<=esmax) then
	      sfa = 1.
	      if(data_idir1(id)==-1)  then
		  ecmi = (data_energies(ip,id)+dataEshift) * ecmrat
		  etai = etarat / SQRT(ecmi)
		  sfa = exp(2d0*PI*etai) *  ecmi
		 endif
    		 write(304,1062) data_energies(ip,id)+dataEshift,d*sfa,e*sfa
		  d = d*sfa; x = data_energies(ip,id)+dataEshift
	  		dplotted=.true.
	    xmin = min(xmin,x)
	    xmax = max(xmax,x)
 	    ymin = min(ymin,d)
	    ymax = max(ymax,d)
!	 	write(6,'(i1,6f10.4)') 1,xmin,xmax,ymin,ymax,x,d
            endif ! energy in range
	   endif  ! type = 1,2,3,4
83	   continue ! ip in dataset id
	   if(.not.dplotted) then
		write(6,*) ' No data plotted, so skip graph ',ig
		write(304,'(a,i3)') '# No data plotted, so skip graph ',ig
		ig=ig-1
                ig = max(ig,0)
		go to 899 !  if data requested, but none actually plotted, skip theory part too!
	 	endif
           endif  ! nodat
	   
	    is=is+1; idp=idp+1 ; 
            icol=mod(idp+1-1,15)+1
            isym = mod(idp,7)+2
	    logy = logy .or. dat_logy(id)
	   if(.not.nodat.or.all_legends) then
            cmlab='CM '; if(data_lab(id)) cmlab='LAB'
	    info = ''
	    ! if(len(trim(dat_info(id)))>0) info = ' @ '//dat_info(id)
	    !! if(dat_info(id).ne.' ') info = ' @ '//dat_info(id)
              t = data_chisq(id)/max(1,datalen(id))
              cmd = dat_file(id)
              dir = ''
              i = index(cmd,'/')
              !# write(6,*) 'cmd:',cmd
              if(i>0) then
                 dir = cmd(1:i)
                 cmd = cmd(i+1:)
                 endif
              !# write(6,*) '/ at',i,'dir:',trim(dir),', cmd:',trim(cmd)
              i = index(cmd,'.data')
              if(i>0) cmd = cmd(1:i-1)
              !# write(6,*) '.data at',i,', cmd:',trim(cmd)
            if(abs(datanorm-1.0)>1e-3) then
            red = '+'
            if(datanorm<1.0) red='-'
             if(abs(dataEshift)<0.002) then
	      write(304,10628) '@'//pre_is,is,post_is,id,trim(cmd)
     x         ,red,abs(datanorm-1.)*100.,t
             else
             shf = '+'
             if(dataEshift<1.0) shf='-'
	      write(304,10629) '@'//pre_is,is,post_is,id,trim(cmd),
     x       red,abs(datanorm-1.)*100.,shf,nint(abs(dataEshift*1000)),t
             endif
            else
	     write(304,10630) '@'//pre_is,is,post_is,id,trim(cmd),t
            endif
10628       format(  '@subtitle size 0.7'/'@legend ON',/,
     x             a15,i0,a,' "',i3,': ',a,
     x             ' ',a1,f5.1,'% @',f5.1,'"')
10629       format(  '@subtitle size 0.7'/'@legend ON',/,
     x             a15,i0,a,' "',i3,': ',a,
     x             ' ',a1,f5.1,'% ',a1,i4,' keV',' @',f5.1,'"')
10630       format(  '@subtitle size 0.7'/'@legend ON',/,
     x             a15,i0,a,' "',i3,': ',a,' @',f5.1,'"')
            if(.not.nodat) then
	    write(304,10631) is,is,is,isym,is,icol ,is,is
10631	    format('@ s',i0,' linestyle 0',/
     x             '@ s',i0,' errorbar length 0.10'/
     x             '@ s',i0,' symbol ',i2/
     x             '@ s',i0,' color ',i2/
     x             '@ s',i0,' symbol fill 1'/
     x             '@ s',i0,' symbol size 0.3')
            if(grace) write(304,10632) is,icol,is,icol,is,icol
10632       format('@ s',i0,' symbol fill color ',i0 /
     x             '@ s',i0,' symbol color ',i0 /
     x             '@ s',i0,' errorbar color  ',i0)
 	    endif ! nodat
	  !endif ! nodat or all_legends
	  write(304,'(a)') '&'   ! end of dataset

	  if (.not.nodat) is=is+1
	  endif ! nodat or all_legends
	  ! endif

            cmlab='  '; 
	    if(data_type(id)<=1) then  ! theoryplot only available in cm..
	       cmlab='CM'
	       if(data_lab(id)) then
	        write(6,*)  '  Note: in search.plot, data is LAB,',
     x             ' continuous theory plot is CM, so not given!'
	       endif
	    endif

	   iwd=1 
	   if(boldplot) then 
        	iwd=2
!	   	icol=mod(icol+4,16)+1
	   	icol=1
             endif
            write(304,1064)   type_tag,cmlab,id,is,icol,is,iwd
1064	    format(A5,  ' xy',/,'#     ',a2,' Dataset ',i3/
     x        '@ s',i0,' color ',i2/'@ s',i0,' linewidth ',i2)
	  if(data_type(id)==-1) then
	     do ia=1,datalen(id)
	      write(304,1061) nint(datangles(ia,id)),theoryvals(ia,id)
	     enddo

	  else if(data_type(id)==0) then
!	    write(1304,*) ' ILEN,id =',ILEN,id
!	    write(1304,*) theoryplot(ia,1,id),theoryplot(ia,ILEN,id)
	   if(.not.data_lab(id).and..not.valsonly) then
             write(304,"('#Theoryplot')")
	     do ia=1,NANGL
	      write(304,1062) thmin+(ia-1)*abs(THINC),
     x				theoryplot(ia,ILEN,id)
!	      write(1304,1062) thmin+(ia-1)*abs(THINC),
!    x				theoryplot(ia,ILEN,id)
	     enddo
             write(304,"('#Theoryvals')")
	     do ia=1,datalen(id)
	      write(304,1062) datangles(ia,id),theoryvals(ia,id)
	     enddo
	   else
	     do ia=1,datalen(id)
	      write(304,1062) datangles(ia,id),theoryvals(ia,id)
	     enddo
	   endif
	   
	  else if(data_type(id)==1) then
	  
	   if(eplot) then ! plot only ien=ILEN 
		write(304,'(a,i3)')  '## NANGL =',NANGL
	   	ien=ILEN
	      do ia=1,NANGL
	        x = thmin+(ia-1)*abs(THINC)
	        d = theoryplot(ia,ien,id)
 	        write(304,1062) x,d !,energy_list(ien)
	    	  xmin = min(xmin,x); xmax = max(xmax,x)
 	    	  ymin = min(ymin,d); ymax = max(ymax,d)
!	 	write(6,'(i1,4f10.4)') 2,xmin,xmax,ymin,ymax
              enddo ! ia	   
	   else if(ILEN>num_energies) then ! ~eplot    ! plot as function of energy, for energies and angles in this dataset
           do ien=1,num_energies ! now select those single-energies appearing in data_energies(:,id)
            if(energy_count(ien)==1) then
		ia = 0
              do ip=1,datalen(id)
               if(abs(data_energies(ip,id)+dataEshift-energy_list(ien))
     x         <1e-6) ia = nint((datangles(ip,id)-THMIN)/abs(THINC))+1  ! nearest theory angle (for now)
              enddo !ip
	      if(ia>0.and.ia<=NANGL)  then
		x = energy_list(ien)
		d = theoryplot(ia,ien,id)
		write(304,1062) x,d
	    	  xmin = min(xmin,x); xmax = max(xmax,x)
 	    	  ymin = min(ymin,d); ymax = max(ymax,d)
!	 	write(6,'(i1,4f10.4)') 3,xmin,xmax,ymin,ymax
		endif
	    endif ! energy_count(ien)>1
           enddo ! ien
	   endif ! eplot or last ILEN
	  
	  
!	   do ien=1,num_energies
!	   if(.not.eplot.or.ien==ILEN) then
!	   do ia=1,NANGL
! 	    write(304,1062) thmin+(ia-1)*abs(THINC),
!     x 				theoryplot(ia,ien,id) !,energy_list(ien)
!         enddo ! ia
!	     if(ien<=num_energies.and..not.eplot) 
!     x            write(304,*) 'Energy',energy_list(ien+1)
!	   endif ! ien==ILEN or ~eplot
!	   enddo ! ien
	   
	   
	   
	  else if(data_type(id)==5) then
	   do ia=1,datalen(id)
	    write(304,10621) bs_no(ia,id),theoryvals(ia,id)
	   enddo !ia

	  else if(data_type(id)==6) then
c		 no nothing here

	  else   ! OTHER type
	   do ia=1,datalen(id)
	    energy = data_energies(ia,id)+dataEshift
 	    if(esmin<=energy .and. energy<=esmax) then
	      sfa = 1.
	      if(data_idir1(id)==-1)  then
		  ecmi = (data_energies(ia,id)+dataEshift) * ecmrat
		  etai = etarat / SQRT(ecmi)

!		  ecmi = (data_energies(ia,id)+dataEshift) *
!     x               MASS(3-LIN,LAB) / (MASS(2,LAB)+MASS(1,LAB))
!		  etai = ETACNS * MASS(2+1,pel) * MASS(2+2,pel)
!     X                          * SQRT(RMASS(pel)/ ecmi)

		  sfa = exp(2d0*PI*etai) *  ecmi
	          endif
	    write(304,1062) data_energies(ia,id)+dataEshift,
     x                      theoryvals(ia,id)*sfa
	    endif ! energy in range
	   enddo ! ia
	  endif
	   write(304,*) '&'
	   endif  ! correct k
	  
!	   if(eplot) then
	 if(abs(xmax-xmin)<1e-5) then
	   	xmin = xmin-1; xmax=xmax+1
	 	endif
	 if(ymin<=0.) logy=.false.
	 if(xmin<xmax) then
	        if(ILEN<=num_energies.and.xmin<xmax) then
      	  	 write(6,841) ig,xmin,xmax,ymin,ymax,energy_list(ILEN),logy
	  	else
	  	 write(6,8411) ig,xmin,xmax,ymin,ymax,logy
	  	endif
	  endif
841	   format(/'  Graph ',i3,': Amin,Amax =',2f10.3,
     x		' ymin,ymax =',2f10.3,' for E =',f8.3, L3)
8411	   format(/'  Graph ',i3,': Emin,Emax =',2f10.3,
     x		' ymin,ymax =',2f10.3,L3)
	 if(xmax>xmin .and. ymax>ymin) then
	    if(logy) then
	     if(ymin>0.) ymin = 10.**int(log10(ymin))
	     if(ymax>0.) ymax = 10.**int(log10(ymax)+1)
	    endif
	   if(xmax-xmin>20) then
   	   	xmin  = 15*int(xmin/15)
	   	xmax  = 15*int(xmax/15)+15
	   	xtmaj = (xmax-xmin)/4
	   	xtmaj  = 15*int(xtmaj/15)
	   else
	   	xtmaj = (xmax-xmin)/4
	   endif
	   	ytmaj = (ymax-ymin)/4
           	xtmin = xtmaj/3
           	ytmin = ytmaj/3
!	 	ymin=ymin-ytmaj; ymax=ymax+ytmaj
!	 	xmin=xmin-xtmaj; xmax=xmax+xtmaj
!	  	xmin=max(xmin,1d-3); 
	  	if(logy) ymin=max(ymin,1d-3);
	    if(logy) then
	    	 ytmaj=1; ytmin=1
	    	 endif
         if(cmd(1:4)/='LINE'.and.cmd(1:4)/='line') then
	 write(304,84) xmin,xmax,ymin,ymax,xtmaj,xtmin,ytmaj,ytmin
84	   format('@    world xmin ',f10.3,/
     x            '@    world xmax ',f10.3,/
     x            '@    world ymin ',f10.3,/
     x            '@    world ymax ',f10.3,:,/
     x            '@    xaxis  tick major ',f10.3,/
     x            '@    xaxis  tick minor ',f10.3,/
     x            '@    yaxis  tick major ',f10.3,/
     x            '@    yaxis  tick minor ',f10.3
     x              )
!	   write(6,84) xmin,xmax,ymin,ymax,xtmaj,xtmin,ytmaj,ytmin
	  	
	    if(logy) write(304,842) ig
842	    format('@G',I1,' type LOGY')
           endif ! not LINE
	  endif ! xmax>xmin .and. ymax>ymin
!	  endif	  
899	  continue   ! id
	  
1068	  continue   ! ILEN
	  close(304)
	  write(6,*) ' xmgr file written: ',trim(plot_file),
     x                    ' with ',ig+1,' graphs'
          if(xmgr) call system(' (xmgr '//trim(plot_file)//' &) &')
	  enddo  ! rank k

        else  if(cmd(1:3)=='MIN'.or.cmd(1:3)=='min') then  ! call minuit
          write(6,'(a)') '>'//trim(cmd)
	  interactive = .true.
!		call fcn(1,srch_error,fval,srch_value,1,0)
!		write(6,*) 'fcn =',fval
	   inquire(105,opened=op)
	   if(.not.op) open(105,file=trim(output_file)//'-trace')
	   inquire(106,opened=op)
	   if(.not.op) open(106,file=trim(output_file)//'-snap')
	  t = 1./ndof
	  write(6,1070) t
1070	  format(/'   Call MINUIT',/
     x            '   with NOGradient, STRategy=0, ERRordef=',f8.4)
	  call MNCOMD(fcn,'set nogradient',ICONDN,0)
!	  call MNCOMD(fcn,'set strat 1',ICONDN,0)
	  call MNEXCM(fcn,'set errordef',T,1,ICONDN,0)
	  call MNINTR(fcn,0)
C                   Find dof for chisq fit
           ndof = dof(datalen,datasets,srch_step,nvars)
         
!      write(6,10425) nparameters,nvarsv,ndatapoints,ndof

	  do ip=1,nvars
!			Get variable values back from Minuit
	   call MNPOUT(ip,srch_name(ip),srch_value(ip),srch_error(ip),
     x		srch_minvalue(ip),srch_maxvalue(ip),ivarbl)
	  enddo
	  noerror = .false.
        else
	  go to 980
	endif
	go to 100
980	  write(6,*) 'Unrecognised command : ',cmd
	go to 100
981	  write(6,*) 'Incomplete input command : ', cmd
	go to 100
982	  write(6,*) 'Parameter number',ip,' outside defined range 1:'
     x                   ,nvars
	go to 100
983	  write(6,*) 'Set number',id,' outside defined range'
	go to 100
999	stop
	end

!############################################################## FCN
	subroutine fcn(npar,grad,fval,xval,iflag,futil)
	use searchpar
	use searchdata
	implicit real*8(a-h,o-z)
	real*8 fval,grad(*),xval(*)
	external futil
	if(iflag==1) then
	  number_calls = 0
!					Initialise
	else if(iflag==2) then
!					Gradients -> grad
	endif
	  number_calls = number_calls+1
	  do ip=1,nvars
	   srch_value(ip) = xval(ip)
	  enddo
	  penalty = 0.
          if(hasEshift) call energylist  ! remake union list of energies

	  call fr

	  totchisq = sum(data_chisq(1:datasets)) + penalty
	  fval = totchisq/ndof
	  ip = min(6,datasets)
	  if(datasets==1) ip=0
	  write(105,2) number_calls,fval, 
     x           (data_chisq(id)/max(1,datalen(id)),id=1,ip)
1	 format(1p,5e12.5)
2	 format(i6,1p,e14.6,1x,0p,6f10.4)
3	 format(2i6,1p,e12.4,1x,1p,10e10.2)
	  ip = min(10,datasets)
	  write(106,3) number_calls,nvars,fval,
     x           (data_chisq(id)/max(1,datalen(id)),id=1,ip)
	  write(106,1) xval(1:nvars)
	  call flush(105)
	  call flush(106)
	
	if(iflag==3) then
!					Wrapup
	endif
	return
	end
	
!############################################################## energylist
	subroutine energylist ! (ndatasets,esmin,esmax)
	use searchdata
	use searchpar, only: nvars,srch_kind,srch_datanorm,srch_value
        use fresco1, only: peli=>pel
	implicit real*8(a-h,o-z)
	integer type,peld
        logical, EXTERNAL:: refer

        n = num_energies
	num_energies=0; energy_count(:) = 0; data_ien(:,:) = 0

	do id=1,datasets
	 type = data_type(id) 
           if(type>=5) cycle

	   ! write(6,*) 'For id',id,' nvars',nvars
	   dataEshift = 0.
          ip = data_shiftvar(id)
          if(ip>0) dataEshift=dataEshift + srch_value(ip)
!           do ip=1,nvars
!           if(srch_kind(ip)==6)  then
!           if(refer(ip,id)) dataEshift=dataEshift + srch_value(ip)
!            endif
!	   enddo

          neng=datalen(id)
          if(type==-1.or.type==0) neng=1
            peld = peli
            if(data_pel(id)>0) peld = data_pel(id)

          do ip=1,neng
            energy=data_energies(ip,id)+dataEshift
        ! write(6,*) id,ip,data_energies(ip,id),dataEshift,esmin,esmax
	  if(energy>0..and.energy>=esmin.and.energy<=esmax) then
	  ! write(6,*) ' Include ',energy
          ien=0
          do ie=1,num_energies
           if(abs(energy-energy_list(ie)-dataEshift)<1e-6 .and.
     x            peld==pel_list(ie)) then
             ien=ie;
!            write(6,*) 'id,type,dtype,energy,ce=',
!    x                id,type,data_type(id),energy,energy_count(ien)
             data_ien(id,ip) = ien
             if(type==-1.or.type==0) data_ien(id,1:datalen(id)) = ien
             if(data_type(id)<3) 
     x              energy_count(ien) = energy_count(ien)+1 ! count energies with multiple angles for eplots
             go to 14           ! found existing energy
           endif
          enddo
          ien = num_energies+1 ! list new energy
          num_energies = ien
          if(ien>maxen) then
            write(6,*) 'Should increase PARAMATER maxen!'
            stop
            endif
          energy_list(ien) = energy
! 	  write(6,*) ' New ',energy
          pel_list(ien) = peld
          energy_count(ien) = 1
          data_ien(id,ip) = ien
           if(type==-1.or.type==0) data_ien(id,1:datalen(id)) = ien
!         if(data_type(id)<3) energy_count(ien) = energy_count(ien) + 1 ! count energies with multiple angles for eplots
   14     continue
	  !else
	    !write(199,*) ' Exclude ',energy
	  endif  ! energy > 0
          enddo  ! ip
	enddo  ! id
	! write(6,*) ' Have ',num_energies
 	! write(6,20) energy_list(1:min(num_energies,20))
20	format('E',20f6.3)
          if(num_energies>n.and.allocated(theoryplot)) then
           deallocate(theoryplot)
           allocate (theoryplot(ntheory_pts,num_energies,datasets))
           endif

	end

        subroutine write_var1(ivar,i,name,kind,kp,pline,col,potential,
     x                        step)
        implicit real*8(a-h,o-z)
        character*15 name
        integer kind,kp,pline,col
        real*8 potential,step
        namelist /variable/ name,kind,kp,pline,col,potential,step,ivar
        write(i,nml=variable)
        end

        subroutine write_var12(ivar,i,name,kind,kp,pline,col,potential,
     x                         step,pline2,col2,ratio2)
        implicit real*8(a-h,o-z)
        character*15 name
        integer kind,kp,pline,col,pline2,col2
        real*8 potential,step,ratio2
        namelist /variable/ name,kind,kp,pline,col,potential,step,ivar,
     x                      pline2,col2,ratio2
        write(i,nml=variable)
        end

        subroutine write_var2(ivar,i,name,kind,nafrac,afrac,step)
        implicit real*8(a-h,o-z)
        character*15 name
        integer kind,nafrac
        real*8 afrac,step
        namelist /variable/ name,kind,nafrac,afrac,step,ivar
        write(i,nml=variable)
        end

        subroutine write_var3(ivar,i,name,kind,term,jtot,par,energy,
     x                        nopot,damp,step)
        implicit real*8(a-h,o-z)
        character*15 name
        integer kind,term,par
        real*8 jtot,energy,step
        logical nopot
        namelist /variable/ name,kind,term,jtot,par,energy,step,
     x           nopot,damp,ivar
        write(i,nml=variable)
        end

        subroutine write_var4(ivar,i,name,kind,term,icch,iach,lch,width,
     x                      sch,jch,channel,step,rwa)
        implicit real*8(a-h,o-z)
        character*15 name
        integer kind,term,channel
        real*8 jch,sch,width
        logical rwa
        namelist /variable/ name,kind,term,icch,iach,lch,width,rwa,
     x                      sch,jch,channel,step,ivar
        write(i,nml=variable)
        end

        subroutine write_var5(ivar,i,name,kind,dataset,reffile,
     x                        datanorm,step)
        implicit real*8(a-h,o-z)
        character*15 name
        character*80 reffile
        integer kind,dataset(2)
        real*8 datanorm,step
        namelist /variable/name,kind,dataset,reffile,datanorm,step,ivar
        write(i,nml=variable)
	end

        subroutine write_var6(ivar,i,name,kind,dataset,reffile,
     x                        dataEshift,step)
        implicit real*8(a-h,o-z)
        character*15 name
        character*80 reffile
        integer kind,dataset(2)
        real*8 dataEshift,step
        namelist /variable/ name,kind,dataset,reffile,
     x                       dataEshift,step,ivar
        write(i,nml=variable)
        end
        
        subroutine write_var7a(ivar,i,name,kind,term,damp,step,
     x                        valmin,valmax)
        implicit real*8(a-h,o-z)
        character*15 name
        integer kind,term
        namelist /variable/ name,kind,term,damp,step,valmin,valmax
        write(i,nml=variable)
        end
        subroutine write_var7b(ivar,i,name,kind,term,damp,step,
     x                        valmin,valmax,leff,energy,brune)
        implicit real*8(a-h,o-z)
        character*15 name
        integer kind,term
        real*8 leff
        logical brune
        namelist /variable/ name,kind,term,damp,step,valmin,valmax,
     x                      leff,energy,brune
        write(i,nml=variable)
        end



        subroutine write_varlink5(i,dataset,reffile,datanorm)
        implicit real*8(a-h,o-z)
        character*80 reffile
        real*8 datanorm
        integer dataset(2)
        namelist /datasetnorm/ reffile, datanorm , dataset
        write(i,nml=datasetnorm)
        end
        subroutine write_varlink6(i,dataset,reffile,dataEshift)
        implicit real*8(a-h,o-z)
        character*80 reffile
        integer dataset(2)
        real*8 dataEshift
        namelist /datasetshift/ reffile, dataEshift , dataset
        write(i,nml=datasetshift)
        end




	subroutine print_data(i,id,type,data_file,points,xmin,delta,lab,
     X   energy,idir,iscale,abserr,ic,ia,ib,data_kind,k,q,angle,jtot,
     x   par,channel,leg,pel,exl,labe,lin,lex,value,error,
     x   detinfo,info,eplots,logy,reffile,term)

	implicit real*8(a-h,o-z)
        integer id,type,points,idir,iscale,k,q,par,channel,
     x          pel,exl,labe,lin,lex,data_kind,term
	logical lab,abserr,detinfo,logy,eplots
	real*8 jtot
        character*14 info
	character*80 data_file,reffile

	if(type==0) then
        if(points>0.or.k>0.or.q>0.or.exl*lin*lex>1.or.labe/=pel) then
         call write_data0(id,i,type,data_file,points,lab,energy,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
	else
         call write_data0q(id,i,type,data_file,points,lab,energy,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
	endif

	else if(type==1) then
        if(points>0.or.k>0.or.q>0.or.exl*lin*lex>1.or.labe/=pel) then
         call write_data1(id,i,type,data_file,points,lab,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
	else
         call write_data1q(id,i,type,data_file,points,lab,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
	endif

	else if(type==2) then
        if(points>0.or.k>0.or.q>0.or.exl*lin*lex>1.or.labe/=pel) then
         call write_data2(id,i,type,data_file,points,lab,angle,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
	else
         call write_data2q(id,i,type,data_file,points,lab,angle,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
	endif

	else if(type==3) then
         call write_data3(id,i,type,data_file,points,lab,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)

	else if(type==4) then
         call write_data4(id,i,type,data_file,points,lab,jtot,par,
     x      channel,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)

	else if(type==5) then
         call write_data5(id,i,type,data_file,points,abserr)

	else if(type==6) then
         if(len(trim(reffile))>2) then
         call write_data6r(id,i,type,reffile,data_kind,value,error,
     x                     abserr)
	 else
         call write_data6p(id,i,type,par,value,error,abserr)
         endif

	else if(type==7.or.type==8) then
         call write_data78(id,i,type,term,value,error,abserr)

	else
	 write(6,*) "Writing data of type",type,"not yet implemented"
	endif
	return
	end

        subroutine write_data0(idat,i,type,data_file,points,lab,energy,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
        implicit real*8(a-h,o-z)
        integer type,points,ic,ia,l,q,pel,exl
        logical lab,abserr
        character*80 data_file
        namelist /data/ type,data_file,points,lab,energy,idat,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex
        write(i,nml=data)
        end

        subroutine write_data1(idat,i,type,data_file,points,lab,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
        implicit real*8(a-h,o-z)
        integer type,points,ic,ia,l,q,pel,exl
        logical lab,abserr
        character*80 data_file
        namelist /data/ type,data_file,points,lab,idat,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex
        write(i,nml=data)
        end

        subroutine write_data2(idat,i,type,data_file,points,lab,angle,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
        implicit real*8(a-h,o-z)
        integer type,points,ic,ia,l,q,pel,exl
        logical lab,abserr
        character*80 data_file
        namelist /data/ type,data_file,points,lab,angle,idat,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex
        write(i,nml=data)
        end

        subroutine write_data0q(idat,i,type,data_file,points,lab,energy,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
        implicit real*8(a-h,o-z)
        integer type,points,ic,ia,l,q,pel,exl
        logical lab,abserr
        character*80 data_file
        namelist /data/ type,data_file,lab,energy,idat,
     x      idir,iscale,abserr,ic,ia,pel
        write(i,nml=data)
        end

        subroutine write_data1q(idat,i,type,data_file,points,lab,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
        implicit real*8(a-h,o-z)
        integer type,points,ic,ia,l,q,pel,exl
        logical lab,abserr
        character*80 data_file
        namelist /data/ type,data_file,lab,idat,
     x      idir,iscale,abserr,ic,ia,pel
        write(i,nml=data)
        end

        subroutine write_data2q(idat,i,type,data_file,points,lab,angle,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
        implicit real*8(a-h,o-z)
        integer type,points,ic,ia,l,q,pel,exl
        logical lab,abserr
        character*80 data_file
        namelist /data/ type,data_file,lab,angle,idat,
     x      idir,iscale,abserr,ic,ia,pel
        write(i,nml=data)
        end

        subroutine write_data3(idat,i,type,data_file,points,lab,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
        implicit real*8(a-h,o-z)
        integer type,points,ic,ia,l,q,pel,exl
        logical lab,abserr
        character*80 data_file
C    	namelist /data/ type,data_file,points,xmin,delta,lab,energy,
C    X    idir,iscale,abserr,ic,ia,k,q,angle,jtot,par,channel,leg,
C    x    pel,exl,labe,lin,lex,value,error,detinfo,info
        namelist /data/ type,data_file,points,lab,idat,
     x      idir,iscale,abserr,ic,ia,k,q,pel,exl,labe,lin,lex
        write(i,nml=data)
        end

        subroutine write_data4(idat,i,type,data_file,points,lab,jtot,
     x      par,channel,abserr,ic,ia,k,q,pel,exl,labe,lin,lex)
        implicit real*8(a-h,o-z)
        integer type,points,ic,ia,l,q,par,pel,exl,channel
        logical lab,abserr
        character*80 data_file
        real*8 jtot
C    	namelist /data/ type,data_file,points,xmin,delta,lab,energy,
C    X    idir,iscale,abserr,ic,ia,k,q,angle,jtot,par,channel,leg,
C    x    pel,exl,labe,lin,lex,value,error,detinfo,info
        namelist /data/ type,data_file,points,lab,jtot,par,idat,
     x      channel,abserr,ic,ia,k,q,pel,exl,labe,lin,lex
        write(i,nml=data)
        end

        subroutine write_data5(idat,i,type,data_file,points,abserr)
        implicit real*8(a-h,o-z)
        integer i,type,points
	logical abserr
        character*80 data_file
        namelist /data/ type,data_file,points,abserr,idat
        write(i,nml=data)
        end

        subroutine write_data6r(idat,i,type,reffile,kind,
     x                          value,error,abserr)
        implicit real*8(a-h,o-z)
        integer type
        logical abserr
        character*80 reffile
        namelist /data/ type,reffile,value,error,abserr,idat,kind
        write(i,nml=data)
        end

        subroutine write_data6p(idat,i,type,par,value,error,abserr)
        implicit real*8(a-h,o-z)
        integer type,par
        logical abserr
        namelist /data/ type,par,value,error,abserr,idat
        write(i,nml=data)
        end

        subroutine write_data78(idat,i,type,term,value,error,abserr)
        implicit real*8(a-h,o-z)
        integer type,term
        logical abserr
        namelist /data/ type,term,value,error,abserr,idat
        write(i,nml=data)
        end
      subroutine shorten(COL) ! remove leading blanks and trailing 0 after . 
      character*80 COL,COL1
!   find first non-blank
      COL1 = COL
      do i=1,80
       if(COL(i:i) /= ' ') go to 20
      enddo
      i = 80
20    if(i>1) COL(1:80-(i-1)) = COL(i:80)
!         write(93,*) ' FNB =',i
!      write(93,'(a)') '<'//trim(COL1)//'> becomes <'//trim(COL)//'>'

! remove trailing 0 after . except after E
      if(index(COL,'.')/=0.and.index(COL,'E')==0) then
      l = lnblnk(COL)
      do i=l,2,-1
       if(COL(i:i)=='.') go to 30
       if(COL(i:i)/='0') go to 30
      enddo
30    l = i
      do i=l+1,80
       COL(i:i) = ' '
       enddo
      endif
      write(93,'(a)') '<'//trim(COL1)//'> becomes <'//trim(COL)//'>'
        return
          end

      integer function dof(datalen,ndatasets,srch_step,nvars)
c
C                   Find dof for chisq fit
c
        integer datalen(ndatasets),ndatasets
        real*8 srch_step(nvars)

        ndof = 0
        do id=1,ndatasets
          ndof = ndof + datalen(id)
        enddo
          ndatapoints = ndof
        nvarsv = 0
        do ip=1,nvars
          if(srch_step(ip)>1e-9) then
             nvarsv = nvarsv + 1
             ndof = ndof - 1
             endif
        enddo
         ndof = max(1,ndof) ! avoid divide by 0
      write(6,10425) nvars,nvarsv,ndatapoints,ndof
10425 format('  ',i4,' parameters (',i4,' variable) and',i7,
     x       ' data points, so dof=',i7)
        dof = ndof
        return
        end
