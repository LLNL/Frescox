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
	program fresco
	use parameters
	use factorials
	use drier
	use io
	use fresco1, only: lmax,jbord,jump,rterms
	use searchpar
	use searchdata
	use parallel
	use trace
#ifdef MPI
        use mpi
#endif /* MPI */
    	character*70 TMP,TMPN
      CHARACTER CHME*3,F3*10
     	common/cfs/ TMP,TMPN  

!				Change stdout recl on some machines
	call stdout(6,140)

	ki = 5
	ko = 3
	koe = 6
*  MACH =
*  1) standard serial F90
*  2) standard serial F90 with BLAS
*  3) MPI on SGI Altix ( at MSU ) which could only read in from the master thread so had to pass all input variables to the other threads
*  4) MPI on all other architectures
*  		(note 3 and 4 both use BLAS too)
*  STDINALL if all threads can read from stdin

	call machine(mach)
        MPIC=.false.
        STDINALL=.false.
        if(MACH==3 .or. MACH==4) MPIC=.true.
        if(MACH==4) STDINALL=.true.
	iame = 0
	mpinodes = 1
	
**MPI*******************************************************************
#ifdef MPI
      if(MPIC) then
       call MPI_init(ierr)
       call MPI_comm_size(MPI_COMM_WORLD,mpinodes,ierr)
       call MPI_comm_rank(MPI_COMM_WORLD,iame,ierr)
!	 write(500+iame,*) ' I am ',iame,' of ',mpinodes
      endif ! MPI
#endif /* MPI */
************************************************************************
c	open(ko,form='formatted',delim='apostrophe',file='fort.3')
	open(ko,form='formatted',delim='apostrophe')

 	datasets = 0
	nvars = 0
	num_energies = 0
	number_calls = -1 	! no searching
	rterms = .false.
	final = .true.
	say = .true.
	jbord = 0; jump=0
	written(:) = .false.
	written(3) = .true.
!			Pre-read input to find SOME parameter limits
        call freadf(ki,ko,koe,TMP,lmax,jbord,jump,.true.)
      
        maxmul = lmax1 + 10
	acc8 = epsilon(acc8)
	fpmax = huge(acc8)**0.8d0
!	fpmax = 1d290
	ki = 3
	ko = 6
	koi= 6
**MPI*******************************************************************
      if(MPIC) then
        if(iame/=0) ko = 51
      endif ! MPI
************************************************************************
	written(ko) = .true.
	rewind ki
       	pi = 4d0 * atan(1d0)
       	r4pi = sqrt(0.25/pi)
	gettheoryplot = .true.
        MAXCH=0 ; MAXICH=0 			! no arrays allocated
C			DO IT!
	call fr
C			DONE IT!
**MPI*******************************************************************
#ifdef MPI
      if(MPIC)then
        call MPI_finalize(ierr)
      endif
	if(numthread>1) write(500+iame,*) ' Stop'
#endif /* MPI */
************************************************************************
 	stop
C
C    FILE ALLOCATIONS
C    ----------------
C
C    File   Fixed/  RECL   Sequential  by     for
C     No.  Variable bytes  /Random
C
C     3      F      80        S        FREAD  local copy of user input
C     4      F      80        S        INTER  external KIND=1,2,9 form factors
C                                      POTENT external potentials
C     5      F      80        S               standard input
C     6      V     133        S               standard output
C     7      F      80        S        DISPX  elastic S-matrix elements
C     8  #   F    N*LC        R        FR,INTER  s/p wfs, channel wfs
C     9  #   F NLO*(NLN-1)*LC R        Q/KERNEL   complex transfer multip
C    10  @   V                S        FR,CRISS  S-matrix elements (cfs)
C    11  #   F NLO*(NLN-1)*LR R        Q/KERNEL   real transfer multipole
C    12  #   V NLL*NLO*LC     S        KERNEL/SOURCEtransfer kernels
C    13      F                S          FR      total cross sections/state
C    14  #   V                S        INTER/CPAIR interaction potentials
C    15      F NLL*NLO*LC     S          FR      local equivalent potent
C    16      F      80        S          CRISS   tables of cross section
C    17      F      80        S          FR      output scattering wave
C    18  #   V                S          FR      wfns of 'best' iterate
C    19      F      80        S        CRISS     Cross sections for plotting
C    20-33				Available for user (eg bound states)
C    34      F      80        S        POTENT    output potentials
C    35      F      80        S        CRISS  input scattering amplitude
C    36                                CRISS  output scattering AMPL amplitudes
C    37                                CRISS  output scattering FAM amplitudes
C    38      F      80        S        DISPX     cross sections for each J/pi
C    39      F      80        S         FR       2 cross sections for each Elab
C    40      F      80        S         FR       all cross sectns. for each Elab
C    41      F      80        S        SOURCE    source terms at each iteration
C    42      F      80        S        SOURCE    bin wavefunctions for each E
C    43      F      80        S        INFORM    bin phase shifts as k functions
C    44      F      80        S        INFORM    bin phase shifts as E functions
C    45      F      80        S        ERWIN     scat phase shift as E functions
C    46      F      80        S        INFORM    bs wave functions & Whit ratios
C    47      F      80        S        MULTIP    reduced matrix elements
C    48  #   V     133        S          FR   concurrency log file
C    49  #   V                S        FR,CRISS  S-m elements concurrent->10
C    50  #   V     133        S           tmp standard output concurrent->6
C    51  #   V     133        S      unwanted standard output for node N>0
C    52  #   V     133        S      unwanted standard output for node N>0
C    53      V                S      spec file for sturmx(x)
C    54      V     bin        S      mel file for sturmx(x)
C    55      F      80        S        INFORM    single-particle wave funtions
C    60      F      80        S        FREAD  local copy of user input
C    71      F      80        S        Elastic phase excitation functions
C    66      V NLL*NLO*LC     S        INTER  KIND=9 nonlocal formfactor
C    99      F      80        S        CRISS     input J/L rescaling factors
C
C   Here, LR == size of REAL*8,   LC == size of COMPLEX*16
C         #  == distributed file: distinct file name on each node
C         @  == shared concurrent file: single file, with writes from each node
C
C NOT USED ANY MORE:
C     1  #   F    NLO*2*8     R        QNN    two-particle multipoles
	end
