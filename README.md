# Frescox
Scattering code Frescox for coupled-channels calculations

FRESCOX   FRXY version 7.2.2 at https://github.com/LLNL/Frescox
LLNL-CODE-811517


This directory contains four sub-directories: source, man, test and util.

The source/ directory contains : Fortran files *.f,
                                 fx*.def files for separate machines

  File nagstub.f in the case you do not have the nag library locally.

The test/   directory contains : at least 6 test jobs xeta, lane20 & f19xfr,
                                             e80f49b, on2 & be11
                                 their various outputs  SUN/*.out 
               (The input files were originally CRAY UNICOS jobs,
                hence the comments at the beginning.)

The man/   directory contains the instruction manual in latex:
                      frescox-input-manual.tex: latex source 
                      frescox-input-manual.pdf: printable output
More documentation is at http://www.fresco.org.uk/documentation.htm

To compile FRESCOX,
 
   1) Enter frxy/source, and then edit the makefile for your target machine,
	by setting the MACH variable as appropriate (either in the makefile
		or in our local shell setup),
	by choosing the file nagstub.f if the NAG library not available
	
	The script 'mk' attempts to guess the correct MACH settings
	for ordinary frescox version AND compile in a corresponding subdirectory.

   2) Edit aliases there,
      to set FRESCOXLIB to point to directory for storing the binary

   3) Copy your aliases to ~/.fresco.aliases
      Edit FRESCOXLIB according to 2) above
      Execute .fresco.aliases e.g. in .cshrc  by including:	
        source ~/.fresco.aliases

If you are to install frescox yourself in a standard bin directory, 
       then steps 3 and 4 may be omitted, and step 2' performed manually.

 If your compiler is gfortran, for example, then:
   2') Compile the subroutines required by:
        
        make MACH=gfortran

   4) Install, to copy `frescox' to the FRESCOXLIB.
        make MACH=gfortran install

   5) Clean up, with:	
        make clean
  
To run FRESCOX,

   1) Enter test/ directory.

   2) The scripts include commands to construct temporary 'data' files. 
	These scripts are run by just saying  e.g.
           frescox < lane20.nin > lane20.out
        See file 'do-all.bat' to run all the test cases

   3) In the test/legacy directory there are executable 'job' run scripts
      To save the output in a file .e.g. `out', run the scripts by
       lane20.job > out &
         or simply (using the 'run' command in the 'aliases' file):;
       run lane20.job 
         to use input file lane20.job and produce output file lane20.out.

   4) To save any other output files from frescox, e.g. fort.16 for 
      cross sections, 
         touch lane20.xsecs
         ln -s lane20.xsecs fort.16
      before running frescox.
      The file fort.16 may have to be called for016.dat on some machines.

Please let me know if you have any questions or problems:
    
   I.Thompson@fresco.org.uk

Cheers, Ian Thompson
November 2022

