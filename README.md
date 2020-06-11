# Frescox
Scattering code Frescox for coupled-channels calculations

FRESCO   FRXY version 7 at https://github.com/LLNL/Frescox
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

To compile FRESCO,
 
   1) Enter frxy/source, and then edit the makefile for your target machine,
	by setting the MACH variable as appropriate (either in the makefile
		or in our local shell setup),
	by choosing the file nagstub.f if the NAG library not available
	
	The scripts 'mk' and 'mkp' attempt to guess the correct MACH settings
	for ordinary fresco and pvm3 frescop versions, respectively,
	AND compile in a corresponding subdirectory.

   2) Edit aliases there,
      to set FRESCOXLIB to point to directory for storing the binary

   3) Copy your aliases to ~/.fresco.aliases
      Edit FRESCOXLIB according to 2) above
      Execute .fresco.aliases e.g. in .cshrc  by including:	
        source ~/.fresco.aliases

   If you are to install fresco yourself in a standard bin directory, 
       then steps 3 and 4 may be omitted, and step 2' performed manually

   2') Compile the subroutines required by:
        make

   4) Install, to copy `fresco' to the FRESCOXLIB.
        make install

   5) Clean up, with:	
        make clean
  
To run FRESCO,

   1) Enter test/ directory.

   2) The scripts include commands to construct temporary 'data' files. 
	These scripts are run by just saying  e.g.
       lane20.job

   3) To save the output in a file .e.g. `out', run the scripts by
       lane20.job > out &
         or simply
       run lane20.job 
         to use input file lane20.job and produce output file lane20.out.

   4) If you have separate `data' or `in' files, the command is
       fresco < lane20.in > lane20.out

   5) To save any other output files from fresco, e.g. fort.16 for 
      cross sections, 
         touch lane20.xsecs
         ln -s lane20.xsecs fort.16
      before running fresco.
      The file fort.16 may have to be called for016.dat on some machines.

Please let me know if you have any questions or problems:
    
   I.Thompson@fresco.org.uk

Cheers, Ian Thompson
November 2018

