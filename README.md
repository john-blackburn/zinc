## Zinc, a multiphysics finite element program

Zinc is a finite element program capable of solving a range of physics problems which can be expressed in the form of partial differential equations (PDEs) of second order. Zinc can solve linear or non-linear problems which are mono- or multiphysics. Furthermore Zinc can solve static (including steady state) problems or transient problems.

Zinc is a mathemetical program meaning it has no built-in knowledge of physics of any problem. Rather the user needs to set the form of the PDEs to be solved by setting various matrices. So to solve a problem in Zinc, you need to write out the equations you want to solve and then plug them into Zinc by reading the user manual.

So far, Zinc has been used to solve: electrostatic, elastic, piezoelectric, multiferroic and fuel cell problems (the last involving diffusion and electrochemical behaviour). Since Zinc can solve almost any 2nd order PDE system, a near infinity of other systems could be solved also.

Zinc consists of a front end mesher/visualisation program (Zmesh), a core solver (Zinc) and a back-end program for visualising the simulation result (Zpp: the Zinc Post Processor). 

To find out more about Zinc and download the manuals, see the website

http://interactive.npl.co.uk/zinc

## Supported operating systems

Zinc is designed to run on Windows but can run on Linux using WINE. The Zinc core and post processor programs (zinc.exe, zpp.exe) should be easy to recompile on Linux as they are standard Fortran. The only thing that would need to be changed is the reference to win32 API LoadLibrary, FreeLibrary and GetProcAddress in matrices.f90. (this loads a DLL file for non-linear problems).

The mesher/visualiser program zmesh uses win32 APIs extensively (including Windows GDI) so needs WINE to run.

## Compilers and prerequisites

Zinc uses the UMFPACK library part of SuiteSparse:

http://faculty.cse.tamu.edu/davis/suitesparse.html

and also the SLATEC library:

http://www.netlib.org/slatec/

In makefile you will see references to these libraries which must be present. Absolute paths are given so you will likely need to change these paths to whereever you have installed the static libs needed.

You will need to compile both the above libraries yourself. Note that UMFPACK depends on BLAS. I used OpenBlas for this purpose.

http://www.openblas.net/

Zinc can be compiled using gcc (on Windows this is called mingw). In particular the programs gcc and gfortran are used (see makefiles). However, the zmesh part must (currently) be compiled using Intel Fortran (ifort) and you can see from zmesh/makefile. This is because Intel Fortran contains bindings from Fortran to the Win32 API. In future I hope to add my own bindings so the whole project can be compiled with gcc.

## Compiling Zinc

In zinc main directory run

<pre>
make
</pre>

This will build zinc.exe, zpp.exe (note, if you install mingw on Windows make is actually provided under the name mingw32-make)

run

<pre>
zfec.bat
zppwinc.bat
</pre>

to build zincwin.exe and zppwin.exe (the GUI equivalents of the above two programs)

To build zmesh.exe and zmeshwin.exe go to /zmesh and run
<pre>
make
</pre>

## Compile the documentation

<pre>
latex zincman.tex
dvips -o zincman.ps zincman.dvi
latex zinctheo.tex
dvips -o zinctheo.ps zinctheo.dvi
</pre>

Now convert both postscript files to PDF using Ghostview.

## Create an installer package

Use NSIS to run zinc.nsi. This will package all parts of Zinc (executables, manuals, examples) into a self-extracting exe file.

## Example simulations

These are present in the examples/ directory. See also fullfuel/ for an advanced fuel cell problem.
