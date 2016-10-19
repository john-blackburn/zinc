## Zinc, a multiphysics finite element program

Zinc is a finite element program capable of solving a range of physics problems which can be expressed in the form of partial differential equations (PDEs) of second order. Zinc can solve linear or non-linear problems which are mono- or multiphysics. Furthermore Zinc can solve static (including steady state) problems or transient problems.

Zinc is a mathemetical program meaning it has no built-in knowledge of physics of any problem. Rather the user needs to set the form of the PDEs to be solved by setting various matrices. So to solve a problem in Zinc, you need to write out the equations you want to solve and then plug them into Zinc by refering to the user manual. You will also need to consider initial and boundary conditions and the manual also describes how to enter these details and many others. (note full user manuals are provided with Zinc releases as PDF files; they are also hosted here on Github as LaTeX files)

So far, Zinc has been used to solve: electrostatic, elastic, piezoelectric, multiferroic and fuel cell problems (the last involving diffusion and electrochemical equations such as Maxwell-Stephan and Nernst-Planck). Since Zinc can solve almost any 2nd order PDE system, a near infinity of other systems could be solved also.

Zinc consists of a front end mesher/visualisation program (Zmesh), a core solver (Zinc) and a back-end program for visualising the simulation result (Zpp: the Zinc Post Processor). 

To find out more about Zinc and download the PDF manuals, see the website

http://interactive.npl.co.uk/zinc

## Supported operating systems

Zinc is designed to run on Windows but can run on Linux using WINE. The Zinc core and post processor programs (zinc.exe, zpp.exe) should be easy to recompile on Linux as they are standard Fortran. The only thing that would need to be changed is the reference to win32 APIs LoadLibrary, FreeLibrary and GetProcAddress in matrices.f90. (these load a user-provided DLL file for non-linear problems).

The mesher/visualiser program zmesh uses Win32 APIs extensively (including Windows GDI) so needs WINE to run. In future it would be nice to convert this to OpenGL so it could be compiled anywhere but that is way off in the future.

## Compilers and prerequisites

Zinc uses the UMFPACK library part of SuiteSparse:

http://faculty.cse.tamu.edu/davis/suitesparse.html

and also the SLATEC library:

http://www.netlib.org/slatec/

In makefile you will see references to these libraries which must be present. Absolute paths are given so you will likely need to change these paths to whereever you have installed the static libs needed.

You will need to compile both the above libraries yourself. Note that UMFPACK depends on BLAS. I used OpenBlas for this purpose.

http://www.openblas.net/

Zinc can be compiled using gcc (on Windows we use mingw from http://mingw.org). In particular the programs gcc and gfortran are used to build Zinc (see makefiles). However, the zmesh part must (currently) be compiled using Intel Fortran (ifort) as you can see from zmesh/makefile. This is because Intel Fortran contains bindings from Fortran to the Win32 API. In future I hope to add my own bindings so the whole project can be compiled with gcc.

## Compiling Zinc

In zinc main directory run

<pre>
make depend
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

Now convert both postscript files to PDF using Ghostview (http://pages.cs.wisc.edu/~ghost/).

## Create an installer package

Use NSIS (https://sourceforge.net/projects/nsis/) to run zinc.nsi. This will package all parts of Zinc (executables, manuals, examples) into a self-extracting exe file.

## Example simulations

These are present in the examples/ directory. See also fullfuel/ for an advanced fuel cell problem.

## Licence

Zinc is provided as free software under the Gnu Public Licence version 3 (GPL v3). See COPYING.txt for details. In particular Zinc comes with NO WARRANTY as stated in the GPL (see COPYING.txt):

<pre>
  15. Disclaimer of Warranty.

  THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

  16. Limitation of Liability.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.
</pre>
