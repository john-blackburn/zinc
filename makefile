###############################################################################
# Makefile for Zinc
# This makefile will do some unnecessary compilation when a module foo.f90
# is changed but its interface stays the same (almost always the case)
# In that case foo.o and foo.mod are recreated and any file that USEs foo
# will be recompiled. In practice when foo.f90 is changed any file above it
# is also recompiled. Thus changing zinc.f90 will recompile only this file
# but changing qpr.f90 forces recompile of most of the project.
# This makefile is more efficient than a batch file but only just!
# (this is a problem with make which cannot deal with auto-generated header files)

# Note that a static library in slatec_zinc is needed for either compiler.
# use "make" or "make gcc" to compile with gfortran
# use "make intel" to compile with intel compiler.
# If changing compilers or changing optimisation flags, do a "make clean" first

# Should run "make depend" every so often as well

###############################################################################
# Options

SHELL := cmd.exe

 # default compiler to use, can be changed on command line
COMPILER = gcc

ifeq ($(COMPILER), intel)
	FC = ifort
	LD = link
	FFLAGS = /O3 /nologo /names:lowercase /libs:static /threads /c /object:$@
	LDFLAGS = /NOLOGO /OUT:$@
	ZINC_LIBS = slatec_zinc/slatec_zinc.lib
	ZPP_LIBS =
else ifeq ($(COMPILER), gcc)
	FC = gfortran
	LD = gfortran
	FFLAGS = -s -c -cpp -O2 -Wall -Wextra -Wno-unused-dummy-argument -Wno-unused-parameter -I d:\suitesparse\umfpack\demo -fno-underscoring -fbounds-check
	LDFLAGS = -o $@
	ZINC_LIBS = d:\slatec\libslatec.a d:\suitesparse\umfpack\lib\libumfpack.a d:\suitesparse\AMD\lib\libamd.a d:\suitesparse\suitesparse_config\libsuitesparseconfig.a d:\openblas\lib\libopenblas.a
	ZPP_LIBS =
endif

###############################################################################
# Files

ZINC_SRC = util.f90 precision.f90 strings.f90 evaluate.f90 common.f90 indexq.f90 solvers.f90 \
    iofile.f90 shape.f90 geom.f90 heap.f90 matrices.f90 interpolate.f90 tree.f90 qpr.f90 cuboid.f90 update.f90 jacobian.f90 zinc.f90
ZINC_OBJS = $(patsubst %.f90, %.o, $(ZINC_SRC))
ZINC_DEPEND = zinc.dep
ZINC = zinc.exe

ZPP_SRC = util.f90 precision.f90 strings.f90 evaluate.f90 common.f90 indexq.f90 iofile.f90 shape.f90 \
          geom.f90 fcnb.f90 matrices.f90 integrate.f90 scan.f90 zpp.f90
ZPP_SRC_FIXED = dnls1edep.f
ZPP_OBJS = $(patsubst %.f90, %.o, $(ZPP_SRC)) $(patsubst %.f, %.o, $(ZPP_SRC_FIXED))
ZPP_DEPEND = zpp.dep
ZPP = zpp.exe

EXTRA_TO_DELETE = *.mod nl.exp nl.lib

###############################################################################
# Targets

.PHONY: all
all : $(ZINC) $(ZPP)

$(ZINC) : $(ZINC_OBJS)
	$(LD) $(LDFLAGS) $(ZINC_OBJS) $(ZINC_LIBS)

$(ZPP) : $(ZPP_OBJS)
	$(LD) $(LDFLAGS) $(ZPP_OBJS) $(ZPP_LIBS)

.PHONY : clean
clean :
	del $(ZINC_OBJS) $(ZPP_OBJS) $(ZINC) $(ZPP) $(EXTRA_TO_DELETE)

.PHONY : depend
depend :
	gfortran -cpp -I d:\suitesparse\umfpack\demo -M $(ZINC_SRC) > $(ZINC_DEPEND)
	gfortran -cpp -M $(ZPP_SRC) $(ZPP_SRC_FIXED) > $(ZPP_DEPEND)

###############################################################################
# Rules. The delete command ensures the .mod file always has the same timestamp
# as the .o file

%.o %.mod : %.f90
	@del $*.mod
	$(FC) $(FFLAGS) $<

%.o : %.f
	$(FC) $(FFLAGS) $<

###############################################################################
# Dependencies

-include $(ZINC_DEPEND)
-include $(ZPP_DEPEND)

dielectric.o: dielectric.f90
