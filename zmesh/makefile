COMPILER = gcc

ifeq ($(COMPILER), intel)
	FC = ifort
	LD = link
	FFLAGS = /O3 /nologo /names:lowercase /libs:static /threads /c /object:$@
	LDFLAGS = /NOLOGO /WINAPP /OUT:$@
	ZMESH_LIBS =
	ZMESHWIN_LIBS =
else ifeq ($(COMPILER), gcc)
	FC = gfortran
	LD = gfortran
	FFLAGS = -s -c -O2 -Wall -Wextra -Wno-unused-dummy-argument -Wno-unused-parameter -fbounds-check
	LDFLAGS = -o $@
	ZMESH_LIBS =
	ZMESHWIN_LIBS =
endif

ZMESHWIN_SRC = snap.f90 heap.f90 const.f90 global.f90 geo.f90 win.f90 createmesh.f90 zmesh.f90
ZMESHWIN_OBJS = $(patsubst %.f90, %.o, $(ZMESHWIN_SRC)) zmesh.res
ZMESHWIN_DEPEND = zmeshwin.dep
ZMESH = zmeshwin.exe

ZMESH_SRC = snap.f90 global.f90 createmesh.f90 zmeshcon.f90
ZMESH_OBJS = $(patsubst %.f90, %.o, $(ZMESH_SRC))
ZMESH_DEPEND = zmesh.dep
ZMESH = zmesh.exe

EXTRA_TO_DELETE = *.mod

.PHONY: all
all : $(ZMESH) $(ZMESHWIN)

$(ZMESH) : $(ZMESH_OBJS)
	$(LD) $(LDFLAGS) $(ZMESH_OBJS) $(ZMESH_LIBS)

$(ZMESHWIN) : $(ZMESHWIN_OBJS)
	$(LD) $(LDFLAGS) $(ZMESHWIN_OBJS) $(ZMESHWIN_LIBS)

.PHONY : clean
clean :
	del $(ZMESH_OBJS) $(ZMESHWIN_OBJS) $(ZMESH) $(ZMESHWIN) $(EXTRA_TO_DELETE)

.PHONY : depend
depend :
	gfortran -cpp -M $(ZMESH_SRC) > $(ZMESH_DEPEND)
	gfortran -cpp -M $(ZMESHWIN_SRC) > $(ZMESHWIN_DEPEND)

%.o %.mod : %.f90
	@del $*.mod
	$(FC) $(FFLAGS) $<

# correct for ifort, needs to be generalised for gfortran
zmesh.res : zmesh.rc
	rc zmesh.rc

-include $(ZMESH_DEPEND)
-include $(ZMESHWIN_DEPEND)
