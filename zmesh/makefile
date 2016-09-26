%.obj : %.f90
	ifort /nologo /c $(flags) $<

zmeshwinobj = snap.obj heap.obj const.obj global.obj geo.obj win.obj createmesh.obj zmesh.obj zmesh.res

zmeshobj = snap.obj global.obj createmesh.obj zmeshcon.obj


zmeshwinflags = /nologo /optimize:3 /winapp /warn
zmeshflags = /nologo /optimize:3 /warn

.PHONY: all
all : zmeshwin.exe zmesh.exe

zmeshwin.exe : $(zmeshwinobj)
	ifort $(zmeshwinflags) /exe:zmeshwin.exe $(zmeshwinobj)

zmesh.exe : $(zmeshobj)
	ifort $(zmeshflags) /exe:zmesh.exe $(zmeshobj)

snap.obj :

heap.obj :

const.obj :

global.obj :

geo.obj : global.obj const.obj

win.obj : global.obj createmesh.obj const.obj geo.obj heap.obj

createmesh.obj : snap.obj global.obj

zmesh.obj : global.obj const.obj win.obj

zmesh.res : zmesh.rc
	  rc zmesh.rc
