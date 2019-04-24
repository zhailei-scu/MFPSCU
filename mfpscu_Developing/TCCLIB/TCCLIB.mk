#compiler
#sorce directories

ifeq ($(origin Ccomp), undefined) 
  Ccomp   := g++
  ifeq ($(operateSystem), CYGWIN)
    Ccomp   := cl
  endif
  
  ifeq ($(operateSystem), LINUX)
    Ccomp   := g++
  endif
endif

ifeq ($(origin ConfigName), undefined) 
  ConfigName := Release
endif

CP := -c
OutF := -o
IncF := -I
ifeq ($(origin Coflags), undefined)
  ifeq ($(operateSystem), CYGWIN)
      ifeq ($(ConfigName), Release)
        Coflags := /Wall /O2
      endif

      ifeq ($(ConfigName), Debug)
        Coflags := /Wall /Zi
      endif

      CP := /c
      OutF := /Fo:
      IncF := /I	
  endif

  ifeq ($(operateSystem), LINUX)
      ifeq ($(ConfigName), Release)
        Coflags := -Wall -Mconcur -O2 -Minform=warn -Minfo=all
      endif

      ifeq ($(ConfigName), Debug)
        Coflags := -Wall -g -O0 -Minform=warn -Minfo=all
      endif

      CP := -c
      OutF := -o
      IncF := -I
  endif

endif

###########################################################
#sorce dir name
objname := TCCLIB

ifeq ($(origin TCCINTERDIRS), undefined) 
  TCCINTERDIRS := $(mcpscusor)$(Segment)TCCLIB$(Segment)sor
endif
sor  := $(TCCINTERDIRS)

#include dir
ifeq ($(origin LIBDIRD), undefined)
LIBDIRD := $(mcworkspace)$(Segment)LIB$(Segment)$(ConfigName)
endif
incdir := $(LIBDIRD)

#target directories
tgt  := $(LIBDIRD)

#target lib name
libname  := lib_$(objname).$(LIB_EXT)

#the Macro

ifeq ($(operateSystem), CYGWIN)
  Macro += /DLIBTCCPATH=\"$(tccpath)$(Segment)win32\" /DLIBTCCHEADERF=\"$(tccpath)$(Segment)win32$(Segment)libtcc$(Segment)libtcc.h\" /DDIFFUSORDEFPATH=\"$(TCCINTERDIRS)\" /D$(operateSystem)
endif

ifeq ($(operateSystem), LINUX)
  Macro += -DLIBTCCPATH=\"$(tccpath)\" -DLIBTCCHEADERF=\"$(tccpath)/libtcc.h\" -DDIFFUSORDEFPATH=\"$(TCCINTERDIRS)\" -D$(operateSystem)
endif

#######################################################          
nlist    :=  DiffusorList	\
             ReactionsList	\
	     TccInterp
             
objects  := $(foreach n, $(nlist), $(tgt)$(Segment)$(n).o)
cfiles   := $(foreach n, $(nlist), $(sor)$(Segment)$(n).cpp)
#######################################################
$(libname) : $(objects)  
	ar -rcs $(libname) $(objects)
	mv $(libname) $(tgt)

$(tgt)$(Segment)DiffusorList.o : $(sor)$(Segment)DiffusorList.cpp
	$(Ccomp) $(CP) $< $(Coflags) $(IncF)$(incdir) $(OutF) $@

$(tgt)$(Segment)ReactionsList.o : $(sor)$(Segment)ReactionsList.cpp
	$(Ccomp) $(CP) $< $(Coflags) $(IncF)$(incdir) $(OutF) $@

$(tgt)$(Segment)TccInterp.o : $(sor)$(Segment)TccInterp.cpp
	$(Ccomp) $(CP) $< $(Coflags) $(Macro) $(IncF)$(incdir) $(OutF) $@

######################################################################
clean:
	-rm $(objects) $(libname) 
