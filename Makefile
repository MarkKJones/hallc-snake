# By Jixie Zhang
# Makefile for snake program
# must be compiled with -m32

#####################################################
PROG   := snake
SRCDIR := src
INCDIR := src
OBJDIR := obj
#####################################################
FF = gfortran
LD = gfortran
FFLAGS = -O2 -ff2c -g -fno-automatic 

################################ Libraries
ifndef SYSLIBS
 SYSLIBS  = -lnsl -lX11
endif

#for SunOS, cernlib 2001 need ibsocket.a
ifeq ($(MYOS),SunOS)
 SYSLIBS  += -lsocket
endif
  
#this program used 2001 cernlib /apps/cernlib/sun4_solaris8/2001 in jlabs2,
CERN_ROOT := /apps/cernlib/x86_64_rhel7/2005

CERNLIBS := -L$(CERN_ROOT)/lib -lpdflib804 -lmathlib -lphtools \
 -lgeant321 -lpawlib -lgraflib -lgrafX11 -lpacklib -lkernlib 


LIBS = $(SYSLIBS) $(CERNLIBS)
###################################
ifndef PROG
  PROG = $(shell grep -i program $(SRCDIR)/*.[fF] | awk '$$2=="program"&&$$4==""{print $$3}' | sed s/\ //g)
endif
ifeq ($(PROG),) 
  PROG = $(shell grep -i program $(SRCDIR)/*.[fF] | awk '$$1=="program"&&$$3==""{print $$2}' | sed s/\ //g)  
endif


PROG := $(shell echo ${PROG} | sed s/\ //g)
EXE   = $(PROG)


SOURCES  = $(wildcard $(SRCDIR)/*.[fF])
#OBJS     = $(shell echo $(SOURCES) | sed s/f\ /o\ /g)
OBJS    := $(addsuffix .o, $(basename $(SOURCES)))
OBJS    := $(patsubst $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))

###################################
$(OBJDIR)/%.o:$(SRCDIR)/%.f
	$(FF) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o:$(SRCDIR)/%.F
	$(FF) $(FFLAGS) -c $< -o $@

%.o:%.f
	$(FF) $(FFLAGS) -c $< -o $@

%.o:%.F
	$(FF) $(FFLAGS) -c $< -o $@

all:  $(OBJDIR) $(EXE)

$(OBJDIR):
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

$(EXE): $(OBJDIR) $(OBJS) 
	$(FF) -o $(EXE) $(OBJS) $(LIBS) $(FFLAGS)
	@echo "Linking $(EXE) ... done!"

#this part is to support CLAS_PACK make mechanism
exe: lib
	$(FF) -o $(EXE) $(OBJS) -L$(OBJDIR) -l$(PROG) $(LIBS) $(FFLAGS)
	@echo "Linking $(EXE) ... done!"

#this part is to support CLAS_PACK make mechanism
lib: $(OBJDIR) $(OBJS)
	@ar -rf $(OBJDIR)/lib${PROG}.a ${OBJS}

clean:
	rm -f $(OBJS) $(OBJDIR)/lib${PROG}.a core

delete:
	rm -f $(EXE)

help:
	@echo PROG = ${PROG}
	@echo SOURCES = ${SOURCES}
	@echo OBJS = ${OBJS}
	@echo EXE = ${EXE}
	@echo LIBS = ${LIBS}
	@echo FFLAGS = ${FFLAGS}
