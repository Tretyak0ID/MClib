#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = ifort
OPTSC   =  -c -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict -module mod
OPTSL   =  -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -module mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules

#compiling rules
#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
all: $(addprefix $(DEXE),$(EXES))
