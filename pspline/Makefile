#GNU Makefile to build pspline Library libpspline.a 
#        and test drivers pspltest and pspltes8
#
#02/09/99 C. Ludescher-Furth

# 
#=================================================================
# Instructions to build shareable image libpshare.so
#=================================================================
#
# I. build shareable image
#
#   If NTCC module:
#
#      1. make so
#         libpspline.so will be in ../<MACHINE>/lib/libpspline.so
#
#   If TRANSP distribution:
#
#      1. cd $CODESYSDIR/pspline
#      2. make -k -f ../pspltest/makefile.inf OBJ=$LOCAL so
#         libpspline.so will be in $LOCAL/lib/libpspline.so
#      3. rm $LOCAL/obj/pspline/*.o
#
# II. Install library
#
#     copy libpspline.so into your "site specific" library directory
#     the library location must be in the LD_LIBRARY_PATH of 
#         - the MDSplus Server
#         - the users 
#     At PPPL the library is in /p/beast/<machine>/lib
#
#-----------------------------------------------------------------

#

ifneq ($(MAKELEVEL),0)
# if OBJ was defined for main Makefile one level down
ifeq ("${OBJ}",".")
	override OBJ=..
endif
endif

# define system dependent flags, etc.
-include ../share/Make.local

SHR = ../share/Make.flags
ifneq ("$(wildcard $(CODESYSDIR)/source/misc/makeflags.mk)","")
        SHR = $(CODESYSDIR)/source/misc/makeflags.mk
endif
include $(SHR)

ifeq ($(MACHINE),HP)
	FTEST := $(FFLAGS) $(STATIC)
else
	FTEST := $(FFLAGS)
endif

ifdef DEBUG
	FTEST := $(DFFLAGS)
	FFLAGS = $(DFFLAGS)
endif

INCFLAGS := -I./ ${INCFLAGS}

ARC = $(OBJ)/lib/libpspline.a
OBJDIR = $(OBJ)/obj/pspline

# don't rebuild library if using installed one in $PREFIX 
ifeq ($(MAKELEVEL),0)
	THISLIB=$(ARC)
endif

srcdir = $(shell pwd)
DATA   = $(shell $(LS) *.ref *.m *output)

CZ_H = czspline_handle_size.h czspline_capi.h
#=============================================================
# f90 modules
 
Mnams = ezspline.mod ezspline_obj.mod ezspline_type.mod pspline_calls.mod \
	czspline_pointer_types.mod
Mobjs = ezspline_mod.o pspline_calls.o czspline_pointer_types.o
MODS = $(foreach file,$(Mobjs),$(ARC)($(file)))

ifeq ($(MODUP),Y)
 MODS0=$(foreach m,$(Mnams),$(shell  echo $(m) | tr 'a-z' 'A-Z'))
 MODULES=$(foreach m,$(MODS0),$(subst .MOD,.$(MODEXT),$(m)))
else
 MODULES = $(foreach m,$(Mnams),$(subst .mod,.$(MODEXT),$(m)))
endif 

# C archive members
C =  $(subst .c,.o,$(wildcard *.c))

# fortran archive members
F  = $(subst .f,.o, $(wildcard *.f)) $(subst .F,.o, $(wildcard *.F)) $(subst .for,.o, $(wildcard *.for))
F90 = $(subst .f90,.o, $(wildcard *.f90)) $(subst .F90,.o, $(wildcard *.F90))
FALL = $(C) $(F) $(F90) 
#remove test programs
FMEM =$(filter-out ezspline_test%, $(filter-out ezspline_perf%, $(filter-out %lookup_test.o, $(filter-out qk_pspline%, $(filter-out ezspline_io_test.o, $(filter-out psplte%,$(FALL)))))))
MEM = $(foreach file,$(FMEM),$(ARC)($(file)))
OBJECTS = $(foreach file, $(FMEM),$(OBJDIR)/$(file))


# Test Programs
PRGS1 = pspltest lookup_test r8lookup_test 
PRGS2 = ezspline_perf_r4 ezspline_perf_r8 ezspline_test_r4 ezspline_test_r8 ezspline_io_test
PRGS3 = qk_pspline
PRGS4 = czspline_test
PRG1 =  $(foreach file,$(PRGS1),$(OBJ)/test/$(file))
PRG2 =  $(foreach file,$(PRGS2),$(OBJ)/test/$(file))
PRG3 =  $(foreach file,$(PRGS3),$(OBJ)/test/$(file))
PRG4 =  $(foreach file,$(PRGS4),$(OBJ)/test/$(file))

LDLIBS = -L$(OBJ)/lib $(LLOC) -lpspline -lportlib

ifndef LOCAL
   SLDLIBS = -L$(OBJ)/lib -lportlib -lezcdf
else
   SLDLIBS = $(OBJ)/lib/portlib.a $(OBJ)/lib/ezcdf.a
endif

.PHONY: clean realclean checkdir all

libs: FORCE $(ARC)
 
all:  libs exec
	@echo done  

$(OLDLIB): timestamp.inf
	@echo "--- DETECTED libpspline.a source update"
	@echo "--- Re-Making $(ARC)"
	@$(MAKE) libs

FORCE:
	@test -d $(OBJ)/lib || mkdir -p $(OBJ)/lib
	@test -d $(MDIR) || mkdir -p $(MDIR)
	@test -d $(OBJ)/obj/pspline || mkdir -p $(OBJ)/obj/pspline


 
# compile c
$(ARC)(%.o): %.c
	$(c_proc)
	$(make_ar)

#pre-process/compile free form fortran 90 code
$(ARC)(%.o): %.F90
	$(F90_proc)
	$(make_ar)

so: FORCE $(MODS) $(OBJ)/lib/libpspline.so

$(OBJ)/lib/libpspline.so: $(OBJECTS)
	$(LD) -o $@ $(X_LDFLAGS) $(OBJECTS) $(SLDLIBS) $(LUSRLIB) $(NETCDF) $(XLIBS) $(SHLIBS) $(XFLIB)


$(ARC)(%.o): %.for
	$(F_proc)

# compile free form fortran 90 code
#	$(FC90) $(FFLAGS) $(MODFLAGS) $(INCFLAGS) $(OUT) $(OBJDIR)/$*.o $<
$(ARC)(%.o): %.f90
	$(F90_proc)
	$(make_ar)


# compile fixed fortran
$(ARC)(%.o): %.f
	$(FC) $(FFLAGS) $(MODFLAGS) $(INCFLAGS) $(OUT) $(OBJDIR)/$*.o $<
	$(make_ar)

# pre-process/compile fixed code
$(ARC)(%.o): %.F
	$(F_proc)
	$(make_ar)

$(ARC): $(MODS) $(MEM)
	@$(RANLIB) $@ > /dev/null
#
#---------------------------------------------------

ckdirs: FORCE
	@test -d $(OBJ)/test || mkdir -p $(OBJ)/test
	@for i in $(DATA); do \
	(cd $(OBJ)/test; test -f $$i || ln -s $(srcdir)/$$i $$i;) done 

# don't compile main with dynamic, it requires too much stack space

do_exec: ckdirs $(PRG1) $(PRG2) $(PRG3) $(PRG4) $(OBJ)/test/pspltes8

exec: 
ifneq ("$(wildcard pspltes8.F)","")
	$(MAKE) do_exec
else
	@echo
endif

$(OBJ)/test/pspltes8: $(OBJDIR)/pspltes8.o $(ARC) 
	$(FC) $(LDFLAGS) -o $@ $< $(LDLIBS) $(XFLIB)

# NOTE: if HP would not require "STATIC"
# one could just use
# $(PROC_F)
#
$(OBJDIR)/pspltes8.o: pspltes8.F
	$(proc_F)


#	@$(PREFPP) $< $*.f $(INCFLAGS) $(DEFS) $(STATIC)
#	$(FC) $(FTEST) $(MODFLAGS) $(OUT) $@ $*.f
#	@if test -z "$(DEBUG)"; then \
#	rm -f $*.f*; fi	


$(OBJDIR)/%.o: %.cc 
	$(CXX) $(CXXFLAGS) -I../include/fpreproc $(DEFS) $(CINCL) -o $@ $<

# for TRANSP
$(OBJDIR)/%.o: %.for
	$(F_proc)

#	$(FC90) $(FTEST) $(MODFLAGS)  $(OUT) $(OBJDIR)/$*.o $<
$(OBJDIR)/%.o: %.f90
	$(F90_proc)

$(OBJDIR)/%.o: %.F90
	$(F90_proc)

$(OBJDIR)/%.o: %.F
	$(F_proc)

$(OBJDIR)/%.o: %.f
	$(FC) $(FTEST) $(OUT) $(OBJDIR)/$*.o $<

$(PRG1): $(OBJ)/test/%: $(OBJDIR)/%.o 
	$(FC) $(LDFLAGS) -o $@  $(OBJDIR)/$*.o $(LDLIBS)  $(XFLIB)

$(PRG2): $(OBJ)/test/%: $(OBJDIR)/%.o 
	$(FC90) $(LDFLAGS) -o $@  $(OBJDIR)/$*.o $(LDLIBS) -lezcdf $(NETCDF) $(XFLIB)

$(PRG3): $(OBJ)/test/%: $(OBJDIR)/%.o 
	$(FC) $(LDFLAGS) -o $@  $(OBJDIR)/$*.o $(LDLIBS) -lezcdf $(NETCDF) $(XFLIB)

$(PRG4): $(OBJ)/test/%: $(OBJDIR)/%.o 
	$(CXX) $(LDFLAGS) -o $@  $(OBJDIR)/$*.o $(LDLIBS) -lezcdf $(NETCDF) $(FORTLIBS) $(XFLIB)

# To link with shared lib:
#	$(CXX) $(LDFLAGS) -o $@  $(OBJDIR)/$*.o $(OBJ)/lib/libpspline.so -L../LINUX/lib 

chkins:
	@test -d $(MODDIR) ||  mkdir -p $(MODDIR)
	@test -d $(INCLDIR)/cpp ||  mkdir -p $(INCLDIR)/cpp
	@test -d $(LIBDIR) ||  mkdir -p $(LIBDIR)
	@test -d $(MANDIR)/man3 || mkdir -p $(MANDIR)/man3

install: chkins
	@if test -f $(ARC); then \
	  umask 133; cp $(ARC) $(LIBDIR)/; \
	else test -f $(LIBDIR)/libpspline.a || echo ">>> can not find $(ARC)";fi;
	@umask 133; cp pspline.3 $(MANDIR)/man3/
	$(foreach file,$(MODULES),$(shell cp $(MDIR)/$(file) $(MODDIR)/))
	$(foreach file,$(CZ_H),$(shell cp ../include/cpp/$(file) $(INCLDIR)/cpp))


uninstall:
	$(foreach m,$(MODULES),$(shell rm -f $(MODDIR)/$(m) ))
	rm -f $(LIBDIR)/libpspline.a
	rm -f $(MANDIR)/man3/pspline.3

clean: 
	@rm -f $(OBJDIR)/*
	@if test -d $(MDIR); then \
	  (cd $(MDIR); rm -f $(MODULES)); fi

realclean: clean
	@rm -f $(ARC)
	@rm -f $(PRG1) $(PRG2)
	@if test -d $(OBJ)/test; then \
	  (cd $(OBJ)/test; rm -f $(DATA)); fi
