#!/usr/bin/make
#----------------------------------------------------------------------------------------------------------------------------------
# make init

# shell
SHELL = /bin/bash
# no verbose
$(VERBOSE).SILENT:
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# User options
OSYSTEM   = uix
PRESET    = no
COMPILER  = intel
DEBUG     = no
F03STD    = no
PROFILING = no
OPTIMIZE  = no
OPENMP    = no
MPI       = no
R16P      = no
NULi      = no
NULj      = no
NULk      = no
WENO      = WENO
VACUUM    = no
SMSW      = SMSWz
HYBRID    = NOHYBRID
PPL       = no
LMA       = no
TECIO     = no
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# compiler specific rules
# GNU
WRN_GNU = -fmax-errors=0 -Wall -Wno-array-temporaries -Warray-bounds -Wcharacter-truncation -Wline-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wunderflow -Wextra -Wuninitialized
CHK_GNU = -fcheck=all
DEB_GNU = -fmodule-private -ffree-line-length-132 -fimplicit-none -ffpe-trap=invalid,overflow -fbacktrace -fdump-core -finit-real=nan #-fno-range-check  ,precision,denormal,underflow
STD_GNU = -std=f2003 -fall-intrinsics
OMP_GNU = -fopenmp
OPT_GNU = -O3
PRF_GNU =
# Intel
WRN_INT = -warn all
CHK_INT = -check all
DEB_INT = -debug all -extend-source 132 -traceback -gen-interfaces#-fpe-all=0 -fp-stack-check -fstack-protector-all -ftrapuv -no-ftz
STD_INT = -std03
OMP_INT = -openmp
OPT_INT = -O3 -ipo -inline all -ipo-jobs4 -vec-report1
PRF_INT = #-p
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
.PHONY : DEFAULTRULE
DEFAULTRULE: $(DEXE)OFF

.PHONY : help
help:
	@echo
	@echo -e '\033[1;31m Make options of OFF codes\033[0m'
	@echo
	@echo -e '\033[1;31m OS running selection: OSYSTEM=$(OSYSTEM)\033[0m\033[1m => default\033[0m'
	@echo -e '\033[1;31m  OSYSTEM=uix      \033[0m\033[1m => OS running is Unix/Linux \033[0m'
	@echo -e '\033[1;31m  OSYSTEM=win      \033[0m\033[1m => OS running is MS Windows \033[0m'
	@echo
	@echo -e '\033[1;31m Preset configurations: PRESET=$(PRESET)\033[0m\033[1m => default\033[0m'
	@echo -e '\033[1;31m  PRESET=debug     \033[0m\033[1m => debug (no optimized) and serial\033[0m'
	@echo -e '\033[1;31m  PRESET=fast      \033[0m\033[1m => optimized (no debug) and serial\033[0m'
	@echo -e '\033[1;31m  PRESET=openmp    \033[0m\033[1m => optimized (no debug) and openmp enabled\033[0m'
	@echo -e '\033[1;31m  PRESET=mpi       \033[0m\033[1m => optimized (no debug) and mpi enabled\033[0m'
	@echo -e '\033[1;31m  PRESET=openmp-mpi\033[0m\033[1m => optimized (no debug) and openmp-mpi enabled\033[0m'
	@echo
	@echo -e '\033[1;31m Compiler choice: COMPILER=$(COMPILER)\033[0m\033[1m => default\033[0m'
	@echo -e '\033[1;31m  COMPILER=gnu  \033[0m\033[1m => GNU gfortran          \033[0m'
	@echo -e '\033[1;31m  COMPILER=intel\033[0m\033[1m => Intel Fortran         \033[0m'
	@echo -e '\033[1;31m  COMPILER=pgi  \033[0m\033[1m => Portland Group Fortran\033[0m'
	@echo -e '\033[1;31m  COMPILER=g95  \033[0m\033[1m => free g95              \033[0m'
	@echo
	@echo -e '\033[1;31m Compiling options\033[0m'
	@echo -e '\033[1;31m  DEBUG=yes(no)    \033[0m\033[1m => on(off) debug                  (default $(DEBUG))\033[0m'
	@echo -e '\033[1;31m  F03STD=yes(no)   \033[0m\033[1m => on(off) check standard fortran (default $(F03STD))\033[0m'
	@echo -e '\033[1;31m  PROFILING=yes(no)\033[0m\033[1m => on(off) profile the code       (default $(PROFILING))\033[0m'
	@echo -e '\033[1;31m  OPTIMIZE=yes(no) \033[0m\033[1m => on(off) optimization           (default $(OPTIMIZE))\033[0m'
	@echo -e '\033[1;31m  OPENMP=yes(no)   \033[0m\033[1m => on(off) OpenMP directives      (default $(OPENMP))\033[0m'
	@echo -e '\033[1;31m  MPI=yes(no)      \033[0m\033[1m => on(off) MPI    directives      (default $(MPI)) \033[0m'
	@echo
	@echo -e '\033[1;31m Preprocessing options\033[0m'
	@echo -e '\033[1;31m  R16P=yes(no)   \033[0m\033[1m => on(off) definition of real with "128 bit" (default $(R16P))\033[0m'
	@echo -e '\033[1;31m  NULi=yes(no)   \033[0m\033[1m => on(off) nullify i direction (1D or 2D)    (default $(NULi))\033[0m'
	@echo -e '\033[1;31m  NULj=yes(no)   \033[0m\033[1m => on(off) nullify j direction (1D or 2D)    (default $(NULj))\033[0m'
	@echo -e '\033[1;31m  NULk=yes(no)   \033[0m\033[1m => on(off) nullify k direction (1D or 2D)    (default $(NULk))\033[0m'
	@echo -e '\033[1;31m  PPL=yes(no) \033[0m\033[1m => on(off) Positivity Preserving Limiter     (default $(PPL))\033[0m'
	@echo -e '\033[1;31m  LMA=yes(no) \033[0m\033[1m => on(off) Low Mach number Adjsutment        (default $(LMA))\033[0m'
	@echo -e '\033[1;31m  WENO=WENO/WENOZ/WENOM\033[0m\033[1m WENO algorithm (default $(WENO))\033[0m'
	@echo -e '\033[1m   WENO  => Original Jiang-Shu\033[0m'
	@echo -e '\033[1m   WENOZ => Improved Borges-Carmona-Costa-Don\033[0m'
	@echo -e '\033[1m   WENOM => Improved Henrick-Aslam-Powers    \033[0m'
	@echo -e '\033[1;31m  SMSW=SMSWz/SMSWliu/SMSWvanleer/SMSWvanalbada/SMSWharten\033[0m\033[1m smoothness switchg algorithm (default $(SMSW))\033[0m'
	@echo -e '\033[1m   SMSWz         => Riemann-solver-like algorithm      \033[0m'
	@echo -e '\033[1m   SMSWliu       => Liu algorithm                      \033[0m'
	@echo -e '\033[1m   SMSWvanleer   => van Leer slope limiter algorithm   \033[0m'
	@echo -e '\033[1m   SMSWvanAlbada => van Albada slope limiter  algorithm\033[0m'
	@echo -e '\033[1m   SMSWharten    => Harten slope limiter algorithm     \033[0m'
	@echo -e '\033[1;31m  HYBRID=NOHYBRID/HYBRID/HYBRIDC\033[0m\033[1m hybrid scheme (default $(HYBRID))\033[0m'
	@echo -e '\033[1m   NOHYBRID => no hybrid scheme                 \033[0m'
	@echo -e '\033[1m   HYBRID   => hybrid weno/weno_optimal scheme  \033[0m'
	@echo -e '\033[1m   HYBRIDC  => hybrid weno/noweno_central scheme\033[0m'
	@echo
	@echo -e '\033[1;31m External libraries\033[0m'
	@echo -e '\033[1;31m  TECIO=yes(no)\033[0m\033[1m => on(off) Tecplot IO library linking (default $(TECIO))\033[0m'
	@echo
	@echo -e '\033[1;31m Provided Rules\033[0m'
	@echo -e '\033[1;31m  Defualt rule =>\033[0m\033[1m OFF\033[0m'
	@echo -e '\033[1;31m  help         =>\033[0m\033[1m printing this help message\033[0m'
	@echo -e '\033[1;31m  OFF          =>\033[0m\033[1m building OFF code\033[0m'
	@echo -e '\033[1;31m  IBM          =>\033[0m\033[1m building IBM code\033[0m'
	@echo -e '\033[1;31m  POG          =>\033[0m\033[1m building POG code\033[0m'
	@echo -e '\033[1;31m  cleanobj     =>\033[0m\033[1m cleaning compiled object\033[0m'
	@echo -e '\033[1;31m  cleanmod     =>\033[0m\033[1m cleaning .mod files\033[0m'
	@echo -e '\033[1;31m  cleanmsg     =>\033[0m\033[1m cleaning make-log massage files\033[0m'
	@echo -e '\033[1;31m  cleanexe     =>\033[0m\033[1m cleaning executable files\033[0m'
	@echo -e '\033[1;31m  clean        =>\033[0m\033[1m running cleanobj, cleanmod and cleanmsg\033[0m'
	@echo -e '\033[1;31m  cleanall     =>\033[0m\033[1m running clean and cleanexe\033[0m'
	@echo -e '\033[1;31m  tar          =>\033[0m\033[1m creating a tar archive of the project\033[0m'
	@echo -e '\033[1;31m  doc          =>\033[0m\033[1m building the documentation\033[0m'
	@echo
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# directory & file
DSRC  = ./src/
DOBJ  = ./obj/
DMOD  = ./mod/
DLIB  = ./lib/
DEXE  = ./
VPATH = $(DSRC) $(DOBJ) $(DMOD) $(DLIB)

MKDIRS = $(DOBJ) $(DMOD) $(DEXE)
LBITS := $(shell getconf LONG_BIT)
ifeq "$(TECIO)" "yes"
  PREPROC = -DTECIO
	ifeq ($(LBITS),64)
  LIBS = $(DLIB)64bit/tecio64.a $(DLIB)64bit/libstdc++64.5.0.7.so
	else
		LIBS = $(DLIB)32bit/tecio.a $(DLIB)32bit/libstdc++.5.0.7.so
	endif
	STATIC =
else
  PREPROC =
  LIBS =
	STATIC = #-static
endif
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# compiling and linking options
# OS running
ifeq "$(OSYSTEM)" "uix"
  PREPROC := $(PREPROC) -D_OSYSTEMuix
endif
ifeq "$(OSYSTEM)" "win"
  PREPROC := $(PREPROC) -D_OSYSTEMwin
endif
# presets
ifeq "$(PRESET)" "debug"
  DEBUG     = yes
  F03STD    = yes
  PROFILING = no
  OPTIMIZE  = no
  OPENMP    = no
  MPI       = no
endif
ifeq "$(PRESET)" "fast"
  DEBUG     = no
  F03STD    = no
  PROFILING = no
  OPTIMIZE  = yes
  OPENMP    = no
  MPI       = no
endif
ifeq "$(PRESET)" "openmp"
  DEBUG     = no
  F03STD    = no
  PROFILING = no
  OPTIMIZE  = yes
  OPENMP    = yes
  MPI       = no
endif
ifeq "$(PRESET)" "mpi"
  DEBUG     = no
  F03STD    = no
  PROFILING = no
  OPTIMIZE  = yes
  OPENMP    = no
  MPI       = yes
endif
ifeq "$(PRESET)" "openmp-mpi"
  DEBUG     = no
  F03STD    = no
  PROFILING = no
  OPTIMIZE  = yes
  OPENMP    = yes
  MPI       = yes
endif
# setting rules according user options
ifeq "$(COMPILER)" "gnu"
  FC = gfortran
  OPTSC = -cpp -c -J$(DMOD) $(STATIC) -fprotect-parens
  OPTSL =                   $(STATIC) -fprotect-parens
  WRN = $(WRN_GNU)
  CHK = $(CHK_GNU)
  DEB = $(DEB_GNU)
  STD = $(STD_GNU)
  OMP = $(OMP_GNU)
  OPT = $(OPT_GNU)
  PRF = $(PRF_GNU)
endif
ifeq "$(COMPILER)" "intel"
  FC = ifort
  OPTSC = -cpp -c -module $(DMOD) $(STATIC) -assume protect_parens -assume realloc_lhs -fp-model source
  OPTSL =                         $(STATIC) -assume protect_parens -assume realloc_lhs -fp-model source
  WRN = $(WRN_INT)
  CHK = $(CHK_INT)
  DEB = $(DEB_INT)
  STD = $(STD_INT)
  OMP = $(OMP_INT)
  OPT = $(OPT_INT)
  PRF = $(PRF_INT)
endif
ifeq "$(DEBUG)" "yes"
  PREPROC := $(PREPROC) -DDEBUG
  OPTSC := $(OPTSC) -O0 -C -g $(WRN) $(CHK) $(DEB)
  OPTSL := $(OPTSL) -O0 -C -g $(WRN) $(CHK) $(DEB)
endif
ifeq "$(F03STD)" "yes"
  OPTSC := $(OPTSC) $(STD)
  OPTSL := $(OPTSL) $(STD)
endif
ifeq "$(PROFILING)" "yes"
  PREPROC := $(PREPROC) -DPROFILING
  OPTSC := $(OPTSC) $(PRF)
  OPTSL := $(OPTSL) $(PRF)
endif
ifeq "$(OPTIMIZE)" "yes"
  OPTSC := $(OPTSC) $(OPT)
  OPTSL := $(OPTSL) $(OPT)
endif
ifeq "$(OPENMP)" "yes"
  PREPROC := $(PREPROC) -DOPENMP
  OPTSC := $(OPTSC) $(OMP)
  OPTSL := $(OPTSL) $(OMP)
endif
ifeq "$(MPI)" "yes"
  PREPROC := $(PREPROC) -D_MPI
  FC = mpif90
endif
# pre-processing options
# R16 precision
R16PCHK = (Unknown R16P switch) Used default R16P=no
ifeq "$(R16P)" "no"
  R16PCHK = (Known R16P switch) Used R16P=$(R16P)
endif
ifeq "$(R16P)" "yes"
  R16PCHK = (Known R16P switch) Used R16P=$(R16P)
  PREPROC := $(PREPROC) -Dr16p
endif
# 1D or 2D solver
NULiCHK = (Unknown NULi switch) Used default NULi=no
ifeq "$(NULi)" "no"
  NULiCHK = (Known NULi switch) Used NULi=$(NULi)
endif
ifeq "$(NULi)" "yes"
  NULiCHK = (Known NULi switch) Used NULi=$(NULi)
  PREPROC := $(PREPROC) -DNULi
endif
NULjCHK = (Unknown NULj switch) Used default NULj=no
ifeq "$(NULj)" "no"
  NULjCHK = (Known NULj switch) Used NULj=$(NULj)
endif
ifeq "$(NULj)" "yes"
  NULjCHK = (Known NULj switch) Used NULj=$(NULj)
  PREPROC := $(PREPROC) -DNULj
endif
NULkCHK = (Unknown NULk switch) Used default NULk=no
ifeq "$(NULk)" "no"
  NULkCHK = (Known NULk switch) Used NULk=$(NULk)
endif
ifeq "$(NULk)" "yes"
  NULkCHK = (Known NULk switch) Used NULk=$(NULk)
  PREPROC := $(PREPROC) -DNULk
endif
# Positivity-Preserving-Limiter
PPLCHK = (Unknown PPL switch) Used default PPL=no
ifeq "$(PPL)" "no"
  PPLCHK = (Known PPL switch) Used PPL=$(PPL)
endif
ifeq "$(PPL)" "yes"
  PPLCHK = (Known PPL switch) Used PPL=$(PPL)
  PREPROC := $(PREPROC) -DPPL
endif
# Low Mach number Adjustment
LMACHK = (Unknown LMA switch) Used default LMA=no
ifeq "$(LMA)" "no"
  LMACHK = (Known LMA switch) Used LMA=$(LMA)
endif
ifeq "$(LMA)" "yes"
  LMACHK = (Known LMA switch) Used LMA=$(LMA)
  PREPROC := $(PREPROC) -DLMA
endif
# WENO algorithms
WENOCHK = (Unknown WENO switch) Used default WENO=WENO
ifeq "$(WENO)" "WENO"
  WENOCHK = (Known WENO switch) Used WENO=$(WENO)
endif
ifeq "$(WENO)" "WENOZ"
  WENOCHK = (Known WENO switch) Used WENO=$(WENO)
  PREPROC := $(PREPROC) -D$(WENO)
endif
ifeq "$(WENO)" "WENOM"
  WENOCHK = (Known WENO switch) Used WENO=$(WENO)
  PREPROC := $(PREPROC) -D$(WENO)
endif
# Smoothness switch algorithm
SMSWCHK = (Unknown SMSW switch) Used default SMSW=SMSWz
ifeq "$(SMSW)" "SMSWz"
  SMSWCHK = (Known SMSW switch) Used SMSW=$(SMSW)
  PREPROC := $(PREPROC) -D$(SMSW)
endif
ifeq "$(SMSW)" "SMSWliu"
  SMSWCHK = (Known SMSW switch) Used SMSW=$(SMSW)
  PREPROC := $(PREPROC) -D$(SMSW)
endif
ifeq "$(SMSW)" "SMSWvanleer"
  SMSWCHK = (Known SMSW switch) Used SMSW=$(SMSW)
  PREPROC := $(PREPROC) -D$(SMSW)
endif
ifeq "$(SMSW)" "SMSWvanalbada"
  SMSWCHK = (Known SMSW switch) Used SMSW=$(SMSW)
  PREPROC := $(PREPROC) -D$(SMSW)
endif
ifeq "$(SMSW)" "SMSWharten"
  SMSWCHK = (Known SMSW switch) Used SMSW=$(SMSW)
  PREPROC := $(PREPROC) -D$(SMSW)
endif
# Hybrid scheme
HYBRIDCHK = (Unknown HYBRID switch) Used default HYBRID=NOHYBRID
ifeq "$(HYBRID)" "NOHYBRID"
  HYBRIDCHK = (Known HYBRID switch) Used HYBRID=$(HYBRID)
endif
ifeq "$(HYBRID)" "HYBRID"
  HYBRIDCHK = (Known HYBRID switch) Used HYBRID=$(HYBRID)
  PREPROC := $(PREPROC) -D$(HYBRID)
endif
ifeq "$(HYBRID)" "HYBRIDC"
  HYBRIDCHK = (Known HYBRID switch) Used HYBRID=$(HYBRID)
  PREPROC := $(PREPROC) -D$(HYBRID)
endif

PREPROC := $(PREPROC) -D_COMPILER='$(shell $(FC) --version | head -n 1)' -D_COMPFLAG='$(OPTSC)' -D_COMPPROC='$(PREPROC)' -D_COMPLIBS='$(LIBS)'

WHICHFC = $(shell which $(FC))
PRINTCHK = "\\033[1;31m Compiler used \\033[0m\\033[1m $(COMPILER) => $(WHICHFC)\\033[0m \n\
            \\033[1;31mSource dir    \\033[0m\\033[1m $(DSRC)\\033[0m \n\
            \\033[1;31mLibraries     \\033[0m\\033[1m $(LIBS)\\033[0m \n \n\
            \\033[1;31m Debug         \\033[0m\\033[1m $(DEBUG)\\033[0m \n\
            \\033[1;31m F-standard    \\033[0m\\033[1m $(F03STD)\\033[0m \n\
            \\033[1;31m Profiling     \\033[0m\\033[1m $(PROFILING)\\033[0m \n\
            \\033[1;31m Optimize      \\033[0m\\033[1m $(OPTIMIZE)\\033[0m \n\
            \\033[1;31m OpenMP        \\033[0m\\033[1m $(OPENMP)\\033[0m \n\
            \\033[1;31m MPI           \\033[0m\\033[1m $(MPI)\\033[0m \n\
            \\033[1;31m R16P          \\033[0m\\033[1m $(R16PCHK)\\033[0m \n\
            \\033[1;31m Nullify i     \\033[0m\\033[1m $(NULiCHK)\\033[0m \n\
            \\033[1;31m Nullify j     \\033[0m\\033[1m $(NULjCHK)\\033[0m \n\
            \\033[1;31m Nullify k     \\033[0m\\033[1m $(NULkCHK)\\033[0m \n\
            \\033[1;31m PP Limiter    \\033[0m\\033[1m $(PPLCHK)\\033[0m \n\
            \\033[1;31m LMA Limiter   \\033[0m\\033[1m $(LMACHK)\\033[0m \n\
            \\033[1;31m WENO          \\033[0m\\033[1m $(WENOCHK)\\033[0m \n\
            \\033[1;31m Riemann Solver\\033[0m\\033[1m $(RSUCHK)\\033[0m \n\
            \\033[1;31m Smooth switch \\033[0m\\033[1m $(SMSWCHK)\\033[0m \n\
            \\033[1;31m Hybrid scheme \\033[0m\\033[1m $(HYBRIDCHK)\\033[0m \n\
            \\033[1;31m TECIO         \\033[0m\\033[1m $(TECIO)\\033[0m"
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
.PHONY : PRINTINFO
.NOTPARALLEL : PRINTINFO
PRINTINFO:
	@echo | tee make.log
	@echo -e $(PRINTCHK) | tee -a make.log
	@echo | tee -a make.log
	@echo -e "\033[1;31m Compiling options\033[0m" | tee -a make.log
	@echo -e "\033[1m [$(OPTSC)]\033[0m" | tee -a make.log
	@echo | tee -a make.log
	@echo -e "\033[1;31m Linking options \033[0m" | tee -a make.log
	@echo -e "\033[1m [$(OPTSL)]\033[0m" | tee -a make.log
	@echo | tee -a make.log
	@echo -e "\033[1;31m Preprocessing options \033[0m" | tee -a make.log
	@echo -e "\033[1m [$(PREPROC)]\033[0m" | tee -a make.log
	@echo | tee -a make.log

.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@

.PHONY : cleanobj
cleanobj:
	@echo -e "\033[1;31m deleting objects \033[0m" | tee make.log
	@rm -fr $(DOBJ)

.PHONY : cleanmod
cleanmod:
	@echo -e "\033[1;31m deleting mods \033[0m" | tee -a make.log
	@rm -fr $(DMOD)

.PHONY : cleanexe
cleanexe:
	@echo -e "\033[1;31m deleting exes \033[0m" | tee -a make.log
	@rm -f $(addprefix $(DEXE),$(EXES))

.PHONY : cleanmsg
cleanmsg:
	@rm -f diagnostic_messages
	@rm -f error_messages

.PHONY : clean
clean: cleanobj cleanmod cleanmsg

.PHONY : cleanall
cleanall: clean cleanexe

.PHONY : tar
tar: cleanall
	@echo -e "\033[1;31m Creating tar archive of the code \033[0m" | tee make.log
	@rm -rf OFF
	@mkdir -p OFF
	@cp -r EXAMPLE.md examples inputs-template lib makefile README.md util .doxygenconfig .gitignore OFF/
	@cp -rL src OFF/
	@tar czf OFF.tgz OFF
	@rm -rf OFF

.PHONY : doc
doc:
	@echo -e "\033[1;31m Building documentation\033[0m" | tee make.log
	@doxygen .doxygenconfig

.PHONY : codes
codes: $(DEXE)IBM $(DEXE)OFF $(DEXE)POG
#----------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------
# rules of linking and compiling
COTEXT  = -e "\033[1;31m Compiling\033[0m\033[1m $(<F)\033[0m"
LITEXT  = -e "\033[1;31m Assembling\033[0m\033[1m $@\033[0m"
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

$(DEXE)IBM : PRINTINFO $(MKDIRS) $(DOBJ)ibm.o
	@rm -f $(filter-out $(DOBJ)ibm.o,$(EXESOBJ))
	@echo | tee -a make.log
	@echo $(LITEXT) | tee -a make.log
	@$(FC) $(OPTSL) $(PREPROC) $(DOBJ)*.o $(LIBS) -o $@ 1>> diagnostic_messages 2>> error_messages
EXES := $(EXES) IBM

$(DEXE)OFF : PRINTINFO $(MKDIRS) $(DOBJ)off.o
	@rm -f $(filter-out $(DOBJ)off.o,$(EXESOBJ))
	@echo | tee -a make.log
	@echo $(LITEXT) | tee -a make.log
	@$(FC) $(OPTSL) $(PREPROC) $(DOBJ)*.o $(LIBS) -o $@ 1>> diagnostic_messages 2>> error_messages
EXES := $(EXES) OFF

$(DEXE)POG : PRINTINFO $(MKDIRS) $(DOBJ)pog.o
	@rm -f $(filter-out $(DOBJ)pog.o,$(EXESOBJ))
	@echo | tee -a make.log
	@echo $(LITEXT) | tee -a make.log
	@$(FC) $(OPTSL) $(PREPROC) $(DOBJ)*.o $(LIBS) -o $@ 1>> diagnostic_messages 2>> error_messages
EXES := $(EXES) POG

$(DOBJ)data_type_adimensional.o : Data_Type_Adimensional.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_amrblock.o : Data_Type_AMRBlock.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_hashid.o \
	$(DOBJ)data_type_hashtcell.o \
	$(DOBJ)data_type_hashtface.o \
	$(DOBJ)data_type_hashtnode.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_bc.o : Data_Type_BC.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_cell_indexes.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_bc_in1.o : Data_Type_BC_in1.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_primitive.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_block_bc.o : Data_Type_Block_BC.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_xml_tag.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_block_dimensions.o : Data_Type_Block_Dimensions.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_xml_tag.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_block_extents.o : Data_Type_Block_Extents.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)data_type_xml_tag.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_cell.o : Data_Type_Cell.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_species.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_variables_conversions.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_cell_indexes.o : Data_Type_Cell_Indexes.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_compiledcode.o : Data_Type_CompiledCode.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_conservative.o : Data_Type_Conservative.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_face.o : Data_Type_Face.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_base.o : Data_Type_File_Base.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_bc.o : Data_Type_File_BC.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_parallel.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_xml_tag.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_blocks_cartesian.o : Data_Type_File_Blocks_Cartesian.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_mesh_dimensions.o \
	$(DOBJ)data_type_region.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_fluid.o : Data_Type_File_Fluid.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_parallel.o \
	$(DOBJ)data_type_species.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_time_step.o \
	$(DOBJ)data_type_xml_tag.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_gnu.o : Data_Type_File_GNU.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_postprocess.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_math.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_ibm_options.o : Data_Type_File_IBM_Options.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_varying_string.o \
	$(DOBJ)data_type_xml_tag.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_icemcfd.o : Data_Type_File_ICEMCFD.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_lock.o : Data_Type_File_Lock.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_time.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_mesh.o : Data_Type_File_Mesh.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_mesh_dimensions.o \
	$(DOBJ)data_type_parallel.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_xml_tag.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_off_options.o : Data_Type_File_OFF_Options.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_xml_tag.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_procmap.o : Data_Type_File_Procmap.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_mesh_dimensions.o \
	$(DOBJ)data_type_parallel.o \
	$(DOBJ)data_type_xml_tag.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_profile.o : Data_Type_File_Profile.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_residuals.o : Data_Type_File_Residuals.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_files.o : Data_Type_Files.f90 \
	$(DOBJ)data_type_file_bc.o \
	$(DOBJ)data_type_file_fluid.o \
	$(DOBJ)data_type_file_blocks_cartesian.o \
	$(DOBJ)data_type_file_gnu.o \
	$(DOBJ)data_type_file_ibm_options.o \
	$(DOBJ)data_type_file_lock.o \
	$(DOBJ)data_type_file_mesh.o \
	$(DOBJ)data_type_file_off_options.o \
	$(DOBJ)data_type_file_procmap.o \
	$(DOBJ)data_type_file_profile.o \
	$(DOBJ)data_type_file_species.o \
	$(DOBJ)data_type_file_solver_options.o \
	$(DOBJ)data_type_file_tec.o \
	$(DOBJ)data_type_file_vtk.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_solver_options.o : Data_Type_File_Solver_Options.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_adimensional.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_space_step.o \
	$(DOBJ)data_type_time_step.o \
	$(DOBJ)data_type_xml_tag.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_species.o : Data_Type_File_Species.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_species.o \
	$(DOBJ)data_type_xml_tag.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_tec.o : Data_Type_File_Tec.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_postprocess.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_file_vtk.o : Data_Type_File_VTK.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_file_base.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_postprocess.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_math.o \
	$(DOBJ)lib_vtk_io.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_global.o : Data_Type_Global.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_adimensional.o \
	$(DOBJ)data_type_bc_in1.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_compiledcode.o \
	$(DOBJ)data_type_face.o \
	$(DOBJ)data_type_file_profile.o \
	$(DOBJ)data_type_mesh_dimensions.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_parallel.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_space_step.o \
	$(DOBJ)data_type_species.o \
	$(DOBJ)data_type_time_step.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_parallel.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_hashid.o : Data_Type_HashID.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)lib_morton.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_hashtcell.o : Data_Type_HashTCell.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_hashid.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_hashtface.o : Data_Type_HashTFace.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_face.o \
	$(DOBJ)data_type_hashid.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_hashtnode.o : Data_Type_HashTNode.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_hashid.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_mesh_dimensions.o : Data_Type_Mesh_Dimensions.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_block_dimensions.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_os.o : Data_Type_OS.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_parallel.o : Data_Type_Parallel.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_postprocess.o : Data_Type_PostProcess.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_sblock.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_primitive1d.o : Data_Type_Primitive1D.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_species.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_primitive.o : Data_Type_Primitive.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_species.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)data_type_xml_tag.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_region.o : Data_Type_Region.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)data_type_xml_tag.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_riemann_conservative1d.o : Data_Type_Riemann_Conservative1D.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_riemann_interstate1d.o : Data_Type_Riemann_InterState1D.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_riemann_primitive1d.o : Data_Type_Riemann_Primitive1D.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_sblock.o : Data_Type_SBlock.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_block_bc.o \
	$(DOBJ)data_type_block_dimensions.o \
	$(DOBJ)data_type_block_extents.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_face.o \
	$(DOBJ)data_type_mesh_dimensions.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_region.o \
	$(DOBJ)data_type_species.o \
	$(DOBJ)data_type_space_step.o \
	$(DOBJ)data_type_time_step.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_fluxes_convective.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_runge_kutta.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o \
	$(DOBJ)lib_variables_conversions.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_sl_list.o : Data_Type_SL_List.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_space_step.o : Data_Type_Space_Step.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_specie.o : Data_Type_Specie.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_xml_tag.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_species.o : Data_Type_Species.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_specie.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_tensor.o : Data_Type_Tensor.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_time.o : Data_Type_Time.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_time_step.o : Data_Type_Time_Step.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_varying_string.o : Data_Type_Varying_String.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_vector.o : Data_Type_Vector.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_xml_tag.o : Data_Type_XML_Tag.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_varying_string.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)ibm.o : IBM.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_files.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)ir_precision.o : IR_Precision.f90
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_base64.o : Lib_Base64.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_fluidynamic.o : Lib_Fluidynamic.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_face.o \
	$(DOBJ)data_type_files.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_fluxes_convective.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_math.o \
	$(DOBJ)lib_profiling.o \
	$(DOBJ)lib_runge_kutta.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o \
	$(DOBJ)lib_parallel.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_fluxes_convective.o : Lib_Fluxes_Convective.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_face.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_primitive1d.o \
	$(DOBJ)data_type_riemann_conservative1d.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)data_type_species.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_riemann_solvers.o \
	$(DOBJ)lib_weno.o \
	$(DOBJ)lib_variables_conversions.o \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_face.o \
	$(DOBJ)data_type_primitive1d.o \
	$(DOBJ)data_type_species.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_fluxes_diffusive.o : Lib_Fluxes_Diffusive.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_tensor.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_io_misc.o : Lib_IO_Misc.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_math.o : Lib_Math.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_sl_list.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_morton.o : Lib_Morton.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_multigrid.o : Lib_Multigrid.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_fluidynamic.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_parallel.o : Lib_Parallel.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_mesh_dimensions.o \
	$(DOBJ)data_type_parallel.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_profiling.o : Lib_Profiling.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann_cvl_solver.o : Lib_Riemann_CVL_Solver.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)data_type_riemann_interstate1d.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann_hllc_solver.o : Lib_Riemann_HLLC_Solver.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_conservative1d.o \
	$(DOBJ)data_type_riemann_interstate1d.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)lib_riemann_cvl_solver.o \
	$(DOBJ)lib_riemann_pvl_solver.o \
	$(DOBJ)lib_riemann_tr_solver.o \
	$(DOBJ)lib_riemann_ts_solver.o \
	$(DOBJ)lib_riemann_z_solver.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)data_type_riemann_interstate1d.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann_lf_solver.o : Lib_Riemann_LF_Solver.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_conservative1d.o \
	$(DOBJ)data_type_riemann_interstate1d.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)lib_riemann_pvl_solver.o \
	$(DOBJ)lib_riemann_z_solver.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)data_type_riemann_interstate1d.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann_pvl_solver.o : Lib_Riemann_PVL_Solver.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)data_type_riemann_interstate1d.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann_solvers.o : Lib_Riemann_Solvers.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_conservative1d.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)lib_riemann_hllc_solver.o \
	$(DOBJ)lib_riemann_lf_solver.o \
	$(DOBJ)data_type_riemann_conservative1d.o \
	$(DOBJ)data_type_riemann_primitive1d.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann_tr_solver.o : Lib_Riemann_TR_Solver.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)data_type_riemann_interstate1d.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann_ts_solver.o : Lib_Riemann_TS_Solver.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)data_type_riemann_interstate1d.o \
	$(DOBJ)lib_riemann_pvl_solver.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann_z_solver.o : Lib_Riemann_Z_Solver.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)data_type_riemann_interstate1d.o \
	$(DOBJ)lib_riemann_pvl_solver.o \
	$(DOBJ)lib_riemann_tr_solver.o \
	$(DOBJ)lib_riemann_ts_solver.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_runge_kutta.o : Lib_Runge_Kutta.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_thermodynamic_laws_ideal.o : Lib_Thermodynamic_Laws_Ideal.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_variables_conversions.o : Lib_Variables_Conversions.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_primitive1d.o \
	$(DOBJ)data_type_riemann_primitive1d.o \
	$(DOBJ)data_type_species.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_vtk_io.o : Lib_VTK_IO.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)lib_base64.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_weno.o : Lib_WENO.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)off.o : OFF.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_files.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_fluxes_convective.o \
	$(DOBJ)lib_riemann_solvers.o \
	$(DOBJ)lib_runge_kutta.o \
	$(DOBJ)lib_weno.o \
	$(DOBJ)lib_parallel.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)pog.o : POG.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_files.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_postprocess.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)%.o : %.f90
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $(PREPROC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#-----------------------------------------------------------------------------------------------------------------------------------
