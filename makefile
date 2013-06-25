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
PRESET    = no
COMPILER  = intel
DEBUG     = no
F03STD    = no
PROFILING = no
OPTIMIZE  = no
OPENMP    = no
MPI       = no
R16P      = no
SYSless   = no
NULi      = no
NULj      = no
NULk      = no
WENO      = WENO
RECV      = RECVC
RSU       = HLLCp
VACUUM    = no
WS        = WSup
SMSW      = SMSWz
HYBRID    = NOHYBRID
PPL       = no
LMA       = no
TECIO     = no

.PHONY : DEFAULTRULE
DEFAULTRULE: $(DEXE)OFF

.PHONY : help
help:
	@echo
	@echo -e '\033[1;31m Make options of OFF codes\033[0m'
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
	@echo -e '\033[1;31m  SYSless=yes(no)\033[0m\033[1m => on(off) system call functions             (default $(SYSless))\033[0m'
	@echo -e '\033[1;31m  NULi=yes(no)   \033[0m\033[1m => on(off) nullify i direction (1D or 2D)    (default $(NULi))\033[0m'
	@echo -e '\033[1;31m  NULj=yes(no)   \033[0m\033[1m => on(off) nullify j direction (1D or 2D)    (default $(NULj))\033[0m'
	@echo -e '\033[1;31m  NULk=yes(no)   \033[0m\033[1m => on(off) nullify k direction (1D or 2D)    (default $(NULk))\033[0m'
	@echo -e '\033[1;31m  PPL=yes(no)    \033[0m\033[1m => on(off) Positivity Preserving Limiter     (default $(PPL))\033[0m'
	@echo -e '\033[1;31m  LMA=yes(no)    \033[0m\033[1m => on(off) Low Mach number Adjsutment        (default $(LMA))\033[0m'
	@echo -e '\033[1;31m  WENO=WENO/WENOZ/WENOM\033[0m\033[1m WENO algorithm (default $(WENO))\033[0m'
	@echo -e '\033[1m   WENO  => Original Jiang-Shu\033[0m'
	@echo -e '\033[1m   WENOZ => Improved Borges-Carmona-Costa-Don\033[0m'
	@echo -e '\033[1m   WENOM => Improved Henrick-Aslam-Powers    \033[0m'
	@echo -e '\033[1;31m  RECV=RECVC/RECVP\033[0m\033[1m reconstruction variables type (default $(RECV))\033[0m'
	@echo -e '\033[1m   RECVC => reconstruction in (local) characteristic variables\033[0m'
	@echo -e '\033[1m   RECVP => reconstruction in primitive variables             \033[0m'
	@echo -e '\033[1;31m  RSU=HLLCb/HLLCc/HLLCp/HLLCt/HLLCz/EXA/PVL/TR/TS/APRS/ALFR/LF/LFz/ROE\033[0m\033[1m Riemann solver algorithm (default $(RSU))\033[0m'
	@echo -e '\033[1m   HLLCb => Approximate HLLC solver using BCLC waves speed estimation         \033[0m'
	@echo -e '\033[1m   HLLCc => Approximate HLLC solver using CVL  waves speed estimation         \033[0m'
	@echo -e '\033[1m   HLLCp => Approximate HLLC solver using PVL  waves speed estimation         \033[0m'
	@echo -e '\033[1m   HLLCt => Approximate HLLC solver using TR   waves speed estimation         \033[0m'
	@echo -e '\033[1m   HLLCz => Approximate HLLC solver using Z    waves speed estimation         \033[0m'
	@echo -e '\033[1m   EXA   => Exact solver                                                      \033[0m'
	@echo -e '\033[1m   PVL   => Approximate PVL solver                                            \033[0m'
	@echo -e '\033[1m   TR    => Approximate TR solver                                             \033[0m'
	@echo -e '\033[1m   TS    => Approximate TS solver                                             \033[0m'
	@echo -e '\033[1m   APRS  => Approximate APRS solver                                           \033[0m'
	@echo -e '\033[1m   ALFR  => Approximate ALFR solver                                           \033[0m'
	@echo -e '\033[1m   LFp   => Approximate Lax-Friedrichs solver using PVL waves speed estimation\033[0m'
	@echo -e '\033[1m   LFz   => Approximate Lax-Friedrichs solver using Z   waves speed estimation\033[0m'
	@echo -e '\033[1m   ROE   => Approximate Roe solver                                            \033[0m'
	@echo -e '\033[1;31m  WS=WSu/WSup\033[0m\033[1m Waves Speed estimation algorithm (default $(WS))\033[0m'
	@echo -e '\033[1m   WSu  => WavesSpeed14u  or WavesSpeed1234u  algorithm\033[0m'
	@echo -e '\033[1m   WSup => WavesSpeed14up or WavesSpeed1234up algorithm\033[0m'
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
else
  PREPROC =
  LIBS =
endif
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# compiling and linking options
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
DEB_INT = -debug all -extend-source 132 -fpe-all=0 -fp-stack-check -fstack-protector-all -ftrapuv -no-ftz -traceback -gen-interfaces
STD_INT = -std03
OMP_INT = -openmp
OPT_INT = -O3 -ipo -inline all -ipo-jobs4 -vec-report1
PRF_INT = #-p
# setting rules according user options
ifeq "$(COMPILER)" "gnu"
  FC = gfortran
  OPTSC = -cpp -c -J$(DMOD) -static -fprotect-parens -fno-realloc-lhs
  OPTSL =
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
  OPTSC = -cpp -c -module $(DMOD) -static -assume protect_parens -assume norealloc_lhs -fp-model source
  OPTSL =
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
  PREPROC := $(PREPROC) -DMPI2
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
# SYSless
SYSlessCHK = (Unknown SYSless switch) Used default SYSless=no
ifeq "$(SYSless)" "no"
  SYSlessCHK = (Known SYSless switch) Used SYSless=$(SYSless)
endif
ifeq "$(SYSless)" "yes"
  SYSlessCHK = (Known SYSless switch) Used SYSless=$(SYSless)
  PREPROC := $(PREPROC) -DSYSTEMless
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
# reconstruction variables type
ifeq "$(RECV)" "RECVC"
  RECVCHK = (Known RECV switch) Used RECV=$(RECV)
  PREPROC := $(PREPROC) -D$(RECV)
else
  ifeq "$(RECV)" "RECVP"
    RECVCHK = (Known RECV switch) Used RECV=$(RECV)
    PREPROC := $(PREPROC) -D$(RECV)
  else
    RECVCHK = (Unknown RECV switch) Used default RECV=RECVC
    PREPROC := $(PREPROC) -DRECVC
  endif
endif
# Riemann solver
RSUCHK = (Unknown RSU switch) Used default RSU=HLLCp
ifeq "$(RSU)" "HLLCb"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "HLLCc"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "HLLCp"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "HLLCt"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "HLLCz"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "EXA"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "PVL"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "TR"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "TS"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "APRS"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "ALFR"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "ROE"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "LFp"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
ifeq "$(RSU)" "LFz"
  RSUCHK = (Known RSU switch) Used RSU=$(RSU)
  RS = -DRS$(RSU)
endif
PREPROC := $(PREPROC) $(RS)
# Waves Speed estimation algorithm
WSCHK = (Unknown WS switch) Used default WS=WSup
ifeq "$(WS)" "WSu"
  WSCHK = (Known WS switch) Used WS=$(WS)
  PREPROC := $(PREPROC) -D$(WS)
endif
ifeq "$(WS)" "WSup"
  WSCHK = (Known WS switch) Used WS=$(WS)
  PREPROC := $(PREPROC) -D$(WS)
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
OPTSC := $(OPTSC) $(PREPROC)
OPTSL := $(OPTSL) $(PREPROC)

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
            \\033[1;31m SYSless       \\033[0m\\033[1m $(SYSlessCHK)\\033[0m \n\
            \\033[1;31m Nullify i     \\033[0m\\033[1m $(NULiCHK)\\033[0m \n\
            \\033[1;31m Nullify j     \\033[0m\\033[1m $(NULjCHK)\\033[0m \n\
            \\033[1;31m Nullify k     \\033[0m\\033[1m $(NULkCHK)\\033[0m \n\
            \\033[1;31m PP Limiter    \\033[0m\\033[1m $(PPLCHK)\\033[0m \n\
            \\033[1;31m LMA Limiter   \\033[0m\\033[1m $(LMACHK)\\033[0m \n\
            \\033[1;31m WENO          \\033[0m\\033[1m $(WENOCHK)\\033[0m \n\
            \\033[1;31m RECV type     \\033[0m\\033[1m $(RECVCHK)\\033[0m \n\
            \\033[1;31m Riemann Solver\\033[0m\\033[1m $(RSUCHK)\\033[0m \n\
            \\033[1;31m WS Algorithm  \\033[0m\\033[1m $(WSCHK)\\033[0m \n\
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
	@cp -r examples OFF/
	@cp -rL inputs-template lib util src makefile README.md .doxygenconfig .gitignore OFF/
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
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@ 1>> diagnostic_messages 2>> error_messages
EXES := $(EXES) IBM

$(DEXE)OFF : PRINTINFO $(MKDIRS) $(DOBJ)off.o
	@rm -f $(filter-out $(DOBJ)off.o,$(EXESOBJ))
	@echo | tee -a make.log
	@echo $(LITEXT) | tee -a make.log
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@ 1>> diagnostic_messages 2>> error_messages
EXES := $(EXES) OFF

$(DEXE)POG : PRINTINFO $(MKDIRS) $(DOBJ)pog.o
	@rm -f $(filter-out $(DOBJ)pog.o,$(EXESOBJ))
	@echo | tee -a make.log
	@echo $(LITEXT) | tee -a make.log
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@ 1>> diagnostic_messages 2>> error_messages
EXES := $(EXES) POG

#$(DOBJ)data_type_amrblock.o : Data_Type_AMRBlock.f90 \
#  $(DOBJ)ir_precision.o \
#  $(DOBJ)data_type_cell.o \
#  $(DOBJ)data_type_global.o \
#  $(DOBJ)data_type_hashid.o \
#  $(DOBJ)data_type_hashtcell.o \
#  $(DOBJ)data_type_hashtface.o \
#  $(DOBJ)data_type_hashtnode.o \
#  $(DOBJ)data_type_vector.o \
#  $(DOBJ)lib_io_misc.o
#  @echo $(COTEXT) | tee -a make.log
#  @$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_bc.o : Data_Type_BC.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_cell.o : Data_Type_Cell.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_conservative.o : Data_Type_Conservative.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_face.o : Data_Type_Face.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_global.o : Data_Type_Global.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#$(DOBJ)data_type_hashid.o : Data_Type_HashID.f90 \
#  $(DOBJ)ir_precision.o \
#  $(DOBJ)lib_morton.o
#  @echo $(COTEXT) | tee -a make.log
#  @$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#$(DOBJ)data_type_hashtcell.o : Data_Type_HashTCell.f90 \
#  $(DOBJ)ir_precision.o \
#  $(DOBJ)data_type_cell.o \
#  $(DOBJ)data_type_hashid.o
#  @echo $(COTEXT) | tee -a make.log
#  @$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#$(DOBJ)data_type_hashtface.o : Data_Type_HashTFace.f90 \
#  $(DOBJ)ir_precision.o \
#  $(DOBJ)data_type_face.o \
#  $(DOBJ)data_type_hashid.o
#  @echo $(COTEXT) | tee -a make.log
#  @$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#$(DOBJ)data_type_hashtnode.o : Data_Type_HashTNode.f90 \
#  $(DOBJ)ir_precision.o \
#  $(DOBJ)data_type_hashid.o \
#  $(DOBJ)data_type_vector.o
#  @echo $(COTEXT) | tee -a make.log
#  @$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_os.o : Data_Type_OS.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_primitive.o : Data_Type_Primitive.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_probe.o : Data_Type_Probe.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_sblock.o : Data_Type_SBlock.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_face.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#$(DOBJ)data_type_sl_list.o : Data_Type_SL_List.f90 \
#  $(DOBJ)ir_precision.o
#  @echo $(COTEXT) | tee -a make.log
#  @$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_tensor.o : Data_Type_Tensor.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_time.o : Data_Type_Time.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_vector.o : Data_Type_Vector.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)ibm.o : IBM.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)ir_precision.o : IR_Precision.f90
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_base64.o : Lib_Base64.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_fluidynamic.o : Lib_Fluidynamic.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_face.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_fluxes_convective.o \
	$(DOBJ)lib_fluxes_diffusive.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_math.o \
	$(DOBJ)lib_parallel.o \
	$(DOBJ)lib_profiling.o \
	$(DOBJ)lib_runge_kutta.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o \
	$(DOBJ)lib_parallel.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_fluxes_convective.o : Lib_Fluxes_Convective.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_face.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_riemann.o \
	$(DOBJ)lib_weno.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_fluxes_diffusive.o : Lib_Fluxes_Diffusive.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_tensor.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_io_misc.o : Lib_IO_Misc.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_os.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#$(DOBJ)lib_math.o : Lib_Math.f90 \
#  $(DOBJ)ir_precision.o \
#  $(DOBJ)data_type_sl_list.o \
#  $(DOBJ)data_type_vector.o
#  @echo $(COTEXT) | tee -a make.log
#  @$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_math.o : Lib_Math.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#$(DOBJ)lib_morton.o : Lib_Morton.f90 \
#  $(DOBJ)ir_precision.o
#  @echo $(COTEXT) | tee -a make.log
#  @$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_multigrid.o : Lib_Multigrid.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_fluidynamic.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_parallel.o : Lib_Parallel.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_postprocessing.o : Lib_PostProcessing.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_vtk_io.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_profiling.o : Lib_Profiling.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann.o : Lib_Riemann.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_runge_kutta.o : Lib_Runge_Kutta.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_conservative.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_thermodynamic_laws_ideal.o : Lib_Thermodynamic_Laws_Ideal.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_vtk_io.o : Lib_VTK_IO.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)lib_base64.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_weno.o : Lib_WENO.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_sblock.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)off.o : OFF.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_probe.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_tensor.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)lib_fluidynamic.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_math.o \
	$(DOBJ)lib_profiling.o \
	$(DOBJ)lib_runge_kutta.o \
	$(DOBJ)lib_weno.o \
	$(DOBJ)lib_parallel.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)pog.o : POG.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_global.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_sblock.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_postprocessing.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)%.o : %.f90
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#-----------------------------------------------------------------------------------------------------------------------------------
