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
COMPILER = gnu
DEBUG    = yes
F03STD   = yes
OPTIMIZE = no
OPENMP   = no
MPI      = no
R16P     = no
NULi     = no
NULj     = no
NULk     = no
WENO     = WENO
RECV     = RECVC
RSU      = HLLCp
VACUUM   = no
WS       = WSup
SMSW     = SMSWz
HYBRID   = NOHYBRID
PPL      = no
TECIO    = no

.PHONY : DEFAULTRULE
DEFAULTRULE: $(DEXE)OFF

.PHONY : help
help:
	@echo
	@echo -e '\033[1;31m Make options of OFF codes\033[0m'
	@echo
	@echo -e '\033[1;31m Compiler choice\033[0m'
	@echo -e '\033[1;31m  COMPILER=gnu  \033[0m\033[1m => GNU gfortran          \033[0m'
	@echo -e '\033[1;31m  COMPILER=intel\033[0m\033[1m => Intel Fortran         \033[0m'
	@echo -e '\033[1;31m  COMPILER=pgi  \033[0m\033[1m => Portland Group Fortran\033[0m'
	@echo -e '\033[1;31m  COMPILER=g95  \033[0m\033[1m => free g95              \033[0m'
	@echo -e '\033[1;31m  COMPILER=$(COMPILER)  \033[0m\033[1m => default         \033[0m'
	@echo
	@echo -e '\033[1;31m Compiling options\033[0m'
	@echo -e '\033[1;31m  DEBUG=yes(no)   \033[0m\033[1m => on(off) debug                  (default $(DEBUG))\033[0m'
	@echo -e '\033[1;31m  F03STD=yes(no)  \033[0m\033[1m => on(off) check standard fortran (default $(F03STD))\033[0m'
	@echo -e '\033[1;31m  OPTIMIZE=yes(no)\033[0m\033[1m => on(off) optimization           (default $(OPTIMIZE))\033[0m'
	@echo -e '\033[1;31m  OPENMP=yes(no)  \033[0m\033[1m => on(off) OpenMP directives      (default $(OPENMP))\033[0m'
	@echo -e '\033[1;31m  MPI=yes(no)     \033[0m\033[1m => on(off) MPI    directives      (default $(MPI)) \033[0m'
	@echo
	@echo -e '\033[1;31m Preprocessing options\033[0m'
	@echo -e '\033[1;31m  R16P=yes(no)\033[0m\033[1m => on(off) definition of real with "128 bit" (default $(R16P))\033[0m'
	@echo -e '\033[1;31m  NULi=yes(no)\033[0m\033[1m => on(off) nullify i direction (1D or 2D)    (default $(NULi))\033[0m'
	@echo -e '\033[1;31m  NULj=yes(no)\033[0m\033[1m => on(off) nullify j direction (1D or 2D)    (default $(NULj))\033[0m'
	@echo -e '\033[1;31m  NULk=yes(no)\033[0m\033[1m => on(off) nullify k direction (1D or 2D)    (default $(NULk))\033[0m'
	@echo -e '\033[1;31m  PPL=yes(no) \033[0m\033[1m => on(off) Positivity Preserving Limiter     (default $(PPL))\033[0m'
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
ifeq "$(COMPILER)" "gnu"
  OPTSC   = -cpp -c -J$(DMOD)
  OPTSL   =
  # debug
  ifeq "$(DEBUG)" "yes"
    PREPROC := $(PREPROC) -DDEBUG
    OPTSC := $(OPTSC) -O0 -Wall -Warray-bounds -fcheck=all -fbacktrace -ffpe-trap=invalid,overflow,underflow,precision,denormal
    OPTSL := $(OPTSL) -O0 -Wall -Warray-bounds -fcheck=all -fbacktrace -ffpe-trap=invalid,overflow,underflow,precision,denormal
#-Warray-temporaries
  endif
  # standard
  ifeq "$(F03STD)" "yes"
    OPTSC := $(OPTSC) -std=f2008 -fall-intrinsics
    OPTSL := $(OPTSL) -std=f2008 -fall-intrinsics
  endif
  # optimization
  ifeq "$(OPTIMIZE)" "yes"
    OPTSC := $(OPTSC) -O3
    OPTSL := $(OPTSL) -O3
  endif
  # openmp
  ifeq "$(OPENMP)" "yes"
    OPTSC := $(OPTSC) -fopenmp
    OPTSL := $(OPTSL) -fopenmp
    PREPROC := $(PREPROC) -DOPENMP
  endif
  # mpi
  ifeq "$(MPI)" "yes"
    PREPROC := $(PREPROC) -DMPI2
    FC = mpif90
  else
    FC = gfortran
  endif
endif
ifeq "$(COMPILER)" "intel"
  OPTSC   = -cpp -c -module $(DMOD)
  OPTSL   =
  # debug
  ifeq "$(DEBUG)" "yes"
    PREPROC := $(PREPROC) -DDEBUG
    CHK = -check all -check noarg_temp_created
    DEB = -debug all
    WRN = -warn all
    OPTSC := $(OPTSC) -O0 -fpe-all=0 -fp-stack-check -traceback $(WRN) $(CHK) $(DEB)
    OPTSL := $(OPTSL) -O0 -fpe-all=0 -fp-stack-check -traceback $(WRN) $(CHK) $(DEB)
  endif
  # standard
  ifeq "$(F03STD)" "yes"
    OPTSC := $(OPTSC) -std03
    OPTSL := $(OPTSL) -std03
  endif
  # optimization
  ifeq "$(OPTIMIZE)" "yes"
    OPTSC := $(OPTSC) -O3 -ipo
    OPTSL := $(OPTSL) -O3 -ipo
  endif
  # openmp
  ifeq "$(OPENMP)" "yes"
    OPTSC := $(OPTSC) -openmp
    OPTSL := $(OPTSL) -openmp
    PREPROC := $(PREPROC) -DOPENMP
  endif
  # mpi
  ifeq "$(MPI)" "yes"
    PREPROC := $(PREPROC) -DMPI2
    FC = mpif90
  else
    FC = ifort
  endif
endif
ifeq "$(COMPILER)" "pgi"
  OPTSC   = -Mpreprocess -c -module $(DMOD)
  OPTSL   =
	PREPROC := $(PREPROC) -Dpgf95
  # debug
  ifeq "$(DEBUG)" "yes"
    PREPROC := $(PREPROC) -DDEBUG
    OPTSC := $(OPTSC) -C -g -Mbounds -Mchkstk
    OPTSL := $(OPTSL) -C -g -Mbounds -Mchkstk
  endif
  # standard
  ifeq "$(F03STD)" "yes"
    OPTSC := $(OPTSC) -Mstandard
    OPTSL := $(OPTSL) -Mstandard
  endif
  # optimization
  ifeq "$(OPTIMIZE)" "yes"
    OPTSC := $(OPTSC) -O1
    OPTSL := $(OPTSL) -O1
  endif
  # openmp
  ifeq "$(OPENMP)" "yes"
    OPTSC := $(OPTSC) -mp
    OPTSL := $(OPTSL) -mp
    PREPROC := $(PREPROC) -DOPENMP
  endif
  # mpi
  ifeq "$(MPI)" "yes"
    PREPROC := $(PREPROC) -DMPI2
    FC = mpif90
  else
    FC = pgf95
  endif
endif
ifeq "$(COMPILER)" "g95"
  OPTSC   = -cpp -c -fmod=$(DMOD)
  OPTSL   =
  # debug
  ifeq "$(DEBUG)" "yes"
    PREPROC := $(PREPROC) -DDEBUG
    OPTSC := $(OPTSC) -O0 -g
    OPTSL := $(OPTSL) -O0 -g
  endif
  # standard
  ifeq "$(F03STD)" "yes"
    OPTSC := $(OPTSC) -std=f2003 -fintrinsic-extensions
    OPTSL := $(OPTSL) -std=f2003 -fintrinsic-extensions
  endif
  # optimization
  ifeq "$(OPTIMIZE)" "yes"
    OPTSC := $(OPTSC) -O3
    OPTSL := $(OPTSL) -O3
  endif
  FC = g95
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
            \\033[1;31m Optimize      \\033[0m\\033[1m $(OPTIMIZE)\\033[0m \n\
            \\033[1;31m OpenMP        \\033[0m\\033[1m $(OPENMP)\\033[0m \n\
            \\033[1;31m MPI           \\033[0m\\033[1m $(MPI)\\033[0m \n\
            \\033[1;31m R16P          \\033[0m\\033[1m $(R16PCHK)\\033[0m \n\
            \\033[1;31m Nullify i     \\033[0m\\033[1m $(NULiCHK)\\033[0m \n\
            \\033[1;31m Nullify j     \\033[0m\\033[1m $(NULjCHK)\\033[0m \n\
            \\033[1;31m Nullify k     \\033[0m\\033[1m $(NULkCHK)\\033[0m \n\
            \\033[1;31m PP Limiter    \\033[0m\\033[1m $(PPLCHK)\\033[0m \n\
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
	@echo -e '\033[1;31m Compiling options\033[0m' | tee -a make.log
	@echo -e '\033[1m [$(OPTSC)]\033[0m' | tee -a make.log
	@echo | tee -a make.log
	@echo -e '\033[1;31m Linking options \033[0m' | tee -a make.log
	@echo -e '\033[1m [$(OPTSL)]\033[0m' | tee -a make.log
	@echo | tee -a make.log

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
	@mkdir -p OFF
	@cp -rL input lib util src makefile OFF/
	@tar czf OFF.tgz OFF
	@rm -rf OFF

.PHONY : doc
doc:
	@echo -e "\033[1;31m Creating documentation\033[0m" | tee make.log
	@mkdir -p doc
	@doxygen .doxygenconfig

.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
#----------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------
# rules of linking and compiling
COTEXT  = -e "\033[1;31m Compiling\033[0m\033[1m $(<F)\033[0m"
LITEXT  = -e "\033[1;31m Assembling\033[0m\033[1m $@\033[0m"
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

$(DEXE)ICG : PRINTINFO $(MKDIRS) $(DOBJ)icg.o
	@rm -f $(filter-out $(DOBJ)icg.o,$(EXESOBJ))
	@echo | tee -a make.log
	@echo $(LITEXT) | tee -a make.log
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@ 1>> diagnostic_messages 2>> error_messages
EXES := $(EXES) ICG

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

$(DOBJ)data_type_bc.o : Data_Type_BC.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_cell.o : Data_Type_Cell.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_conservative.o : Data_Type_Conservative.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)data_type_globals.o : Data_Type_Globals.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

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

$(DOBJ)icg.o : ICG.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_cell.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_globals.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)ir_precision.o : IR_Precision.f90
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_fluidynamic.o : Lib_Fluidynamic.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_conservative.o \
	$(DOBJ)data_type_globals.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_math.o \
	$(DOBJ)lib_parallel.o \
	$(DOBJ)lib_riemann.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o \
	$(DOBJ)lib_weno.o \
	$(DOBJ)lib_parallel.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_io_misc.o : Lib_IO_Misc.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_os.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_math.o : Lib_Math.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_mesh.o : Lib_Mesh.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_globals.o \
	$(DOBJ)data_type_vector.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

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
	$(DOBJ)data_type_globals.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)lib_io_misc.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_postprocessing.o : Lib_PostProcessing.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_globals.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_vtk_io.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_riemann.o : Lib_Riemann.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)lib_thermodynamic_laws_ideal.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_thermodynamic_laws_ideal.o : Lib_Thermodynamic_Laws_Ideal.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_vtk_io.o : Lib_VTK_IO.f90 \
	$(DOBJ)ir_precision.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)lib_weno.o : Lib_WENO.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_globals.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)off.o : OFF.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_bc.o \
	$(DOBJ)data_type_globals.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_primitive.o \
	$(DOBJ)data_type_probe.o \
	$(DOBJ)data_type_tensor.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)lib_fluidynamic.o \
	$(DOBJ)lib_math.o \
	$(DOBJ)lib_mesh.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_weno.o \
	$(DOBJ)lib_parallel.o \
	$(DOBJ)lib_parallel.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

$(DOBJ)pog.o : POG.f90 \
	$(DOBJ)ir_precision.o \
	$(DOBJ)data_type_globals.o \
	$(DOBJ)data_type_os.o \
	$(DOBJ)data_type_time.o \
	$(DOBJ)data_type_vector.o \
	$(DOBJ)lib_io_misc.o \
	$(DOBJ)lib_mesh.o \
	$(DOBJ)lib_postprocessing.o
	@echo $(COTEXT) | tee -a make.log
	@$(FC) $(OPTSC) $< -o $@ 1>> diagnostic_messages 2>> error_messages

#-----------------------------------------------------------------------------------------------------------------------------------
