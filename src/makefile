#==================================================================================================#
#  MAKEFILE
#==================================================================================================#
#
#--------------------------------------------------------------------------------------------------#
INFO_COMPILING=\e[1;34m[COMPILING]\e[0m
INFO_LINKING=\e[1;34m[LINKING]\e[0m
#--------------------------------------------------------------------------------------------------#
#
#--------------------------------------------------------------------------------------------------#
COMPILER=INTEL
#--------------------------------------------------------------------------------------------------#
ifeq ($(findstring GNU,$(COMPILER)),GNU)
  FC=gfortran
  INSTRUCTION=-march=native
  FCFLAGS  = -xf95-cpp-input -cpp -fdefault-real-8 -fdefault-double-8 -fbackslash
  FLFLAGS  = -fdefault-real-8 -fbackslash
  ifeq ($(findstring DEBUG,$(COMPILER)),DEBUG)
    FCFLAGS  += -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace
    FLFLAGS  += -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace
  else
    FCFLAGS  += -O3 $(INSTRUCTION)
    FLFLAGS  += -O3 $(INSTRUCTION)
  endif
endif
#--------------------------------------------------------------------------------------------------#
ifeq ($(findstring INTEL,$(COMPILER)),INTEL)
  FC=ifort
  INSTRUCTION=-O3 -xCORE-AVX2 -ip
  FCFLAGS  =-fpp -allow nofpp_comments -assume bscc -r8 -i4 -traceback -warn all
  FLFLAGS  =-r8 -i4 -traceback -assume bscc
  ifeq ($(findstring DEBUG,$(COMPILER)),DEBUG)
    FCFLAGS  +=-g -O0 -fpe0 -traceback
    FCFLAGS  +=-check all,noarg_temp_created,noformat,nooutput_conversion,pointer,bounds,uninit
    FCFLAGS  +=-init=snan -init=arrays
    FLFLAGS  +=-g -O0
  else
    FCFLAGS  +=$(INSTRUCTION)
    FLFLAGS  +=$(INSTRUCTION)
  endif
endif
#--------------------------------------------------------------------------------------------------#
#
#--------------------------------------------------------------------------------------------------#
AR       =ar
AR_FLAGS =clr
#--------------------------------------------------------------------------------------------------#
BLAS_LAPACK = -llapack -lblas
#--------------------------------------------------------------------------------------------------#
MAIN_LIB = libmain.a
#--------------------------------------------------------------------------------------------------#
FCFLAGS += -I.
FLFLAGS +=
#--------------------------------------------------------------------------------------------------#
SRC = $(shell cat sources.mk)
OBJ = $(SRC:.f90=.o)
#--------------------------------------------------------------------------------------------------#
#
#--------------------------------------------------------------------------------------------------#
all:
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "COMPILING MAIN CODE...\n"
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@$(MAKE) -s $(MAIN_LIB)
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "$(INFO_COMPILING)   $(FC) -c main.f90 -o main.o\n";
	@$(FC) $(FCFLAGS) -c main.f90 -o main.o
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "LINKING MAIN CODE...\n"
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "$(INFO_LINKING)     $(FC) main.o ../lib/$(MAIN_LIB) -o ../bin/main\n";
	@$(FC) $(FLFLAGS) $(BLAS_LAPACK) main.o ../lib/$(MAIN_LIB) -o ../bin/main
	@printf "%`tput cols`s\n"|tr ' ' '-'
	@printf "COMPILING AND LINKING MAIN CODE...DONE!\n"
	@printf "%`tput cols`s\n"|tr ' ' '-'
#--------------------------------------------------------------------------------------------------#
$(MAIN_LIB): $(OBJ)
	$(AR) $(AR_FLAGS) ../lib/$(MAIN_LIB) $(OBJ)
#--------------------------------------------------------------------------------------------------#
%.o: %.f90
	@printf "$(INFO_COMPILING)   $(FC) $<\n";
	@$(FC) $(FCFLAGS) -c $< -o $@
#--------------------------------------------------------------------------------------------------#
clean:
	rm -f $(OBJ) *.mod *.i main.o ../bin/main ../lib/$(MAIN_LIB)
#--------------------------------------------------------------------------------------------------#
#
#--------------------------------------------------------------------------------------------------#
# DEPENDENCIES
#--------------------------------------------------------------------------------------------------#
include deplist.mk
#--------------------------------------------------------------------------------------------------#
#
#--------------------------------------------------------------------------------------------------#
# Utility Targets
#--------------------------------------------------------------------------------------------------#
.PHONY: all clean 
#--------------------------------------------------------------------------------------------------#
