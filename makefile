SRC_DIR := src
OBJS_DIR := objs

SRC_F90 := $(notdir $(shell find $(SRC_DIR) -maxdepth 1 -name '*.f90' | sed "s|^\./||"))
SRC_F := $(notdir $(shell find $(SRC_DIR) -maxdepth 1 -name '*.f' | sed "s|^\./||"))

$(info $(SRC_F90))
$(info $(SRC_F))

OBJS_F90 := $(subst .f90,.o,$(SRC_F90))
OBJS_F := $(subst .f,.o,$(SRC_F))
OBJS_LINK := $(addprefix $(OBJS_DIR)/,$(OBJS_F90)) $(addprefix $(OBJS_DIR)/,$(OBJS_F))

$(info $(OBJS_F90))
$(info $(OBJS_F))
$(info $(OBJS_LINK))

INCS := -I$(OBJS_DIR) -J$(OBJS_DIR)

#ifeq ($(NERSC_HOST),cori)
# makefile for NERSC Edison and Cori
# You must first load the cray-netcdf module and python module:
#   module load cray-netcdf python
# It is convenient to run
#   module unload cray-libsci
# to avoid warning messages about libsci during compiling.

#  FC := ftn
  ## NERSC documentation recommends against specifying -O2
  ## -mkl MUST APPEAR AT THE END!!
#  EXTRA_COMPILE_FLAGS += -qopenmp -mkl
#  EXTRA_LINK_FLAGS +=  -qopenmp -mkl -Wl,-ydgemm_
  # Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
#else
FC := mpifort
FCFLAGS := -O3 -march=skylake-avx512 -fopenmp -I/opt/gcc/netcdf-fortran-4.4.5/include -ffree-line-length-none
LDFLAGS := -O3 -march=skylake-avx512 -fopenmp -lnetcdff -lnetcdf -lopenblas
#endif
# End of system-dependent variable assignments

# Promote real to double, as gs2 does:
FCFLAGS += -fdefault-real-8 -fdefault-double-8

# The variable LIBSTELL_DIR should either be "mini_libstell", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a directory containing libstell .mod files
# elsewhere on your system.
LIBSTELL_DIR := mini_libstell

# The variable LIBSTELL_FOR_SFINCS should either be "mini_libstell/mini_libstell.a", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a libstell.a library elsewhere on your system.
LIBSTELL := $(LIBSTELL_DIR)/mini_libstell.a

EXEC := vmec2sfl
all: $(EXEC)
libstell: $(LIBSTELL)

export

include makefile.depend

$(OBJS_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) $(INCS) -I $(LIBSTELL_DIR) -c -o $@ $<

$(OBJS_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FCFLAGS) $(INCS) -I $(LIBSTELL_DIR) -c -o $@ $<

$(EXEC): $(OBJS_LINK)
	$(FC) -o $@ $^ $(LIBSTELL) $(LDFLAGS)

$(LIBSTELL):
	$(MAKE) -C mini_libstell

.PHONY: all clean cleanexec

clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.MOD $(LIBSTELL_DIR)/*.o $(LIBSTELL_DIR)/*.mod $(LIBSTELL_DIR)/*.MOD *~

cleanexec:
	vmec2sfl
	
#cd mini_libstell; rm -f *.o *.mod *.MOD *.a

test_make:
	@echo HOSTNAME is $(HOST)
	@echo FC is $(FC)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
