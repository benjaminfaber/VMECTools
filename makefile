SRC_DIR := src
OBJS_DIR := objs
LIB_DIR := lib
$(shell mkdir -p $(OBJS_DIR))

SRC_F90 := $(notdir $(shell find $(SRC_DIR) -maxdepth 1 -name '*.f90' | sed "s|^\./||"))
SRC_F := $(notdir $(shell find $(SRC_DIR) -maxdepth 1 -name '*.f' | sed "s|^\./||"))
SRC_CXX := $(notdir $(shell find $(SRC_DIR) -maxdepth 1 -name '*.cpp' | sed "s|^\./||"))

OBJS_F90 := $(subst .f90,.o,$(SRC_F90))
OBJS_F := $(subst .f,.o,$(SRC_F))
OBJS_CXX := $(subst .cpp,.o,$(SRC_CXX))
OBJS_LINK_F90 := $(addprefix $(OBJS_DIR)/,$(OBJS_F90))
OBJS_LINK_F := $(addprefix $(OBJS_DIR)/,$(OBJS_F))
OBJS_LINK_CXX := $(addprefix $(OBJS_DIR)/,$(OBJS_CXX))

FC_INCS := -I$(OBJS_DIR) -J$(OBJS_DIR) -I$(LIB_DIR)
CXX_INCS := -I$(OBJS_DIR) -I$(LIB_DIR)

# Define the NetCDF libraries here, if not already in the compiler search path
NETCDF_F_DIR := /opt/gcc/netcdf-fortran-4.4.5
NETCDF_C_DIR := /opt/gcc/netcdf-c-4.6.3
BLAS_DIR := /opt/gcc/openblas-0.3.6


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
CXX := mpicxx
FCFLAGS := -O3 -march=skylake-avx512 -fopenmp -I$(NETCDF_F_DIR)/include -ffree-line-length-none
CXXFLAGS := -O3 -march=skylake-avx512 -fopenmp
FLDFLAGS := -O3 -march=skylake-avx512 -fopenmp -L$(NETCDF_F_DIR)/lib -L$(NETCDF_C_DIR)/lib -L$(BLAS_DIR)/lib -lnetcdff -lnetcdf -lopenblas
CXXLDFLAGS := $(FLDFLAGS) -lgfortran
#endif
# End of system-dependent variable assignments

# Promote real to double, as gs2 does:
FCFLAGS += -fdefault-real-8 -fdefault-double-8

# The variable LIBSTELL_DIR should either be "mini_libstell", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a directory containing libstell .mod files
# elsewhere on your system.
MINI_LIBSTELL_DIR := mini_libstell

# The variable LIBSTELL_FOR_SFINCS should either be "mini_libstell/mini_libstell.a", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a libstell.a library elsewhere on your system.
MINI_LIBSTELL := $(addprefix $(LIB_DIR)/,mini_libstell.a)

# If libstell.a does not exist, set it to the $(MINI_LIBSTELL) options
LIBSTELL_DIR := /home/bfaber/projects/stellopt/src/LIBSTELL/Release
LIBSTELL := $(addprefix $(LIBSTELL_DIR)/,libstell.a)

VMEC2PEST := vmec2pest
C_EXEC := c_exec
V2PLIB := $(addprefix $(LIB_DIR)/,libvmec2pest.a)
all: v2p
v2p: $(VMEC2PEST)
lib: $(V2PLIB)
ctest: $(V2PLIB) $(C_EXEC)
mini_libstell: $(MINI_LIBSTELL)

export

include makefile.depend

$(OBJS_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) $(FC_INCS) -I $(LIBSTELL_DIR) -c -o $@ $<

$(OBJS_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FCFLAGS) $(FC_INCS) -I $(LIBSTELL_DIR) -c -o $@ $<

$(OBJS_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXX_INCS) -c -o $@ $<

$(VMEC2PEST): $(OBJS_LINK_F90) $(OBJS_LINK_F)
	$(FC) -o $@ $^ $(LIBSTELL) $(FLDFLAGS)

$(C_EXEC): $(OBJS_LINK_CXX)
	$(CXX) -o $@ $^ $(V2PLIB) $(LIBSTELL) $(CXXLDFLAGS)

$(V2PLIB): $(OBJS_LINK_F90) $(OBJS_LINK_F)
	ar crs $@ $^

$(MINI_LIBSTELL):
	$(MAKE) -C mini_libstell

$(LIBSTELL):

.PHONY: all allclean cleanexec libclean objclean libstellclean

allclean:
	rm -f $(OBJS_DIR)/*.o $(OBJS_DIR)/*.mod $(OBJS_DIR)/*.MOD $(LIBSTELL_DIR)/*.o $(LIBSTELL_DIR)/*.mod $(LIBSTELL_DIR)/*.MOD $(LIB_DIR)/*.a *~

cleanexec:
	rm -f $(VMEC2PEST) $(C_EXEC)

objclean:
	rm -f $(OBJS_DIR)/*.o $(OBJS_DIR)/*.mod $(OBJS_DIR)/*.MOD

libclean:
	rm -f $(LIB_DIR)/*.a

libstellclean:
	rm -f $(LIBSTELL_DIR)/*.o $(LIBSTELL_DIR)/*.mod $(LIBSTELL_DIR)/*.MOD

test_make:
	@echo HOSTNAME is $(HOST)
	@echo FC is $(FC)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
