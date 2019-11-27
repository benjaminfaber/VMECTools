SRC_DIR := src
OBJS_DIR := objs
LIB_DIR := lib
$(shell mkdir -p $(OBJS_DIR))
$(shell mkdir -p $(LIB_DIR))

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
DOCKER_ENV_FILE := $(shell find / -maxdepth 1 -name '.dockerenv')
ifneq ($(DOCKER_ENV_FILE),/.dockerenv)
NETCDF_F_INC_DIR = /opt/gcc/openmpi/netcdf-fortran-4.4.5/include
NETCDF_F_LIB_DIR = /opt/gcc/openmpi/netcdf-fortran-4.4.5/include
NETCDF_C_INC_DIR = /opt/gcc/netcdf-c-4.7.0/lib
NETCDF_C_LIB_DIR = /opt/gcc/netcdf-c-4.7.0/lib
BLAS_LIB_DIR = /opt/gcc/openblas-0.3.7/lib
else
NETCDF_F_INC_DIR = /usr/include
NETCDF_F_LIB_DIR = /usr/lib
NETCDF_C_INC_DIR = /usr/include
NETCDF_C_LIB_DIR = /usr/lib
BLAS_LIB_DIR = /usr/lib
endif

FC := mpifort
CXX := mpicxx
ifneq ($(DOCKER_ENV_FILE),/.dockerenv)
MARCH_FLAG = native
else
MARCH_FLAG = x86_64
endif
COMPILE_FLAGS := -O3 -march=$(MARCH_FLAG) -fopenmp
FC_FLAGS := $(COMPILE_FLAGS) -I$(NETCDF_F_INC_DIR) -ffree-line-length-none
CXX_FLAGS := $(COMPILE_FLAGS)
F_LD_FLAGS := $(COMPILE_FLAGS) -L$(NETCDF_F_LIB_DIR) -L$(NETCDF_C_LIB_DIR) -L$(BLAS_LIB_DIR) -lnetcdff -lnetcdf -lopenblas
CXX_LD_FLAGS := $(F_LD_FLAGS) -lgfortran
# End of system-dependent variable assignments

# Promote real to double, as gs2 does:
FC_FLAGS += -fdefault-real-8 -fdefault-double-8

USE_MINI_LIBSTELL = 0
# The variable LIBSTELL_DIR should either be "mini_libstell", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a directory containing libstell .mod files
# elsewhere on your system.
MINI_LIBSTELL_DIR = mini_libstell

# The variable LIBSTELL_FOR_SFINCS should either be "mini_libstell/mini_libstell.a", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a libstell.a library elsewhere on your system.
MINI_LIBSTELL = $(addprefix $(LIB_DIR)/,mini_libstell.a)

# If libstell.a does not exist, set it to the $(MINI_LIBSTELL) options
LIBSTELL_DIR = /home/bfaber/projects/stellopt/src/LIBSTELL/Release

ifeq ($(USE_MINI_LIBSTELL),0)
ifeq ($(wildcard $(LIBSTELL_DIR)/*.a),$(addprefix $(LIBSTELL_DIR)/,libstell.a))
LIBSTELL = $(addprefix $(LIBSTELL_DIR)/,libstell.a)
LIBSTELL_INC = $(LIBSTELL_DIR)
else
LIBSTELL = $(MINI_LIBSTELL)
LIBSTELL_INC = $(MINI_LIBSTELL_DIR)
endif
else
LIBSTELL = $(MINI_LIBSTELL)
LIBSTELL_INC = $(MINI_LIBSTELL_DIR)
endif


VMEC2PEST := vmec2pest
C_EXEC := c_exec
V2PLIB := $(addprefix $(LIB_DIR)/,libvmectools.a)
all: v2p
v2p: $(LIBSTELL) $(VMEC2PEST)
lib: $(V2PLIB)
ctest: $(V2PLIB) $(C_EXEC)
mini_libstell: $(MINI_LIBSTELL)

export

include makefile.depend

$(OBJS_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FC_FLAGS) $(FC_INCS) -I $(LIBSTELL_INC) -c -o $@ $<

$(OBJS_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FC_FLAGS) $(FC_INCS) -I $(LIBSTELL_INC) -c -o $@ $<

$(OBJS_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXX_FLAGS) $(CXX_INCS) -c -o $@ $<

$(VMEC2PEST): $(OBJS_LINK_F90) $(OBJS_LINK_F)
	$(FC) -o $@ $^ $(LIBSTELL) $(F_LD_FLAGS)

$(C_EXEC): $(OBJS_LINK_CXX)
	$(CXX) -o $@ $^ $(V2PLIB) $(LIBSTELL) $(CXX_LD_FLAGS)

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
