# poisson2d Makefile
# Builds: libpoisson2d.a, test_GRF_disk_r64, test_GRF_disk_r128, poisson2d_mex.<ext>
#
# Targets:
#   make         -- build core + mex + run test64 + run test128
#   make core    -- build libpoisson2d_core.a (no mex wrapper, for tests)
#   make lib     -- build libpoisson2d.a (includes mex wrapper)
#   make test64  -- build + run double-precision disk test
#   make test128 -- build + run quad-precision disk test
#   make mex     -- build MATLAB MEX (requires mwrap)
#   make clean   -- remove build artifacts

SHELL := /bin/bash

ROOT       := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ROOT       := $(patsubst %/,%,$(ROOT))
SRC_DIR    := $(ROOT)/src
TEST_DIR   := $(ROOT)/test
MATLAB_DIR := $(ROOT)/matlab
BLD_DIR    := $(ROOT)/build

FC := gfortran-14
CC := gcc-14
CXX := g++-14
MW := ~/mwrap/mwrap

GMSH_SRC_DIR  := $(ROOT)/external/gmsh
GMSH_BLD_DIR  := $(ROOT)/external/gmsh_build
GMSH_INC      := $(ROOT)/external/gmsh/api
GMSH_LIB      := $(ROOT)/external/libgmsh.a
HB_PREFIX     := /opt/homebrew

GMSH_CMAKE_OPTS := \
  -DCMAKE_C_COMPILER=gcc-14 \
  -DCMAKE_CXX_COMPILER=g++-14 \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DENABLE_OPENMP=OFF \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_BUILD_LIB=ON \
  -DENABLE_BUILD_SHARED=OFF \
  -DENABLE_BUILD_DYNAMIC=OFF \
  -DENABLE_FLTK=OFF \
  -DENABLE_GRAPHICS=OFF \
  -DENABLE_POST=ON \
  -DENABLE_PLUGINS=ON \
  -DENABLE_OCC=OFF \
  -DENABLE_PARSER=ON \
  -DENABLE_MESH=ON \
  -DENABLE_CGNS=OFF \
  -DENABLE_MED=OFF \
  -DENABLE_HXT=OFF

MWFLAGS := -c99complex -i8 -mex

UNAME := $(shell uname)
ARCH  := $(shell uname -m)

ifeq ($(UNAME), Darwin)
  MATLAB_ROOT := $(shell ls -d /Applications/MATLAB_R*.app 2>/dev/null | sort | tail -n1)
  ifeq ($(ARCH), arm64)
    MATLAB_ARCH := maca64
    MEX_EXT     := mexmaca64
  else
    MATLAB_ARCH := maci64
    MEX_EXT     := mexmaci64
  endif
  MATLAB_INC  := -I$(MATLAB_ROOT)/extern/include
  MATLAB_LIBS := $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmx.dylib \
                 $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmex.dylib \
                 $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmat.dylib -lm
  MEX_LDFLAGS := -bundle -Wl,-undefined,dynamic_lookup
else
  MATLAB_ROOT := $(shell ls -d /usr/local/MATLAB/R* 2>/dev/null | sort | tail -n1)
  MATLAB_ARCH := glnxa64
  MEX_EXT     := mexa64
  MATLAB_INC  := -I$(MATLAB_ROOT)/extern/include
  MATLAB_LIBS := -L$(MATLAB_ROOT)/bin/$(MATLAB_ARCH) -lmx -lmex -lmat -lm
  MEX_LDFLAGS := -shared
endif

FFLAGS := -g -O2 -fPIC \
           -fdefault-integer-8 \
           -frecursive \
           -std=legacy -w \
           -fallow-argument-mismatch \
           -J$(BLD_DIR) -I$(BLD_DIR) -I$(SRC_DIR)

# Core sources (compile order matters: dependencies first)
CORE_SOURCES := $(SRC_DIR)/utils_mod.f90 \
                $(SRC_DIR)/geometry_mod.f90 \
                $(SRC_DIR)/lap2d_mod.f90 \
                $(SRC_DIR)/specialquad_mod.f90 \
                $(SRC_DIR)/triquad_mod.f90
CORE_OBJECTS := $(patsubst $(SRC_DIR)/%.f90, $(BLD_DIR)/%.o, $(CORE_SOURCES))
CORE_LIB     := $(BLD_DIR)/libpoisson2d_core.a

# Full library (core + mex wrapper)
LIB_SOURCES  := $(CORE_SOURCES) $(SRC_DIR)/poisson2d_mex.f90
LIB_OBJECTS  := $(patsubst $(SRC_DIR)/%.f90, $(BLD_DIR)/%.o, $(LIB_SOURCES))
LIB          := $(BLD_DIR)/libpoisson2d.a

TEST64_SRC   := $(TEST_DIR)/test_GRF_disk_r64.f90
TEST64_BIN   := $(BLD_DIR)/test_GRF_disk_r64
TEST128_SRC  := $(TEST_DIR)/test_GRF_disk_r128.f90
TEST128_BIN  := $(BLD_DIR)/test_GRF_disk_r128

MW_SRC  := $(MATLAB_DIR)/poisson2d.mw
MEX_C   := $(MATLAB_DIR)/poisson2d_mex.c
MEX_OUT := $(MATLAB_DIR)/poisson2d_mex.$(MEX_EXT)

.PHONY: all core lib test64 test128 gmsh mex clean clean_gmsh

all: core mex test64 test128

core: $(CORE_LIB)

lib: $(LIB)

test64: $(TEST64_BIN)
	$(TEST64_BIN)

test128: $(TEST128_BIN)
	$(TEST128_BIN)

mex: $(MEX_OUT)

gmsh: $(GMSH_LIB)

$(GMSH_LIB):
	mkdir -p $(GMSH_BLD_DIR)
	cd $(GMSH_BLD_DIR) && cmake $(GMSH_CMAKE_OPTS) $(GMSH_SRC_DIR)
	cd $(GMSH_BLD_DIR) && $(MAKE) -j8 lib
	@if [ -f $(GMSH_BLD_DIR)/libgmsh.a ]; then \
		cp $(GMSH_BLD_DIR)/libgmsh.a $(GMSH_LIB); \
		echo "Success: $(GMSH_LIB) created."; \
	else \
		echo "Error: libgmsh.a not found in $(GMSH_BLD_DIR)"; exit 1; \
	fi

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BLD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

GMSH_WRAP_OBJ := $(BLD_DIR)/gmsh_wrap.o

$(GMSH_WRAP_OBJ): $(SRC_DIR)/gmsh_wrap.cpp | $(BLD_DIR)
	$(CXX) -O2 -std=gnu++17 -fPIC -I$(GMSH_INC) -c $< -o $@

$(CORE_LIB): $(CORE_OBJECTS)
	ar rcs $@ $^

$(LIB): $(LIB_OBJECTS)
	ar rcs $@ $^

$(TEST64_BIN): $(CORE_LIB) $(TEST64_SRC) | $(BLD_DIR)
	$(FC) $(FFLAGS) $(TEST64_SRC) -L$(BLD_DIR) -lpoisson2d_core \
	  -lgfortran -lm -o $(TEST64_BIN)

$(TEST128_BIN): $(CORE_LIB) $(TEST128_SRC) | $(BLD_DIR)
	$(FC) $(FFLAGS) $(TEST128_SRC) -L$(BLD_DIR) -lpoisson2d_core \
	  -lgfortran -lm -o $(TEST128_BIN)

$(MEX_C): $(MW_SRC) | $(BLD_DIR)
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) poisson2d_mex -mb -list poisson2d.mw
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) poisson2d_mex -c poisson2d_mex.c poisson2d.mw
	perl -pi -e 's/_{2,}/_/g' $(MEX_C)

MEX_C_OBJ := $(BLD_DIR)/poisson2d_mex_c.o

$(MEX_C_OBJ): $(MEX_C) | $(BLD_DIR)
	$(CC) -c -fPIC \
	  -DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2018a -DMX_COMPAT_32=0 \
	  $(MATLAB_INC) -include $(SRC_DIR)/gmsh_wrap.h \
	  $(MEX_C) -o $(MEX_C_OBJ)

$(MEX_OUT): $(LIB) $(MEX_C_OBJ) $(GMSH_WRAP_OBJ) $(GMSH_LIB)
	$(CXX) $(MEX_LDFLAGS) -fPIC \
	  $(MEX_C_OBJ) $(GMSH_WRAP_OBJ) \
	  -L$(BLD_DIR) -lpoisson2d \
	  $(GMSH_LIB) \
	  -L$(HB_PREFIX)/lib -lgmp -lfreetype \
	  $(MATLAB_LIBS) \
	  -lgfortran -lm \
	  -o $(MEX_OUT)

clean:
	rm -rf $(BLD_DIR) $(MEX_OUT)

clean_gmsh:
	rm -rf $(GMSH_BLD_DIR) $(GMSH_LIB)
