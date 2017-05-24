FXDR_INC_DIR=/home/jjl/util-library/fortran/bfield_library_jdl/include/
FXDR_DIR=/home/jjl/util-library/fortran/bfield_library_jdl/lib/

# Location of m3dc1 libraries
LIBM3DC1_DIR=$(FIO_INSTALL_DIR)/lib
LM3D = -lfusionio -lm3dc1 -lz -lhdf5 -lstdc++ -llapack

FF = -O3

xdr_to_ascii: dump_xdr_to_ascii.f90 $(FXDR_DIR)libbjdl.a
	gfortran -o $@ $(FF) $^ -I${FXDR_INC_DIR}  -L${FXDR_DIR} -lbjdl -L$(LIBM3DC1_DIR) $(LM3D)

clean:
	-rm *.o *.mod xdr_to_ascii
	@echo "Cleaned object and mod files from "$(PWD)"."
