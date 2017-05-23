FXDR_INC_DIR=/home/jjl/w7/xdr/fxdr/build
FXDR_DIR=/home/jjl/w7/xdr/fxdr/build

dump_it: dump_xdr_to_ascii.f90 
	gfortran44 -o $@ dump_xdr_to_ascii.f90 -I${FXDR_INC_DIR}  -L${FXDR_DIR} -lfxdr_gfortran44

clean:
	-rm *.o *.mod dump_it
	@echo "Cleaned object and mod files from "$(PWD)"."