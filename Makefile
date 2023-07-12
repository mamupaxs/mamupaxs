all: amp.cpython-38-x86_64-linux-gnu.so

# debug: debug.f95
# 	gfortran --free-line-length-512 -O3 $^ -o debug.out

amp.cpython-38-x86_64-linux-gnu.so: amp.f95 
	f2py3 -m amp -c $^ --f90flags=--free-line-length-512 --opt=-O3 skip: kin_3p kin_4p kin_mt4p 

.PHONY: clean
clean:
	rm -fr *.o *.so *.mod
