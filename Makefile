DIR_MAIN = ./
DIR_INCLUDE = $(DIR_MAIN)amplitude/
DIR_SRC = $(DIR_MAIN)src/

FLAGS = --f90flags=--free-line-length-512 --opt=-O3

OLIB=amp.cpython-38-x86_64-linux-gnu.so
SRC=$(DIR_SRC)amp.f95 
SKIP=skip: kin_3p kin_4p kin_mt4p 


all: $(OLIB)

$(OLIB): $(SRC)
	f2py -m amp -c $^ $(FLAGS) $(SKIP)

.PHONY: clean
clean:
	rm -fr $(DIR_MAIN)*.o $(DIR_MAIN)*.so $(DIR_MAIN)*.mod
