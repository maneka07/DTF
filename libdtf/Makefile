DNAME=libdtf.so
SNAME=libdtf.a
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

CFLAGS=-Wall -g -Wno-comment -pedantic
MPICC=mpicc

#K
#CFLAGS=-Wall -Xg -Wno-comment
#MPICC=mpifccpx


all: $(SNAME) $(DNAME)

$(SNAME): $(OBJ)
	ar rcs $(SNAME) $(OBJ)
	ranlib $(SNAME)

$(DNAME): $(OBJ)
	$(MPICC) -shared -o $@ $(OBJ) -lm

%.o: %.c
	$(MPICC) -c -fPIC $< $(CFLAGS)

clean:
	rm -f $(OBJ) $(SNAME) $(DNAME)
