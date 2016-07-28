DNAME=libpfarb.so
SNAME=libpfarb.a
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

CFLAGS=-Wall -g -Wno-comment
MPICC=mpicc

all: $(SNAME) $(DNAME)

$(SNAME): $(OBJ)
	ar r $(SNAME) $(OBJ)
	ranlib $(SNAME)

$(DNAME): $(OBJ)
	$(MPICC) -shared -o $@ $(OBJ)

%.o: %.c
	$(MPICC) -c -fPIC $< $(CFLAGS)

clean:
	rm -f $(OBJ) $(SNAME) $(DNAME)
