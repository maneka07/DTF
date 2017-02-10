DTF_DIR=../../libdtf
PNETCDF_DIR=$(HOME)/apps/pnetcdf
CFLAGS=-I$(DTF_DIR) -I$(PNETCDF_DIR)/include
LDFLAGS=-L$(PNETCDF_DIR)/lib -lpnetcdf -L$(DTF_DIR) -ldtf
CC=mpicc

all: reader writer

reader.o: reader-1.c 
	$(CC) -c -o $@ $< $(CFLAGS)

writer.o: writer-1.c 
	$(CC) -c -o $@ $< $(CFLAGS)

reader: reader.o
	$(CC) -o $@ $^ $(LDFLAGS)

writer: writer.o
	$(CC) -o $@ $^ $(LDFLAGS)	

clean:
	rm -f *.o reader writer
