PFARB_DIR=../../libfarb
PNETCDF_DIR=$(HOME)/apps/pnetcdf
CFLAGS=-I$(PFARB_DIR) -I$(PNETCDF_DIR)/include
LDFLAGS=-L$(PNETCDF_DIR)/lib -lpnetcdf -L$(PFARB_DIR) -lpfarb
CC=mpicc

all: reader writer

reader.o: reader.c 
	$(CC) -c -o $@ $< $(CFLAGS)

writer.o: writer.c 
	$(CC) -c -o $@ $< $(CFLAGS)

reader: reader.o
	$(CC) -o $@ $^ $(LDFLAGS)

writer: writer.o
	$(CC) -o $@ $^ $(LDFLAGS)	

clean:
	rm -f *.o reader writer
