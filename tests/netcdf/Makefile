FARB_DIR=../comm_lib/libfarb
NETCDF_DIR=
CC=
CFLAGS=-I$(FARB_DIR) -I$(NETCDF_DIR)/include
LDFLAGS=-L$(FARB_DIR) -lfarb -L$(NETCDF_DIR)/lib -lnetcdf


all: reader writer

simple_xy_rd.o: simple_xy_rd.c 
	$(CC) -c -o $@ $< $(CFLAGS)

simple_xy_wr.o: simple_xy_wr.c 
	$(CC) -c -o $@ $< $(CFLAGS)

reader: simple_xy_rd.o
	$(CC) -o $@ $^ $(LDFLAGS)

writer: simple_xy_wr.o
	$(CC) -o $@ $^ $(LDFLAGS)	

clean:
	rm -f *.o reader writer
