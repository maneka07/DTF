#TARGET = client server da.exe sim.exe
TARGET = client server

#CC=mpifccpx
CC=mpicc
FC=mpif90
CFLAGS=-O0 -Wall -g -fopenmp -I../../libfarb
O = bin

LIBCOMM_LDFLAGS=-L../../libfarb -Wl,-rpath=../../libfarb -lfarb -ldl

LD_FLAGS=$(LD_FLAGS) $(LIBCOMM_LDFLAGS)
all: $(TARGET)

client:	client.o
	$(CC) -o client client.o $(LIBCOMM_LDFLAGS)
server: server.o
	$(CC) -o server server.o $(LIBCOMM_LDFLAGS)

sim.exe: sim.o
	$(FC) -o sim.exe sim.o $(LIBCOMM_LDFLAGS)

da.exe: da.o
	$(FC) -o da.exe da.o  $(LIBCOMM_LDFLAGS)

da.o: da.f90
	$(FC) -c da.f90 -o da.o $(CFLAGS)
sim.o: sim.f90
	$(FC) -c sim.f90 -o sim.o $(CFLAGS)

clean:
	rm -f $(TARGET) *.o
