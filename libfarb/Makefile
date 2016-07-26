BIN=libfarb.so
SOURCES = $(wildcard *.c)
OBJECTS = $(SOURCES:.c=.o)

CFLAGS=-Wall -g -Wno-comment
MPICC=

all: $(BIN)

libfarb.so: $(OBJECTS)
	$(MPICC) -shared -o $@ $(OBJECTS)

%.o: %.c
	$(MPICC) -c -fPIC $< $(CFLAGS)

clean:
	rm -f $(OBJECTS) $(BIN)
