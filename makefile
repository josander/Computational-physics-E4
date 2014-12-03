
CC = gcc
CFLAGS = -O3
LIBS = -lm -lgsl -lgslcblas


HEADERS = E4_func.h
OBJECTS = E4_func.o E4.o 
PROGRAM = E4

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c


