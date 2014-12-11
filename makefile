
CC = gcc
CFLAGS = -O3
LIBS = -lm -lgsl -lgslcblas


HEADERS = fft_func.h E4_func.h
OBJECTS = fft_func.o E4_func.o E4.o 
PROGRAM = E4

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c


