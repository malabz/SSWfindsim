CC = gcc

CFLAGS= -g -O3 -fopenmp -Wall -pipe

PREFIX = /usr/local
BINDIR = $(PREFIX)/bin

LOBJS = sswfindsim.o utils.o ../lib/ssw/src/ssw.o

PROG = sswfindsim


.PHONY: all default clean

ifdef __arm__ # (M1)
$(PROG) $(EXAMPLE): $(LOBJS)
	$(CC) -o $@ $(filter-out %.h,$^) $(CFLAGS) $(LDFLAGS) -march=armv8-a+fp+simd+crypto+crc
else # x86(Intel)
$(PROG) $(EXAMPLE): $(LOBJS)
	$(CC) -o $@ $(filter-out %.h,$^) $(CFLAGS) $(LDFLAGS)
endif

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROG)
	@echo done.

clean:
	rm -r $(LOBJS) $(PROG)
