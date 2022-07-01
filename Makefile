CC = gcc
CFLAGS = -Wall -Werror
LIBS = -L. -llise -lm
AR = ar
ARFLAGS = rcs
LOBJECTS = ipro.o iproatom.o ilig.o iligatom.o iprostats.o stats.o grid.o results.o

lise: main.o liblise.a
	$(CC) main.o $(LIBS) -o $@

main.o: main.c src/lise.h
	$(CC) $(CFLAGS) -c $*.c $(LIBS)

lib: liblise.a

liblise.a: $(LOBJECTS)
	$(AR) $(ARFLAGS) $@ $?

$(LOBJECTS): src/lise.h
	$(CC) $(CFLAGS) -c ./src/$*.c

remove: clean
	rm -f liblise.a lise

clean:
	rm -f ./*.o
