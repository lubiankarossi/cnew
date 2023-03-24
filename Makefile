CC       = g++
LDFLAGS  = 
CPPFLAGS = -mpc32 -std=c++11

CPPSOURCES = storage.c structures.c functions.c main.cpp

all: cinew

cinew:    $(CPPSOURCES:.c=.o)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CC) -c $< $(CPPFLAGS) -o $@

clean:
	-rm -f *.o ftest *~

