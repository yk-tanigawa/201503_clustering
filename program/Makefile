CC = g++
LD = g++
CFLAGS =  -O2 -Warning
LDFLAGS =  -lpthread
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.dep)
EXEC = $(SRCS:.cpp=)
RM = rm -f


all: rand_data_generator Lloyd Hamerly

rand_data_generator: rand_data_generator.o
	$(LD)  -o $@ $^ $(LDFLAGS)

Lloyd.o: my_vector.hpp time_bench.hpp
Lloyd: Lloyd.o 
	$(LD)  -o $@ $^ $(LDFLAGS)

Hamerly.o: my_vector.hpp time_bench.hpp
Hamerly: Hamerly.o
	$(LD)  -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJS) $(EXEC) *~

check-syntax:
	$(CC) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)

.PHONY: clean 
	all clean
