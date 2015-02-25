CC = g++
LD = g++
CXXFLAGS := -std=c++0x -O2 -Wall -Wextra -pedantic-errors 
LDFLAGS = -std=c++0x -lpthread
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.dep)
EXEC = $(SRCS:.cpp=)
RM = rm -f


all: rand_data_generator Lloyd Hamerly kmeans

rand_data_generator: rand_data_generator.o
	$(LD)  -o $@ $^ $(LDFLAGS)

Lloyd.o: my_vector.hpp time_bench.hpp
Lloyd: Lloyd.o 
	$(LD)  -o $@ $^ $(LDFLAGS)

Hamerly.o: my_vector.hpp time_bench.hpp
Hamerly: Hamerly.o
	$(LD)  -o $@ $^ $(LDFLAGS)

kmeans.o: time_bench.hpp
kmeans: kmeans.o
	$(LD)  -o $@ $^ $(LDFLAGS)


clean:
	$(RM) $(OBJS) $(EXEC) *~

check-syntax:
	$(CC) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)

.PHONY: clean 
	all clean
