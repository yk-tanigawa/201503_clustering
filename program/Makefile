CC = g++
LD = g++
CXXFLAGS := -std=c++0x -O2 -Wall -Wextra -pedantic-errors 
LDFLAGS = -std=c++0x -lpthread
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.dep)
EXEC = $(SRCS:.cpp=)
RM = rm -f


all: rand_data_generator kmeans

rand_data_generator: rand_data_generator.o
	$(LD)  -o $@ $^ $(LDFLAGS)

Lloyd.hpp: kmeans.hpp my_ublas.hpp
Hamerly.hpp: kmeans.hpp my_ublas.hpp
kmeans.o: time_bench.hpp Lloyd.hpp Hamerly.hpp
kmeans: kmeans.o
	$(LD)  -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJS) $(EXEC) *.txt *~

check-syntax:
	$(CC) -std=c++0x -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)

.PHONY: clean 
	all clean
