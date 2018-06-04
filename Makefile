all: main

CXX = g++
CXXFLAGS = -std=c++1z -Wall -Werror -fno-exceptions   # -fo-rtti wasn't working on my mac

main: main.o io_plink.o methods.o test.o
	$(CXX) $(CXXFLAGS) -o main main.o io_plink.o methods.o test.o

objects = main.o io_plink.o methods.o test.o

main.o: main.cc
	$(CXX) $(CXXFLAGS) -c main.cc

io_plink.o: io_plink.cc
	$(CXX) $(CXXFLAGS) -c io_plink.cc

methods.o: methods.cc
	$(CXX) $(CXXFLAGS) -c methods.cc

test.o: test.cc
	$(CXX) $(CXXFLAGS) -c test.cc

.PHONY : clean
clean :
	-rm edit $(objects)
