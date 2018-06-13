all: main

CXX = g++
CXXFLAGS = -std=c++1z -Wall -Werror -fno-exceptions   # -fo-rtti wasn't working on my mac

main: main.o multiarray.o io_plink.o methods.o li_stephens.o test.o
	$(CXX) $(CXXFLAGS) -o main main.o multiarray.o io_plink.o methods.o li_stephens.o test.o

objects = main.o multiarray.o io_plink.o methods.o li_stephens.o test.o

main.o: main.cc
	$(CXX) $(CXXFLAGS) -c main.cc

multiarray.o: multiarray.cc
	$(CXX) $(CXXFLAGS) -c multiarray.cc

io_plink.o: io_plink.cc
	$(CXX) $(CXXFLAGS) -c io_plink.cc

methods.o: methods.cc
	$(CXX) $(CXXFLAGS) -c methods.cc
	
li_stephens.o: li_stephens.cc
	$(CXX) $(CXXFLAGS) -c li_stephens.cc	

test.o: test.cc
	$(CXX) $(CXXFLAGS) -c test.cc

.PHONY : clean
clean :
	-rm edit $(objects)
