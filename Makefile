all: main

main: main.cc
	g++ -std=c++1z -Wall -Werror -fno-exceptions main.cc -o main # -fo-rtti wasn't working on my mac
