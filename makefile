# compiler
CC = g++

CFLAGS = -g -Wall -I./../eigen


house: house.o

house.o: house.cpp house.hpp
