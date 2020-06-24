# make file for suvendu's code, calculating 
#-------------------------------------
CC = g++
CFLAGS = -Wall -g
#-------------------------------------
myPro: chain.o p_vec_ten.o main.o 
	${CC} ${CFLAGS} chain.o p_vec_ten.o main.o -o myPro

p_vec_ten.o: p_vec_ten.cc p_vec_ten.hh
	${CC} ${CFLAGS} -c p_vec_ten.cc

chain.o: chain.cpp chain.h
	${CC} ${CFLAGS} -c chain.cpp

main.o: main.cpp p_vec_ten.hh p_vec_ten.cc
	${CC} ${CFLAGS} -c main.cpp







