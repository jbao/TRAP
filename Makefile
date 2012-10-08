# $Id: Makefile 339 2012-09-25 09:16:14Z jbao $

all:
	g++ TRAP.cpp -o TRAP -I/opt/local/include/ -g

clean:
	rm TRAP
