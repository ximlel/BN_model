GRPLAG.exe: GRPLAG.o VIPLimiter.o
	g++ -g -o GRPLAG.exe GRPLAG.o VIPLimiter.o -lm
GRPLAG.o: GRPLAG.cpp
	g++ -g -c GRPLAG.cpp
VIPLimiter.o: VIPLimiter.cpp
	g++ -g -c VIPLimiter.cpp

clean:
	rm -f GRPLAG.exe GRPLAG.o VIPLimiter.o
.PHONY:all clean
