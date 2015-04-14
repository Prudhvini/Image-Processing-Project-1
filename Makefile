all: CImg.h p1.cpp
	g++ p1.cpp -o p1 -lX11 -lpthread

clean:
	rm p1
