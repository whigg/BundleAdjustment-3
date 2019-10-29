
LLIB_OPENCV = -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_calib3d -lopencv_imgcodecs
CC = g++

CFLAGS = -std=c++11 -Wall

ba: ba.o
	${CC} ba.o -o ba 
ba.o: ba.cpp ba.hpp
	${CC} -c ba.cpp ${CFLAGS} 

.PHONY : clean
clean:
	rm -rf *.o
