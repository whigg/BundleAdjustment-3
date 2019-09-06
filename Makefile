
LLIB_OPENCV = -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_calib3d
CC = g++

ba: ba.o
	${CC} ba.o -o ba ${LLIB_OPENCV}
ba.o: ba.cpp ba.hpp
	${CC} -c ba.cpp

.PHONY : clean
clean:
	rm -rf *.o
