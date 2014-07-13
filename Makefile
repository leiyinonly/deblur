#INCLUDE = $(shell pkg-config --cflags opencv)
#LIBS = $(shell pkg-config --libs opencv)
INCLUDE = `pkg-config --cflags opencv`
LIBS = `pkg-config --libs opencv`

OBJECTS = deblur.o basicOperation.o restore.o assessment.o estimatePar.o
SOURCE = deblur.cpp basicOperation.cpp restore.cpp assessment.cpp estimatePar.cpp
BIN = deblur

$(BIN):$(OBJECTS)
	g++ -g -o $(BIN) $(OBJECTS) $(INCLUDE) $(LIBS)
deblur.o:deblur.cpp basicoperation.h
	g++ -c  deblur.cpp
basicOperation.o:basicOperation.cpp
	g++ -c basicOperation.cpp
restore.o:restore.cpp
	g++ -c restore.cpp
assessment.o:assessment.cpp
	g++ -c assessment.cpp
estmatePar.o:estimatePar.cpp
	g++ -c estimatePar.cpp
clean:
	rm $(OBJECTS) $(BIN)
