CC = g++
SRC = ../src
TARGET = ../build
SRC_H = $(SRC)/h
SRC_CPP = $(SRC)/cpp
All = main.o CMesh.o BaseMesh.o IDT.o Utils.o

IDT: $(All)
	$(CC) -o $(TARGET)/IDT $(All)
main.o: $(SRC)/main.cpp $(SRC_H)/CMesh.h
	$(CC) -c -std=c++11 $(SRC)/main.cpp -g
CMesh.o: $(SRC_H)/CMesh.h $(SRC_H)/Geometry.h $(SRC_H)/BaseMesh.h $(SRC_CPP)/CMesh.cpp
	$(CC) -c -std=c++11 $(SRC_CPP)/CMesh.cpp -g
BaseMesh.o: $(SRC_H)/BaseMesh.h $(SRC_H)/Geometry.h $(SRC_H)/Utils.h $(SRC_CPP)/BaseMesh.cpp
	$(CC) -c -std=c++11 $(SRC_CPP)/BaseMesh.cpp -g
IDT.o: $(SRC_H)/IDT.h $(SRC_CPP)/IDT.cpp
	$(CC) -c -std=c++11 $(SRC_CPP)/IDT.cpp -g
Utils.o: $(SRC_H)/Utils.h $(SRC_CPP)/Utils.cpp
	$(CC) -c -std=c++11 $(SRC_CPP)/Utils.cpp -g
clean:
	rm $(TARGET)/IDT $(All)