CC=g++
TARGET=HomologyLocalization_Code
IncludeDir=./Third_Party/
RootDir=./$(TARGET)
SRC=$(wildcard $(TARGET)/*.cpp)
CFLAGS=

all:
	$(CC) -g -std=c++11 $(SRC) -o HomologyLocalization -I$(IncludeDir) -I$(RootDir)/ -I$(RootDir)/Algorithms -I$(RootDir)/DataReaders -I$(RootDir)/External -I$(RootDir)/Filtration
