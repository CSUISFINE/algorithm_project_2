CXX = g++
CXXFLAGS = -std=c++23 -O3 -Wall -g

LDFLAGS = -lm

SRC = my_algo.cc

TARGET = my_algo

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o

.PHONY: all clean