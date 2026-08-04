CXX      ?= g++
CXXFLAGS ?= -std=c++11 -O3 -Wall
TARGET   := voidVolume3D
SRC      := voidVolume3D.cpp

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(SRC) $(CXXFLAGS) -o $(TARGET)

# Example: make run FILE=N3000GaussianDistribution.dat PROBE=0.0
run: $(TARGET)
	./$(TARGET) $(FILE) $(PROBE)

clean:
	rm -f $(TARGET)
