CXX = g++
CXXFLAGS = -std=c++11 --debug
INCLUDE	= -I/usr/local/include -I./include
LDFLAGS = -std=c++11
LDLIBS	= -lpthread -lm
EXECUTABLE= bin/main
SOURCES = $(wildcard src/*.cpp)
HEADERS = $(wildcard includes/*.h)
OBJECTS = $(patsubst src/%.cpp, obj/%.o, $(SOURCES))

BASE = $(USER)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

obj/%.o: src/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

clean:
	rm $(OBJECTS)