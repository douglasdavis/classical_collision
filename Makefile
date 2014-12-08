CXX      := c++ 
CXXFLAGS := -O2 -std=c++11 -I.
LDFLAGS  := -L.
LDLIBS   :=
HEADERS  := $(wildcard *.h)

all: main

main: main.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LDLIBS)

main.o: main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	-@$(RM) main main.o
