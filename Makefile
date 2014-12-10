CXX      := c++
CXXFLAGS := -O2 -std=c++11 -I. `root-config --cflags` -I/usr/include/boost
LDFLAGS  := -O2 -L. `root-config --ldflags --glibs` -L/usr/lib
LDLIBS   := -lboost_program_options
TARGET   := main
HEADERS  := $(wildcard *.h)

$(TARGET): $(TARGET).o
	$(CXX) $(LDFLAGS) -o $@ $< $(LDLIBS)

$(TARGET).o: $(TARGET).cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	-@$(RM) $(TARGET) $(TARGET).o
