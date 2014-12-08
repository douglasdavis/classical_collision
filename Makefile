CXX      := c++
CXXFLAGS := -O2 -std=c++11 -I.
LDFLAGS  := -L.
LDLIBS   :=
TARGET   := main
HEADERS  := $(wildcard *.h)

$(TARGET): $(TARGET).o
	$(CXX) $(LDFLAGS) -o $@ $< $(LDLIBS)

$(TARGET).o: $(TARGET).cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

clean:
	-@$(RM) $(TARGET) $(TARGET).o
