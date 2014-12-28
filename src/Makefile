CXX=g++
CXXFLAGS = -ggdb -std=c++11
LDFLAGS=-lstdc++ -lm

TARGET	= indelible
SRCS	= control.cpp  indelible.cpp  models.cpp  paml.cpp  randoms.cpp # M5-13.cpp
OBJS   	= $(SRCS:.cpp=.o)
DEPS   	= $(SRCS:.cpp=.depends)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

# GNU makefile automatic dependency tracking cribbed from 
# http://stackoverflow.com/a/1510648/1135316

%.depends: %.cpp
	$(CXX) -MM $(CXXFLAGS) $< > $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)
