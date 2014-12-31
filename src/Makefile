CXX=g++
CXXFLAGS = -ggdb -std=c++11 -Wunused-variable
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

test : test_nucleotide test_indel 

test_nucleotide:
	rm -rf new old
	mkdir new old
	cd new && ../indelible -q ../../control.txt
	cd old && ../../src.master/indelible ../../control.txt
	diff old/outputname.fas new/outputname.fas
	diff old/outputname_TRUE.phy new/outputname_TRUE.phy

test_indel:
	rm -rf new old
	mkdir new old
	cd new && ../indelible -q ../../tutorial/indels.txt
	cd old && ../../src.master/indelible ../../tutorial/indels.txt
	diff old/outputname1.fas new/outputname1.fas
	diff old/outputname2.fas new/outputname2.fas
	diff old/outputname3.fas new/outputname3.fas
	diff old/outputname4.fas new/outputname4.fas
	diff old/outputname1_TRUE.phy new/outputname1_TRUE.phy
	diff old/outputname2_TRUE.phy new/outputname2_TRUE.phy
	diff old/outputname3_TRUE.phy new/outputname3_TRUE.phy
	diff old/outputname4_TRUE.phy new/outputname4_TRUE.phy
