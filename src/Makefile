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

test : test_nucleotide test_indel test_aminoacid

indelible.orig:
	wget https://github.com/fhcrc/indelible/releases/download/baseline/indelible -O $@
	chmod +x $@

test_nucleotide: indelible.orig
	rm -rf new old
	mkdir new old
	cd new && ../indelible -q ../../testing/basicnucleotide/basicnucleotide.txt
	cd old && ../indelible.orig ../../testing/basicnucleotide/basicnucleotide.txt "./"
	diff -x LOG.txt old/ new/

test_indel: indelible.orig
	rm -rf new old
	mkdir new old
	cd new && ../indelible -q ../../testing/indels/indels.txt
	cd old && ../indelible.orig ../../testing/indels/indels.txt "./"
	diff -x LOG.txt old/ new/

test_aminoacid: indelible.orig
	rm -rf new old
	mkdir new old
	cd new && ../indelible -q ../../testing/AMINOACID/AMINOACID.txt
	cd old && ../indelible.orig ../../testing/AMINOACID/AMINOACID.txt "./"
	diff -x LOG.txt old/ new/

test_codon: indelible.orig
	rm -rf new old
	mkdir new old
	cd new && ../indelible -q ../../testing/CODON/CODON.txt
	cd old && ../indelible.orig ../../testing/CODON/CODON.txt "./"
	diff -x LOG.txt old/ new/
