CC=g++
CXXFLAGS = -Dcheckindeldistributions=1
LIBS=-lstdc++ -lm
sources = control.cpp  indelible.cpp  models.cpp  paml.cpp  randoms.cpp # M5-13.cpp
objects := $(patsubst %.cpp,%.o,$(sources))

indelible: $(objects)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm *.o

clobber: clean
	rm indelible
