CXX   ?= c++
RM    ?= rm -rf
MKDIR ?= mkdir

.PHONY: all clean

all: lib/libsass.so

lib/libsass.so: libsass.cpp libsass.h
	$(MKDIR) -p lib
	$(CXX) libsass.cpp -I. -fPIC -shared -o lib/libsass.so

clean:
	rm -rf lib
