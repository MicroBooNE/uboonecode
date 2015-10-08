#
# Makefile for gmbrowser.
#
#
# Location of kerberos may have to be changed.
#
KRB5 =  /usr/kerberos

LIBS =  -L`root-config --libdir` -L/usr/local/lib \
	`root-config --glibs`
#        -lCore -lPhysics -lCint -lPostscript -lTree -lHist -lMatrix \
#       -lGraf -lGX11 -lGpad -lGui \
#	-L $(KRB5)/lib -lkrb5 -lk5crypto -lcom_err 
 
OPTCOMP =   -g
CXXFLAGS =  `root-config --cflags` -I.. -I$(KRB5)/include
OBJECTS = src/GM.o src/gmbrowser.o src/GMdict.o 

GCC = g++ 

all: bin/gmplotter bin/gmbrowser


bin/gmbrowser: $(OBJECTS)
	$(GCC) $(OPTCOMP) $(CXXFLAGS) -o bin/gmbrowser -g $(OBJECTS) $(LIBS)

bin/gmplotter: src/GMPlot.o src/gmplotter.o
	$(GCC) $(OPTCOMP) $(CXXFLAGS) -o bin/gmplotter -g \
	src/GMPlot.o src/gmplotter.o \
	$(LIBS)

src/GM.o: src/GM.cpp
	cd src; $(GCC) $(OPTCOMP) $(CXXFLAGS) -c GM.cpp; cd ..

src/GMPlot.o: src/GMPlot.cpp gmbrowser/GMPlot.hpp
	cd src; $(GCC) $(OPTCOMP) $(CXXFLAGS) -c GMPlot.cpp; cd ..

src/gmbrowser.o: src/gmbrowser.cpp
	cd src; $(GCC) $(OPTCOMP) $(CXXFLAGS) -c gmbrowser.cpp; cd ..

src/gmplotter.o: src/gmplotter.cpp
	cd src; $(GCC) $(OPTCOMP) $(CXXFLAGS) -c gmplotter.cpp; cd ..

src/GMdict.o: src/GMdict.cpp
	cd src; $(GCC) $(OPTCOMP) $(CXXFLAGS) -c GMdict.cpp; cd ..

src/GMdict.cpp: gmbrowser/GMBrowser.hpp
	cd src; rootcint -f GMdict.cpp -c ../gmbrowser/GMBrowser.hpp ../gmbrowser/LinkDef.h


clean:
	rm -f src/*.o
	rm -f src/GMdict.*
	rm -f bin/gmbrowser
	rm -f *~;
