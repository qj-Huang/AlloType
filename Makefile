CXX = g++
CXXFLAGS = -I $${LOCAL}/include/eigen -I $${LOCAL}/include/boost -std=c++11

objects = main.o read_pdb.o method.o Pro.o ProAnalysis.o
exe = AlloType.exe

$(exe) : $(objects)
	$(CXX) $(CXXFLAGS) -o $(exe) $(objects) $${LOCAL}/include/boost/stage/lib/*.a

main.o : Pro.h ProAnalysis.h
read_pdb.o : read_pdb.h
method.o : method.h
Pro.o : Pro.h read_pdb.h
ProAnalysis.o : ProAnalysis.h method.h

.PHONY : clean
clean :
	rm $(exe) $(objects)
