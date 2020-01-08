OBJ = src/gate.o src/qs_ancillas.o src/qs_errors.o
OBJ += src/qs_make.o src/qs_evol.o src/qs_measure.o
OBJ += src/qs_utility.o src/qsystem.o src/utility.o
OBJ += src/microtar.o
HEADER = $(wildcard header/*.h)

OUT = qsystem/_qsystem.so

PYTHON = /usr/include/python3.8/

CXXFLAGS = -Wall -O2 -fPIC -std=c++17 -I$(PYTHON)
CLINK = -shared -Xlinker -export-dynamic -lboost_serialization 

all: $(OBJ)
	$(CXX) $(OBJ) -o $(OUT) $(CXXFLAGS) $(CLINK)

%.o: %.cpp $(HEADER)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

%.o: %.c 
	$(CC) -c $< -o $@ -O2 -fPIC

%.cpp: %.i $(HEADER)
	swig -c++ -python -o $@ $<  
	mkdir -p qsystem
	mv src/qsystem.py qsystem/__init__.py
	patch qsystem/__init__.py src/qsystem.py.patch

dist: src/qsystem.cpp qsystem/__init__.py armadillo-code 
	python setup.py sdist 

armadillo-code:
	git clone https://gitlab.com/conradsnicta/armadillo-code.git --branch 9.800.x

install: dist
	pip install dist/* --user 

doc:
	doxygen doxygen/Doxyfile

clean:
	rm -rf __pycache__ qsystem doxygen/documentation
	rm -rf src/{qsystem.cpp,qsystem.py,*.o}
	rm -rf build dist qsystem QSystem.egg-info 
	
clean-all: clean
	rm -rf armadillo-code
