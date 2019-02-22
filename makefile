OBJ = src/gate.o src/microtar.o src/qs_ancillas.o src/qs_errors.o
OBJ += src/qs_make.o src/qs_evol.o src/qs_measure.o
OBJ += src/qs_utility.o src/qsystem.o
HEADER = $(wildcard header/*.h)

OUT = _qsystem.so

PYTHON = /usr/include/python3.7m/

CFLAGS = -Wall -O2 -fPIC
CXXFLAGS = $(CFLAGS) -std=c++17 -I$(PYTHON)
CLINK = -shared -Xlinker -export-dynamic

all: $(OBJ) qsystem.py
	$(CXX) $(OBJ) -o $(OUT) $(CXXFLAGS) $(CLINK)

%.o: %.cpp $(HEADER)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

%.cpp: %.i $(HEADER)
	swig -c++ -python -o $@ $< 

qsystem.py:
	ln -s src/qsystem.py $@

dist: src/qsystem.cpp qsystem/__init__.py armadillo-code
	python setup.py sdist 

qsystem/__init__.py:
	mkdir -p qsystem
	ln -s ../src/qsystem.py $@

armadillo-code:
	git clone https://gitlab.com/conradsnicta/armadillo-code.git

install: dist
	pip install dist/*

clean:
	rm -rf $(OUT) __pycache__ qsystem.py
	rm -rf src/{qsystem.cpp,qsystem.py,*.o}
	rm -rf build dist qsystem QSystem.egg-info armadillo-code

