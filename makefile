SRC = $(wildcard src/*.cpp) src/qsystem.cpp
OBJ = $(SRC:.cpp=.o) src/microtar.o
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

dist: src/qsystem.cpp qsystem/__init__.py
	python setup.py sdist bdist_wheel

qsystem/__init__.py:
	mkdir -p qsystem
	ln -s ../src/qsystem.py $@

install: dist
	pip install dist/*

clean:
	rm -rf $(OUT) __pycache__ qsystem.py
	rm -rf src/{qsystem.cpp,qsystem.py,*.o}
	rm -rf build dist qsystem QSystem.egg-info

