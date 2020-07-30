all: src/qsystem.cpp

%.cpp: %.i $(HEADER)
	swig -c++ -python -o $@ $<  
	mkdir -p qsystem
	mv src/qsystem.py qsystem/__init__.py
	patch qsystem/__init__.py src/qsystem.py.patch

clean:
	rm -rf __pycache__ qsystem doxygen/documentation
	rm -rf src/{qsystem.cpp,qsystem.py,*.o}
	rm -rf build dist qsystem QSystem.egg-info 
	