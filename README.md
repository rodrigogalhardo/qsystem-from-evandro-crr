# QSystem
![PyPI](https://img.shields.io/pypi/v/qsystem.svg)
![PyPI - License](https://img.shields.io/pypi/l/qsystem.svg?color=brightgree)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/qsystem.svg?color=red)

A quantum computing simulator for Python.

------------------------
The QSystem simulator is inspired in the quantum circuit model, so it's easy
to convert any quantum circuit to  Python, like the follow example:

![circ](https://gitlab.com/evandro-crr/qsystem/raw/master/circ.svg?inline=false)

```python
from qsystem import Gates, QSystem
from cmath import exp, pi
gates = Gates()

q = QSystem(3, gates, 24)
q.evol(gate='H', qbit=0, count=3)
q.add_ancillas(4)

q.evol(gate='X', qbit=6)
q.cnot(target=4, control=[2])    
q.cnot(5, [2])    
q.cnot(5, [3])    
q.cnot(3, [1, 5]) 
q.cnot(5, [3])    
q.cnot(4, [6])    
q.cnot(6, [1, 4]) 
q.cnot(4, [6])    

q.measure(qbit=3, count=4)
print('ancillas measurement =', q.bits()[3:])
q.rm_ancillas()

q.evol('H', 0)
q.cphase(phase=1j, target=1, control=[0])
q.evol('H', 1)
q.cphase(exp(pi*1j/4), 2, [0])
q.cphase(1j, 2, [1])
q.evol('H', 2)
q.swap(0, 2)

q.measure(0, 3)
print('final measurement =', q.bits())
```

Seed the [wiki](https://gitlab.com/evandro-crr/qsystem/wikis/home) for
documentation.

---------------------------
This software is suported by 
<img src="http://www.fapesc.sc.gov.br/wp-content/uploads/2014/09/logo-Fapesc-fundo-transparente.png"
alt="FEPESC" width="54">
