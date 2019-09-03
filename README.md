[![PyPI](https://img.shields.io/pypi/v/qsystem.svg)](https://pypi.org/project/QSystem/)
[![PyPI - License](https://img.shields.io/pypi/l/qsystem.svg?color=brightgree)](https://gitlab.com/evandro-crr/qsystem/blob/master/LICENSE)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/qsystem.svg?color=red)](https://www.python.org/)
[![Wiki](https://img.shields.io/badge/wiki-available-sucess.svg)](https://gitlab.com/evandro-crr/qsystem/wikis/home)
[![Doc](https://img.shields.io/badge/doc-available-succes.svg)](https://evandro-crr.gitlab.io/qsystem/index.html)

# QSystem
A quantum computing simulator for Python.

------------------------
The QSystem simulator is inspired in the quantum circuit model, so it's easy to
convert any quantum circuit to  Python.

Like the follow example:

![circ](https://gitlab.com/evandro-crr/qsystem/raw/1.1.0/circ.svg)

```python
from qsystem import QSystem
from cmath import exp, pi
q = QSystem(3, 24)                 # init q0, q1, q2

q.evol(gate='H', qbit=0, count=3)         # H q0; H q1; H q2
q.add_ancillas(4)                         # init a0, a1, a2, a3
          
q.evol(gate='X', qbit=6)                  # X a3
q.cnot(target=4, control=[2])             # CNOT a1, q2
q.cnot(5, [2])                            # CNOT a2, q2
q.cnot(5, [3])                            # CNOT a2, a0
q.cnot(3, [1, 5])                         # Toffoli a1, q1, a2
q.cnot(5, [3])                            # CNOT a2, a0
q.cnot(4, [6])                            # CNOT a1, a3
q.cnot(6, [1, 4])                         # Toffoli a3, q1, a1
q.cnot(4, [6])                            # CNOT a1, a3

q.measure(qbit=3, count=4)                # measure a0, a1, a2, a3
print('ancillas measurement =', q.bits()[3:])
# ancillas measurement = [0, 1, 0, 0]
q.rm_ancillas()                           # rm a0, a1, a2, a3

q.evol('H', 0)                            # H q0                ┐
q.cphase(phase=1j, target=1, control=[0]) # Controlled S q1, q0 │
q.evol('H', 1)                            # H q1                │
q.cphase(exp(pi*1j/4), 2, [0])            # Controlled T q2, q0 │ = q.qft(0, 3)
q.cphase(1j, 2, [1])                      # Controlled S q2, q1 │
q.evol('H', 2)                            # H q1                │
q.swap(0, 2)                              # SWAP q0, q2         ┘

q.measure(0, 3)                           # measure q0, q1, q2
print('final measurement =', q.bits())
# final measurement = [1, 0, 0]
```

* [Shor's factoring algorithm using QSystem](https://nbviewer.jupyter.org/urls/gitlab.com/evandro-crr/qsystem/raw/1.2.x/example/factoring.ipynb)

# Installation
QSystem depends on [Boost C++ Libraries](https://www.boost.org) and requires a 
C/C++ compiler. 

To install use the follow command:
```
pip install QSystem
```
---------------------------
Seed the [wiki](https://gitlab.com/evandro-crr/qsystem/wikis/home) and the
[documentation](https://evandro-crr.gitlab.io/qsystem/index.html).

---------------------------
This software is supported by 
<img src="http://www.fapesc.sc.gov.br/wp-content/uploads/2014/09/logo-Fapesc-fundo-transparente.png"
alt="FEPESC" width="54">
