from quantum_simulator_py2 import *

#Write code here

quantum_state = create_state(2,[0])
quantum_state = apply_unitary(Y,[0,1],quantum_state)
