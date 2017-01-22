from quantum_simulator import *

mynewstate = create_state(2, [0,complex(1,1),complex(1,1),complex(1,1)])

print_me(mynewstate)

for i in range(10):
    mynewstate = grover_iteration(mynewstate, [0])
    print_me(mynewstate, 'probabilities')
    print is_normalised(mynewstate)
g
