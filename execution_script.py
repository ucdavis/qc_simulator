from quantum_simulator import *
'''
mynewstate = create_state(2, [0,complex(1,1),complex(1,1),complex(1,1)])

print_me(mynewstate)

for i in range(10):
    mynewstate = grover_iteration(mynewstate, [0])
    print_me(mynewstate, 'probabilities')
    print is_normalised(mynewstate)

'''
#################### Execution

quantum_state = create_state(3,[0])
print Rx(-math.pi/4)
quantum_state = apply_unitary(Y, [1,2], quantum_state)
#project_on_blochsphere(quantum_state)

'''
quantum_state = create_state(4,[15])
print_me(quantum_state, 'full')

quantum_state = apply_total_unitary(X, [0,1,2], quantum_state)

print_me(quantum_state, 'full')
#project_on_blochsphere(state)
measure(quantum_state, 4)
'''
