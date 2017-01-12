#!/usr/bin/env python

import numpy as np # for numerics
import random # for random number generation
from math import sqrt,pi,e # some common function
import qutip as qp # qutip library for Bloch sphere visualisations
import cmath # library for complex numbers



#################### State functions
def create_state(num_qubits, lis):
    state = np.zeros(num_qubits)
    # if the list has less entries than amplitudes (2**n),
    # interpret the entries as positions
    if len(lis) < 2**num_qubits:
        # check if all entries are valid positions
        if any(isinstance(item, complex) or item > 2**num_qubits \
                   or not isinstance(item, int) for item in lis) :
            raise StandardError('Cannot interpret input of State() creator.'\
                                    ' Please enter a list of valid positions.')
        # initialise state
        state = np.array([1./sqrt(len(lis)) if i in lis else 0 \
                                 for i in range(2**num_qubits)])
    # else if the list has as many entries as amplitudes (2**n),
    # interpret the entries as amplitudes
    elif len(lis) == 2**num_qubits:
        state = np.array(lis)
        if not is_normalised(state):
            state = renormalise(state)
            print 'Note thate the state you generated was normalised ' \
                    'automatically '

    else:
        raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid amplitudes or positions.')
    return state



def renormalise(state): # Renormalise the amplitude vector to unit length
        normalis_factor = float(np.sqrt(np.vdot(state,state)))
        state = state/normalis_factor
        return state

def is_normalised(state): #Check if a state is normalised
    if np.isclose(float(np.vdot(state,state)), float(1.0)):
        return True
    else:
        return False

def measure(state, runs = 1, output = 'outcomes'): #perform measurements on the state
# simulates a repeated generation and projective measurement of the state
# Options: 'outcomes' prints every result while 'stats' prints an overview
    results = np.random.choice(len(state), runs, p=[abs(el)**2 for el in state])
    if output == 'outcomes':
        print "Measurement Results"
        print "Index  Basis state "
        print "-----   ----------- "
        for el_res in results:
            print "{0:04}".format(el_res),'   |', "".join(( "{0:0", \
                                str(int(np.log2(len(state)))),"b}")).format(el_res),'>'
    #if output == 'stats':
    #still to do
    return None

def print_me(state, style = None): # print out current state.
# Options:
# None/empty - simple output of state information
# 'slim' - As in 'None' but without zero amplitude entries
# 'amplitudes' - only array of amplitudes
    np.set_printoptions(precision=3, suppress = True)
    print
    if style == None: # print all nonzero amplitudes
        print "Index Probability Amplitude Basis state "
        print "----- ----------- --------- ----------- "
        for i in range(len(state)):
            if not np.isclose(state[i],0.0):
                basis_string = "".join(( "{0:0", \
                                str(int(np.log2(len(state)))),"b}"))
                print '', "{0:04}".format(i), '    ', \
                        "{0:.3f}".format(abs(state[i])**2), '   ', \
                        "{0:.3f}".format(state[i]), \
                        "".join(('  |',basis_string.format(i) , '>' ))
    if style == 'full': # print all amplitudes
        print "Index Probability Amplitude Basis state "
        print "----- ----------- --------- ----------- "
        for i in range(len(state)):
            basis_string = "".join(( "{0:0", str(int(np.log2(len(state)))),"b}"))
            print '', "{0:04}".format(i), '    ', \
                        "{0:.3f}".format(abs(state[i])**2), '   ', \
                        "{0:.3f}".format(state[i]), \
                        "".join(('  |',basis_string.format(i) , '>' ))

    if style == 'amplitudes':
        print "Amplitudes: ", state

    print
    return None

#################### Gate functions



#################### Execution

state = create_state(2,[1,1,1,1])
print_me(state, 'full')
measure(state, 4)
