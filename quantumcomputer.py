#!/usr/bin/env python

import numpy as np #for numerics
import random # for random number generation
from tabulate import tabulate # for nice printing
from math import sqrt,pi,e # some common function
import qutip as qp


####
## States
####
    
class State:
# call with State(number_of_qubits, list_of_basisstatepositions_OR_amplitudes)

    # ground state as default instance of a State 
    def __init__(self, lis = [0]):

        # if the list has less entries than amplitudes 2**n,
        # interpret the entries as positions
        if len(lis) < 2**num_qubits:
            # check if all entries are valid positions
            if any(isinstance(item, complex) \
                   or item > 2**num_qubits \
                   or not isinstance(item, int) \
                   for item in lis) :
                raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid positions.')
            # initialise state
            self.state = np.array([1./sqrt(len(lis)) if i in lis else 0 \
                                   for i in range(2**num_qubits)])
        # else if the list has as many entries as amplitudes 2**n,
        # interpret the entries as amplitudes
        elif len(lis) == 2**num_qubits:
            #TODO: check if normalised
            self.state = np.array(lis)
            
        else:
            raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid amplitudes or positions.')

                
    def print_amplitudes(self):
        np.set_printoptions(precision=3, suppress = True)
        print self.state

    def print_basisstates(self):
        np.set_printoptions(precision=3, suppress = True)
        print self.state

    def print_full(self):
        np.set_printoptions(precision=3, suppress = True)
        for i in range(len(self.state)):
            basis_string = "".join(( "{0:0", str(int(np.log2(len(self.state)))),"b}"))
            print round(self.state[i],2), "".join(('  |',basis_string.format(i) , '>' ))

    def print_full_pretty(self):
        np.set_printoptions(precision=3, suppress = True)
        data = []
        for i in range(len(self.state)):
            basis_string = "".join(( "{0:0", str(int(np.log2(len(self.state)))),"b}"))
            row = [ self.state[i],  "".join(('  |',basis_string.format(i) , '>' )) ]
            data.append(row)
        print
        print tabulate(data, headers=['Amplitude',' Basis state'] )
        print


####
## Gates
####

class Gate:

    # identity as default instance of a Gate 
    def __init__(self, unitary, qubit_pos = -1):
        #if not is_unitary(unitary):
        #    raise StandardError('Cannot create new Gate().'\
        #                        'Input matrix must be unitary.')
        if qubit_pos == -1:
            if cmp(unitary.shape, (2**num_qubits,2**num_qubits)) != 0:
                raise StandardError('Cannot create new Gate().'\
                   'Input matrix must be 2^n x 2^n.')
            self.gate = np.mat(unitary)

        if 0 <= qubit_pos < num_qubits:
            if cmp(unitary.shape, (2,2)) != 0:
                raise StandardError('Cannot create new Gate().'\
                                'Input matrix must be 2x2.')
            unitary_list = [unitary if qubit_pos == i else Gate.eye for i in range(num_qubits)]

            u_new = np.kron(unitary_list[0], unitary_list[1])
            for k in range(num_qubits-2):
                 u_new = np.kron(u_new, unitary_list[k+2])
            self.gate = u_new
                
    
    # define elementary gates
    i_ = np.complex(0,1)
    H = 1./sqrt(2)*np.matrix('1 1; 1 -1')
    X = np.matrix('0 1; 1 0')
    Y = np.matrix([[0, -i_],[i_, 0]])
    Z = np.matrix([[1,0],[0,-1]])
    eye = np.eye(2,2)
    S=np.matrix([[1,0],[0,i_]])
    Sdagger=np.matrix([[1,0],[0,-i_]])
    T=np.matrix([[1,0],[0, e**(i_*pi/4.)]])
    Tdagger=np.matrix([[1,0],[0, e**(-i_*pi/4.)]])
    #TODO: CNOT
# hallo


#SCRIPT
global num_qubits;
num_qubits = 3



newGate = Gate(Gate.X, 2)
myState1 = State( [0,1])
myState1.print_full()

#TESTMARIA
#comment master
