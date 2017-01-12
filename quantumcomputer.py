#!/usr/bin/env python

import numpy as np #for numerics
import random # for random number generation
from math import sqrt,pi,e # some common function
import qutip as qp


####
## States
####
    
class State:
# call with State(list_of_basisstatepositions)
# or State(list_of_amplitudes)
# needs global variable num_qubits


    def __init__(self, lis = [0]): # Initialise a State object
        # if the list has less entries than amplitudes (2**n),
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
        # else if the list has as many entries as amplitudes (2**n),
        # interpret the entries as amplitudes
        elif len(lis) == 2**num_qubits:
            self.state = np.array(lis)
            if not self.is_normalised():
                self.renormalise()
                print 'Note thate the state you generated was normalised ' \
                    'automatically '
        else:
            raise StandardError('Cannot interpret input of State() creator.'\
                                ' Please enter a list of valid amplitudes or positions.')
        
    def renormalise(self): # Renormalise the amplitude vector to unit length
        normalis_factor = float(np.sqrt(np.vdot(self.state,self.state)))
        self.state = self.state/normalis_factor
        return None

    def is_normalised(self): #Check if a state is normalised
        if np.isclose(float(np.vdot(self.state,self.state)), float(1.0)):
            return True
        else:
            return False

    def measure(self, runs = 1, output = 'outcomes'): #perform measurements on the state
    # simulates a repeated generation and projective measurement of the state
    # Options: 'outcomes' prints every result while 'stats' prints an overview
        results = np.random.choice(len(self.state), runs, p=[abs(el)**2 for el in self.state])
        if output == 'outcomes':
            print "Measurement Results"
            print "Index  Basis state "
            print "-----   ----------- "
            for el_res in results:
                print "{0:04}".format(el_res),'   |', "".join(( "{0:0", \
                                str(int(np.log2(len(self.state)))),"b}")).format(el_res),'>'
        if output == 'stats':
            #still to do
            
                 
    def print_me(self, style = None): # print out current state.
    # Options:
    # None/empty - simple output of state information
    # 'slim' - As in 'None' but without zero amplitude entries
    # 'amplitudes' - only array of amplitudes
        np.set_printoptions(precision=3, suppress = True)
        print
        if style == None: # print all nonzero amplitudes
            print "Index Probability Amplitude Basis state "
            print "----- ----------- --------- ----------- "
            for i in range(len(self.state)):
                if not np.isclose(self.state[i],0.0):
                    basis_string = "".join(( "{0:0", \
                                str(int(np.log2(len(self.state)))),"b}"))
                    print '', "{0:04}".format(i), '    ', \
                        "{0:.3f}".format(abs(self.state[i])**2), '   ', \
                        "{0:.3f}".format(self.state[i]), \
                        "".join(('  |',basis_string.format(i) , '>' ))
        if style == 'full': # print all amplitudes
            print "Index Probability Amplitude Basis state "
            print "----- ----------- --------- ----------- "
            for i in range(len(self.state)):
                    basis_string = "".join(( "{0:0", str(int(np.log2(len(self.state)))),"b}"))
                    print '', "{0:04}".format(i), '    ', \
                        "{0:.3f}".format(abs(self.state[i])**2), '   ', \
                        "{0:.3f}".format(self.state[i]), \
                        "".join(('  |',basis_string.format(i) , '>' ))
            
        if style == 'amplitudes':
            print "Amplitudes: ", self.state

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


#SCRIPT
global num_qubits;
num_qubits = 10



# Test Maria's printing
num_qubits = 5
myState1 = State( [1,2,20] )
myState1.print_me()
myState1.print_me('full')
myState1.print_me('amplitudes')
# Test Maria's measuring
myState1.measure(50)


