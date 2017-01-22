# qc_simulator
A simple educational quantum computing simulator in Python (for Quantum Machine Learning Summer School Jan 2017)

*******************
CREATE A NEW STATE:

create_state(num_qubits, lis)

INPUTS:
-> num_qubits: 
				number of qubits of the quantum system
-> lis: 
				if lis is a list of integers with length < num_qubits, it is 
				interpreted as a list of positions of nonzero amplitudes in a 
				uniform superposition 

				if lis is a list of numbers with length = num_qubits, it is 
				interpreted as an amplitude vector

OUTPUTS:
-> numpy array of amplitudes representing a quantum state


-----Examples: 

mynewstate = create_state(2, [0,3])
print mynewstate

>>> [ 0.707  0.     0.     0.707]

mynewstate = create_state(2, [1/sqrt(2),1/sqrt(2), 0,0])
print mynewstate

>>> [ 0.707  0.707  0.     0.   ]

mynewstate = create_state(2, [1, 2, 3, 4])
print mynewstate

>>> Note thate the state you generated was normalised automatically 
>>> [ 0.183  0.365  0.548  0.73 ]

mynewstate = create_state(2, [1,2,3,4,5])
print mynewstate

>>> StandardError: Cannot interpret input of State() creator. Please enter a list of valid amplitudes or positions.

******************
RENORMALISE VECTOR:

renormalise(state)

INPUTS: 
			-> state: a numpy array representing the amplitude vector of
								a non-normalised quantum state
OUTPUTS:
			-> the same amplitude vector but normalised to unit length


-----Examples:
 
my_normalised_state = normalise(state)
print my_normalised_state

>>> [ 0.707  0.     0.     0.707]

******************
FIND OUT IF VECTOR IS NORMALISED:

is_normalised(state)

INPUTS: 
			-> state: a numpy array representing the amplitude vector of
								a non-normalised quantum state
OUTPUTS:
			-> True or False (boolean), depending on whether the state is normalised
					(error tolerance is set to 1e-03 = 0.001)

-----Examples: 

print is_normalised([ 0.707  0.     0.     0.707])

>>> True

print is_normalised([ 1, 2, 3, 4])

>>> False

***********************
MEASURE A STATE

measure(state, runs, output)

INPUTS:
			-> state: a numpy array representing the amplitude vector of
								a quantum state
			-> runs:  number of times the measurement has to be repeated 
								(simulating a repeated preparation of the state before 
								each measurement)
			-> output: Type of printout
								'outcomes' - prints a table of the (stochastic) measurement
														outcomes
								'stats' - prints a histogram of the (stochastic) measurement
													outcomes

OUTPUTS:
			-> prints the measurement results


-----Examples: 

measure(create_state(4, [0,2,4,6]),2, 'outcomes') 

>>> Measurement Results: 
>>>   Run    Index  Basis state
>>> -----  -------  -------------
>>>     1        4  |0100>
>>>     2        0  |0000>



measure(create_state(4, [0,2,4,6]),4, 'stats') 

>>> Measurement Statistics: 
>>>   rel. occ.    Index  Basis state
>>> -----------  -------  -------------
>>>        0.25        2  |0010>
>>>        0.5         4  |0100>
>>>        0.25        6  |0110>





****************************
PRINT FUNCTION

print_me(state, style)

INPUT: 
			-> state: a numpy array representing the amplitude vector of
								a quantum state
			-> style: How to print
								-- 'None' or no entry: prints a table of only
										the nonzero basis states
								-- 'full': prints all states
								-- 'amplitudes': prints the amplitude vector only
OUTPUT:
			-> printout of measurements

----Examples:

print print_me(create_state(2, [0]), None)

>>> Quantum State:
>>>   Index    Probability    Amplitude  Basis state
>>> -------  -------------  -----------  -------------
>>>       0            0.5        0.707  |00>
>>>       2            0.5        0.707  |10>



print_me(create_state(2, [0,2]), 'full')

>>>  Quantum State:
>>>   Index    Probability    Amplitude  Basis state
>>> -------  -------------  -----------  -------------
>>>       0            0.5        0.707  |00>
>>>       1            0          0      |01>
>>>       2            0.5        0.707  |10>
>>>       3            0          0      |11>




*******************************
GROVER ITERATION:

grover_iteration(state, marked_pos)

INPUT:
			-> state
			-> marked_pos: list of integers of the marked positions



