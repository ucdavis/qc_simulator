# Quantum Computing Simulator (qc_simulator)

[![PyPI](https://img.shields.io/pypi/pyversions/Django.svg)]()

A simple educational **quantum computing simulator** written in Python. The library was specifically designed for the **28th Chris Engelbrecht Summer School of the South African National Institute for Theoretical Physics**.

This repository contains:

1. [The quantum computing simulator library for Python 2.7](quantum_simulator_py2.py).
2. [The quantum computing simulator library for Python 3.4, 3.5 and 3.6](quantum_simulator_py3.py).
3. [An execution script for Python 2.7](execution_script_py2.py) to use the quantum computing simulator library.
4. For people using Python 3.4, 3.5 or 3.6, please use [the execution script for these versions](execution_script_py3.py).
5. [The README](README.md) - containing the full documentation for the library.


## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
	- [Initialising quantum states](#initialising-quantum-states)
	- [Amplitude normalisation](#amplitude-normalisation)
	- [Printing quantum states](#printing-quantum-states)
	- [Quantum gates](#quantum-gates)
	- [Special functions](#special-functions)
	- [Measurement & statistics](#measurement-and-statistics)
- [Contribute](#contribute)
- [License](#license)

## Background

This library was designed and written by [@mariaschuld](https://github.com/mariaschuld) and [@markf94](https://github.com/markf94) from the [Centre for Quantum Technology](http://quantum.ukzn.ac.za/) at the [University of KwaZulu-Natal](www.ukzn.ac.za). It is meant as an educational tool to simulate quantum computations. The main aim was to come up with a syntax that is easy to read & understand for beginners. Additionally, the developed library is scalable to any number of qubits only limited by the amount of available RAM in your computer. The library includes initialisation of arbitrary quantum states and amplitude distributions, Bloch sphere projections, a variety of single-qubit gates, all-to-all connected controlled-U gates, where U can be any unitary 2x2 matrix and lastly, all-to-all connected Toffoli (CCNOT) gates. With easy syntax and without visible for-loops the simulation can be repeated many times with clear output of the resulting statistics.

## Install

This project uses [Python](http://python.org/) and is compatible with versions 2.7, 3.4, 3.5 and 3.6. The following packages are required: [numpy](http://www.numpy.org), [QuTiP](http://qutip.org), [random](https://docs.python.org/2/library/random.html), [collections](https://docs.python.org/2/library/collections) and [tabulate](https://pypi.python.org/pypi/tabulate). **Go check them out if you don't have them locally installed.**

Just download this repository as a ZIP file or clone it via your favourite shell. Unzip the folder and you are ready to go. No additional installation required.

## Usage

This section is a full documentation of the possibilities, functions and
syntax of this Python quantum computing simulator library. Use the
[table of contents](#table-of-contents) to navigate your way through the large
number of examples and explanations.

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------


### Initialising quantum states

General syntax:

```
create_state(num_qubits, list_type, lis)
```

INPUTS:

-> `num_qubits`: number of qubits of the quantum system.

-> `list_type`: A string with two options:

>-- `"indices"` if you want to initialise a specific qubit state with probability
1 or a uniform superposition over particular qubit states

>-- `"amplitudes"` if you want to specify the individual amplitudes of the state

-> `lis`: Two possibilities:

>-- if `list_type` = `indices`, `lis` is interpreted as a list of indices
of qubit states which will be put in a uniform superposition.

>--If `list_type` = `amplitudes`,`lis` must be a list of numbers with
length = 2^num_qubits. The quantum state will then be initialized with the specified amplitudes.

OUTPUT:

-> numpy array of amplitudes representing a quantum state.

--------------------------------------------------------------------------------

**Example 1:**

Suppose you want to initialise a uniform superposition over the 0th and 3rd
two-qubit states |00> and |11>. Therefore, we interpret the second input parameter
`lis` as a list of positions of nonzero amplitudes in a
uniform superposition. The syntax is as follows:

```
mynewstate = create_state(2, [0,3])
print(mynewstate)
```

This outputs the desired amplitude vector
`[ 0.70710678+0.j  0.00000000+0.j  0.00000000+0.j  0.70710678+0.j]
`.

--------------------------------------------------------------------------------

**Example 2:**

Suppose you want a superposition over the first two two-qubit states |00> and |01>.
Instead of using the input variable `lis` as a list of qubit positions we can fill
it with four amplitudes corresponding to the four two-qubit states.
To initialise a uniform superpositon over |00> and |01>, the following syntax can
be used:
```
mynewstate = create_state(2, [1/sqrt(2),1/sqrt(2), 0,0])
print(mynewstate)
```
This outputs the desired amplitude vector
`[ 0.70710678+0.j  0.70710678+0.j  0.00000000+0.j  0.00000000+0.j]`.

--------------------------------------------------------------------------------

**Example 3:**

If you choose two qubits with four amplitudes and if you specify four numbers,
the create_state() function interprets those as amplitudes. If your chosen amplitudes
are not normalised then the function will automatically do it for you. For example:

```
mynewstate = create_state(2, [1, 2, 3, 4])
print(mynewstate)
```

This will normalise the amplitude vector and inform you about it. The console
 output reads:
```
Note thate the state you generated was normalised automatically
[ 0.18257419+0.j  0.36514837+0.j  0.54772256+0.j  0.73029674+0.j]
```

--------------------------------------------------------------------------------

**Example 4:**

If you specify more than 2^n amplitudes for n qubits the function will return
an error. For example:

```
mynewstate = create_state(2, [1,2,3,4,5])
print(mynewstate)
```
leads to the console output:
```
StandardError: Cannot interpret input of State() creator.
Please enter a list of valid amplitudes or positions.
```

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------

### Amplitude normalisation

**Renormalising an amplitude vector**

If you want to renormalise an amplitude vector at any point, the general syntax
is:

```
renormalise(state)
```

INPUTS:

-> `state`: a numpy array representing the amplitude vector of a non-normalised quantum state

OUTPUT:

-> the same amplitude vector but normalised to unit length

--------------------------------------------------------------------------------

**Example 1:**

```
state = np.array([1,2,3,4])
my_normalised_state = renormalise(state)
print(my_normalised_state)
```
This outputs the normalised amplitude vector
`[ 0.18257419  0.36514837  0.54772256  0.73029674]`.

--------------------------------------------------------------------------------

**Finding out if an amplitude vector is normalised**

If you want to find out if an amplitude vector is normalised the syntax is:

```
is_normalised(state)
```

INPUTS:

-> `state`: a numpy array representing the amplitude vector of
								a non-normalised quantum state

OUTPUT:

-> `True` or `False` (boolean), depending on whether the state is normalised
					(error tolerance is set to 1e-03 = 0.001)

--------------------------------------------------------------------------------

**Example 1:**

```
print(is_normalised([ 0.707,  0.,     0.,     0.707]))
```
This outputs `True`.


```
print(is_normalised([ 1, 2, 3, 4]))
```
This outputs `False`.

--------------------------------------------------------------------------------

--------------------------------------------------------------------------------

### Printing quantum states

To print the amplitudes and probabilities of a quantum state, the
general syntax is:

```
print_me(state, style)
```

INPUT:

-> `state`: a numpy array representing the amplitude vector of a quantum state

-> `style`: How to print. Option values are:

>-- `None` or no entry: prints a table of only the nonzero basis states

>-- `full`: prints all states

>-- `amplitudes`: prints the amplitude vector only

OUTPUT:

-> printout of measurements

--------------------------------------------------------------------------------

**Example 1:**

This example initialises two qubits in uniform superposition over the states
 |00> and |01> and prints with the option `None` which only prints the
 qubit states with non-zero values. Syntax:
```
print_me(create_state(2, [0, 1]), None)
```
The printed output will look like this:
```
Quantum State:
 Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
		 0            0.5  0.707+0.000j  |00>
		 1            0.5  0.707+0.000j  |01>

```

--------------------------------------------------------------------------------

**Example 2:**

This example again initialises two qubits in uniform superposition over the states
 |00> and |10> and prints it with the option value `full`. Syntax:
```
print_me(create_state(2, [0,2]), 'full')
```

This will also print qubit states with zero amplitudes:
```
Quantum State:
 Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
		 0            0.5  0.707+0.000j  |00>
		 1            0    0.000+0.000j  |01>
		 2            0.5  0.707+0.000j  |10>
		 3            0    0.000+0.000j  |11>
```
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------

### Quantum gates

Analogously to a classical computer, a quantum computer processes information
and performs quantum computations by manipulating qubits using quantum logic
gates, often just called quantum gates. Usually, a sequence of such quantum
gates is required to perform a certain task or solve a particular problem on
a quantum computer. Such a sequence of quantum gates is called
a quantum algorithm. In this library the following single-qubit quantum logic
gates are defined:

The identity or idle gate:
> eye

The three gates based on Pauli matrices:
> X, Y, Z

The Hadamard gate:
> H

Various z-rotation gates:
> S, Sdagger, T, Tdagger

Three generalized rotation gates that allow for rotations by `angle` radians
around the x, y and z axes on the Bloch sphere:
> Rx(angle), Ry(angle), Rz(angle)

--------------------------------------------------------------------------------

To apply a single-qubit quantum logic gate to a specific qubit in a quantum
state you will need the following function:

```
apply_unitary(gate, qubit_pos, state)
```

INPUT:

-> `gate`: A 2x2 unitary matrix representing a single-qubit quantum logic gate.

-> `qubit_pos`: A list of qubit positions. There are three different possibilities:

>-- [t] If only one qubit position `t` is specified, `gate` is applied to the
specified target qubit. The identity gate is applied to all other qubits in
`quantum_state`.

>-- [c1,t] If two qubit positions `c` and `t` are specified, `gate` is applied
to the qubit at position `t` if and only if the qubit `c` is in state |1>.
This, essentially, is a controlled-U operation where U can be any unitary 2x2
matrix.

>-- [c1,c2,t] If three qubit positions `c1`, `c2` and `t` are specified,
`gate` is applied to the qubit at position `t` if and only if the qubits at
positions `c1` and `c2` are in state |1>. This represents a controlled controlled
U gate. However, the current version of this library ONLY allows for a controlled
controlled X gate! This is often called Toffoli or CCNOT gate. Future versions
might include controlled controlled U operations for any 2x2 unitary matrix U.

-> `state`: a numpy array representing the amplitude vector of a quantum state


OUTPUT:

-> numpy array of amplitudes representing the unitarily transformed quantum state.

--------------------------------------------------------------------------------

**Example 1:**

We first initialise a single qubit in state |0>, then we apply an `X` gate to this
qubit and finally print the full amplitude and probability distribution of the
`quantum_state`. The syntax looks as follows:

```
quantum_state = create_state(1, [0])
apply_unitary(X, [0], quantum_state)
print_me(quantum_state, 'full')
```

The printed output will look like this:

```
Quantum State:
  Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
      0              1  1.000+0.000j  |0>
      1              0  0.000+0.000j  |1>



 Quantum State:
  Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
      0              0  0.000+0.000j  |0>
      1              1  1.000+0.000j  |1>

```
--------------------------------------------------------------------------------

**Example 2:**

We first initialise a single qubit in state |0>, then we apply a `H` gate to this
qubit and finally print the full amplitude and probability distribution of the
`quantum_state`. The syntax looks as follows:

```
quantum_state = create_state(1, [0])
print_me(quantum_state, 'full')
quantum_state = apply_unitary(H, [0], quantum_state)
print_me(quantum_state, 'full')
```

Since the H gate puts the qubit into uniform superposition over the states |0>
and |1> the printed output will look like this:

```
Quantum State:
 Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
		 0              1  1.000+0.000j  |0>
		 1              0  0.000+0.000j  |1>



Quantum State:
 Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
		 0            0.5  0.707+0.000j  |0>
		 1            0.5  0.707+0.000j  |1>
```
--------------------------------------------------------------------------------

**Example 3:**

In this example we initialise two qubits in the state |11> and print it to the
console. Next, we apply a CNOT (controlled X) gate to the `quantum_state` with
the 0th qubit being the control and the 1st qubit being the target and print
the full amplitude and probability distribution into the console. The syntax
looks like:

```
quantum_state = create_state(2, [3])
print_me(quantum_state, 'full')
quantum_state = apply_unitary(X, [0,1], quantum_state)
print_me(quantum_state, 'full')
```

Since the 0th qubit is in state |1>, the 1st qubit is flipped to the |0> state
resulting in the final `quantum_state` |01>.
In the printed output below you can see the how the CNOT gate changed the qubit
state |11> into |01>:

```
Quantum State:
 Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
		 0              0  0.000+0.000j  |00>
		 1              0  0.000+0.000j  |01>
		 2              0  0.000+0.000j  |10>
		 3              1  1.000+0.000j  |11>



Quantum State:
 Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
		 0              0  0.000+0.000j  |00>
		 1              0  0.000+0.000j  |01>
		 2              1  1.000+0.000j  |10>
		 3              0  0.000+0.000j  |11>

```
--------------------------------------------------------------------------------

**Example 4:**

In this example we again initialise two qubits in the state |11> and print it to
the console. Next, we apply a controlled x-rotation gate with `angle=-math.pi/4`
to the `quantum_state` with the 0th qubit being the control and the 1st qubit
being the target. Lastly, we print the full amplitude and probability
distribution into the console. The syntax then looks like:

```
quantum_state = create_state(2, [3])
print_me(quantum_state, 'full')
quantum_state = apply_unitary(Rx(-math.pi/4), [0,1], quantum_state)
print_me(quantum_state, 'full')
```

The printed output below shows the evolution from the |11>
state to a superposition of the |10> and |11> state since we rotated the 1st
qubit around the x-axis of the Bloch sphere.
```
Quantum State:
 Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
		 0              0  0.000+0.000j  |00>
		 1              0  0.000+0.000j  |01>
		 2              0  0.000+0.000j  |10>
		 3              1  1.000+0.000j  |11>



Quantum State:
 Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
		 0          0      0.000+0.000j  |00>
		 1          0      0.000+0.000j  |01>
		 2          0.146  0.000+0.383j  |10>
		 3          0.854  0.924+0.000j  |11>
```

--------------------------------------------------------------------------------

**Example 5:**

In this example we initialise three qubits in the state |111> and print it to
the console. Next, we apply a Toffoli (controlled controlled X or CCNOT) gate to
the `quantum_state` with the 0th and 1st qubits being the control qubits and the
2nd qubit being the target. Finally, we print the full amplitude and probability
distribution into the console. The syntax then looks like:
```
quantum_state = create_state(3, [7])
print_me(quantum_state, 'full')
quantum_state = apply_unitary(X, [0,1,2], quantum_state)
print_me(quantum_state, 'full')
```
Since the 0th and 1st qubit are both in state |1> the 2nd qubit is flipped into
state |0>. The printed output below shows the evolution of the |111> state into
the |110> state:
```
Quantum State:
 Index    Probability  Amplitude     Basis state
-------  -------------  ------------  -------------
		 0              0  0.000+0.000j  |000>
		 1              0  0.000+0.000j  |001>
		 2              0  0.000+0.000j  |010>
		 3              0  0.000+0.000j  |011>
		 4              0  0.000+0.000j  |100>
		 5              0  0.000+0.000j  |101>
		 6              0  0.000+0.000j  |110>
		 7              1  1.000+0.000j  |111>



Quantum State:
 Index    Probability  Amplitude      Basis state
-------  -------------  -------------  -------------
		 0              0  0.000-0.000j   |000>
		 1              0  0.000-0.000j   |001>
		 2              0  -0.000-0.000j  |010>
		 3              0  0.000-0.000j   |011>
		 4              0  -0.000+0.000j  |100>
		 5              0  0.000-0.000j   |101>
		 6              1  1.000-0.000j   |110>
		 7              0  -0.000+0.000j  |111>
```
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------

### Special functions

#### Grover iteration

```
grover_iteration(state, marked_pos)
```

INPUT:

-> `state`: a numpy array representing the amplitude vector of a quantum state

-> `marked_pos`: list of integers of the marked positions

--------------------------------------------------------------------------------

**Example 1:**

to be filled...
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------

### Measurement and statistics

In order to retrieve the result of any quantum computation one needs to measure
the qubits in the system. In our simulation library all qubits are measured at
the same time. In this version, measuring only a few qubits within a larger
multi-qubit system is not supported. Qubit measurements are done with the
following syntax:

```
measure(state, runs, output)
```

INPUTS:

-> `state`: a numpy array representing the amplitude vector of a quantum state

-> `runs`:  number of times the measurement has to be repeated
(simulating a repeated preparation of the state before each measurement)

-> `output`: Type of printout. Option values are:

>-- `outcomes` - prints a table with all individual (stochastic) measurement outcomes

>-- `stats` - prints the stochastic measurement outcomes for the individual qubit states

OUTPUT:

			-> prints the measurement results

--------------------------------------------------------------------------------

**Example 1:**

We initialise a uniform superposition of the states |0000>, |0010>, |0100> and
|0110> and measure the quantum state. This is repeated two times (`runs = 2`)
and by selecting the option value `outcomes` the two measurement results are
printed to the console. The syntax looks as follows:

```
measure(create_state(4, [0,2,4,6]),2, 'outcomes')
```
This will print the results of the individual qubit measurements.

```
Measurement Results:
  Run    Index  Basis state
-----  -------  -------------
    1        4  |0100>
    2        0  |0000>
```
--------------------------------------------------------------------------------

**Example 2:**

In this example, we again initialise a uniform superposition over the
states |0000>, |0010>, |0100> and |0110> and measure the state.
This is repeated 100 times (`runs = 100`) and by
selecting the option value `stats` the probability distribution over the
measurement results is printed to the console. The syntax looks as follows:

```
measure(create_state(4, [0,2,4,6]),100, 'stats')
```
This will print the retrieved statistics from the 100 simulation runs:
```
Measurement Statistics:
 rel. occ.    Index  Basis state
-----------  -------  -------------
			0.16        0  |0000>
			0.38        2  |0010>
			0.24        4  |0100>
			0.22        6  |0110>
```

## Contribute

Feel free to dive in! [Open an issue](https://github.com/mariaschuld/qc_simulator/issues/new)

If you have any suggestions on how we can simplify the documentation or improve
the simulator library, please send [Maria](mailto:mariaschuld@gmail.com) or
[Mark](mailto:markfingerhuth@protonmail.com) an email.

## License

GNU (c) Maria Schuld & Mark Fingerhuth
