'''
The main code for the project.

It became clear that for the extension part of the project it would make sense to use an OOP approach

This file defines new classes such as the quantum register, the register trajectory and the gate class.
'''




from __future__ import division
import numpy as np
from scipy import linalg as lin
import random
from tqdm import tqdm
from matplotlib import widgets
from matplotlib import pyplot as pp
import matplotlib

matplotlib.rcParams['font.serif'] = "Times New Roman"
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams["mathtext.fontset"] = "dejavuserif"

def norm(state_vector):
    b = np.conjugate(state_vector.transpose())
    a = np.dot(b, state_vector)
    return np.sqrt(np.absolute(a))[0][0]


def normalize(state_vector):
    if norm(state_vector) != 0:
        const = np.sqrt(1 / norm(state_vector))
    else:
        const = 1
        print('Divide by zero encountered in normalization, skipping..')
    return const * state_vector


def evolution_operator(hamiltonian, timestep):
    '''Neglects the reduced planck constant (at this stage). H must be a numpy array.'''
    M = -1j * timestep * hamiltonian
    return lin.expm(M)


def residuals(ydata, ycompare):
    '''Returns two arrays, the residuals and the normalised residuals'''
    resids = []
    nresids = []
    for i in range(len(ydata)):
        diff = ycompare[i] - ydata[i]
        resids.append(diff)
        nresids.append(diff / ydata[i])
    return resids, nresids


def evolve_MCWF(initial_wavefunction, evolution_op, decay_lifetime, omega, output_target=None, timestep=0.01,
                n_timesteps=500, iterations=1, decay_off=False):
    steps_per_cycle = omega / timestep
    threshold = decay_lifetime / steps_per_cycle
    if decay_off:
        threshold = 0
    for k in tqdm(range(iterations)):
        psi = initial_wavefunction
        pdata = []
        for i in range(n_timesteps):
            epsilon = random.random()
            if epsilon < threshold:
                psi[0][0] = 1 + 0j
                psi[1][0] = 0 + 0j
            else:
                psi = np.dot(evolution_op, psi)
            psi = normalize(psi)
            pdata.append(np.real(psi[1][0]) ** 2 + np.imag(psi[1][0]) ** 2)
        if output_target == None:
            pass
        else:
            output_target.append(pdata)
    return pdata


def cos_MCWF(omega, decay_rate, output_target=[], timestep=0.01, n_timesteps=500, iterations=1):
    '''Function will append output_target with an array for every iteration and return the final array.'''
    steps_per_cycle = omega / timestep
    threshold = decay_rate / steps_per_cycle
    for k in tqdm(range(iterations)):
        pdata = []
        counter = 0
        for i in range(n_timesteps):
            epsilon = random.random()
            if epsilon < threshold:
                prob = 0
                counter = 0
                # jump condition, collapse to ground state
            else:
                prob = -0.5 * np.cos(counter * timestep * omega) + 0.5
                # this is the no jump condition, non hermitian Hamiltionian followed be renormalisation
            pdata.append(prob)
            counter += 1
        output_target.append(pdata)
    return pdata


def average_WF(input_file_path, output_file_path):
    file = np.load(input_file_path)
    pdata = []
    iterations = len(file)
    time_length = len(file[0])
    for i in tqdm(range(time_length)):
        running_total = 0
        for k in range(iterations):
            running_total += file[k][i]
        value = running_total / iterations
        pdata.append(value)
    np.save(output_file_path, pdata)


def average_WF2(data):
    file = data
    pdata = []
    iterations = len(file)
    time_length = len(file[0])
    for i in tqdm(range(time_length)):
        running_total = 0
        for k in range(iterations):
            running_total += file[k][i]
        value = running_total / iterations
        pdata.append(value)
    return pdata


def stitch(lists):
    '''It appears this is entirely redundant as python does this automatically when you use + on two lists... oops'''
    result = []
    for i in lists:
        for k in i:
            result.append(k)
    return result


def unwrap(object):
    if len(object.shape) > 1:
        try:
            object = object[0]
            object = unwrap(object)
        except:
            print('Object could not be unwrapped from nested array')
    else:
        pass
    return object

def mat_to_a_power(mat, no_times):
    assert no_times >= 0
    if no_times == 0:
        return np.identity(mat.shape[0], dtype='complex128')
    elif no_times == 1:
        return mat
    else:
        result = mat
        for i in range(no_times-1):
            result = np.dot(result, mat)
        return result


def f_expm(mat, accuracy):
    mat = np.array(mat)
    total = np.zeros(mat.shape, dtype='complex128')
    for i in range(accuracy):
        prefactor = 1 / np.math.factorial(i)
        term = mat_to_a_power(mat, i)
        #print('Total:')
        #print(total)
        #print('Prefactor:')
        #print(prefactor)
        #print('New term:')
        #print(term)

        total += prefactor * term
        #print(str(i) + 'th term in expansion is:')
        #print(term)
    return total


class QRegister:
    def __init__(self, n_qubits, rabi_freq=0.5, decay_rate = 0.1, bitflip_rate = 0.1, phaseflip_rate = 0.1, time_unit=0.01, decay_off=False):
        self.n_qubits = n_qubits
        H = 0.5 * np.array([[0, rabi_freq], [rabi_freq, 0]])
        self.hamiltonian = H
        self.rabi_freq = rabi_freq
        self.time_unit = time_unit
        for i in range(n_qubits - 1):
            self.hamiltonian = np.kron(self.hamiltonian, H)
        self.zero = np.array([[1.0 + 0.0j],
                              [0.0 + 0.0j]])
        self.register = self.zero
        for i in range(n_qubits - 1):
            self.register = self.extend_reg(self.zero)
        self.bin_labels = []
        for i in range(len(self.register)):
            self.bin_labels.append(format(i, '0' + str(self.n_qubits) + 'b'))
        self.decay_rate = decay_rate
        self.bitflip_rate = bitflip_rate
        self.phaseflip_rate = phaseflip_rate
        self.decay_off = decay_off

    def extend_reg(self, new_qubit):
        return np.kron(self.register, new_qubit)

    def get_register(self):
        return self.register

    def get_labels(self):
        return self.bin_labels

    def output(self):
        for i in range(len(self.register)):
            print(self.bin_labels[i] + str(self.register[i]))

    def measure(self, position=None):
        'if position left blank, the entire reigster is measured'
        if position == None:
            probs = []
            normalize(self.register)
            for i in self.register:
                probs.append((np.real(i) ** 2 + np.imag(i) ** 2)[0])
            #print('Probability = ' + str(sum(probs)))
            state = np.random.choice(self.bin_labels, p=probs)
            for i in range(len(self.bin_labels)):
                if self.bin_labels[i] == state:
                    self.register[i] = 1
                else:
                    self.register[i] = 0
            return state

        else:
            one = 0
            zero = 0
            # add up probability from different states
            for i in range(len(self.register)):
                if self.bin_labels[i][position] == '1':
                    one += (np.real(self.register[i]) ** 2 + np.imag(self.register[i]) ** 2)
                elif self.bin_labels[i][position] == '0':
                    zero += (np.real(self.register[i]) ** 2 + np.imag(self.register[i]) ** 2)
                else:
                    print('Something very strange happened at index ' + str(i))
            # choose a state
            epsilon = random.random()
            # set states we know we're not in to zero
            #print('epsilon:' + str(epsilon) + ' prob1:' + str(one) + ' prob0:' + str(zero))
            if epsilon < one:
                result = 1
                for i in range(len(self.register)):
                    if self.bin_labels[i][position] == '0':
                        self.register[i] = 0
            elif epsilon >= one:
                result = 0
                for i in range(len(self.register)):
                    if self.bin_labels[i][position] == '1':
                        self.register[i] = 0
            # normalise the register WF
            normalize(self.register)
            return result

    def evolve(self):
        '''A function designed to find the register state one timestep in the future'''
        # define evolution operator from hamiltonian
        evop = evolution_operator(self.hamiltonian, self.time_unit)
        # calculate a probability threshold
        steps_per_cycle = self.rabi_freq / self.time_unit
        threshold = self.decay_rate / steps_per_cycle
        bitflip_thresh = self.bitflip_rate / steps_per_cycle
        phaseflip_thresh = self.phaseflip_rate / steps_per_cycle
        if self.decay_off:
            threshold = 0
            bitflip_thresh = 0
            phaseflip_thresh = 0
        decayed = False
        for i in range(self.n_qubits):
            # for each qubit check if it has decayed using RNG
            epsilon = random.random()
            if epsilon < threshold:
                # if it has decayed call decay process and record that it's happened
                self.decay(i)
                decayed = True
                #print('Decayed!')
            epsilon = random.random()
            if epsilon < bitflip_thresh:
                # check for and implement bitflip error
                self.bitflip(i)
                decayed = True
                #print('Flipped a bit!')
            epsilon = random.random()
            if epsilon < phaseflip_thresh:
                # check for and implement phaseflip error
                self.phaseflip(i)
                decayed = True
                #print('Flipped a phase!')
        if not decayed:
            # system evolved normally ONLY if no decay processes have happened
            self.register = np.dot(evop, self.register)
            self.register = normalize(self.register)
        return self.register

    ### IS THIS IMPLIMENTATION CORRECT??? SHOULDN'T THE EVOLUTION STACK WITH THE DECAY IN THE SAME TIMESTEP RATHER THAN BE REPLACED??

    def decay(self, position):
        '''A function that causes a certain qubit in the register to decay'''
        for i in range(len(self.bin_labels)):
            # for every state in the register check if the qubit in position of interest is excited
            if self.bin_labels[i][position] == '1':
                # if the qubit of interest is found to be excited, the probability is set to zero
                self.register[i] = 0.0 + 0.0j
        # register is normalized to ensure total probability of 1
        self.register = normalize(self.register)

    def bitflip(self, position):
        '''Executes a bitflip operation on a qubit in a given position'''
        #print('flipped a bit from :')
        #self.output()
        flipper = Gate('Pauli_x', position=position, n_qubits=self.n_qubits, lead_time=0, out_time=0, duration=self.time_unit)
        flipper.apply(self)
        #print('to: ')
        #self.output()

    def phaseflip(self, position):
        '''Executes a phaseflip operation on a qubit in a given position'''
        flipper = Gate('Pauli_z', position=position, n_qubits=self.n_qubits, lead_time=0, out_time=0, duration=self.time_unit)
        flipper.apply(self)




class Gate:
    def __init__(self, identifier, position, n_qubits, duration=0.2, lead_time=0.1, out_time=0.1, gate_param=0,
                 timestep=0.01):
        I_2 = np.array([[1, 0], [0, 1]])

        self.identifier = identifier
        self.position = position
        self.label = identifier
        self.lead_time = lead_time
        self.out_time = out_time
        self.timestep = timestep
        gatestring = ''  # indicates which qubit the gate should be applied to
        for i in range(n_qubits):
            if i == position:
                gatestring += '1'
            else:
                gatestring += '0'
        self.gatestring = gatestring
        self.duration = duration
        self.n_qubits=n_qubits
        applied_iterations = int(duration / timestep)
        self.gparam = gate_param / applied_iterations

        if identifier == 'Pauli_x':
            # print('Pauli x selected')
            base_mat = np.array([[0, 1], [1, 0]]) / applied_iterations
        elif identifier == 'Pauli_y':
            # print('Pauli y selected')
            base_mat = np.array([[0, -1j], [1j, 0]]) / applied_iterations
        elif identifier == 'Pauli_z':
            # print('Pauli z selected')
            base_mat = np.array([[1, 0], [0, -1]]) / applied_iterations
        elif identifier == 'Hadamard':
            # print('Hadamard selected')
            base_mat = np.array([[1, 1], [1, -1]]) / np.sqrt(2) * applied_iterations
        elif identifier == 'Rx':
            # print('Rx selected')
            base_mat = np.array([[np.cos(gate_param / (2 * applied_iterations)),
                                  -1j * np.sin(gate_param / (2 * applied_iterations))],
                                 [-1j * np.sin(gate_param / (2 * applied_iterations)),
                                  np.cos(gate_param / (2 * applied_iterations))]])
        elif identifier == 'Ry':
            # print('Ry selected')
            base_mat = np.array([[np.cos(gate_param / (2 * applied_iterations)),
                                  -1 * np.sin(gate_param / (2 * applied_iterations))],
                                 [np.sin(gate_param / (2 * applied_iterations)),
                                  np.cos(gate_param / (2 * applied_iterations))]])
        elif identifier == 'Rz':
            # print('Rz selected')
            base_mat = np.array([[np.exp(-0.5j * gate_param), 0],
                                 [0, np.exp(0.5j * gate_param)]])
        elif identifier == 'Identity':
            base_mat = np.array([[1, 0], [0, 1]])

            ##can add more types of base_mat in elifs here

        if gatestring[0] == '1':
            self.matrix = base_mat
        elif gatestring[0] == '0':
            self.matrix = I_2
        for i in gatestring[1:]:
            if i == '0':
                self.matrix = np.kron(self.matrix, I_2)
            elif i == '1':
                self.matrix = np.kron(self.matrix, base_mat)

    def __repr__(self):
        return self.label

    def get_gatestring(self):
        return self.gatestring

    def apply(self, register, disable_inout=False):
        '''Returns a trajectory slice which results from the action of the gate on a register'''
        trajectory = []
        without = []
        if self.lead_time == 0 or disable_inout:
            pass
        else:
            for i in range(int(self.lead_time / self.timestep)):
                register.evolve()
                trajectory.append(register.register)
        for i in range(int(self.duration / self.timestep)):
            register.register = np.dot(self.matrix, register.register)
            register.register = normalize(register.register)
            register.evolve()
            trajectory.append(register.register)
            without.append(register.register)
        if self.out_time == 0 or disable_inout:
            pass
        else:
            for i in range(int(self.out_time / self.timestep)):
                register.evolve()
                trajectory.append(register.register)
        return trajectory, without

    def flash(self, register):
        '''Returns the effect of a gate on a register without changing the register'''
        a = np.dot(self.matrix, register.register)
        a = normalize(a)
        return a





class BCGate(Gate):
    def __init__(self, identifier, positions, n_qubits, duration=0.2, lead_time=0.1, out_time=0.1, gate_param=0,
                 timestep=0.01):
        I_2 = np.array([[1, 0], [0, 1]])

        self.identifier = identifier
        self.positions = positions
        self.label = identifier
        self.lead_time = lead_time
        self.out_time = out_time
        self.timestep = timestep
        gatestring = ''  # indicates which qubit(s) the gate should be applied to
        for i in range(n_qubits):
            if i in positions:
                gatestring += '1'
            else:
                gatestring += '0'
        self.gatestring = gatestring
        self.duration = duration
        self.n_qubits=n_qubits
        applied_iterations = int(duration / timestep)
        self.gparam = gate_param / applied_iterations

        # define gate types (Remember t parameter in my code needs to include the pi/2 not included in the notes)
        if identifier == 'Ising':
            # print('Ising selected')
            self.matrix = self.generate_canonical_binary(self.gparam, 0, 0)
        elif identifier == 'CNOT':
            # print('CNOT selected')
            self.matrix = self.generate_canonical_binary(0.25*np.pi, 0, 0)
        elif identifier == 'SWAP':
            self.matrix = self.generate_canonical_binary(0.25*np.pi, 0.25*np.pi, 0.25*np.pi)
        elif identifier == 'iSWAP':
            self.matrix = self.generate_canonical_binary(0.25*np.pi, 0.25*np.pi, 0)
        elif identifier == 'rtSWAP':
            self.matrix = self.generate_canonical_binary(0.125*np.pi, 0.125*np.pi, 0.125*np.pi)
        elif identifier == 'QFT':
            self.matrix = self.generate_canonical_binary(0.25*np.pi, 0.25*np.pi, 0.125*np.pi)

    def generate_canonical_binary(self, tx, ty, tz):
        """Produces the register specific block matrix required for any binary gate acting on arbitrary qubits"""
        # create the basic matrices needed for a canonical gate
        pauli_x = np.array([[0j, 1+0j], [1+0j, 0j]])
        pauli_y = np.array([[0j, -1j], [1j, 0j]])
        pauli_z = np.array([[1+0j, 0j], [0j, -1+0j]])
        I_2 = np.array([[1, 0j], [0j, 1]])
        swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])

        # extract the positions of the relevant qubits
        first_qubit_pos = self.positions[0]
        second_qubit_pos = self.positions[1]
        print('Position Info:')
        print(self.gatestring)
        print(first_qubit_pos)
        print(second_qubit_pos)

        # check if adjacent
        if second_qubit_pos - first_qubit_pos == 0:
            places_shifted = 0
        else:
            places_shifted = second_qubit_pos - first_qubit_pos

        # move the second qubit to be adjacent to the first
        adjacent_maker = np.identity(2 ** self.n_qubits)
        for i in range(places_shifted-1):
            swapper = np.identity(2 * first_qubit_pos + 2)
            for m in range(places_shifted - 2 - i):
                swapper = np.kron(swapper, I_2)
            swapper = np.kron(swapper, swap)
            for k in range(self.n_qubits - second_qubit_pos - 2 + i):
                swapper = np.kron(swapper, I_2)
            adjacent_maker = np.dot(adjacent_maker, swapper)

        # apply canonical gate on the adjacent qubits
        applied = False
        if self.gatestring[0] == '1':
            matrix = tx * np.kron(pauli_x, pauli_x) + ty * np.kron(pauli_y, pauli_y) + tz * np.kron(pauli_z, pauli_z)
            for i in matrix:
                for k in i:
                    k *= -1j
            gmatrix = f_expm(matrix, 15)
            applied = True
        elif self.gatestring[0] == '0':
            gmatrix = I_2
        for i in range(1, len(self.gatestring)):
            if self.gatestring[i] == '1' and not applied:
                matrix = tx * np.kron(pauli_x, pauli_x) + ty * np.kron(pauli_y, pauli_y) + tz * np.kron(pauli_z,
                                                                                                         pauli_z)
                for i in matrix:
                    for k in i:
                        k *= -1j
                matrix = f_expm(matrix, 15)
                applied = True
                gmatrix = np.kron(gmatrix, matrix)
            elif self.gatestring[i] == '0':
                gmatrix = np.kron(gmatrix, I_2)

        # swap the second qubit back to its original position
        adjacent_unmaker = np.identity(2 ** self.n_qubits)
        for i in range(places_shifted-1):
            swapper = np.identity(2 * first_qubit_pos + 2)
            for m in range(i):
                swapper = np.kron(swapper, I_2)
            swapper = np.kron(swapper, swap)
            for k in range(self.n_qubits - second_qubit_pos - 2 - i):
                swapper = np.kron(swapper, I_2)
            print('i:' + str(i))
            print(swapper.shape)
            print(adjacent_unmaker.shape)
            adjacent_unmaker = np.dot(adjacent_unmaker, swapper)

        # combine these matrices
        result = np.dot(adjacent_maker, gmatrix)
        result = np.dot(result, adjacent_unmaker)

        return result

    def CB_matrix_ver2(self, tx, ty, tz):
        return None

class Adjacent_Binary_Gate(Gate):
    def __init__(self, identifier, positions, n_qubits, duration=0.2, lead_time=0.1, out_time=0.1, gate_param=0,
                 timestep=0.01):
        I_2 = np.array([[1, 0], [0, 1]])

        self.identifier = identifier
        self.positions = positions
        self.label = identifier
        self.lead_time = lead_time
        self.out_time = out_time
        self.timestep = timestep
        gatestring = ''  # indicates which qubit(s) the gate should be applied to
        for i in range(n_qubits):
            if i in positions:
                gatestring += '1'
            else:
                gatestring += '0'
        self.gatestring = gatestring
        self.duration = duration
        self.n_qubits=n_qubits
        applied_iterations = int(duration / timestep)
        self.gparam = gate_param / applied_iterations

        # define gate types (Remember t parameter in my code needs to include the pi/2 not included in the notes)
        if identifier == 'Ising':
            # print('Ising selected')
            gate_mat = np.array([[np.cos(self.gparam), 0, 0, -1j*np.sin(self.gparam)],
                                 [0, np.cos(self.gparam), -1j*np.sin(self.gparam), 0],
                                 [0, -1j*np.sin(self.gparam), np.cos(self.gparam), 0],
                                 [-1j*np.sin(self.gparam), 0, 0, np.cos(self.gparam)]])
        elif identifier == 'CNOT':
            gate_mat = np.array([[1, 0, 0, 0],
                                 [0, 1, 0, 0],
                                 [0, 0, 0, 1],
                                 [0, 0, 1, 0]])
        elif identifier == 'SWAP':
            gate_mat = np.array([[1, 0, 0, 0],
                                 [0, 0, 1, 0],
                                 [0, 1, 0, 0],
                                 [0, 0, 0, 1]])

        if gatestring[0] == '1':
            self.matrix = gate_mat
            applied = True
        elif gatestring[0] == '0':
            self.matrix = I_2
            applied = False
        for i in range(1, len(gatestring)):
            if gatestring[i] == '1' and not applied:
                self.matrix = np.kron(self.matrix, gate_mat)
                applied = True
            elif gatestring[i] == '0':
                self.matrix = np.kron(self.matrix, I_2)



class Trajectory:
    def __init__(self, gate_list, register, success_condition=None):
        '''If checking for success give a function to the success condition agrument that
        returns True if the operation was performed correctly.'''
        lists = []
        self.register = register
        self.gate_list = gate_list
        self.labels = []
        self.timestep = gate_list[0].timestep
        self.change_times = [[0.0, register.register]]
        self.region_boundaries = []
        self.failed = True
        time = 0.0
        for i in gate_list:
            assert i.timestep == self.timestep, 'Timesteps must be equal for all stages of evolution'
            full_trajectory, gate_only = i.apply(register)
            lists.append(full_trajectory)
            duration = i.lead_time + i.duration + i.out_time
            time += duration
            gate_t = time - i.out_time
            self.region_boundaries.append([[gate_t - i.duration, gate_t], gate_only])
            self.change_times.append([time, register.register])
            self.labels.append([gate_t - i.duration / 2, i.label])
        self.trajectory = np.array(stitch(lists)) ##### THIS IS A POTENTIAL AREA FOR EFFICIENCY IMPROVEMENT
        self.p_trajectory = self.probabilities()
        self.pauli_matrices = []
        for i in range(self.register.n_qubits):
            self.pauli_matrices.append([Gate('Pauli_x', position=i, n_qubits=self.register.n_qubits, duration=0.01),
                                        Gate('Pauli_y', position=i, n_qubits=self.register.n_qubits, duration=0.01),
                                        Gate('Pauli_z', position=i, n_qubits=self.register.n_qubits, duration=0.01)])
        if success_condition is None:
            pass
        elif success_condition(self.register):
            self.failed = False

        # self.p_trajectory = self.restructure(self.p_trajectory)

    def probabilities(self):
        probs = []
        for i in self.trajectory:
            probs.append(np.real(i) ** 2 + np.imag(i) ** 2)
        return probs

    def restructure(self, list_of_lists):
        new_list = []
        for i in list_of_lists[0]:
            new_list.append([])
        for i in list_of_lists:
            for k in range(len(i)):
                new_list[k].append(i[k])
        return new_list

    def plot(self, enable_axspan=True, enable_boundaries=True, rotate_labels=False, fig_width=3.4):
        labels = self.register.get_labels()
        iterations = len(self.p_trajectory)
        timedata = self.timestep * np.array(range(iterations))
        fig = pp.figure(figsize=(fig_width, 0.6*self.register.n_qubits**2), dpi=300)
        fig.subplots(1, len(self.p_trajectory[0]), sharex=True)
        grid = pp.GridSpec(len(self.p_trajectory[0]), 1, hspace=0)
        self.p_trajectory = np.array(self.p_trajectory)
        for i in range(len(self.p_trajectory[0])):
            ax = pp.subplot(grid[i, 0])
            if enable_axspan:
                for k in self.region_boundaries:
                    pp.axvspan(k[0][0]-0.01, k[0][1]-0.01, alpha=0.5, color=(0.8, 0.8, 0.8))
            if enable_boundaries:
                for k in self.change_times:
                    pp.plot([k[0], k[0]], [-0.2, 1.2], linestyle='--', color=(0.8, 0.8, 0.8))
            ax.plot(timedata, self.p_trajectory[:, i])
            pp.ylim(-0.15, 1.15)
            pp.ylabel(labels[i], rotation=90, labelpad=-0.2)
            ax.get_xaxis().set_visible(False)
        swap_stacked = False
        for i in self.labels:
            coord = ax.transData.transform([i[0], (2 ** (self.register.n_qubits - 2)) * 5.2])
            inv = fig.transFigure.inverted()
            coord = inv.transform(coord)
            if i[1] == 'SWAP':
                print('we got a swap')
                if swap_stacked:
                    print('we passed')
                    pass
                else:
                    if rotate_labels:
                        pp.gcf().text(coord[0], coord[1] + 0.005, 'SWAP(s)', horizontalalignment='center', rotation=90)
                    else:
                        pp.gcf().text(coord[0], coord[1], 'SWAP(s)', horizontalalignment='center')
                    print('we executed this block')
                    swap_stacked = True
                    print(swap_stacked)
            else:
                if rotate_labels:
                    pp.gcf().text(coord[0], coord[1] + 0.005, i[1], horizontalalignment='center', rotation=90)
                else:
                    pp.gcf().text(coord[0], coord[1], i[1], horizontalalignment='center')
                print('This also got executed')
                swap_stacked = False
        ax.get_xaxis().set_visible(True)
        pp.xlabel('Time / 100 time-steps')
        pp.gcf().text(0.01, 0.5, 'Probability of Each State', rotation=90, verticalalignment='center')
        pp.show()

    def get_bloch_coords(self, reg_state):
        coord_set = []
        for k in self.pauli_matrices:
            a = reg_state
            b = np.conjugate(a)
            coord = []
            for i in k:
                thing = np.dot(np.transpose(b), np.dot(i.matrix, a))
                coord.append(np.sign(thing) * norm(thing))
            coord_set.append(coord)
        return coord_set

    def bloch_plot(self, enable_gates=False):

        r = 1
        pi = np.pi
        cos = np.cos
        sin = np.sin
        phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0 * pi:100j]
        x = r * sin(phi) * cos(theta)
        y = r * sin(phi) * sin(theta)
        z = r * cos(phi)
        cy1 = 0 * sin(theta)
        cz1 = r * sin(theta)
        cx1 = cos(theta)
        cx3 = 0 * sin(theta)
        cy3 = cos(theta)
        cx2 = r * cos(2 * phi)
        cy2 = r * sin(2 * phi)

        if enable_gates:
            grid = pp.GridSpec(self.register.n_qubits, len(self.region_boundaries), hspace=0, wspace=0)
            fig = pp.figure(figsize=(8 * len(self.region_boundaries), 4 * self.register.n_qubits))
            for i in range(len(self.region_boundaries)):
                points = []
                for l in range(len(self.region_boundaries[i][1])):
                    points.append(self.get_bloch_coords(self.region_boundaries[i][1][l]))
                for k in range(self.register.n_qubits):
                    ax = pp.subplot(grid[k, i], projection='3d')
                    ax.plot_surface(x, y, z, alpha=0.1)
                    ax.plot(xs=cx1.flatten(), ys=cy1.flatten(), zs=cz1.flatten(), color=(0.8, 0.8, 0.8), linewidth=0.1,
                            alpha=1)
                    ax.plot(xs=cx3.flatten(), ys=cy3.flatten(), zs=cz1.flatten(), color=(0.8, 0.8, 0.8), linewidth=0.1,
                            alpha=1)
                    ax.plot(cx2, cy2, color=(0.8, 0.8, 0.8), linewidth=0.1, alpha=1)
                    ax.plot([0, 0], [0, 0], [-1.2, 1.2], color=(0.9, 0.9, 0.9))
                    ax.quiver(0, 0, -1.2, 0, 0, 2.6, color=(0.9, 0.9, 0.9), arrow_length_ratio=0.1)
                    ax.plot([0, 0], [-1.2, 1.2], [0, 0], color=(0.9, 0.9, 0.9))
                    ax.plot([-1.2, 1.2], [0, 0], [0, 0], color=(0.9, 0.9, 0.9))
                    pp.title(self.labels[i][1])
                    for l in range(len(points)):
                        ax.quiver(0, 0, 0, points[l][k][0], points[l][k][1], points[l][k][2], arrow_length_ratio=0.1,
                                  alpha=((0.8 * l) / (len(points)) + 0.1))
                    ax.set_axis_off()
        else:
            grid = pp.GridSpec(self.register.n_qubits, len(self.change_times), hspace=0, wspace=0)
            fig = pp.figure(figsize=(14, 4))
            for i in range(len(self.change_times)):
                points = self.get_bloch_coords(self.change_times[i][1])
                for k in range(len(points)):
                    ax = pp.subplot(grid[k, i], projection='3d')

                    ax.plot_surface(x, y, z, alpha=0.1)
                    ax.plot([0, 0], [0, 0], [-1.2, 1.2], color=(0.9, 0.9, 0.9))
                    ax.quiver(0, 0, -1.2, 0, 0, 2.6, color=(0.9, 0.9, 0.9), arrow_length_ratio=0.1)
                    ax.plot([0, 0], [-1.2, 1.2], [0, 0], color=(0.9, 0.9, 0.9))
                    ax.plot([-1.2, 1.2], [0, 0], [0, 0], color=(0.9, 0.9, 0.9))
                    ax.quiver(0, 0, 0, points[k][0][0], points[k][0][1], points[k][1][2], arrow_length_ratio=0.1)
                    ax.set_axis_off()

                    '''
                    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
                    ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
                    ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
                    ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
                    '''
        pp.show()


class AvTrajectory(Trajectory):
    def __init__(self, trajectories):
        self.register = trajectories[0].register
        self.timestep = trajectories[0].timestep
        self.trajectory = []
        self.region_boundaries = trajectories[0].region_boundaries
        self.change_times = trajectories[0].change_times
        self.labels = trajectories[0].labels
        #print(len(trajectories[0].trajectory))
        #print(trajectories[0].trajectory[:][0].shape)
        #print(trajectories[0].trajectory[:][0])
        #print(trajectories[0].trajectory[:][1])
        #print(trajectories[0].trajectory[:][218])
        big_bastard = []
        self.number_fails = 0
        self.number_successes = 0
        for k in range(len(trajectories)):
            big_bastard.append(trajectories[k].probabilities())
            if trajectories[k].failed:
                self.number_fails += 1
            else:
                self.number_successes += 1
        prob_traj = []
        error_traj = []
        for k in range(len(big_bastard[0])):
            prob_step = []
            error_step = []
            for l in range(len(big_bastard[0][0])):
                running_tot = 0
                datapoint_vals = []
                for i in range(len(big_bastard)):
                    running_tot += big_bastard[i][k][l]
                    datapoint_vals.append(big_bastard[i][k][l])
                running_tot /= len(big_bastard)
                error = np.std(datapoint_vals)
                error /= np.sqrt(len(datapoint_vals))
                #print(str(datapoint_vals) + str(error))
                prob_step.append(running_tot)
                error_step.append(error)
            #print(error_step)
            prob_traj.append(prob_step)
            error_traj.append(error_step)
        #print(error_traj)
        self.p_trajectory = np.array(prob_traj)
        self.errors = error_traj
        self.fidelity_error = 0
        #print(self.p_trajectory)

    def mcwf_plot(self):
        iterations = len(self.p_trajectory)
        timedata = self.timestep * np.array(range(iterations))
        pp.figure(figsize=(3.4, 2.5), dpi=300)
        self.p_trajectory = np.array(self.p_trajectory)
        upper = []
        lower = []
        for i in range(len(self.errors[1])):
            upper.append(self.p_trajectory[1][i] + self.errors[1][i])
            lower.append(self.p_trajectory[1][i] + self.errors[1][i])
        for i in range(len(upper)):
            pp.fill_between(timedata, upper[i], lower[i], alpha=0.2)
        pp.plot(timedata, self.p_trajectory[:, 1])
        pp.ylim(-0.05, 1.05)
        pp.xlabel('Time / 100 time-steps')
        pp.ylabel('Probability of Excited State')
        pp.show()

    def fidelity_stats(self):
        total_events = self.number_successes + self.number_fails
        fidelity = self.number_successes / total_events
        failure_rate = self.number_fails / total_events
        things = []
        for i in range(self.number_fails):
            things.append(0)
        for i in range(self.number_successes):
            things.append(1)
        std = np.std(things)
        error = std / np.sqrt(len(things))
        self.fidelity_error = error
        print('Fidelity statistics (' + str(self.register.decay_rate) + '-' + str(self.register.bitflip_rate) + '-' + str(self.register.phaseflip_rate) + '):')
        print('Success rate: ' + str(fidelity) + '+/-' + str(error))
        print('Failure rate: ' + str(failure_rate))
        return fidelity, error

    def get_errors(self):
        return self.errors
