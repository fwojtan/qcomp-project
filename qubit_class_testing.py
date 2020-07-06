from futils import*
from mpl_toolkits.mplot3d import Axes3D

np.set_printoptions(threshold=np.inf)
psi = np.array([[1.0 +0.0j],
                [0.0+0.0j]])
psi2 = np.array([[0.0 +0.0j],
                 [1.0+0.0j]])
#I_2 = np.array([[1, 0], [0, 1]])
NOT = np.array([[0, 1], [1, 0]])


n_qubits = 2
H = 0.5 * np.array([[0, 0.5], [0.5, 0]])
# define gates

flip1 = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=np.pi, out_time=0.1, lead_time=0.1)
flip2 = Gate('Ry', position=1, n_qubits=n_qubits, gate_param=np.pi, out_time=0.3, lead_time=0.1)
Ry = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=(np.pi/2), out_time=0.1, lead_time=0.3)
XX = BCGate('Ising', positions=[0,1], n_qubits=n_qubits, gate_param=(np.pi/4), out_time=0.1, lead_time=0.1)
negRx = Gate('Rx', position=1, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
negRy = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
negRz = Gate('Rz', position=0, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
Rynegpi1 = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=(-np.pi), out_time=0.1, lead_time=0.3)
Ry2 = Gate('Ry', position=1, n_qubits=n_qubits, gate_param=(np.pi/2), out_time=0.1, lead_time=0.3)
Rynegpi2 = Gate('Ry', position=1, n_qubits=n_qubits, gate_param=(-np.pi), out_time=0.1, lead_time=0.3)

UGate = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=(-np.pi), out_time=0.1, lead_time=0.3, duration=1)

UGate.matrix = np.array([[1, 1, 0, 0],
                         [0, 0, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])



H1 = [Rynegpi1, Ry]
H2 = [Rynegpi2, Ry2]
U = [UGate]
CNOT12 = [Ry, XX, negRy, negRx, negRz]
Deutsch = H1 + H2 + U + H1

def CNOT_success(register):
    if register.measure(1) == 1:
        return True

test_range = np.linspace(0, 0.5, 10)
n_iterations = 500

dec_fidelities = []
for i in test_range:
    trajectories = []
    for k in range(n_iterations):
        reg = QRegister(n_qubits=n_qubits, rabi_freq=0.5, decay_rate=i, bitflip_rate=0.0, phaseflip_rate=0.0)
        flip1.apply(reg)
        cnot_traj = Trajectory(CNOT12, reg, CNOT_success)
        trajectories.append(cnot_traj)
    average = AvTrajectory(trajectories)
    fidelity = average.fidelity_stats()
    dec_fidelities.append(fidelity)

bit_fidelities = []
for i in test_range:
    trajectories = []
    for k in range(n_iterations):
        reg = QRegister(n_qubits=n_qubits, rabi_freq=0.5, decay_rate=0.0, bitflip_rate=i, phaseflip_rate=0.0)
        flip1.apply(reg)
        cnot_traj = Trajectory(CNOT12, reg, CNOT_success)
        trajectories.append(cnot_traj)
    average = AvTrajectory(trajectories)
    fidelity = average.fidelity_stats()
    bit_fidelities.append(fidelity)

pha_fidelities = []
for i in test_range:
    trajectories = []
    for k in range(n_iterations):
        reg = QRegister(n_qubits=n_qubits, rabi_freq=0.5, decay_rate=0.0, bitflip_rate=0.0, phaseflip_rate=i)
        flip1.apply(reg)
        cnot_traj = Trajectory(CNOT12, reg, CNOT_success)
        trajectories.append(cnot_traj)
    average = AvTrajectory(trajectories)
    fidelity = average.fidelity_stats()
    pha_fidelities.append(fidelity)

pp.figure()
pp.plot(test_range, dec_fidelities, label='Decay')
pp.plot(test_range, bit_fidelities, label='Bitflip')
pp.plot(test_range, pha_fidelities, label='Phaseflip')
pp.legend()
pp.xlabel('Decay Rate')
pp.ylabel('Fidelity')
pp.show()

np.save('Data/fidelity_v2.npy', [dec_fidelities, bit_fidelities, pha_fidelities])


# initialize register
#flip1.apply(register=myReg)
#flip2.apply(register=myReg)


#results.bloch_plot(enable_gates=True)

# only remaining issue to resolve - how to get two qubit gates to operate on non-adjacent qubits DONE
# then only improvements to make gate definitions concise and easy to use IN PROG
# THEN after implementing algorithms look at reintroducing monte carlo decays to see effect on gate fidelity







'''
H = 0.5 * np.array([[0, 0.5], [0.5, 0]])
H2 = 0.5 * np.array([[1, 0.5], [0.5, 1]])

myReg = QRegister(2, H)
compReg = QRegister(2,H)

# test_gate = Gate('Pauli_y', positions=[0], n_qubits=2, duration=1)
# gate2 = Gate('Pauli_x', positions=[0], n_qubits=2, duration=1)

# Define Some Useful Gates
Ry = Gate('Ry', positions=[0,1], n_qubits=2, duration=15, gate_param=(np.pi / 2))
negRy = Gate('Ry', positions=[0], n_qubits=2, duration=15, gate_param=(-np.pi / 2))
negRx = Gate('Rx', positions=[1], n_qubits=2, duration=15, gate_param=(-np.pi / 2))
Rz = Gate('Rz', positions=[0], n_qubits=2, duration=15, gate_param=(np.pi / 2))
Rx = Gate('Rx', positions=[1], n_qubits=2, duration=15, gate_param=(np.pi / 2))
negRz = Gate('Rz', positions=[0], n_qubits=2, duration=15, gate_param=(np.pi/2))
negRx.matrix = np.linalg.inv(Rx.matrix)
negRy.matrix = np.linalg.inv(Ry.matrix)
negRz.matrix = np.linalg.inv(Rz.matrix)
negPiRx = Gate('Rx', positions=[0,1], n_qubits=2, duration=15, gate_param=(-np.pi))
#negRy_with_negRx = Gate() ###NEED TO SOLVE MULTI SIMUL GATES
Ising = Gate('Ising', positions=[0], n_qubits=2, duration=1, gate_param=(np.pi / 4))
flip1 = Gate('Rx', positions=[0], n_qubits=2, duration=15, gate_param=(np.pi))
flip2 = Gate('Rx', positions=[1], n_qubits=2, duration=15, gate_param=np.pi)
CNOT = Gate('CNOT', positions=[0], n_qubits=2, duration=1)

tinit = myReg.evolve_with_gate(flip1, gate_times=[0.5], iterations=50, timestep=0.01)
t0 = myReg.evolve_with_gate(Ry, gate_times=[0.5], iterations=50, timestep=0.01)
t1 = myReg.evolve_with_gate(Ising, gate_times=[0.5], iterations = 50, timestep=0.01)
t2 = myReg.evolve_with_gate(negRy, gate_times=[0.5], iterations = 50, timestep=0.01)
t3 = myReg.evolve_with_gate(negRx, gate_times=[0.5], iterations = 50, timestep=0.01)

t4 = myReg.evolve_with_gate(negRz, gate_times=[0.5], iterations = 50, timestep=0.01)

tc1 = compReg.evolve_with_gate(flip1, gate_times=[0.2], iterations=50, timestep=0.01)
#tc3 = compReg.evolve_with_gate(flip2, gate_times=[0.2], iterations=50, timestep=0.01)
tc2 = compReg.evolve_with_gate(CNOT, gate_times=[0.2], iterations=50, timestep=0.01)

### SO ISSUE SEEMS TO BE Ry and negRy are not inverse operations

t_comparison = stitch([tc1, tc2])

t_combo = stitch([tinit, t0, t1, t2, t3, t4])

tdata = np.load('Data/timedata3000.npy')


def plot_trajectory(trajectory, timedata):
    labels = myReg.get_labels()
    iterations = len(trajectory)
    pp.subplots(1, len(trajectory[0]), sharex=True)
    grid = pp.GridSpec(len(trajectory[0]), 1, hspace=0)
    trajectory = np.array(trajectory)
    for i in range(len(trajectory[0])):
        pp.subplot(grid[i, 0])
        pp.plot(timedata[:iterations], trajectory[:,i])
        pp.ylim(-0.15, 1.15)
        pp.ylabel(labels[i])
    pp.xlabel('Time')
    pp.show()


plot_trajectory(t_comparison, tdata)
plot_trajectory(t_combo, tdata)

print(myReg.output())

# CNOT Output:
# 00 -> 00
# 01 -> 01
# 10 -> 11
# 11 -> 10

# Ising Output:
# 00 -> rt2/2 |00> + (i-1)/2 |11>
# 01 -> rt2/2 |01> -i rt2/2 |10>
# 10 -> rt2/2 |10> -i rt2/2 |01>
# 11 -> rt2/2 |11> -(i+1)/2 |00>

### Is it possible to reverse engineer the hamiltonians from gate matrices? How do we go from instantaneous gates to continuous ones??
'''