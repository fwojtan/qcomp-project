from futils import*
import multiprocessing as mp

max_reg_size = 6

n_qubits = 3

# initialize a bunch of gates
Ry = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=(np.pi/2), out_time=0.1, lead_time=0.1)
XX = Adjacent_Binary_Gate('Ising', positions=[0,1], n_qubits=n_qubits, gate_param=(np.pi/4), out_time=0.1, lead_time=0.1)
negRx = Gate('Rx', position=1, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
negRy = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
negRz = Gate('Rz', position=0, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)

back_Ry2 = Gate('Ry', position=1, n_qubits=n_qubits, gate_param=(np.pi / 2), out_time=0.1, lead_time=0.1)
back_negRx2 = Gate('Rx', position=0, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)
back_negRy2 = Gate('Ry', position=1, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)
back_negRz2 = Gate('Rz', position=1, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)

swap12 = Adjacent_Binary_Gate('SWAP', positions=[0, 1], n_qubits=n_qubits, duration=0.01, lead_time=0.01, out_time=0.01)

CNOT12 = [Ry, XX, negRy, negRx, negRz]
CNOT21 = [back_Ry2, XX, back_negRy2, back_negRx2, back_negRz2]


Ry2 = Gate('Ry', position=1, n_qubits=n_qubits, gate_param=(np.pi/2), out_time=0.1, lead_time=0.1)
XX2 = Adjacent_Binary_Gate('Ising', positions=[1,2], n_qubits=n_qubits, gate_param=(np.pi/4), out_time=0.1, lead_time=0.1)
negRx2 = Gate('Rx', position=2, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
negRy2 = Gate('Ry', position=1, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
negRz2 = Gate('Rz', position=1, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)

back_Ry3 = Gate('Ry', position=2, n_qubits=n_qubits, gate_param=(np.pi / 2), out_time=0.1, lead_time=0.1)
back_negRx3 = Gate('Rx', position=1, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)
back_negRy3 = Gate('Ry', position=2, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)
back_negRz3 = Gate('Rz', position=2, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)

swap23 = Adjacent_Binary_Gate('SWAP', positions=[1, 2], n_qubits=n_qubits, duration=0.01, lead_time=0.01, out_time=0.01)

CNOT23 = [Ry2, XX2, negRy2, negRx2, negRz2]
CNOT13_v2 = [swap12, Ry2, XX2, negRy2, negRx2, negRz2, swap12]
CNOT31 = [swap12, back_Ry3, XX2, back_negRy3, back_negRx3, back_negRz3, swap12]
CNOT32 = [back_Ry3, XX2, back_negRy3, back_negRx3, back_negRz3]


flip1 = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=np.pi, out_time=0.1, lead_time=0.1)
flip2 = Gate('Ry', position=1, n_qubits=n_qubits, gate_param=np.pi, out_time=0.1, lead_time=0.1)


gates3 = CNOT12 + CNOT23 + CNOT31 + CNOT32


# write a success function
def success(reg):
    result = reg.measure(position=2)
    result2 = reg.measure(position=1)
    result3 = reg.measure(position=0)
    if result == 1 and result2 == 0 and result3 == 0:
        return True



# loop over many trajectories
iterations = 500

rate_vals = [0.0, 0.00001, 0.000025, 0.00005, 0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5]

output = mp.Queue()

def find_fidelity(decay_rate, output):
    strip = []
    for k in rate_vals:
        print('Started set for k='+str(k))
        trajectories = []
        for i in tqdm(range(iterations)):
            #print('Did this for '+ str(decay_rate))
            reg = QRegister(n_qubits=n_qubits, rabi_freq=0.5, decay_rate=decay_rate, bitflip_rate=0.0, phaseflip_rate=k)
            flip1.apply(reg)
            traj = Trajectory(gates3, reg, success)
            trajectories.append(traj)
        average = AvTrajectory(trajectories)
        fidelity, ferror = average.fidelity_stats()
        strip.append([fidelity, ferror])
        print('Completed set for k=' + str(k))
    output.put((decay_rate, strip))


processes = [mp.Process(target=find_fidelity, args=(m, output)) for m in rate_vals]

for p in processes:
    p.start()
for p in processes:
    p.join()
print('Threads joined')
print('Retrieving output')
rate_fidelities = [output.get() for p in processes]
print('Ouput retrieved')
print('Sorting output')
rate_fidelities.sort()
print('Output sorted')
print('Trimming array')
rate_fidelities = [r[1:] for r in rate_fidelities]
print('Array trimmed')



print('Saving data...')
np.save('Data/Final_Dataset/decay_and_phaseflip.npy', rate_fidelities)
print('Saved..!')
