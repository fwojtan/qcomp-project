from futils import*

n_qubits = 2
qreg = QRegister(n_qubits, decay_off=True)

Ry = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=(np.pi / 2), out_time=0.1, lead_time=0.1)
XX = Adjacent_Binary_Gate('Ising', positions=[0, 1], n_qubits=n_qubits, gate_param=(np.pi / 4), out_time=0.1, lead_time=0.1)
negRx = Gate('Rx', position=1, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)
negRy = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)
negRz = Gate('Rz', position=0, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)

flip1 = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=np.pi, out_time=0.1, lead_time=0.1)

CNOT12 = [Ry, XX, negRy, negRx, negRz]

# NEEDS FLIP HERE
flip1.apply(qreg)
traj = Trajectory(CNOT12, qreg)

# NEEDS AXIS LABELS AND FONT CHANGES
#traj.plot()



IdGate = Gate('Identity', 0, 1, 1, lead_time=0, out_time=75)
#trajectories = []
#for i in tqdm(range(500)):
single = QRegister(1, decay_rate=0.025, bitflip_rate=0, phaseflip_rate=0)
traj2 = Trajectory([IdGate], single)
#trajectories.append(traj2)

data = np.load('Data/Final_Dataset/mcwf_data.npy')
err_data = np.load('Data/Final_Dataset/mcwf_errors.npy')

#average = AvTrajectory(trajectories)
iterations = len(traj2.p_trajectory)
timedata = traj2.timestep * np.array(range(iterations))
pp.figure(figsize=(3.4, 2.5), dpi=300)
#average.p_trajectory = np.array(average.p_trajectory)
upper = err_data[0]
lower = err_data[1]
#for i in range(len(average.errors)):
#    up = average.p_trajectory[i][1] + average.errors[i][1]
#    down = average.p_trajectory[i][1] - average.errors[i][1]
#    upper.append(up[0])
#    lower.append(down[0])
#print(average.p_trajectory)
#print(average.errors)
#print(upper)
#print(lower)
pp.fill_between(timedata, upper, lower, alpha=0.2)
#savedata = average.p_trajectory[:, 1]
#err_data = [upper, lower]
#np.save('Data/Final_Dataset/mcwf_data.npy', savedata)
#np.save('Data/Final_Dataset/mcwf_errors.npy', err_data)
pp.plot(timedata, data, linewidth=0.5)
pp.ylim(-0.05, 1.05)
pp.xlabel('Time / 100 time-steps')
pp.ylabel('Probability of Excited State')
pp.show()

#average.mcwf_plot()
