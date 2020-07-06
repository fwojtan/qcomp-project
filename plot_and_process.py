from futils import*
import multiprocessing

max_reg_size = 6

fidelities = []

for k in tqdm(range(2, max_reg_size)):
    n_qubits = k

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

    if k == 3 or 4 or 5:
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

    if k == 4 or 5:
        Ry3 = Gate('Ry', position=2, n_qubits=n_qubits, gate_param=(np.pi/2), out_time=0.1, lead_time=0.1)
        XX3 = Adjacent_Binary_Gate('Ising', positions=[2,3], n_qubits=n_qubits, gate_param=(np.pi/4), out_time=0.1, lead_time=0.1)
        negRx3 = Gate('Rx', position=3, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
        negRy3 = Gate('Ry', position=2, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
        negRz3 = Gate('Rz', position=2, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)

        back_Ry4 = Gate('Ry', position=3, n_qubits=n_qubits, gate_param=(np.pi / 2), out_time=0.1, lead_time=0.1)
        back_negRx4 = Gate('Rx', position=2, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
        back_negRy4 = Gate('Ry', position=3, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)
        back_negRz4 = Gate('Rz', position=3, n_qubits=n_qubits, gate_param=(-np.pi / 2), out_time=0.1, lead_time=0.1)

        swap34 = Adjacent_Binary_Gate('SWAP', positions=[2, 3], n_qubits=n_qubits, duration=0.01, lead_time=0.01, out_time=0.01)

        CNOT34 = [Ry3, XX3, negRy3, negRx3, negRz3]
        CNOT41 = [swap12, swap23, back_Ry4, XX3, back_negRy4, back_negRx4, back_negRz4, swap23, swap12]
        CNOT42 = [swap23, back_Ry4, XX3, back_negRy4, back_negRx4, back_negRz4, swap23]
        CNOT43 = [back_Ry4, XX3, back_negRy4, back_negRx4, back_negRz4]

    if k == 5:
        Ry4 = Gate('Ry', position=3, n_qubits=n_qubits, gate_param=(np.pi/2), out_time=0.1, lead_time=0.1)
        XX4 = Adjacent_Binary_Gate('Ising', positions=[3,4], n_qubits=n_qubits, gate_param=(np.pi/4), out_time=0.1, lead_time=0.1)
        negRx4 = Gate('Rx', position=4, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
        negRy4 = Gate('Ry', position=3, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
        negRz4 = Gate('Rz', position=3, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)

        Ry5 = Gate('Ry', position=4, n_qubits=n_qubits, gate_param=(np.pi/2), out_time=0.1, lead_time=0.1)
        #XX5 = Adjacent_Binary_Gate('Ising', positions=[4,5], n_qubits=n_qubits, gate_param=(np.pi/4), out_time=0.1, lead_time=0.1)
        negRx5 = Gate('Rx', position=3, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
        negRy5 = Gate('Ry', position=4, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)
        negRz5 = Gate('Rz', position=4, n_qubits=n_qubits, gate_param=(-np.pi/2), out_time=0.1, lead_time=0.1)

        CNOT45 = [Ry4, XX4, negRy4, negRx4, negRz4]
        CNOT51 = [swap12, swap23, swap34, Ry5, XX4, negRy5, negRx5, negRz5, swap34, swap23, swap12]
        CNOT52 = [swap23, swap34, Ry5, XX4, negRy5, negRx5, negRz5, swap34, swap23]
        CNOT53 = [swap34, Ry5, XX4, negRy5, negRx5, negRz5, swap34]
        CNOT54 = [Ry5, XX4, negRy5, negRx5, negRz5]


    flip1 = Gate('Ry', position=0, n_qubits=n_qubits, gate_param=np.pi, out_time=0.1, lead_time=0.1)
    flip2 = Gate('Ry', position=1, n_qubits=n_qubits, gate_param=np.pi, out_time=0.1, lead_time=0.1)


    if k == 2:
        gates3 = CNOT12 + CNOT21 + CNOT21 + CNOT21 + CNOT21 + CNOT21 + CNOT21 + CNOT21
    elif k == 3:
        gates3 = CNOT12 + CNOT23 + CNOT31 + CNOT32 + CNOT31 + CNOT31 + CNOT31 + CNOT31
    elif k == 4:
        gates3 = CNOT12 + CNOT23 + CNOT34 + CNOT41 + CNOT42 + CNOT43 + CNOT41 + CNOT41
    elif k == 5:
        gates3 = CNOT12 + CNOT23 + CNOT34 + CNOT45 + CNOT51 + CNOT52 + CNOT53 + CNOT54








    #trueCNOT = Adjacent_Binary_Gate('CNOT', positions=[0,1], n_qubits=n_qubits, duration=0.01)

    #print(CNOT12[1].matrix)

    #reg = QRegister(n_qubits=n_qubits, rabi_freq=0.5, decay_rate=0.0, bitflip_rate=0.0, phaseflip_rate=0.0)
    #flip1.apply(reg)
    #print(trueCNOT.flash(reg))
    #traj = Trajectory(gates3, reg)

    #traj.plot(enable_axspan=True, enable_boundaries=False, rotate_labels=True, fig_width=12)


    # write a success function
    def success(reg):
        result = reg.measure(position=k-1)
        if result == 1:
            return True

    # loop over many trajectories
    iterations = 500

    rate_vals = [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5]

    def find_fidelity(decay_rate, target)
        trajectories = []
        for i in range(iterations):
            reg = QRegister(n_qubits=n_qubits, rabi_freq=0.5, decay_rate=decay_rate, bitflip_rate=0.0, phaseflip_rate=0.0)
            flip1.apply(reg)
            traj = Trajectory(gates3, reg, success)
            trajectories.append(traj)
        average = AvTrajectory(trajectories)
        fidelity = average.fidelity_stats()
        target.append(fidelity)


    rate_fidelities = []
    for m in tqdm(rate_vals):
        find_fidelity(m, rate_fidelities)

    fidelities.append(rate_fidelities)



np.save('Data/Final_Dataset/rate_and_number_decay_only_fidelities.npy', fidelities)

