from futils import*
import matplotlib
matplotlib.rcParams['font.serif'] = "Times New Roman"
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams["mathtext.fontset"] = "dejavuserif"

data2 = np.load('Data/Final_Dataset/errors_included_rate_and_number_decay_only_fidelities.npy')
xdata2 = [0.00001, 0.000025, 0.00005, 0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5]

print(data2)
data3 = []
errors3 = []
errors4 = []
for y in data2:
    data = [x[0] for x in y]
    upper_errors = [x[0] + x[1] for x in y]
    lower_errors = [x[0] - x[1] for x in y]
    data3.append(data)
    errors3.append(upper_errors)
    errors4.append(lower_errors)

n_qubits = [x for x in range(2,6)]
markers = ['o', '^', 's', 'p']

pp.figure(figsize=(3.4, 2.5), dpi=300)
for i in range(len(errors3)):
    pp.fill_between(xdata2, errors3[i], errors4[i], alpha=0.2)
for i in range(len(data3)):
    pp.plot(xdata2, data3[i], linewidth=0.4, marker=markers[i], markersize=3, fillstyle='none',
            label=str(n_qubits[i]), markeredgewidth=0.5)
pp.semilogx()
pp.legend(title='$n_{qubits}$', fontsize='x-small', title_fontsize='x-small')
pp.xlabel('Decay Rate / $\dfrac{\Omega}{2}$')
pp.ylabel('Overall Fidelity')

pp.show()
