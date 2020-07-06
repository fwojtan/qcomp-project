from futils import*
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

matplotlib.rcParams['font.serif'] = "Times New Roman"
matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams["mathtext.fontset"] = "dejavuserif"

data = np.load('Data/Final_Dataset/decay_and_bitflip.npy')


data2 = np.array([data[x][0] for x in range(len(data))])
xydata = [0.0, 0.00001, 0.000025, 0.00005, 0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5]
data3 = []
data4 = []
chunk_data = []

xdata = []
ydata = []

for i in range(len(xydata)):
    xslice = [xydata[i] for x in xydata]
    yslice = xydata
    xdata += xslice
    ydata += yslice

for i in range(len(data2)):
    slice = []
    for k in range(len(data2[i])):
        slice.append(data2[i][k][:1][0])
        chonk = [xydata[i], xydata[k], data2[i][k][:1][0]]
        chunk_data += chonk
    data3.append(slice)
    data4 += slice

data3 = np.array(data3)
print(data3)

np.save('Data/Final_Dataset/bitdata1.npy', data4)

n_qubits = [x for x in range(2,6)]
markers = ['o', '^', 's', 'p']

fig = pp.figure(figsize=(3.4, 2.5), dpi=300)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.log10(xdata), np.log10(ydata), data4)
ax.view_init(45, 45)
#pp.semilogx()
pp.xlabel('Decay Rate / Arbitrary Units', fontdict={'fontsize':8})
pp.ylabel('Bit-flip Rate / Arbitrary Units', fontdict={'fontsize':8})
ax.set_zlabel('Fidelity', fontdict={'fontsize':8})
#pp.xticks([-1, -2. -3, -4, -5], ['$10^{-1}$', '$10^{-2}$', '$10^{-3}$', '$10^{-4}$', '$10^{-5}$'])
ax.set_xticklabels(['0', '10$^{-5}$', '10$^{-4}$', '10$^{-3}$', '10$^{-2}$', '10$^{-1}$'], fontdict={'fontsize':8})
pp.tick_params(axis='x', pad=-3, labelrotation=-5)
ax.set_yticklabels(['0', '10$^{-5}$', '10$^{-4}$', '10$^{-3}$', '10$^{-2}$', '10$^{-1}$'], fontdict={'fontsize':8})
pp.tick_params(axis='y', pad=-3, labelrotation=5)
ax.set_zticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'])#, fontdict={'fontsize':8})
pp.tick_params(axis='z', pad=-2.5)
pp.show()
