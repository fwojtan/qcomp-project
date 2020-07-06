# a script to look at the 'detuning' effect on the project milestone

from futils import*

timestep = 0.1

detuning = 0.001
omega = 1

H = 0.5 * np.array([[0, omega], [omega, 0]])
H2 = 0.5 * np.array([[detuning, omega], [omega, detuning]])
psi = np.array([[1.0 +0.0j],
               [0.0+0.0j]])
U = evolution_operator(H, timestep)
U2 = evolution_operator(H2, timestep)
decay_lifetime = 0.1
masterdata = []
thing1 = evolve_MCWF(psi, U, decay_lifetime, 1, masterdata, timestep, 1000, 1, decay_off= True)
trajectory = average_WF2(masterdata)
m2data = []
thing2 = evolve_MCWF(psi, U2, decay_lifetime, 1, m2data, timestep, 1000, 1, decay_off= True)
t2 = average_WF2(m2data)

a = [1, 2]
b = [2, 3]


# before I realised you could simply join lists together using +
def stitch(lists):
    result = []
    for i in lists:
        result.append(i)
    return result

'''
tdata = np.load('Data/timedata3000.npy')
pp.plot(tdata[:1000], thing1)
pp.plot(tdata[:1000], thing2)
pp.show()
'''
