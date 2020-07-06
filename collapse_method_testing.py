# a script testing a procedure to collapse the wavefunction of a system

import random
from matplotlib import pyplot as pp
import numpy as np
from futils import*

timestep = 0.01
decay_lifetime = 0.1

H = 0.5 * np.array([[0, 1+0j], [1+0j, 0]])
psi = np.array([[1.0 +0.0j],
                [0.0+0.0j]])
U = evolution_operator(H, timestep)


tdata = np.load('Data/timedata3000.npy')

pp.plot(tdata, cos_MCWF(1, decay_lifetime, n_timesteps= 3000))
#pp.plot(tdata, evolve_MCWF(psi, U, decay_lifetime, 1, timestep = timestep, n_timesteps = 3000))
pp.title('Two state system')
pp.xlabel('Time')
pp.ylabel('Probability amplitude')
#pp.legend()
pp.show()

