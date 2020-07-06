from futils import*
import threading
from multiprocessing import Process, Pool

if __name__ == '__main__':
    timestep = 0.01
    H = 0.5 * np.array([[0, 1+0j], [1+0j, 0]])
    psi = np.array([[1.0 +0.0j],
                   [0.0+0.0j]])
    U = evolution_operator(H, timestep)
    decay_lifetime = 0.1  # decay rate in rabi period

    temp_path = 'Data/temp_path.npy'
    #step_space = np.linspace(60, 80, 11)

    masterdata = []

    #cos_MCWF(1, decay_lifetime, masterdata, timestep, n_timesteps = 3000, iterations= 10000)
    #np.save(temp_path, masterdata)
    #average_WF(temp_path, 'Data/Milestone_Data/cosine_method_10k.npy')

    evolve_MCWF(psi, U, decay_lifetime, 1, masterdata, timestep, 3000, 10000)
    np.save(temp_path, masterdata)
    average_WF(temp_path, 'Data/Milestone_Data/evolution_method_10k.npy')

    #pool = Pool(processes=8)
'''
    threads = Pool()
    what = threads.apply_async(evolve_MCWF, (psi, U, decay_lifetime, masterdata, timestep, 500, 10000))
    threads.close()
    threads.join()
    print(what)
'''



  #  print('100%')

