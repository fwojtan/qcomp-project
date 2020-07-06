from futils import*
import threading

from multiprocessing import Process, Pool

if __name__ == '__main__':
    timestep = 0.01
    H = 0.5 * np.array([[0, 5+0j], [5+0j, 0]])
    psi = np.array([[0.0 +0.0j],
                   [1.0+0.0j]])
    U = evolution_operator(H, timestep)
    decay_lifetime = 60  # decay rate in timesteps

    masterdata = []
    #step_space = np.linspace(60, 80, 11)

    #for i in step_space:
    evolve_MCWF(psi, U, decay_lifetime, masterdata, timestep, 500, 1000000)

    np.save('Data/Decay_Rate_Tests/many_steps_' + str(i) + '.npy', masterdata)

    #pool = Pool(processes=8)
'''
    threads = Pool()
    what = threads.apply_async(evolve_MCWF, (psi, U, decay_lifetime, masterdata, timestep, 500, 10000))
    threads.close()
    threads.join()
    print(what)
'''



  #  print('100%')

