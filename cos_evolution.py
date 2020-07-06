# a script where the behaviour of the system is hardcoded as the expected solution

from futils import*

timestep = 0.01
decay_lifetime = 0.1
n_timesteps = 3000

masterdata= []
for k in tqdm(range(10000)):

    pdata = []
    prob = 1
    c
    for i in range(n_timesteps):

        threshold = decay_lifetime / n_timesteps

        epsilon = random.random()

        if epsilon < threshold:
            prob = 0
            counter = 0
            #jump condition, collapse to ground state
        else:
            prob = -0.5*np.cos(i*timestep*1) + 0.5
            #this is the no jump condition, non hermitian Hamiltionian followed be renormalisation
        pdata.append(prob)
        counter += 1
    masterdata.append(pdata)

np.save('Data/Post_Advice_06_cos.npy', masterdata)
