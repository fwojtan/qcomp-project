# a script that generates solutions to the optical bloch equations for comparison with the 
# numerical approach for the project milestone

from futils import*

omega = 1.0
gamma = 0.1
gamma = 0.1 * 4 / 3

tdata = np.load('Data/timedata3000.npy')
pdata = []
p2data = []
p3data = []

def probability(omega, gamma, t):
    A = 1 / (2 + (gamma ** 2) / (omega ** 2))
    B = (3 * A * gamma) / (4 * np.sqrt(omega ** 2 - (gamma ** 2) / 16))

    term1 = A * np.cos(t * np.sqrt(omega**2 - (gamma**2)/16))
    term2 = B * np.sin(t * np.sqrt(omega**2 - (gamma**2)/16))
    term3 = A
    return (np.exp((-0.75)*gamma*t) * (term1 + term2) + term3) * (-1) + 1

def probability2(omega, gamma, t):
    A = 1 / (2 + (gamma ** 2) / (omega ** 2))
    B = (3 * A * gamma) / (4 * np.sqrt(omega ** 2 - (gamma ** 2) / 16))

    term1 = A * np.cos(t*np.sqrt(omega**2 - (gamma**2)/16))
    term2 = B * np.sin(t*np.sqrt(omega**2 - (gamma**2)/16))
    term3 = A
    return (1 - (np.exp((-0.75)*gamma*t) * (term1 + term2) + term3))/2

def alt_prob(omega, gamma, t):
    A = 1 / (2 + (gamma ** 2) / (omega ** 2))
    B = (3 * A * gamma) / (4 * np.sqrt(omega ** 2 - (gamma ** 2) / 16))

    term1 = A * np.cos(t * np.sqrt(omega**2 - (gamma**2)/16))
    term2 = B * np.sin(t * np.sqrt(omega**2 - (gamma**2)/16))
    term3 = 1 / (1 + 2*(omega**2 / gamma**2))
    pdiff =  (np.exp((-0.75)*gamma*t) * (term1 + term2) + term3)
    return (0.5 - pdiff)

for t in tdata:
    pdata.append(alt_prob(omega, gamma, t))
    p2data.append(probability(omega, gamma, t))
    p3data.append(probability2(omega, gamma, t))
pp.plot(tdata, pdata, label = 'diff method', alpha = 0.4)
pp.plot(tdata, p2data, label = 'first method', alpha = 0.4)
pp.plot(tdata, p3data, label = 'trial method', alpha = 0.4)
pp.legend()
pp.show()



np.save('Data/OBE_Solution.npy', pdata)

