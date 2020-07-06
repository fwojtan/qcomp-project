# a simple script to reformat data

from futils import*
tdata = np.load('Data/timedata1500.npy')
values = [15, 20, 25, 30]

file = np.load('Data/Post_Advice_06_cos.npy')
pdata = []
iterations = len(file)
time_length = len(file[0])
for i in tqdm(range(time_length)):
    running_total = 0
    for k in range(iterations):
        running_total += file[k][i]
    value = running_total / iterations
    pdata.append(value)
np.save('Data/Post_Advice_06_cos_.npy', pdata)
#pp.show()

#print(file)