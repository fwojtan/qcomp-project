# a script exploring the possibility of animating plots live while the data is being generated
# there was no practical reason to do this other than the fact I thought it might look cool

from futils import*

def generate_animation_data(data):
    masterdata = np.zeros((len(data), len(data[0])))
    for i in tqdm(range(len(data))):
        a = np.array(data[i])
        if i == 0:
            masterdata[i] = a
        else:
            masterdata[i] = (a + i * masterdata[i-1]) / (i+1)
    return masterdata

def plot_animation(pdata, xdata, obe_soln, frame_interval, length):
    pp.ion()
    grid = pp.GridSpec(2, 1, hspace = 0, height_ratios=[3, 1])
    fig = pp.figure()
    ax = fig.add_subplot(grid[0, 0])
    line3, = ax.plot(xdata[:length], obe_soln[:length], color=(1, 0, 0.25), linestyle='--')
    line, = ax.plot(xdata[:length], pdata[0][:length])
    pp.xlabel('Time')
    pp.ylabel('Excited state population')
    pp.title('Wavefunction Trajectory as N increases')
    ax2 = fig.add_subplot(grid[1, 0])
    resids, nresids = residuals(pdata[0][:length], obe_soln[:length])
    ax2.set_ylim(-.1, .1)
    line4, = ax2.plot((xdata[0], xdata[length]), (0.01, 0.01), color=(0.75, 0.25, 0.25))
    line5, = ax2.plot((xdata[0], xdata[length]), (-0.01, -0.01), color=(0.75, 0.25, 0.25))
    line2, = ax2.plot(xdata[:length], resids)

    pp.ylabel('Residuals')

    pp.show()

    for i in range(len(pdata)):
        if i % frame_interval == 0:
            resids, nresids = residuals(pdata[i][:length], obe_soln[:length])
            line2.set_ydata(resids)
            line.set_ydata(pdata[i][:length])
            #ax.set_ylim((pdata.min() / (i + 1)), (pdata.max() / (i + 1)))
            fig.canvas.draw()
            fig.canvas.flush_events()
        else:
            pass

if __name__ == '__main__':

    plot_animation(np.load('Data/animation_test_data.npy'), np.load('Data/timedata3000.npy'), np.load('Data/OBE_Solution.npy'), 2, 1000)

