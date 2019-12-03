'''
Author: Wenbin, FAN
Update: Dec. 2, 2019

[Intro]
Plot the evolution of xi and its histogram.

[Usage]
The directory should contain a series of file including the coordinates.
'''
import os

import matplotlib.pyplot as plt
import numpy as np


def plotPara():
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"  # The font closet to Times New Roman
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['ytick.right'] = True
    plt.rcParams['ytick.left'] = True
    plt.rcParams['xtick.top'] = True
    plt.minorticks_on()  # Turn on minor ticks


def main(path=None):
    if path == None:
        path = input('Please input the trajectory folder: \n')

    plotPara()
    fig = plt.figure(figsize=(4, 3))
    grid = plt.GridSpec(1, 2, wspace=0, width_ratios=[3, 1])
    a0 = fig.add_subplot(grid[0])
    a1 = fig.add_subplot(grid[1])  # , sharey=a0)

    fileList = os.listdir(path)

    # Descending sort
    name = []
    value = []
    for file in fileList:
        if str(file[12:18]) != 'centro':
            name.append(str(file[12:16]))
            value.append(np.float(file[12:16]))
    z = zip(value, name)
    z = sorted(z, reverse=True)
    value, name = zip(*z)

    # compute mean and variance for current kforce
    mean = []
    var = []

    all = []
    for i in range(len(name)):
        f = open(path + 'equilibrate_' + name[i] + '.xyz', 'r')
        xi_current = value[i]

        # Read data in each Xi file
        xi = []
        # lineCount = 1
        atomCount = 0
        Natoms = 0
        for line in f:
            if atomCount == 1:
                float = np.float(line.split()[-1])  # - xi_current
                xi.append(float)
                all.append(float)
            if atomCount == 0:
                Natoms = np.int(line)
            atomCount += 1
            if atomCount == Natoms + 2:
                atomCount = 0
            # lineCount += 1
        # xi_count.append(len(xi))
        f.close()

        mean_current = np.mean(xi)
        var_current = np.var(xi)
        mean.append(mean_current - xi_current)
        var.append(np.var(xi))

        # Get histogram (right panel)
        a, b = np.histogram(xi, bins=50)
        a1.plot(a, b[:-1], lw=0.4, alpha=0.2, c='black')

        print("Current xi = {:.2f}, {} points loaded. Mean xi is {:.3f}, and variance is {:.5f}".
              format(xi_current, len(xi), mean_current, var_current))

        # Plot the last xi for the spawn trajectories.
        if np.abs(xi[-1] - xi_current) > 2E-3:  # threshold for deviation # 偏差阈值，默认 2E-3 有点严格，可改大
            a0.axhline(y=xi[-1], linestyle='-', lw=0.2, c='red')
        else:
            a0.axhline(y=xi[-1], linestyle='--', lw=0.1, c='black', alpha=0.5)

    print('\n')
    print("Mean of the mean xi, mean of variance: ")
    print(np.mean(mean), '\t', np.mean(var))

    # Get total histogram
    c, d = np.histogram(all, bins=600)
    a1.plot(c / 2, d[:-1], lw=0.5, c='black')

    xiLimitBelow = -0.1
    xiLimitAbove = 1.1

    a1.set_xticks([])
    a1.set_yticklabels([])
    a1.set_ylim([xiLimitBelow, xiLimitAbove])
    a1.set_xlim([0, len(all) / len(name) / 5.0])

    # Evolution
    bins = 480
    N = 10  # for time interval: N * 0.1 fs
    Y = np.linspace(xiLimitBelow, xiLimitAbove, bins)

    allDataLen = len(all)
    interval = int(allDataLen / N)  # the number of xi in each interval

    dens = np.zeros((interval * bins))

    for j in range(interval):
        current_data = all[N * j:N * (j + 1)]
        a, b = np.histogram(current_data, bins=bins, range=(xiLimitBelow, xiLimitAbove))
        dens[bins * j:bins * (j + 1)] = a

    Z = dens.reshape(interval, bins).transpose()
    Z /= np.max(Z)
    X = np.linspace(0, allDataLen * 0.1 * 0.001, num=interval)

    a0.pcolormesh(X, Y, Z, cmap='Greys', vmax=0.5)
    # a0.colorbar()

    a0.set(xlabel='time / ps', ylabel='Reaction Coordinate')
    # a0.set_yticks(Y)

    # fig.show()
    fig.tight_layout()
    fig.savefig('UmbrellaConfEvolution.png', format='png', dpi=600)  # .png is recommended!
    fig.clf()
    # fig.close()


main()
