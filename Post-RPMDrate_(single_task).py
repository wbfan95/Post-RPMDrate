'''
Author: Wenbin FAN
First Release: Oct. 30, 2019
Modified: Dec. 4, 2019
Verision: 1.2

[Intro]
Plot
# (1) the overlap of each xi window,
# (2) the variance in each trajectory,
# (3) potential of mean force,
# (4) recrossing factor,
# (5) the evolution of normalized population of xi,
for a single task.

[Usage]
run `python <this file name>` then you will be ask to input a path containing a result for a single task.
Then you'll get all four figures above if your task ended normally.

Attention, please! The former figures will be deleted when the program started running.

[Contact]
Mail: fanwenbin@shu.edu.cn, langzihuigu@qq.com
WeChat: 178-6567-1650

Thanks for using this python program, with the respect to my supervisor Prof. Yongle!
Great thanks to Yang Hui and Fang Junhua.

[Bug Fixing]
V1.2:
1) PMF: Plot range modified.
'''

import os

import matplotlib.pyplot as plt
import numpy as np

color = ['#00447c', '#ae0d16', '#47872c', '#800964']


# SHU Blue, Weichang Red, Willow Green, SHU Purple
# This color scheme can be easily obtained on the official website `vi.shu.edu.cn`.

def input_path():
    path = input('Please input the result folder: \n')
    return path


def plot_parameters(title):
    print('Plotting {}! '.format(title))
    plt.figure(figsize=(4, 3))
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"  # The font closet to Times New Roman
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['ytick.right'] = True
    plt.rcParams['ytick.left'] = True
    plt.rcParams['xtick.top'] = True
    plt.minorticks_on()  # Turn on minor ticks

    file = '{}.png'.format(title)
    if os.path.exists(file):
        os.remove(file)


def plot_save(name):
    plt.tight_layout()
    plt.savefig(name + '.png', format='png', dpi=600)  # .svg is recommended!
    plt.clf()
    plt.close()


def get_tail(fname):
    fp = open(fname, 'r')
    last_line = fp.readlines()
    # print(last_line)
    fp.close()
    return last_line


def get_xav(fname):
    last_line = get_tail(fname)
    line = last_line[-1]
    line_element = line.split()
    Nl = len(line_element)
    xav = 0.0
    xav2 = 0.0
    if Nl < 3:
        raise ValueError("The last line doesn't contain 'xav' and 'xav2'. ")
    else:
        if line_element[3] == "===============" or line_element[3] == "*":
            print("[ERROR] This task didn't cycle 1 time. ")
            print(fname)
        else:
            xav = np.float(line_element[3])
            xav2 = np.float(line_element[4])
    return xav, xav2


def my_gaussian(x, xav, xav2):
    y = (1.0 / (np.sqrt(2.0 * np.pi * xav2))) * np.exp(-(x - xav) ** 2 / (2.0 * xav2))
    return y


def get_xilist(path):
    fileList = os.listdir(path)
    xiList = []
    for file in fileList:
        if file[0:17] == 'umbrella_sampling':
            xiList.append(np.float(file[18:25].rstrip('.')))
    return xiList


def plot_overlap(path, xiList):
    title = 'Overlap'
    plot_parameters(title)

    resolution = 2000
    extend = 0.02  # 3E-2

    xiMin = np.min(xiList)
    xiMax = np.max(xiList)
    length = len(xiList)

    x_new = np.linspace(xiMin - extend, xiMax + extend, resolution)
    y_sum = np.zeros((resolution))  # Total density line

    maxPop = 0  # maximum of the summation of all population
    for i in range(length):
        fname = path + "/umbrella_sampling_{0:.4f}.dat".format(xiList[i])
        if os.path.isfile(fname):
            # Gaussian smearing
            xav, xav2 = get_xav(fname)
            y_new = my_gaussian(x_new, xav, xav2)

            # Find biggest population
            if max(y_new) > maxPop:
                maxPop = max(y_new)

            y_sum += y_new  # sum all population
            if xav2 > 5.0E-5:
                print("[WARNING] May be too various in xi = {}! ".format(xiList[i]))
                plt.plot(x_new, y_new, lw=1, c=color[1], alpha=0.8)
            else:
                plt.plot(x_new, y_new, lw=0.5, c=color[0], alpha=.3)

    # Plot summation and difference
    plt.plot(x_new, y_sum, lw=1, c=color[0], label='Summation of all populations')  # SHU Blue

    # plt.xlabel('Reaction Coordinate / Å')
    plt.xlabel('Reaction Coordinate')
    plt.ylabel('Population')
    # plt.legend(loc="best") # No legend

    plt.xlim(xiMin - extend, xiMax + extend)
    # plt.ylim(0, maxPop*1.1)
    # plt.ylim(0, max(y_sum)*1.1)

    plt.yticks([])  # No ticks and labels in y axis

    plot_save(title)


def plot_variance(path, xiList):
    title = 'Variance'
    plot_parameters(title)

    xiMin = np.min(xiList)
    xiMax = np.max(xiList)
    length = len(xiList)

    for i in range(length):
        fname = path + "/umbrella_sampling_{0:.4f}.dat".format(xiList[i])
        # print(fname)
        f = open(fname, 'r')
        fLines = f.readlines()
        f.close()

        timeEvolution = []
        xivar = []
        for line in fLines[15:-1]:
            lineSplit = line.split()
            timeEvolution.append(np.float(lineSplit[2]))
            xivar.append(np.float(lineSplit[-1]))

        timeStep = np.float(fLines[6].split()[3])
        timeEvolution = [x * timeStep * 10E-4 for x in timeEvolution]

        # xivarDelta = []
        # for i in range(len(xivar)-1):
        #     xivarDelta.append(np.abs(xivar[i+1] - xivar[i]))
        # x = range(len(xivar))
        # plt.yscale('log')

        # # Shifted
        # for i in range(len(xivar)):
        #     xivar[i] = xivar[i] - xivar[0]

        plt.xlabel('$t$ / ns')
        plt.ylabel('Variance')

        # Scientific notation for y axis
        # # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        from matplotlib import ticker
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        plt.gca().yaxis.set_major_formatter(formatter)

        # color = (np.random.rand(), np.random.rand(), np.random.rand())
        plt.plot(timeEvolution, xivar, lw=0.2, c=color[0], alpha=0.5)

    plot_save(title)


def plot_pmf(path):
    title = 'PMF'
    plot_parameters(title)

    try:
        f = open(path + '/potential_of_mean_force.dat', 'r')
    except FileNotFoundError:
        print('[ERROR] {} file not found! '.format(title))
    else:
        fLines = f.readlines()
        f.close()

        xi = []
        pmf = []
        for i in fLines[12:-1]:
            xi.append(np.float(i.split()[0]))
            pmf.append(np.float(i.split()[1]))

        # Let W(xi=0) = 0!
        xiAbs = np.abs(xi)
        xiZeroIndex = list(xiAbs).index(min(np.abs(xi)))
        pmf = [x - pmf[xiZeroIndex] for x in pmf]

        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$W(\xi)$ / eV')

        # Choose the fit range of this plot
        plt.xlim(xi[0], xi[-1])
        yRange = max(pmf) - min(pmf)
        plt.ylim(min(pmf) - yRange * 0.1,
                 max(pmf) + yRange * 0.1)  # adds 0.1*yRange to the top and bottom

        plt.plot(xi, pmf, c=color[0])

        # # plot a zoomed subfigure
        # xiMaxIndex = pmf.index(max(pmf)) # the position of maximum
        # extend = 80 # Extra points to plot
        # ximax = xi[xiMaxIndex - extend : xiMaxIndex + extend]
        # pmfmax = pmf[xiMaxIndex - extend : xiMaxIndex + extend]
        #
        # # # Find maximum of xi
        # f = itp(ximax, pmfmax, k=4)
        # cr_pts = f.derivative().roots()
        # cr_pts = np.append(cr_pts, (ximax[0], ximax[-1]))  # also check the endpoints of the interval
        # cr_vals = f(cr_pts)
        # #! min_index = np.argmin(cr_vals)
        # max_index = np.argmax(cr_vals)
        # pmfMax, xiMax = cr_vals[max_index], cr_pts[max_index]
        #
        # subfig = plt.axes([.3, .5, .5, .4])
        # subfig.plot(ximax, pmfmax, c=color[0])
        #
        # subfig.axvline(x=xiMax, c=color[0], lw=0.5, linestyle='--')
        # subfig.axhline(y=pmfMax, c=color[0], lw=0.5, linestyle='--')
        #
        # plt.setp(subfig, xlim=[min(ximax), max(ximax)])

        plot_save(title)


def plot_rexFactor(path):
    title = 'Transmission_Coefficient'
    plot_parameters(title)

    # Find the file first!
    fileList = os.listdir(path)
    rexFileName = ''
    for file in fileList:
        if file[:18] == 'recrossing_factor_':
            rexFileName = file

    try:
        f = open(path + '/' + rexFileName, 'r')
    except FileNotFoundError:
        print('[ERROR] {} file not found! '.format(title))
    except PermissionError:
        print('[ERROR] {} file not found! '.format(title))
    else:
        fLines = f.readlines()
        f.close()
        time = []
        kappa = []
        for i in fLines[17:-1]:
            ele = i.split()
            time.append(np.float(ele[0]))
            kappa.append(np.float(ele[-1]))

        plt.xlabel('$t$ / fs')
        plt.ylabel('$\kappa(t)$')

        # plt.xscale('log')

        plt.xlim(time[0], time[-1])
        # endRF = np.mean(kappa[-5:])
        plt.axhline(y=kappa[-1], c=color[0], lw=0.5, linestyle='--')

        plt.plot(time, kappa, c=color[0])
        plot_save(title)


def plot_overlap_density(path, xiList):
    title = 'Overlap_Density'
    plot_parameters(title)

    resolution = 2000
    extend = 0.03  # 3E-2

    xiList = sorted(xiList)  # Sort them

    xiMin = np.min(xiList)
    xiMax = np.max(xiList)
    sizeH = len(xiList)

    x_new = np.linspace(xiMin - extend, xiMax + extend, resolution)

    # Count the total lines of xi and xvar
    sizeV = 5000  # a number enough large
    for i in range(sizeH):
        with open(path + "/umbrella_sampling_{0:.4f}.dat".format(xiList[i]), 'r') as tempFile:
            for i, l in enumerate(tempFile):
                pass
        temp = i + 1 - 15  # 15 info lines # +1 means the number of lines
        if temp < sizeV:
            sizeV = temp

    # Read time unit
    tempFile = open(path + "/umbrella_sampling_{0:.4f}.dat".format(xiList[0]), 'r')
    lines = tempFile.readlines()
    timeSep = np.float(lines[9].split()[4]) / 1000.0  # to ns

    # Read all xi and var
    xiArray = np.zeros((sizeV, sizeH))
    varArray = np.zeros((sizeV, sizeH))

    # Read in all data
    for i in range(sizeH):
        fname = path + "/umbrella_sampling_{0:.4f}.dat".format(xiList[i])
        f = open(fname, 'r').readlines()

        for j in range(sizeV):
            line = f[15 + j].split()
            # print(i, line[-2], line[-1])
            xiArray[j, i] = line[-2]
            varArray[j, i] = line[-1]

    z = np.zeros((sizeV * resolution))
    y = np.linspace(xiMin - extend, xiMax + extend, resolution)

    # Gaussian summation
    for j in range(sizeV):  # xi
        y_sum = np.zeros((resolution))
        for i in range(sizeH):  # var
            y_new = my_gaussian(x_new, xiArray[j, i], varArray[j, i])
            y_sum += y_new
            z[j * resolution:(j + 1) * (resolution)] = y_sum

    # plt.show()
    z = z.reshape(sizeV, resolution).transpose()
    z /= np.max(z)
    x = np.linspace(0, timeSep * sizeV, sizeV)
    plt.pcolormesh(x, y, z, cmap='Greens', vmax=1.0)  # pcolormesh # contourf # Greys_r

    plt.xlabel('time / ns')
    plt.ylabel('Reaction Coordinate')

    plt.colorbar()
    # plt.title('The Evolution of Normalized Population')
    plot_save(title)


def main(path=None):
    if path == None:
        path = input_path()
    xiList = get_xilist(path)

    plot_overlap(path, xiList)
    plot_variance(path, xiList)
    plot_pmf(path)
    plot_rexFactor(path)
    plot_overlap_density(path, xiList)


main()
# main()
