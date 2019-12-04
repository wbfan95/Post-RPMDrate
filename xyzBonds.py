'''
Program: Plot bond length from .xyz file
Author: Wenbin, FAN

V 1.1 # 20191204 13:20:56 Wenbin, FAN @ SHU
1) Plot OHCH4 system only.
'''

import matplotlib.pyplot as plt
import numpy as np

color = ['#00447c', '#ae0d16', '#47872c', '#800964']
# SHU Blue, Weichang Red, Willow Green, SHU Purple
# This color scheme can be easily obtained on the official website `vi.shu.edu.cn`.

def plotPara():
    plt.figure(figsize=(4, 3))
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"  # The font closet to Times New Roman
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['ytick.right'] = True
    plt.rcParams['ytick.left'] = True
    plt.rcParams['xtick.top'] = True
    plt.minorticks_on()  # Turn on minor ticks


def plotSave(title):
    plt.tight_layout()
    foo_fig = plt.gcf()  # 'get current figure'
    foo_fig.savefig(title + '.png', format='png', dpi=600)
    plt.show()
    # plt.clf()


def readCoord(file=None, Nbeads=None):
    if file is None:
        file = input('Please input the trajectory path: \n')
    print('[INFO] Reading File... ')
    f = open(file, 'r')
    flines = f.readlines()
    # f.close()

    if Nbeads is None:
        Nbeads = int(input('Please input the number of beads (default `1`): \n') or '1')
    Nbeads = np.int(Nbeads)
    Nab = np.int(flines[0])
    Natoms = np.int(Nab / Nbeads)
    Nframe = np.int(len(flines) / (Nab + 2))  # `2` is amount line and title line.

    global atomList
    atomList = []
    for i in range(Natoms):
        atomList.append(flines[i * Nbeads + 2].split()[0])

    assert Nframe * (Nbeads * Natoms + 2) == len(flines)

    # print('Information about current trajectory: ')
    print(r'[INFO] Atoms: {}, Beads: {}, frames: {}'.format(Natoms, Nbeads, Nframe))
    global coord, titles
    coord = np.zeros([3, Natoms, Nbeads, Nframe])
    titles = np.zeros(Nframe)

    for frame in range(Nframe):
        titles[frame] = flines[frame * (Nab + 2) + 1]
        for atom in range(Natoms):
            for bead in range(Nbeads):
                basicFrame = frame * (Nab + 2)
                lines = flines[basicFrame + atom * Nbeads + bead + 2]
                line = lines.split()
                for i in range(3):
                    coord[i, atom, bead, frame] = line[i + 1]

    return coord, atomList


def getCentroid(coord):
    global Natoms, Nbeads, Nframe
    xyz, Natoms, Nbeads, Nframe = np.shape(coord)
    global centr
    centr = np.zeros([xyz, Natoms, Nframe])

    if Nbeads > 1:
        print('[INFO] Computing centroid... ')
        for frame in range(Nframe):
            # for bead in range(Nbeads):
            for atom in range(Natoms):
                for i in range(xyz):
                    centr[i, atom, frame] = np.average(coord[i, atom, :, frame])
    elif Nbeads == 1:
        centr[:, :, :] = coord[:, :, 0, :]

    return centr


def getBondLen(a, b):
    # assert np.shape(a) == np.shape(b)
    print('[INFO] Computing bond length {}({}) -- {}({})...'.format(a, atomList[a], b, atomList[b]))

    dist = np.zeros(Nframe)
    for frame in range(Nframe):
        dist[frame] = getDist(centr[:, a, frame], centr[:, b, frame])

    return dist


def getDist(a, b):
    assert len(a) == len(b) == 3

    squaredSum = 0E0
    for i in range(3):
        squaredSum += (a[i] - b[i]) ** 2E0

    return np.sqrt(squaredSum)


def getBeadBondLen(a, b):
    print('[INFO] Computing bond length {}({}) -- {}({})...'.format(a, atomList[a], b, atomList[b]))

    beadDist = np.zeros([Nframe, Nbeads])
    # minDist = np.zeros(Nframe)
    # maxDist = np.zeros(Nframe)
    for frame in range(Nframe):
        for bead in range(Nbeads):
            beadDist[frame, bead] = getDist(coord[:, a, bead, frame], coord[:, b, bead, frame])
        # minDist[frame] = np.min(beadDist[frame, :])
        # maxDist[frame] = np.max(beadDist[frame, :])

    return beadDist  # , maxDist, minDist


def plotBond(atomA, atomB):
    bondc = getBondLen(atomA, atomB)
    plotPara()
    plt.xlim(0, Nframe - 1)

    label = '{} and {}'.format(atomList[atomA], atomList[atomB])

    if Nbeads > 1:
        bond = getBeadBondLen(atomA, atomB)
        maxDist = np.max(bond, axis=1)
        minDist = np.min(bond, axis=1)
        meanDist = np.mean(bond, axis=1)

        plt.fill_between(list(range(Nframe)), maxDist, minDist, color=color[0], alpha=0.2, lw=0)
        plt.plot(bondc, c=color[0], label=r'$r_{\mathrm{c}}$')
        plt.plot(meanDist, ':', c=color[1], label=r'$\bar{r}$')

        plt.legend()
        plt.xlabel('Frame')
        plt.ylabel('Bond Length of {} / Å'.format(label))
        plotSave('bondLength')
    elif Nbeads == 1:
        plt.plot(bondc, c=color[0])

        plt.xlabel('Frame')
        plt.ylabel('Bond Length of {} / Å'.format(label))
        plotSave('bondLength')


def plotBeads():
    end = ''
    while True:
        if end != 'end' or '0' or 'exit':
            end = input('Please input the TWO order of atom (from `0`, e.g. `1 2`): \n')
            atom = end.split()
            assert len(atom) == 2
            plotBond(int(atom[0]), int(atom[1]))
        else:
            print('[INFO] Exit...')
            break


# Unfinished # 20191204 12:18:51 Wenbin, FAN @ SHU
def plotMulti(list):
    global Nbeads
    Nbeads = 1

    for i in list:  # each group was
        assert len(list[i]) == 2

    length = len(list)


def plotOHCH4():
    bondCH = (getBondLen(5, 1) + getBondLen(5, 2) + getBondLen(5, 3) + getBondLen(5, 4)) / 4.0
    bondOH = getBondLen(0, 6)
    bondCO = getBondLen(5, 6)

    plotPara()

    ticks = np.array(titles, dtype=float)

    plt.plot(ticks, bondOH, label='O - H', c=color[1])
    plt.plot(ticks, bondCO, label='C - O', c=color[0])
    plt.plot(ticks, bondCH, label='C - other H', c=color[2])

    plt.legend()
    plt.xlabel(r'$\xi$')
    plt.ylabel('Bond Length / Å')

    plt.xlim(min(ticks), max(ticks))

    plotSave('bondLength')

    plt.show()

def main():
    path = input('Please input the path of `.xyz`: \n')

    coord, atomList = readCoord(path)
    centr = getCentroid(coord)

    plotOHCH4()



main()

atomicMass = {
    "Mu": 0.113,
    "H": 1.00782503207,
    "D": 2.0141017778,
    "T": 3.0160492777,
    "He3": 3.0160293191,
    "He": 4.00260325415,
    "Li6": 6.015122795,
    "Li": 7.01600455,
    "Be": 9.0121822,
    "B5": 10.0129370,
    "B": 11.0093054,
    "C": 12.0,
    "C13": 13.0033548378,
    "N": 14.0030740048,
    "N15": 15.0001088982,
    "O": 15.99491461956,
    "O18": 17.9991610,
    "F": 18.99840322,
    "Ne": 19.9924401754,
    "Ne22": 21.991385114,
    "Na": 22.9897692809,
    "Mg": 23.985041700,
    "Mg25": 24.98583692,
    "Mg26": 25.982592929,
    "Al": 26.98153863,
    "Si": 27.9769265325,
    "Si29": 28.976494700,
    "Si30": 29.97377017,
    "P": 30.97376163,
    "S": 31.97207100,
    "S34": 33.96786690,
    "Cl": 34.96885268,
    "Cl37": 36.96590259,
    "Ar": 39.9623831225,
    "K": 38.96370668,
    "K41": 40.96182576,
    "Ca": 39.96259098,
    "Sc": 44.9559119,
    "Ti": 47.9479463,
    "V": 50.9439595,
    "Cr50": 49.9460442,
    "Cr": 51.9405075,
    "Cr53": 52.9406494,
    "Cr54": 53.9388804,
    "Mn": 54.9380451,
    "Fe54": 53.9396105,
    "Fe": 55.9349375,
    "Fe57": 56.9353940,
    "Co": 58.9331950,
    "Ni": 57.9353429,
    "Cu": 62.9295975,
    "Cu65": 64.9277895,
    "Zn64": 63.9291422,
    "Zn66": 65.9260334,
    "Zn67": 66.9271273,
    "Zn68": 67.9248442,
    "Ga": 68.9255736,
    "Ga71": 70.9247013,
    "Ge70": 69.9242474,
    "Ge72": 71.9220758,
    "Ge73": 72.9234589,
    "Ge74": 73.9211778,
    "Ge76": 75.9214026,
    "As": 74.9215965,
    "Se74": 73.9224764,
    "Se76": 75.9192136,
    "Se77": 76.9199140,
    "Se78": 77.9173091,
    "Se80": 79.9165213,
    "Se82": 81.9166994,
    "Br": 78.9183371,
    "Br81": 80.9162906,
    "Kr": 83.911507,
    "Rb": 84.911789738,
    "Rb87": 86.909180527,
    "Sr": 87.9056121,
    "Y": 88.9058483,
    "Zr": 89.9047044,
    "Zr91": 90.9056458,
    "Zr92": 91.9050408,
    "Zr94": 93.9063152,
    "Zr96": 95.9082734,
}
