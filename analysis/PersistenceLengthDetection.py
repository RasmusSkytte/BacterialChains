#!/usr/bin/env python3

import os
from math import sqrt
import numpy as np
from scipy.optimize import curve_fit
from multiprocessing import Pool

def R2_fun( l, p ) :
    return 2 * p * l * (1 - p / l * (1 - np.exp(- l / p)))

def analyzer( path ) :

    # Loop through cell data to find last data point
    with open(path + 'CellData.txt', 'r') as f_cell :

        # Read the first line
        line = f_cell.readline()
        datum = [float(l.strip()) for l in line.split('\t')]

        # Current time
        T = datum[0]

        while not len(line) == 0 :
            datum = [float(l.strip()) for l in line.split('\t')]

            # Current time
            T = datum[0]

            # Read next line
            line = f_cell.readline()

    # Load the bacteria into memory
    with open(path + 'CellData.txt', 'r') as f_cell :

        # Read the lines until datum[0] = T
        line = f_cell.readline()
        datum = [float(l.strip()) for l in line.split('\t')]

        # Read until datum[0] == T
        while (not len(line) == 0) and (not datum[0] == T) :

            # Read next line
            line = f_cell.readline()

            # Store data point
            datum = [float(l.strip()) for l in line.split('\t')]

        # Store the last read value
        cells = [datum[1:]]

        # Read the next line
        line = f_cell.readline()

        # Loop over remaining lines in the file
        while not len(line) == 0 :

            datum = [float(l.strip()) for l in line.split('\t')]
            cells.append(datum[1:])

            # Read next line
            line = f_cell.readline()

    # Convert to numpy array
    cells = np.array(cells)

    # Skip N = 1
    if cells.ndim == 1 :
        return

    # Extract coordinates from cells
    P = cells[:, :3]
    Q = cells[:, 3:6]

    # Interlace P and Q
    xx = [val for pair in zip(P, Q) for val in pair]

    # Convert to numpy array
    xx = np.array(xx)

    # Compute the distance between neighboring points
    d = np.sqrt( np.sum( np.square(xx[0:-1] - xx[1:]), axis=1 ) )

    # Compute the arc length
    arcLen = np.cumsum(d)
    arcLen = np.insert(arcLen, 0, 0)

    # Determine the fluctuation in R^2 ############################################################
    distance = []
    R2 = []
    for k in range(len(xx) - 1) :
        r2 = np.sum(np.square(xx[k+1::] - xx[k]), axis=1)
        ell = arcLen[k+1::] - arcLen[k]

        # Write to file
        for i in range(len(r2)) :
            distance.append(ell[i])
            R2.append(r2[i])


    # Fit the relation:
    # 2 * p * l * (1 - p / l * (1 - np.exp(- l / p)))
    p_R2, p_cov_R2 = curve_fit(R2_fun, distance, R2, p0=1, bounds=(1e-8, np.inf))

    # Determine the radius of gyration (simple) ###################################################
    # Determine center of mass
    com = np.mean(xx, axis = 0)

    # Determine squared deviation
    s2 = np.sum((xx - com)**2, axis = 1)

    # Determine radius of gyration (rms)
    Rg = np.sqrt(np.mean(s2))

    with open(path + 'PersistenceLength.txt', 'wt') as f_out :
        f_out.write('%.5f\t%5.5f\n' % (p_R2,     np.sqrt(p_cov_R2)))

    with open(path + 'GyrationRadius.txt', 'wt') as f_out :
        f_out.write('%.5f\n' % Rg)

def iterator( paths, testpath='none' ):

    # Get directory information
    cwd = os.getcwd()

    # Initially set current working directory as the path
    if testpath == 'none':
        testpath = os.getcwd() + '/'
        testpath = testpath + 'data/Chain/'

    # Get list of files in the current folder
    contents = os.listdir(testpath)

    # Does the folder contain the right files?
    if ('CellData.txt' in contents) and (not (('PersistenceLength.txt' in contents) and ('GyrationRadius.txt' in contents))) :
        if not ((testpath.find('N_1/') > -1) or (testpath.find('lysis') > -1)) :
            paths.append(testpath)

    # Continue deeper into subfolders
    for item in contents:
        nextPath = testpath + item + '/'
        if os.path.isdir(nextPath):
            newPaths = []
            iterator(newPaths, testpath=nextPath)
            for path in newPaths :
                paths.append(path)

    # Go back to home path
    os.chdir(cwd)

    # Return the paths
    return paths


def process( path ) :

    # Start the parallel pool
    p = Pool(4)

    # Get the paths
    paths = iterator([], testpath=path)

    # Start the processing
    p.map(analyzer, paths)

# Start the iterations
process(os.getcwd() + '/data/Chain/')