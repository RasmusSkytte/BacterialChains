#!/usr/bin/env python3

import os
from subprocess import call
from math import sqrt
import gzip
import numpy as np
from scipy.optimize import curve_fit
from multiprocessing import Pool

def angle_fun( l, p ) :
    return np.exp( - l / p)

def R2_fun( l, p ) :
    return 2 * p * l * (1 - p / l * (1 - np.exp(- l / p)))

def delta2_fun( L, p ) :
    return L**3 / (24 * p)

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

    # Compute the distance to initial
    distances = np.sqrt( np.sum( np.square(xx[0] - xx[1::]), axis = 1 ) )

    with open(path + 'RandomWalkDistance.txt', 'wt') as f_out :
        for d in distances :
            f_out.write('%.5f\n' % d)

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
    if ('CellData.txt' in contents) and (not 'RandomWalkDistance.txt' in contents) :
        if not ((testpath.find('N_1/') > -1) or (testpath.find('theta_0.000') > -1) or (testpath.find('lysis') > -1)) :
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