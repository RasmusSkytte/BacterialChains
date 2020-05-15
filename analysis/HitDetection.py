#!/usr/bin/env python3

import os
from subprocess import call
from math import sqrt
import gzip
from multiprocessing import Pool

def HitDetector ( cells, phage, R ):

    # Store phage coordinates in S vector
    xS = phage[0]
    yS = phage[1]
    zS = phage[2]

    # Check for overlap
    c = 1
    for cell in cells :

        # Store cell coordinates in Q vector and P vector
        xQ = cell[0]
        yQ = cell[1]
        zQ = cell[2]
        xP = cell[3]
        yP = cell[4]
        zP = cell[5]

        # Vector between Q and P
        xQP = xQ - xP
        yQP = yQ - yP
        zQP = zQ - zP

        # Vector between P and S
        xSP = xS - xP
        ySP = yS - yP
        zSP = zS - zP

        # Compute the projection number u
        u = (xQP*xSP + yQP*ySP + zQP*zSP) / (xQP*xQP + yQP*yQP + zQP*zQP)

        # If u lands outside segment, convert back to line segment
        u = min(1, max(u, 0))

        # Calculate the shortest distance vector
        xU = xSP - u*xQP
        yU = ySP - u*yQP
        zU = zSP - u*zQP

        # Calculate the distance
        d = sqrt( xU*xU + yU*yU + zU*zU )

        # If hit is detected, return which cell was hit
        if d < R :
            return c

        # Increment c
        c = c + 1

    # If no hit is detected, return -1
    return -1


def analyzer( path ) :

    if path.find('lysis') > -1 :
        filter = False
    else :
        filter = True

    # Open log file to get simulation space size
    log = {}
    with open(path + 'log.txt', 'r') as f_log:

        for line in f_log:
            if not line[0:2] == '>>':
                (key, val) = line.strip().split(' = ')
                log[key] = val

    # Define R for faster computation
    R = float(log['R'])

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
    with open(path + 'CellData.txt', 'r') as f_cell:

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

    # Determine the center of the colony
    center = [0, 0, 0]
    for cell in cells :
        center[0] += (cell[0] + cell[3]) / (2 * len(cells))
        center[1] += (cell[1] + cell[4]) / (2 * len(cells))
        center[2] += (cell[2] + cell[5]) / (2 * len(cells))

    # Analyze the data
    with open(path + 'PhageData.txt', 'rt') as f_phage :
        with open(path + 'PhageLocation.txt', 'wt') as f_out :

            # Loop over lines in PhageData.txt
            line = f_phage.readline()
            while not len(line) == 0 :

                # Get data from line
                datum = [float(l.strip()) for l in line.split('\t')]

                # Detect hits
                c = HitDetector( cells, datum[1:], R )

                if filter and c > 0 :
                    # Write location of filtered phage to file
                    f_out.write('%.8f\t%.8f\t%.8f\t%d\n' % (datum[1], datum[2], datum[3], c))
                elif not filter :
                    # Write hit location of all phage to file
                    f_out.write('%d\n' % c)

                # Read next line
                line = f_phage.readline()


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
    if ('Completed.txt' in contents) and ('log.txt' in contents) and ('PhageData.txt' in contents) and (not 'PhageLocation.txt' in contents) :
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
process(os.getcwd() + '/data/WellMixed/')
process(os.getcwd() + '/data/Chain/')
process(os.getcwd() + '/data/SphericalColony/')