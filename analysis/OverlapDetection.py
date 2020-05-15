#!/usr/bin/env python3

import os
from subprocess import call
from math import sqrt
import gzip
from multiprocessing import Pool

def OverlapDetector ( i, cells, R ):

    # Store cell coordinates in R and S vectors
    xR = cells[i][0]
    yR = cells[i][1]
    zR = cells[i][2]

    xS = cells[i][3]
    yS = cells[i][4]
    zS = cells[i][5]

    # Check for overlap
    overlap = 0
    for j in range(len(cells)) :

        if j == i :
            continue

        # Store cell coordinates in Q vector and P vector
        xP = cells[j][0]
        yP = cells[j][1]
        zP = cells[j][2]
        xQ = cells[j][3]
        yQ = cells[j][4]
        zQ = cells[j][5]

        # Vector between Q and P
        xQP = xQ - xP
        yQP = yQ - yP
        zQP = zQ - zP

        # Vector between S and R
        xSR = xS - xR
        ySR = yS - yR
        zSR = zS - zR

        # Vector between P and R
        xPR = xP - xR
        yPR = yP - yR
        zPR = zP - zR

        # Compute the relevant dot products
        dQP_QP = xQP*xQP + yQP*yQP + zQP*zQP
        dSR_SR = xSR*xSR + ySR*ySR + zSR*zSR
        dQP_SR = xQP*xSR + yQP*ySR + zQP*zSR
        dQP_PR = xQP*xPR + yQP*yPR + zQP*zPR
        dSR_PR = xSR*xPR + ySR*yPR + zSR*zPR

        # Compute denominator
        h = dQP_QP*dSR_SR - dQP_SR*dQP_SR
        hu = h
        hv = h

        # Compute u and v
        u = dQP_SR*dSR_PR - dSR_SR*dQP_PR
        v = dQP_QP*dSR_PR - dQP_SR*dQP_PR

        # If u or v lands outside segment, convert back to line segment
        if (u < 0.0) :
            u  = 0.0
            v  = dSR_PR
            hv = dSR_SR
        elif (u > h) :
            u  = h
            v  = dSR_PR + dQP_SR
            hv = dSR_SR

        if (v < 0.0) :
            v = 0.0

            if (dQP_PR > 0.0) :
                u = 0.0
            elif (-dQP_PR > dQP_QP) :
                u = h
            else :
                u  = -dQP_PR
                hu = dQP_QP
        elif (v > hv) :
            v = hv

            if ((dQP_SR - dQP_PR) < 0.0) :
                u = 0
            elif ((dQP_SR - dQP_PR) > dQP_QP) :
                u = h
            else :
                u  = dQP_SR - dQP_PR
                hu = dQP_QP

        # Compute u and v
        u = u / hu
        v = v / hv

        # Compute the shortest distance vector
        xU = xPR + u * xQP - v * xSR
        yU = yPR + u * yQP - v * ySR
        zU = zPR + u * zQP - v * zSR

        # Calculate the distance
        d = sqrt( xU*xU + yU*yU + zU*zU )

        # Store overlap if it is the largest
        if 2 * R - d > overlap :
            overlap = 2 * R - d

    # Return the measured overlap
    return overlap

def OverlapDetector_EndToEnd ( cell, cells, R ):

    # Store cell coordinates in S vector
    xS = cell[3]
    yS = cell[4]
    zS = cell[5]

    # Check for overlap
    overlap = 0
    for c in cells :

        # Store cell coordinates in P vector
        xP = c[0]
        yP = c[1]
        zP = c[2]

        # Vector between P and S
        xU = xP - xS
        yU = yP - yS
        zU = zP - zS

        # Calculate the distance
        d = sqrt( xU*xU + yU*yU + zU*zU )

        # Store overlap if it is the largest
        if 2 * R - d > overlap :
            overlap = 2 * R - d

    # Return the measured overlap
    return overlap


def analyzer( path ) :

    # Open log file to get simulation space size
    log = {}
    with open(path + 'log.txt', 'r') as f_log:

        for line in f_log:
            if not line[0:2] == '>>':
                (key, val) = line.strip().split(' = ')
                log[key] = val

    # Define R for faster computation
    theta = float(log['theta'])
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


    # Determine overlaps
    with open(path + 'Overlaps.txt', 'wt') as f_out :
        maxOverlap = 0
        for i in range(len(cells)):
            if theta == 0:
                overlap = OverlapDetector_EndToEnd(i, cells, R)
            else :
                overlap = OverlapDetector(i, cells, R)

            if overlap > maxOverlap :
                maxOverlap = overlap

            f_out.write('%.8f\n' % overlap)

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
    if ('Completed.txt' in contents) and ('log.txt' in contents) and ('CellData.txt' in contents) and (not 'Overlaps.txt' in contents) :
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