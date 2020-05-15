#!/usr/bin/env python3

import os
from subprocess import call, STDOUT
from math import floor, sqrt, ceil
from multiprocessing import Pool

def prepareOutputFile ( f_POV, extent=50, L = 0 ) :

    # Write header to outputfile
    f_POV.write('#include "colors.inc"\n')
    f_POV.write('#include "textures.inc"\n')
    f_POV.write('#include "shapes.inc"\n\n')

    # Add camera to outputfile
    f_POV.write('camera {\n')
    f_POV.write('\tlocation <%d, %d, %d>\n' % (extent,extent,extent))
    f_POV.write('\tlook_at  <0, 0, 0>\n')
    f_POV.write('}\n\n')

    # Add light source to outputfile
    f_POV.write('light_source {\n')
    f_POV.write('\t<0, %d, %d>\n' % (extent,extent))
    f_POV.write('\tcolor White\n')
    f_POV.write('\tshadowless\n')
    f_POV.write('}\n\n')

    # Add white background
    f_POV.write('background {\n')
    f_POV.write('\tcolor rgb <1, 1, 1>')
    f_POV.write('}\n\n')

def writeDataPointToFile( f_POV, data, center, R ) :

    x1 = data[1] - center[0]
    y1 = data[2] - center[1]
    z1 = data[3] - center[2]
    x2 = data[4] - center[0]
    y2 = data[5] - center[1]
    z2 = data[6] - center[2]

    d = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

    Ax = str(x1 + 0.5*(x1-x2)/d)
    Ay = str(y1 + 0.5*(y1-y2)/d)
    Az = str(z1 + 0.5*(z1-z2)/d)

    Bx = str(x2 + 0.5*(x2-x1)/d)
    By = str(y2 + 0.5*(y2-y1)/d)
    Bz = str(z2 + 0.5*(z2-z1)/d)

    f_POV.write('object {\n')
    f_POV.write('\tRound_Cylinder(<'+Ax+','+Ay+','+Az+',>,<'+Bx+','+By+','+Bz+'>,'+str(R)+','+str(0.99*R)+',0)\n')
    f_POV.write('\ttexture {\n')
    f_POV.write('\t\tpigment { color Blue }\n')
    f_POV.write('\t}\n')
    f_POV.write('}\n')

def writeDataToFile( data, R, it = 0) :

    # Compute center
    center = [0, 0, 0]
    for d in data :
        center[0] += (d[1]+d[4])/(2*len(data))
        center[1] += (d[2]+d[5])/(2*len(data))
        center[2] += (d[3]+d[6])/(2*len(data))

    # Compute extent
    extent = 0
    for d in data :
        extent = max(extent, (d[1] - center[0])**2 + (d[2] - center[1])**2 + (d[3] - center[2])**2)
        extent = max(extent, (d[4] - center[0])**2 + (d[5] - center[1])**2 + (d[6] - center[2])**2)

    # Compute actual extent
    extent = max(5,1.5*ceil(sqrt(extent)))

    # Open output path
    f_POV = open('POVray/'+str(it)+'.pov','w')

    prepareOutputFile( f_POV, extent )
    # Write to the file
    for d in data :
        writeDataPointToFile( f_POV, d, center, R )

    # Close the output file
    f_POV.close()

def generatePOVRay( path ) :

    # Get directory information
    cwd = os.getcwd()
    os.chdir(path)

    # Seconds per simulation hour
    T_s = -1   # seconds	// -1 to plot all frames

    # Create output folder
    if not os.path.isdir('POVray') :
        os.mkdir('POVray')

    # Create output folder
    if not os.path.isdir('Visualizations') :
        os.mkdir('Visualizations')

    # Open log file to get simulation info
    log = {}
    with open(path + 'log.txt','r') as f_log:

        for line in f_log:
            if line[1] == '>':
                break
            (key, val) = line.strip().split(' = ')
            log[key] = val

    # Define R for faster computation
    R = float(log['R'])

    # Import the data
    with open('CellData.txt','r') as f_data :

        # Time step size
        dT = 0.0

        # Read the first line
        line = f_data.readline().split('\t')
        datum = [float(l.strip()) for l in line]

        # Current time
        cT = datum[0]

        # Get the path to write to
        it = 0

        # If path exists, read until path does not exist
        while os.path.isfile('Visualizations/'+str(it)+'.png') :

            # Read until next time step
            while datum[0] == cT :
                line = f_data.readline()

                if len(line) == 0 :
                    return
                datum = [float(l.strip()) for l in line.split('\t')]

            # Increase counter
            it = it + 1

            # Update current time
            cT = datum[0]


        # Output path does not exist, time to create it
        # Loop over remaining lines in the file
        data = [datum]
        while not len(line) == 0 :

            # Read next line
            line = f_data.readline()
            if len(line) == 0 :
                break
            datum = [float(l.strip()) for l in line.split('\t')]

            # If the time is the current time step
            if (datum[0] == cT) :
                data.append(datum)

            # if the time is the next time step
            elif (datum[0] > cT) :

                # Write current data to file
                writeDataToFile( data, R, it )

                # Determine the time step
                if dT == 0.0 :
                    dT = datum[0] - cT

                # Compute the next time data should be outputted
                if (T_s > 0) :
                    cT = round(10000*(cT + floor(1/(dT*30*T_s))*dT))/10000
                else :
                    cT = datum[0]

                while datum[0] < cT :
                    line = f_data.readline()
                    if len(line) == 0 :
                        break
                    datum = [float(l.strip()) for l in line.split('\t')]

                # Store the last read value in data
                data = [datum]

                # Increase counter
                it = it + 1

        # Write the last data to file
        writeDataToFile( data, R, it )

    # Render and delete the files
    it = 0
    Inputpath = 'POVray/'+str(it)+'.pov'
    while os.path.isfile(Inputpath) :
        with open(os.devnull, 'w') as f_NULL :

            # Render
            call(['povray','-H1800','-W1800','+A','-D','-ORender_'+str(it)+'.png', 'POVray/'+str(it)+'.pov'], stdout=f_NULL, stderr=STDOUT)

            # Delete
            os.remove('POVray/'+str(it)+'.pov')

            # Increase counter
            it = it + 1

            # Check next path
            Inputpath = 'POVray/'+str(it)+'.pov'

    os.rmdir('POVray')
    os.chdir(cwd)

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
    if ('Completed.txt' in contents) and ('CellData.txt' in contents) and (not 'Visualization' in contents) :
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
    p.map(generatePOVRay, paths)

# Start the iterations
process(os.getcwd() + '/data/SphericalColony/')