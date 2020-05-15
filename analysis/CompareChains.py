#!/usr/bin/env python3

import os
from subprocess import call, STDOUT
from math import floor, sqrt, ceil

def prepareOutputFile ( f_POV, extent=45 ) :

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

def writeDataToFile( data, R, extent = 50) :

    # Compute center
    center = [0, 0, 0]
    for d in data :
        center[0] += (d[1]+d[4])/(2*len(data))
        center[1] += (d[2]+d[5])/(2*len(data))
        center[2] += (d[3]+d[6])/(2*len(data))

    # Open output path
    f_POV = open('Render.pov', 'w')

    prepareOutputFile( f_POV, extent )
    # Write to the file
    for d in data :
        writeDataPointToFile( f_POV, d, center, R )

    # Close the output file
    f_POV.close()

def GeneratePOVRay( path, extent) :

    # Get directory information
    cwd = os.getcwd()
    os.chdir(path)

    # Seconds per simulation hour
    T_s = 2   # seconds	// -1 to plot all frames

    # Open log file to get simulation info
    log = {}
    with open(path + 'log.txt','r') as f_log:

        for line in f_log:
            (key, val) = line.strip().split(' = ')
            log[key] = val

    # Define R for faster computation
    R = float(log['R'])

    # Import the data
    with open('CellData.txt','r') as f_data :

        # Read the first line
        line = f_data.readline().split('\t')
        datum = [float(l.strip()) for l in line]

        # Loop over remaining lines in the file
        data = [datum]
        while not len(line) == 0 :

            # Read next line
            line = f_data.readline()
            if len(line) == 0 :
                break
            datum = [float(l.strip()) for l in line.split('\t')]

            # If the time is the current time step
            data.append(datum)

        # Write the last data to file
        writeDataToFile( data, R, extent )

    # Render and delete the POV file
    with open(os.devnull, 'w') as f_NULL :

        # Render
        call(['povray','-H1800','-W1800','+A','-D','-ORender.png', 'Render.pov'], stdout=f_NULL, stderr=STDOUT)

        # Delete
        os.remove('Render.pov')

    os.chdir(cwd)

def iterator(path='none', extent=50):

    # Initially set current working directory as the path
    if path == 'none':
        print('No path given! Exiting...')
        return

    # Check if path exists
    if not os.path.exists(path) :
        return

    # Get list of files in the current folder
    contents = os.listdir(path)

    # Does the folder contain the right content?
    if ('Completed.txt' in contents) and ('CellData.txt' in contents) and (not 'Render.png' in contents) :
        GeneratePOVRay(path, extent)

# Start the iterations
iterator(os.getcwd() + '/data/BendingAngleValidation/theta_0.300_pi/repeat_6/', 50)
iterator(os.getcwd() + '/data/BendingAngleValidation/theta_0.600_pi/repeat_6/', 50)
iterator(os.getcwd() + '/data/BendingAngleValidation/theta_0.900_pi/repeat_2/', 50)
iterator(os.getcwd() + '/data/SphericalColony/N_100/repeat_0/', 30)