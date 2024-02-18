import numpy as np

def get_center_of_mass(data):
    # calculate the center of mass of the data
    # data is a 2d array of shape (n, 3)
    # where n is the number of points and 3 is the
    # number of coordinates
    # return the center of mass as a 1d array of shape (3,)
    com = np.mean(data, axis=0)
    return com

def get_closest_point(data, com):
    # find the point in data that is closest to the center of mass
    # data is a 2d array of shape (n, 3)
    # where n is the number of points and 3 is the
    # number of coordinates
    # com is the center of mass of the data
    # return the index of the closest point
    dist = np.linalg.norm(data - com, axis=1)
    closest = np.argmin(dist)
    return closest

def get_shift(point1, point2):
    # calculate the shift between two points
    # point1 and point2 are 1d arrays of shape (3,)
    # return the shift as a 1d array of shape (3,)
    shift = point1 - point2
    return shift
    
def get_shifted_data(data, shift):
    # shift the data by the shift
    # data is a 2d array of shape (n, 3)
    # where n is the number of points and 3 is the
    # number of coordinates
    # shift is a 1d array of shape (3,)
    # return the shifted data
    newdata = data + shift
    return newdata

def writeCOMpositions(positions,filename):
    # write the positions to a file
    # positions is a 2d array of shape (n, 3)
    # where n is the number of points and 3 is the
    # number of coordinates
    # filename is the name of the file to write
    if filename.endswith('.txt'):
        np.savetxt(filename, positions)
    elif filename.endswith('.pdb'):
        # createPDBfile(filename,positions)
        write_pdb(positions, filename)

def createPDBfile(filename,positions):
    templt = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
    with open (filename,'w') as f:
        for idx in range(len(positions)):
            #f.write(f'MDL    {idx+1}\n')
            line = templt.format(
                    'ATOM',
                    idx+1,
                    'CA',
                    '',
                    'GLY',
                    'A',
                    1,
                    '',
                    positions[idx][0],
                    positions[idx][1],
                    positions[idx][2],
                    0,
                    0,
                    'CA',
                    '')
            f.write(line)
            #f.write('ENDMDL\n')
    f.close()

def write_pdb(coords, filename):
    """Write a pdb file from a numpy array of coordinates"""
    with open(filename, 'w') as f:
        for i, coord in enumerate(coords):
            print(coord)
            f.write(f"ATOM {i+1:5d}  CA  ALA A   1    {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00\n")
