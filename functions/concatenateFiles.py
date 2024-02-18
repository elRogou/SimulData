import numpy as np
def concatenateFiles(listoffiles,cols):
    if cols:
        allDataArray = np.empty([0,cols])
    else:
        allDataArray = np.empty([0])
        
    for file in listoffiles:
        data = np.loadtxt(file)
        allDataArray = np.concatenate((allDataArray, data))
    return allDataArray