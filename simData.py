# %%
import numpy as np
import glob
from top.Topology_CG import *
from traj.Trajectory_CG import *
from functions.CG_traj_analysis import *

class simData:
    
    def __init__(self, folder):

        self.location = folder
        self.input = folder + '/input/'
        self.output = folder + '/output/'
        self.traj = folder + '/output/Traj/'
    
        self.__createObjects()

    def __createObjects(self):
        # topology
        self.__gettopologyObject()
        # trajectory
        self.__gettrajectoryList()

    def __gettopologyObject(self):
        dat_files = set(glob.glob(self.input +'*dat'))\
             - set(glob.glob(self.input +'*settings*'))\
             - set(glob.glob(self.input +'*random*'))\
             - set(glob.glob(self.input +'*Final*'))
        try:
            assert len(dat_files) == 1
            top_file = topologyFile(list(dat_files)[0])
        except:
            print(f'FAILED TO CREATE THE TOPOLOGY\n{dat_files}')

        self.top = top_file.createTopologyObject()

    def __gettrajectoryList(self):
        self.traj_files = glob.glob(self.traj + '*.dat')
        if len(self.traj_files) == 0:
            print('NO TRAJECTORIES WERE FOUND')
        else:
            print(f'{len(self.traj_files)} trajectories were found.\n')

    def runPDAnalysis(self):
        """
        protein-DNA analysis
        includes:
        * protein-DNA interactions in time
        * dbGs in time, if any are found
        * binding analysis - 
        """
        pass

# %%