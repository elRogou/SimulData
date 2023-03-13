# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'

# %%
from traj.Trajectory_CG import *
from top.Toplogy_CG import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from random import randint
import random
import sys
import re
from time import time
import glob


# %%
def random_color():
    rgb = np.random.rand(3,)
    return rgb


# %%
def get_ij_distance_in_time(traj_object, idx):
    key_i = idx[0] -1
    key_j = idx[1] -1
    n_steps = int(traj_object.timesteps.n_steps/1000)
    distances = np.zeros((n_steps))
#     print(distances.shape)
    
    for i in range(n_steps):
#         print(i)
        distance = traj_object.timesteps.steps[i].coordinates[key_i].get_distance(
        traj_object.timesteps.steps[i].coordinates[key_j])
        
        distances[i] = distance
        
    return distances

# %%
def get_ij_distance_in_single_timestep(tmstp,idx):
    key_i = idx[0] - 1
    key_j = idx[1] -1 
    distance = tmstp.coordinates[key_i].get_distance(
        tmstp.coordinates[key_j])

    return distance
# %%
def plot_dist1_dist2_jointplot(x,y,title='',xlabel='',ylabel='', save=False, loc=''):
    sns.set_style("darkgrid", {"axes.facecolor": ".9"})
    g = sns.jointplot(x,y, cmap='Blues',
                 kind='kde',
                 rug=True)
    g.fig.suptitle(title)
    g.set_axis_labels(xlabel, ylabel, fontsize=25)
#     plt.suptitle(title)
    g.ax_marg_x.set_xlim(0,100)
    g.ax_marg_x.set_ylim(0,100)
    if save and loc!='':
        g.fig.savefig(loc+'joint_{}.png'.format(title))
    
    plt.show()


# %%
def plot_dist1_dist2_kdeplot(x,y,title='', save=False, loc='',idx1='idx1',idx2='idx2'):
    
    
    t = sns.kdeplot(x)
    t.figure.suptitle(title)
    if save and loc!='':
        t.figure.savefig(loc+'{}-{}_{}.png'.format(idx1[0],idx1[1],title))

    s = sns.kdeplot(y)
    s.figure.suptitle(title)
    if save and loc!='':
        s.figure.savefig(loc+'{}-{}_{}.png'.format(idx2[0],idx2[1],title))
    
    plt.show()


# %%
def plot_distance_in_trajectory(traj_objects,
                                idx,
                                baseline1 = 0,
                                baseline2 = 0,
                                legend = '',
                                title='',
                                save=False,
                                filepath=''):
    
    fig, axes = plt.subplots(figsize = (15, 7))
    
    for k in  range(len(traj_objects)):
        traj_object = traj_objects[k]
        rgb = random_color()
        array = np.zeros(traj_object.timesteps.n_steps // 1000)
        i = 0
        for ts in traj_object.timesteps.steps:
            array[i] = ts.coordinates[idx[0]-1].getDistance(ts.coordinates[idx[1]-1])
            i += 1
    

#         print(legend[k])
        axes.plot(list(range(traj_object.timesteps.n_steps // 1000)),array,
#                   label = f'{idx[0]} - {idx[1]}\n{legend}',
                  color = '#c6d7eb')
    #     axes.plot(list(range(traj_object.timesteps.n_steps // 1000)), energy)
        if baseline1:
            base = np.ones(traj_object.timesteps.n_steps // 1000) * baseline1
            axes.plot(list(range(traj_object.timesteps.n_steps // 1000)),
             base,
            linewidth = 4,
#              label = 'S',
             color = '#ff6e40',
             alpha=0.8)

        if baseline2:
            base = np.ones(traj_object.timesteps.n_steps // 1000) * baseline2
            axes.plot(list(range(traj_object.timesteps.n_steps // 1000)),
             base,
            linewidth=4,
#              label = 'R',
             color = '#1e847f',
             alpha=0.8)
            
        axes.set_ylabel(r'$\AA$', fontsize = 25)
        axes.set_xlabel('step', fontsize = 25)
        plt.title(f'Distance {idx}\n'+title, fontsize = 20)
        axes.legend(fontsize=17)
        k += 1
        if save and filepath!='':
            fig.savefig(filepath+f'{idx[0]}_{idx[1]}.png')
        else:
            pass
#             print('INCORRECT INPUT TO SAVE IMAGE')
#             sys.exit()


# %%
def plot_dbG_dihedrals_section_in_trajectory(topology, traj_object):
    for key, value in topology.dbG_dihedrals.items():
        baseline1 = value.mu1
        baseline2 = value.mu2
        
        plot_dihedral_in_trajectory(traj_object, key, baseline1, baseline2)


# %%
def get_dihed_in_trajectory(stp, key):
    dihed = calc_dihedrals(stp.coordinates[key[0]-1].pos,
                          stp.coordinates[key[1]-1].pos,
                          stp.coordinates[key[2]-1].pos,
                          stp.coordinates[key[3]-1].pos)
    if dihed < 0:
        dihed += (2 * pi)
    i += 1
    return dihed


# %%
def plot_fraction_dihedral_in_trajectory(top, traj_object, baseline1 = 0, baseline2 =0):
    frac_array_mu1 = np.zeros(len(top.dbG_dihedrals.keys()))
    frac_array_mu2 = np.zeros(len(top.dbG_dihedrals.keys()))
    
    for key in top.dbG_dihedrals.keys():
        for stp in traj_object.timesteps.steps:
            dihed = get_dihed_in_trajectory(stp, key)


# %%
def plot_dihedral_in_trajectory(traj_object, key, baseline1 = 0, baseline2 =0):
    dihed_array = np.zeros(traj_object.timesteps.steps.shape[0])
    i = 0
    for stp in traj_object.timesteps.steps:
        dihed = calc_dihedrals(stp.coordinates[key[0]-1].pos,
                              stp.coordinates[key[1]-1].pos,
                              stp.coordinates[key[2]-1].pos,
                              stp.coordinates[key[3]-1].pos)
        if dihed < 0:
            dihed += (2 * pi)
        dihed_array[i] = dihed
        i += 1

    fig, axes = plt.subplots(figsize=(15, 7))
    if baseline1:
        base = np.ones(traj_object.timesteps.n_steps // 1000) * baseline1
        axes.plot(list(range(traj_object.timesteps.n_steps // 1000)),
        base,
        label = 'expected dihed1.',
        color = 'blue',
        linewidth = 4,
        alpha = 0.5)
        
    if baseline2:
        base = np.ones(traj_object.timesteps.n_steps // 1000) * baseline2
        axes.plot(list(range(traj_object.timesteps.n_steps // 1000)),
        base,
        label = 'expected dihed2.',
        color = 'red',
        linewidth = 4,
        alpha = 0.5)
        
    axes.plot(range(dihed_array.shape[0]), dihed_array,
             label = f'ijkl = {key}', color = 'steelblue')
#     axes.plot(list(range(traj_object.timesteps.n_steps // 1000)), energy)
    axes.set_ylabel(r'Radians', fontsize = 17)
    axes.set_xlabel('step', fontsize = 17)
    plt.title('DIHEDRAL ANGLE ijkl', fontsize = 20)
    axes.legend(fontsize=17)
    plt.show()
#     r = randint(0,100)
#     plt.savefig(f'{r}-dihedral.png')
#     print(dihed_array)


# %%
def print_pymol_color_regions_command(value, color='warmpink'):
    if isinstance(value, tuple):
        msg = f'select {value[0]}, index {value[0]};        select {value[1]}, index {value[1]};        color {color}, {value[0]}; color {color}, {value[1]};'
    print(msg)


# %%
def print_pymol_show_command(value,color='skyblue'):
    
    if isinstance(value, dbG_contact):
        msg = f'select {value.key[0]}, index {value.key[0]};        select {value.key[1]}, index {value.key[1]};        distance {value.key[0]}-{value.key[1]}, {value.key[0]},{value.key[1]};'
        
    if isinstance(value, angleTrio):
        msg = f'select {value.key[0]}, index {value.key[0]};        select {value.key[1]}, index {value.key[1]};        select {value.key[2]}, index {value.key[2]};        distance {value.key[0]}-{value.key[1]}, {value.key[0]},{value.key[1]};        distance {value.key[1]}-{value.key[2]}, {value.key[1]},{value.key[2]};'
        
    elif isinstance(value, dbG_dihedralQuart):
        msg = f'select {value.key[0]}, index {value.key[0]};        select {value.key[1]}, index {value.key[1]};        select {value.key[2]}, index {value.key[2]};        select {value.key[3]}, index {value.key[3]};        color firebrick, {value.key[0]};        color firebrick, {value.key[1]};        color firebrick, {value.key[2]};        color firebrick, {value.key[3]};'
    elif isinstance(value, dihedralQuart):
        msg = f'select {value.key[0]}, index {value.key[0]};        select {value.key[1]}, index {value.key[1]};        select {value.key[2]}, index {value.key[2]};        select {value.key[3]}, index {value.key[3]};        color firebrick, {value.key[0]}, {value.key[1]};        color firebrick, {value.key[2]}, {value.key[3]};        distance {value.key[0]}-{value.key[1]}, {value.key[0]},{value.key[1]};        distance {value.key[1]}-{value.key[2]}, {value.key[1]},{value.key[2]};        distance {value.key[2]}-{value.key[3]}, {value.key[2]},{value.key[3]};'
    elif isinstance(value, tuple):
        msg = f'select {value[0]}, index {value[0]};        select {value[1]}, index {value[1]};        distance {value[0]}-{value[1]}, {value[0]},{value[1]};        color {color}, {value[0]}-{value[1]};'
    print(msg)


# %%
def dat2pdb(file,save=False):
    """Writes single .dat file with one model to .pdb file"""
    print('OPEN'.center(100,'-'),file,sep='\n')
    datLines = open(file,'r').readlines()
    pdbLines = []
    for line in datLines:
        #if 'atom positions' in line and 'is the size of chain' not in datLines[datLines.index(line)+1]:
            #print(line)
        if 'is the size of chain' in line and 'is the size of chain' not in datLines[datLines.index(line)+1]:
            coords = datLines.index(line)+1
#     print(datLines[coords])
#     print(datLines[-1])
    lines = np.arange(coords,len(datLines))
#     print(lines)
    for item in lines:
#         print(datLines[item])
#         regex = '\s+(\d+)\s+(\d+)\s+(\w+)\s+(\w+)(.*\d+\.*\s+)\s(.*\d+\.*\s)\s+(.*\d+\.*\d+\s)(.*\d+\.*\d+)'
#         match = re.findall(regex, datLines[item])
#         print(match)
#         spltLines = match[0]
#         spltLines = datLines[item].split()#
# #         if spltLines[2]=='P':
# #             break
#         i = spltLines.pop(0)#
#         spltLines.insert(3,i)#
#         print(spltLines)#
        spltLines = ParseDatCoordinates(datLines[item])
        pdbstring = 'ATOM{:>7}{:>4}{:>5}{:>6}{:>12}{:>8}{:>8}{:>6}  0.00\n'.format(*spltLines)
#         print(item,pdbstring)
        pdbLines.append(pdbstring)
    if save==True:
        pdbFile = open(file+'.pdb','w')
        for line in pdbLines:
            pdbFile.write(line)
        pdbFile.close()
        print('SAVED'.center(100,'='),file+'.pdb',sep='\n')


# %%
def ParseDatCoordinates(s):
    print(s)
    n = [5,4,4,3,8,8,8,8]
    idx_str, index_str, name, residue, x_str, y_str, z_str, mass_str     = [s[sum(n[:i]):sum(n[:i+1])].strip() for i in range(len(n))]
#     print(idx_str, index_str, name, residue, x_str, y_str, z_str, mass_str)
    idx, index, x, y, z, mass     = int(idx_str), int(index_str), float(x_str), float(y_str), float(z_str), float(mass_str)
    return idx, name, residue, index, x, y, z, mass


# %%
def tictoc(start):
    
    elapsed = (time() - start)

    if elapsed < 60:
        print(f'{elapsed:.2f} seconds')
    elif elapsed < 120*60:
        print(f'{(elapsed/60):.2f} minutes')
    else:
        print(f'{(elapsed/3600):.2f} hours')


def getLJenergy(traj_object, top, targetA, targetB):
    """
    ---------------------------------------------------------------
    |  Returns a numpy array of the total LJ interaction 
    |  energy between targetA range and targetB range.
    |  The taget ranges should be passed as tuples.
    |  
    |  EXAMPLE:
    |  targetA = (28, 52)
    |  targetB = (10,18)
    |  
    |  The totalEnergy is calculated by goin over the 
    |  contacts in the topology (top) object and if the contact
    |  satisfies the ranges, the interaction energy is calculated.
    ---------------------------------------------------------------
    """
    totalEnergy = np.zeros(len(traj_object.timesteps.steps))
    stepsCount = 0
    for tmstp in traj_object.timesteps.steps:
        tmstpEnergy = getLJenergy_for_single_step(tmstp, top, targetA, targetB)
        totalEnergy[stepsCount] = tmstpEnergy
        stepsCount += 1
    return totalEnergy

def getLJenergy_for_single_step(tmstp, top, targetA, targetB):
    Ecurrent = 0
    for key, value in top.contacts.items():
        i = key[0] -1
        j = key[1] -1

        rangeA = range(targetA[0], targetA[1]+1)
        rangeB = range(targetB[0], targetB[1]+1)

        if (i in rangeA and j in rangeB) or\
           (i in rangeB and j in rangeA):
        #    calculate the distance
            r = tmstp.coordinates[i].get_distance(tmstp.coordinates[j])
            r2 = r**2
            tempE = calculateLJenergy(r2, value)
            Ecurrent += tempE
    return Ecurrent

def calculateLJenergy(r2, value):
    sigma = value.r2
    epsC = value.epsilon
    rm2 = 1.0/r2
    rm2 = rm2*sigma
    rm10 = rm2**5

    tempE = epsC*rm10*(5*rm2-6)

    return tempE

def getTopologyFromInpFolder(inp):
    dat_files = set(glob.glob(inp+'*dat')) - set(glob.glob(inp+'*settings*')) - set(glob.glob(inp+'*random*'))

    if len(dat_files) == 1:
        top_file = topologyFile(list(dat_files)[0])
        top = top_file.createTopologyObject()
        return top
    else:
        print(dat_files,'\nZERO OR MORE THAN ONE TOPOLOGY FILE WAS FOUND\n')
        return
