# %%
from math import sqrt
from math import pi
import numpy as np
import re
import sys
from MDAnalysis.lib.distances import calc_angles
from MDAnalysis.lib.distances import calc_dihedrals
import copy
from top.HPSparams import *

# %% ---------------------------------------------------------------------------
def ParseCoordinates(s):
    n = [5,4,4,3,8,8,8,8]
    idx_str, index_str, name, residue, x_str, y_str, z_str, mass_str \
    = [s[sum(n[:i]):sum(n[:i+1])].strip() for i in range(len(n))]

    idx, index, x, y, z, mass \
    = int(idx_str), int(index_str), float(x_str), float(y_str), float(z_str), float(mass_str)
    return idx, index, name, residue, x, y, z, mass

# %% ---------------------------------------------------------------------------
class bondedPair:
    """
    
    """
#     dynamic DNA parameter values:
    default_epsilon = 100.00
    perez_epsilon = 20.00
    
    def __init__(self, idx, i, j, length, epsilon):
        self.idx     = int(idx)
        self.i       = int(i)
        self.j       = int(j)
        self.length  = float(length)
        self.epsilon = float(epsilon)
        self.key     = (self.i, self.j)
    
    def __str__(self):
        self.msg     = '{:>5}{:>5}{:>5}{:>8.3f}{:>8.3f}\n'.format(
                        self.idx, self.i, self.j,
                        self.length, self.epsilon)
        return self.msg

# %% ---------------------------------------------------------------------------
class dihedralQuart:
    """
   
    """
    default_kd = 1.0
    perez_kd = 0.3
    def __init__(self, idx,
                i, j, k, l, angle,
                kd, arg1= 0.0, arg2 = 0.0):
        self.idx   = int(idx)
        self.ijkl  = tuple(map(int, [i, j, k, l]))
        self.angle = float(angle)
        self.kd    = float(kd)
        self.arg1  = float(arg1)
        self.arg2  = float(arg2)
        self.key   = self.ijkl
        
    def __str__(self):
        self.msg   = '{:>5}{:>5}{:>5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(
                        self.idx, *self.ijkl, self.angle, self.kd, self.arg1, self.arg2)
        return self.msg

# %% ---------------------------------------------------------------------------
class angleTrio:
    """
    
    """
    default_epsilon = 100.0
    perez_epsilon = 53.0
    idr_epsilon = 10.0
    
    def __init__(self, idx, i, j, k, angle, epsilon):
        self.idx      = int(idx)
        self.ijk      = tuple(map(int, [i, j, k]))
        self.angle    = float(angle)
        self.epsilon  = float(epsilon)
        self.key      = self.ijk
    
    def __str__(self):
        self.msg      = '{:>5}{:>5}{:>5}{:>5}{:>8.3f}{:>8.3f}\n'.format(
                        self.idx, *self.ijk,
                        self.angle, self.epsilon)
        return self.msg

# %% ---------------------------------------------------------------------------
class dbG_dihedralQuart:
    """
    
    """
    def __init__(self, idx, i, j,k,l,
                eps1, mu1, sig1,
                eps2, mu2, sig2, a, phase):
        self.idx   = int(idx)
        self.i     = int(i)
        self.j     = int(j)
        self.k     = int(k)
        self.l     = int(l)
        self.eps1  = float(eps1)
        self.mu1   = float(mu1)
        self.sig1  = float(sig1)
        self.eps2  = float(eps2)
        self.mu2   = float(mu2)
        self.sig2  = float(sig2)
        self.a     = float(a) # is not used in the dihedral potential calculation
        self.key   = (self.i, self.j, self.k, self.l)
        self.phase = phase
        
    def __str__(self):
        self.msg = '{:>5}{:>5}{:>5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(
                    self.idx, self.i, self.j, self.k, self.l, self.eps1,
                    self.mu1, self.sig1, self.eps2, self.mu2,
                    self.sig2, self.a, self.phase)
        return self.msg

# %% ---------------------------------------------------------------------------
class contact:
    """
    
    """
    default_epsilon = 1.0
    r_threshold = 100.0
    stacking_r_threshold = 100.0
    def __init__(self, idx, i, j, r2, epsilon):
        self.idx     = int(idx)
        self.ij      = (int(i), int(j))
        self.r2      = float(r2)
        self.r       = sqrt(self.r2)
        self.epsilon = float(epsilon)
        self.key     = self.ij
        
    def __str__(self):
        self.msg = '{:>5}{:>5}{:>5}{:>10.3f}{:>9.6f}\n'.format(
                    self.idx, *self.ij, self.r2, self.epsilon)
        return self.msg

# %% ---------------------------------------------------------------------------
class dbG_contact:
    """
    
    """
    def __init__(self, idx, i, j, 
                eps1, mu1, sig1,
                eps2, mu2, sig2, a):
        self.idx  = int(idx)
        self.i    = int(i)
        self.j    = int(j)
        self.eps1 = float(eps1)
        self.mu1  = float(mu1)
        self.sig1 = float(sig1)
        self.eps2 = float(eps2)
        self.mu2  = float(mu2)
        self.sig2 = float(sig2)
        self.a    = float(a)
        self.key  = (self.i, self.j)
        
    def __str__(self):
        self.msg = '{:>5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(
                    self.idx, self.i, self.j, self.eps1,
                    self.mu1, self.sig1, self.eps2, self.mu2,
                    self.sig2, self.a)
        return self.msg

# %% ---------------------------------------------------------------------------
class  repulsivePair:
    """
    
    """
    default_distance = 16.0
    protDNA_distance = 32.490
    default_epsilon = 1.0
    non_repulsive_chunk = 3
    def __init__(self, idx, i, j, r2, epsilon):
        self.idx     = int(idx)
        self.i       = int(i)
        self.j       = int(j)
        self.r2      = float(r2)
        self.epsilon = float(epsilon)
        self.key     = (self.i, self.j)
        
    def __str__(self):
        self.msg = '{:>8}{:>5}{:>5}{:>8.3f}{:>8.3f}\n'.format(
                    self.idx, self.i, self.j,
                    self.r2, self.epsilon)
        return self.msg

# %% ---------------------------------------------------------------------------
class electrostaticResidue:
    """
    
    """
    def __init__(self, idx, residx, charge):
        self.idx    = int(idx)
        self.residx = int(residx)
        self.charge = float(charge)
        self.key    = int(idx)
#         print('DEBUG esRes object created')
        
    def __str__(self):
        self.msg = '{:>8}{:>5}{:>10.3f}\n'.format(
                    self.idx, self.residx, self.charge)
        return self.msg

# %% ---------------------------------------------------------------------------
class chiralConstraint:
    """
    (5I5,2F8.3)
    """
    def __init__(self, idx, i, j, k, l, angle, kd):
        self.idx   = int(idx)
        self.i     = int(i)
        self.j     = int(j)
        self.k     = int(k)
        self.l     = int(l)
        self.angle = float(angle)
        self.kd    = float(kd)
        self.key   = (self.i,self.j,self.k,self.l,)
        
    def __str__(self):
        self.msg = '{:>5}{:>5}{:>5}{:>5}{:>5}{:>8.3f}{:>8.3f}\n'.format(
                    self.idx,
                    self.i,
                    self.j,
                    self.k,
                    self.l,
                    self.angle,
                    self.kd)
        return self.msg

# %% ---------------------------------------------------------------------------
class atomPosition:
    """
    (I5,I4,A4,A3,4F8.3)
    """
    def __init__(self, idx, index, name, amino_acid,
                 xpos, ypos, zpos, mass):
        self.idx   = idx
        self.index = int(index)
        self.name  = name
        self.aa    = amino_acid
        self.xpos  = float(xpos)
        self.ypos  = float(ypos)
        self.zpos  = float(zpos)
        self.mass  = float(mass)
        self.key   = (self.idx)
        self.xyz   = np.array([self.xpos, self.ypos, self.zpos])
        
    def __str__(self):
        self.msg = '{:>5}{:>4}{:>3}{:>4}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(
                    self.idx, self.index, self.name, self.aa,
                    self.xpos, self.ypos, self.zpos,
                    self.mass)
        return self.msg
    
    def get_distance(self, position):
        distance = np.linalg.norm(self.xyz - position.xyz) # in nm
        return distance
    
    def get_angle(self, position1, position2):
        angle = calc_angles(self.xyz,
                           position1.xyz,
                           position2.xyz)
        if angle < 0:
            angle += (2 * pi)
            
        return angle
    
    def get_dihedral(self,pos1, pos2, pos3):
        dihedral = calc_dihedrals(self.xyz,
                                 pos1.xyz,
                                 pos2.xyz,
                                 pos3.xyz)
        if dihedral < 0:
            dihedral += (2 * pi)
            
        return dihedral

# %% ---------------------------------------------------------------------------
class topologyFile:
    """
    
    """
    def __init__(self, filename, dbgc = 0, dbgd = 0):
        self.filename = filename
        self.dbgc_flag = dbgc
        self.dbgd_flag = dbgd
        
    def createTopologyObject(self):
        bonds, angles,\
        dihedrals, contacts,\
        repulsive, positions,\
        dbG_dihedrals, dbG_contacts,\
        electrostatic_resid, chiral = self.__getObjects()
        
        Top = Topology(bonds = bonds, angles = angles,
                dihedrals = dihedrals, contacts = contacts,
                repulsions = repulsive, positions = positions,
                dbG_dihedrals = dbG_dihedrals, dbG_contacts = dbG_contacts, 
                electrostatic_resid = electrostatic_resid,
                chains = self.n_positions_in_chains, chiral = chiral)
        return Top
        
    def __getObjects(self):
        with open(self.filename) as top:
            bonds         = self.__readSection(top, 'bond')
            angles        = self.__readSection(top, 'angle')
            dihedrals     = self.__readSection(top, 'dihedral')
            try:
                chiral   = self.__readSection(top, 'chiral')
                if not chiral:
                    chiral = []
                    # print(f'did not find chiral')
            except:
                chiral = []
                print(f'did not find chiral')
            dbG_dihedrals = self.__readSection(top, 'dbG_dihedral')
            contacts      = self.__readSection(top, 'contacts')
            dbG_contacts  = self.__readSection(top, 'dbG_contacts')
            repulsive     = self.__readSection(top, 'repulsive')
            electrostatic_resid = self.__readSection(top, 'electrostatic')
            n_chains, n_positions = self.__assertChains(top)
            positions     = self.__readPositions(top, n_chains, n_positions)
            
        return  bonds, angles,\
                dihedrals,\
                contacts,repulsive,\
                positions,\
                dbG_dihedrals, dbG_contacts,\
                electrostatic_resid, chiral
    
    def __assertChains(self, top):
        line = top.readline()
#         print('positions debug\n', line)
#         print(line)
        assert "positions" in line
        n_positions = int(line.split()[0])
        
        line = top.readline()
        assert "chains" in line
        n_chains = int(line.split()[0])
        self.n_positions_in_chains = np.zeros(n_chains, dtype=np.int32)
#         print('DEBUG: ', self.n_positions_in_chains)
        return n_chains, n_positions
        
    def __readPositions(self, top, n_chains, n_positions):
        positions = {}
        for i in range(n_chains):
#             print(i)
            line = top.readline()
            split = line.split()
            n_pos = int(split[0])
            self.n_positions_in_chains[i] = n_pos
#             print('DEBUG: ', self.n_positions_in_chains)
            

        for i_chains in range(n_chains):
            positions[i_chains] = {}
            for i_positions in range(self.n_positions_in_chains[i_chains]):
                line = top.readline()
                split = ParseCoordinates(line)
#                 split = re.findall('.{1,8}', line)
                if not split:
                    continue
#                 print(split)
                pos = atomPosition(*split)
                positions[i_chains][pos.idx] = pos
#         print(positions)
        return positions
    
#         last_pos = top.tell()
#         top.readline()
#         line = top.readline()
#         split = line.split()
#         chains = split[0]
#         print(line)
#         for i in range(1,int(chains)+1):
#             pos = self.__readSection(top, 'chain')
#             positions[i] = pos
#         return positions
        
    def __readSection(self, top, section):
        last_pos = top.tell()
        line = top.readline()
#         print(f' debug readSections\n{line}') #DEBUG
        try:
            assert section in line
            dictionary = {}
            split = line.split()
            iterate = split[0]
#             print(f'debug iterate: {iterate}') # DEBUG
            try:
                assert iterate != 0
            except:
                print(f'{iterate} in {section}')
#                 top.seek(last_pos)
                return {}
            
        except:
            print(f'{section} section is not found')
            top.seek(last_pos)
            return {}
            
#         print('readSections debug ', iterate, section) # DEBUG
        for i in range(int(iterate)):
            line = top.readline()
            split = line.split()
#             print(line)
#             print(section)
            if section == 'bond':
                obj = bondedPair(*split)
            elif section == 'angle':
                obj = angleTrio(*split)
            elif section == 'dihedral': 
                obj = dihedralQuart(*split)
            elif section == 'dbG_dihedral':
                obj = dbG_dihedralQuart(*split)
            elif section == 'contacts':
                obj = contact(*split)
            elif section == 'dbG_contacts':
                obj = dbG_contact(*split)
            elif section == 'repulsive':
                obj = repulsivePair(*split)
            elif section == 'electrostatic':
                obj = electrostaticResidue(*split)
            elif section == 'chiral':
                obj = chiralConstraint(*split)
            else:
                raise ValueError('{} is not a recognized section'.format(section))
                
            dictionary[obj.key] = obj
        return dictionary

# %%
class Topology:
    """
    holds the following sections:
    bonds, angles, dihedrals, contacts,
    dbG_contacts, repulsions, coordinates
    """
    beads = {'protein':['CA','CB'], 'dna':['P', 'S', 'B']}
    dna_bead_mass = 1.0
    
    def __init__(self, bonds, angles,
                 dihedrals, contacts,
                 repulsions,positions,
                 chiral,
                 dbG_dihedrals = {}, dbG_contacts = {},
                 chains = 1, electrostatic_resid = {}):
        
        self.bonds         = copy.deepcopy(bonds)
        self.angles        = copy.deepcopy(angles)
        self.dihedrals     = copy.deepcopy(dihedrals)
        self.chiral        = copy.deepcopy(chiral)
        self.dbG_dihedrals = copy.deepcopy(dbG_dihedrals)
        self.contacts      = copy.deepcopy(contacts)
        self.dbG_contacts  = copy.deepcopy(dbG_contacts)
        self.repulsive     = copy.deepcopy(repulsions)
        self.electrostatic = copy.deepcopy(electrostatic_resid)
        self.n_positions_in_chains        = copy.deepcopy(chains)
        self.positions     = copy.deepcopy(positions)
        self.dna_positions = []
        
    def print(self):
        print('bonds\n')
        for value in self.bonds.values():
            print(value)
        print('angles\n')
        for value in self.angles.values():
            print(value)
        for value in self.dihedrals.values():
            print(value)
        for value in self.contacts.values():
            print(value)
        for value in self.dbG_contacts.values():
            print(value)
        for value in self.repulsive.values():
            print(value)
        for value in self.positions.values():
            print(value)

    def defineIDRs(self,ranges=[],cModel='None',structured=[],epsilon=False):
        """
        ===========================================================
        INPUT: list of *ranges*; each range is a tuple
        containing the first and the last IDR residue
        *cModel* type; can be hps or kh, anything else will not be
        accepted and the model will be contactless
        OUTPUT: no explicit output, the topology object is modified
        """
        
        """
        !! NOTE: The dihedral angles are DELETED when
        applying this function. No option to revert that back
        unless creating a new object from scratch.

        The number of ranges is not limited.
        Can be used to create IDPs as well.
        ===========================================================
        """
        if not ranges:
            print('TO DEFINE IDRs PASS BEAD RANGES'.center(60,'='))
        else:
            for item in ranges:
                for i in range(item[0],item[1]+1):
                    self.__relaxAngle(i)
                    self.__removeDihedral(i)
                    self.__removeFromContacts(i)
            
            self.updateSection(self.contacts,self.contacts)
            self.updateSection(self.dihedrals,self.dihedrals)
            

            # try:
            self.__defineHPScontacts(ranges,structured,cModel,epsilon)
            # except:
            #     print(' THE IDRs LACK CONTACTS '.center(50,'-'))
            if structured:
                self.__defineRepulsiveRanges(structured,ranges)
                pass
            self.__removeRepulsionPair()
            self.updateSection(self.repulsive,self.repulsive)

    def __relaxAngle(self,i):
        angle_key = (i,i+1,i+2)
        if angle_key not in self.angles.keys():
            print(f'Angle {angle_key} was not found\
            \nPossibly {i} is the end of the chain')
            return
        else:
            self.angles[angle_key].epsilon = angleTrio.idr_epsilon
        
    def __removeDihedral(self,i):
        # dihedral angle
        dihedral_key = (i,i+1,i+2,i+3)
        if dihedral_key not in self.dihedrals.keys():
            print(f'Possibly {i} is the end of the chain')
            return
        else:
            del self.dihedrals[dihedral_key]

    def __removeFromContacts(self,i):
        # contact
        for key in list(self.contacts.keys()):
            if i in key:
                del self.contacts[key]

    def __defineRepulsiveRanges(self,structured,ranges):
        for item in structured:
            for i in range(item[0],item[1]+1):
                if abs(i-item[1]>4):
                    for j in range(i+4,item[1]+1):
                        if (i,j) not in self.contacts.keys() and\
                            (j,i) not in self.contacts.keys():
                                repulsion_object = self.__createRepulsionsObject(i, j)
                                self.__appendObject(repulsion_object, self.repulsive)

        # # for j in range(1,len(self.bonds.keys())+1):
        # #     for item in ranges:
        # #         for i in range(item[0],item[1]+1):
        # #             if abs(i-j) > repulsivePair.non_repulsive_chunk:
        # #                 if (i,j) not in self.contacts.keys() and\
        # #                    (j,i) not in self.contacts.keys():
        # #                     repulsion_object = self.__createRepulsionsObject(i, j)
        # #                     self.__appendObject(repulsion_object, self.repulsive)
        # for item in ranges:
        #     for i in range(item[0],item[1]+1):
        #         for struct_item in structured:
        #             for j in range(struct_item[0],struct_item[1]):
        #                 if abs(i-j) > repulsivePair.non_repulsive_chunk:
        #                     if (i,j) not in self.contacts.keys() and\
        #                     (j,i) not in self.contacts.keys():
        #                         repulsion_object = self.__createRepulsionsObject(i, j)
        #                         self.__appendObject(repulsion_object, self.repulsive)


        self.updateSection(self.repulsive,self.repulsive)
    
    def __createHPScontact(self,i,j,cModel,epsilon):
        params = self.__getHPSparams(i,j,cModel,epsilon)
        contact_object = contact(*params)
        self.__appendObject(contact_object,self.contacts)
        # self.__removeRepulsionPair(i,j)

    def __removeRepulsionPair(self):
        for key in list(self.contacts.keys()):
            self.removeOneRepulsion(key[0],key[1])

    def removeOneRepulsion(self,i,j):
        if (i,j) in self.repulsive.keys():
            del self.repulsive[(i,j)]
            print(f'Repulsion pair {i} {j} was removed')
        if (j,i) in self.repulsive.keys():
            del self.repulsive[(j,i)]
            print(f'Repulsion pair {j} {i} was removed')
        else:
            print(f'Repulsion pair {i} {j} was not found')

    def __defineHPScontacts(self,ranges,structured,cModel,epsilon):
        # hps contacts of the one chain with itself
        for item in ranges:
            for i in range(item[0],item[1]+1):
                for j in range(i,item[1]+1):
                    if abs(i-j) > repulsivePair.non_repulsive_chunk:
                        self.__createHPScontact(i,j,cModel,epsilon)

        # add intersegments hps contacts
        for i_count in range(len(ranges)):
            item = ranges[i_count]
            if i_count < len(ranges)-1:
                for i in range(item[0],item[1]+1):
                    intersegment = ranges[i_count+1]
                    for j in range(intersegment[0],intersegment[1]+1):
                        if abs(i-j) > repulsivePair.non_repulsive_chunk:
                            self.__createHPScontact(i,j,cModel,epsilon)
        
        # add hps contact with structured target
        if structured:
            for item in ranges:
                for i in range(item[0],item[1]+1):
                    for struct in structured:
                        for j in range(struct[0],struct[1]+1):
                            if abs(i-j) > repulsivePair.non_repulsive_chunk:
                                self.__createHPScontact(i,j,cModel,epsilon)

    def createContact(self,i,j,epsilon,r2=0):
        idx = len(self.contacts.keys())+1
        if not r2:
            r2 = self.__getDistance(i,j)**2
        contact_object = contact(idx,i,j,r2,epsilon)
        self.__appendObject(contact_object,self.contacts)
        self.removeOneRepulsion(i,j)
        self.updateSection(self.contacts,self.contacts)
        self.updateSection(self.repulsive,self.repulsive)

    def __getDistance(self,i,j):
        for chain in self.positions.keys():
            if i in self.positions[chain].keys():
                chaini = chain
            if j in self.positions[chain].keys():
                chainj = chain
        return self.positions[chaini][i].get_distance(self.positions[chainj][j])

    def __getHPSparams(self,i,j,cModel,epsilon):
        # read the params file
        # get residue name from ccordinates
        iname = self.positions[0][i].aa
        jname = self.positions[0][j].aa
        # calculate:
        # distance sigmaij
        distance = (HPSparams[cModel][iname]['sigma'] + HPSparams[cModel][jname]['sigma']) /2
        r2 = distance*distance
        # epsilon ij
        if not epsilon:
            eps = (HPSparams[cModel][iname]['epsilon'] + HPSparams[cModel][jname]['epsilon']) /2
        if epsilon:
            eps = epsilon
        lmbd = (HPSparams[cModel][iname]['lambda'] + HPSparams[cModel][jname]['lambda']) /2
        idx = len(self.contacts.keys())+1
        weighted_epsilon = eps * lmbd
        return idx,i,j,r2,weighted_epsilon
    

    def __defineKHcontacts(self,ranges):
        pass


    
    def perezDNA(self):
        dna_chains = self.findDNAchains()
        self.__createDynamicDNAObjects(dna_chains)
        self.__createPerezObjects(dna_chains)
        self.__createDNArepulsions(dna_chains)
        self.__reduceDNAmass(dna_chains)
        
    def dynamicDNA(self):
        # check if dna in positions
        dna_chains = self.findDNAchains()
        self.__createDynamicDNAObjects(dna_chains)
#         self.__createDNAbondedPairs(dna_chains)
#         self.__getangleTrio(dna_chains)
        self.__createDNAcontacts(dna_chains)
        self.__createDNArepulsions(dna_chains)
        self.__reduceDNAmass(dna_chains)
        self.__findDNApositions()
        # create angles
        # create dihedrals
        # create contacts : stacking and base pairs
        # create repulsions
        
    def __reduceDNAmass(self, dna_chains):
        for chain in dna_chains:
            for value in self.positions[chain].values():
                if value.mass != self.dna_bead_mass:
#                     print(f'changing the dna bead mass to {self.dna_bead_mass}')
                    value.mass = self.dna_bead_mass
    
    def __createPerezObjects(self, dna_chains):
#         try:
        assert len(dna_chains) == 2, f'The number of dna chains is {len(dna_chains)}'
        self.__createPerezBonds(dna_chains)
#             self.__createPerezWCBonds(dna_chains)
#             self.__createPerezBBBonds(dna_chains)
        self.__createPerezAngles(dna_chains)
        self.__createPerezDihedrals(dna_chains)
#         except:
#             sys.exit()
        
    def __createPerezBonds(self, dna_chains):
        ################################################
        # creating the following bonds:
        # BB-WC base pairing
        # 5'-BB-3' stacking interaction
        ################################################
#             making sure there are two DNA chains
        length_chain1 = len(self.positions[dna_chains[0]].keys())
        length_chain2 = len(self.positions[dna_chains[1]].keys())
        assert length_chain1 == length_chain2,\
        'perez DNA cannot be created, chains are not equal in length'
            
        self.__createPerezWCBonds(dna_chains)
        self.__createPerezBBBonds(dna_chains)
        print('created Bonds')
        
    def __getBPairsIdx(self, count, dna_chains):
        chain1 = dna_chains[0]
        chain2 = dna_chains[1]
        
        i_start = self.positions[chain1][list(self.positions[chain1].keys())[0]].idx+2
        i_end = self.positions[chain1][list(self.positions[chain1].keys())[-1]].idx

        j_start = self.positions[chain2][list(self.positions[chain2].keys())[0]].idx+2
        j_end = self.positions[chain2][list(self.positions[chain2].keys())[-1]].idx

        i = np.arange(i_start, i_end+3,3)
        j = np.arange(j_end, j_start-3,-3)

        if count == len(i):
            return
        else:  
            key = (i[count], j[count])
        return key

    def __createPerezBBBonds(self, dna_chains):
        for chain in dna_chains:
#             The starting index is the frist base bead (shoud it be +2?)
            start_idx = list(self.positions[chain].keys())[0] +2
#             The end index is the last base bead in the chain
            end_idx = self.positions[chain][list(self.positions[chain].keys())[-1]].idx
            print(start_idx, f'in chain {chain}')
            # for the 5'-BB-3':
            for i in range(start_idx,
                           end_idx - 2,
                           3):
                assert self.positions[chain][i].name == 'B',\
                f'The {i} bead is not a base!'
                assert self.positions[chain][i+3].name == 'B',\
                f'The {i}+3 bead is not a base!'
                self.__createPerezBondedPair(self.positions[chain][i],
                                        self.positions[chain][i+3])
                    
    def __createPerezWCBonds(self,dna_chains):
#       for the BB-WC:
        chain1 = dna_chains[0]
        chain2 = dna_chains[1]
        
        start_idx = self.positions[chain1][list(self.positions[chain1].keys())[2]].idx
        end_idx = self.positions[chain1][list(self.positions[chain1].keys())[-1]].idx
        count = 0
#         first chain dihedrals
        for i in range(start_idx,
                      end_idx+1,
                      3):
            bb_pair = self.__getBPairsIdx(count, dna_chains)
#         for i in range(list(self.positions[dna_chains[0]].keys())[0] +2,\
#                       len(self.positions[dna_chains[0]].keys())+3,3):
#             i_pair = 2*len(self.positions[dna_chains[1]].keys())+7 -i
            assert self.positions[dna_chains[0]][i].name == 'B',\
            f'{i} and {bb_pair[1]} cannot be bonded | Check if base beads'
            assert self.positions[dna_chains[1]][bb_pair[1]].name == 'B',\
            f'{i} and {bb_pair[1]} cannot be bonded | Check if base beads'

            self.__createPerezBondedPair(self.positions[dna_chains[0]][i],
                                    self.positions[dna_chains[1]][bb_pair[1]])
            count += 1
            
    def __createPerezAngles(self, dna_chains):
        ################################################
        # creating the following angle trios:
        # SBB an angle along the BB-WC bond
        ################################################
        chain1 = dna_chains[0]
        chain2 = dna_chains[1]
        
        start_idx = self.positions[chain1][list(self.positions[chain1].keys())[1]].idx
        end_idx = self.positions[chain1][list(self.positions[chain1].keys())[-2]].idx
        count = 0
#         first chain dihedrals
        for i in range(start_idx,
                      end_idx+1,
                      3):
            bb_pair = self.__getBPairsIdx(count, dna_chains)
            self.__getPerezangleTrio(chain1,
                                     chain2,
                                     i,
                                     bb_pair[1])
            count += 1
            
        start_idx = self.positions[chain2][list(self.positions[chain2].keys())[1]].idx
        end_idx = self.positions[chain2][list(self.positions[chain2].keys())[-2]].idx
        count = -1
#         second chain dihedrals
        for i in range(start_idx,
                      end_idx+1,
                      3):
            bb_pair = self.__getBPairsIdx(count, dna_chains)
            self.__getPerezangleTrio(chain2,
                                     chain1,
                                     i,
                                     bb_pair[0])
            count -= 1
        
    def __createPerezDihedrals(self, dna_chains):
        ################################################
        # creating the following dihedral quartet:
        # 3'-PSBB-5' along the BB-WC bond
        # 5'-PSBB-3' along the BB-WC bond
        ################################################
        chain1 = dna_chains[0]
        chain2 = dna_chains[1]
        
        start_idx = self.positions[chain1][list(self.positions[chain1].keys())[0]].idx
        end_idx = self.positions[chain1][list(self.positions[chain1].keys())[-1]].idx
        print('START: ', start_idx, 'END:', end_idx)
        count = 0
#         first chain dihedrals
        for i in range(start_idx,
                      end_idx+1,
                      3):
            bb_pair = self.__getBPairsIdx(count, dna_chains)
            self.__getPerezdihedralQuartet(chain1,
                                          chain2,
                                          i,
                                          bb_pair[1])
            count += 1
            
        start_idx = self.positions[chain2][list(self.positions[chain2].keys())[0]].idx
        end_idx = self.positions[chain2][list(self.positions[chain2].keys())[-1]].idx
        count = -1
#         second chain dihedrals
        for i in range(start_idx,
                      end_idx+1,
                      3):
            bb_pair = self.__getBPairsIdx(count, dna_chains)
            self.__getPerezdihedralQuartet(chain2,
                                          chain1,
                                          i,
                                          bb_pair[0])
            count -= 1
            
    def __getPerezdihedralQuartet(self, chain1, chain2, i, i_pair):
        try:
            assert self.positions[chain1][i].name == 'P', f'{i} is not a P bead'
            assert self.positions[chain1][i-1].name == 'B', f'{i-1} is not a B bead'
            assert self.positions[chain1][i+2].name == 'B', f'{i+2} is not a B bead'
            assert self.positions[chain2][i_pair].name == 'B', f'{i_pair} is not a B bead'
        except:
            print(f'{i} is the begining/end of the chain')
        try:
#           5'-PSBB-3' along the BB-WC bond:
            p_s_b_b53_dihedral = self.__createPerezdihedralQuart(
                                                 self.positions[chain1][i],
                                                 self.positions[chain1][i+1],
                                                 self.positions[chain1][i+2],
                                                 self.positions[chain2][i_pair])
    #           3'-PSBB-5' along the BB-WC bond:
            b_s_b_b35_dihedral = self.__createPerezdihedralQuart(
                                                 self.positions[chain1][i],
                                                 self.positions[chain1][i-2],
                                                 self.positions[chain1][i-1],
                                                 self.positions[chain2][i_pair+3])
        except:
            print(f'{i}P is either start, or an end of the chain')
                
    def __createDynamicDNAObjects(self, dna_chains):
        for chain in dna_chains:
            start_idx = list(self.positions[chain].keys())[0]
            # print(start_idx)
            for i in range(start_idx,
                           start_idx + int(len(self.positions[chain].keys())),
                           3):
                try:
                    assert self.positions[chain][i].name == 'P'
                    assert self.positions[chain][i+1].name == 'S'
                    assert self.positions[chain][i+2].name == 'B'
#                     assert self.positions[chain][i+3].name == 'P'
#                     assert self.positions[chain][i+4].name == 'S'
                    
                    self.__createDNAbondedPairs(chain, i)
                    self.__getangleTrio(chain, i)
                    self.__getdihedralQuart(chain, i)
                    
                except:
                    print(f'{self.positions[chain][i].name} did not meet requirements')
    def findChain(self,i):
        for chain in self.positions.keys():
            for key,value in self.positions[chain].items():
                if value.index == i:
                              return chain, key
                else:
                    pass

    def findDNAchains(self):
        dna_chains = []
        for chain in self.positions.keys():
            for value in self.positions[chain].values():
                if value.name in self.beads['dna']:
                    dna_chains.append(chain)
                    break
#         print(f'chains {dna_chains} is dna')
        return dna_chains

    def __findDNApositions(self):
        self.Basebeads_positions = []
        for chain in self.positions.keys():
            for value in self.positions[chain].values():
                if value.name in self.beads['dna']:
                    self.dna_positions.append(value.index)
                if value.name == 'B':
                    self.Basebeads_positions.append(value.index)
#                     print(f'{value.name} is not a DNA bead. SKIPPING')
    
    def __createDNAbondedPairs(self, chain, i):
        p_s_bond = self.__createBondedPair(self.positions[chain][i],
                                          self.positions[chain][i+1])
        
        s_b_bond = self.__createBondedPair(self.positions[chain][i+1],
                                          self.positions[chain][i+2])
        
        try:
            s_p_bond = self.__createBondedPair(self.positions[chain][i+1],
                                              self.positions[chain][i+3])
        except:
            print(f'reached the end of the chain at {self.positions[chain][i+3]}')
        
    def __createBondedPair(self, position1, position2):
        idx = len(self.bonds.keys()) + 1
        distance = position1.get_distance(position2)
        bonded_pair_object = bondedPair(idx,
                                       position1.index,
                                       position2.index,
                                       distance,
                                       bondedPair.default_epsilon)
        self.__appendObject(bonded_pair_object, self.bonds)
        
    def __createPerezBondedPair(self, position1, position2):
        idx = len(self.bonds.keys()) + 1
        distance = position1.get_distance(position2)
        bonded_pair_object = bondedPair(idx,
                                       position1.index,
                                       position2.index,
                                       distance,
                                       bondedPair.perez_epsilon)
        self.__appendObject(bonded_pair_object, self.bonds)
        
    def __createangleTrioObject(self, position1, position2, position3):
        idx = len(self.angles.keys()) + 1
        angle = position1.get_angle(position2, position3)
        angle_trio_object = angleTrio(idx,
                                      position1.index,
                                      position2.index,
                                      position3.index,
                                      angle,
                                      angleTrio.default_epsilon)
        return angle_trio_object
    
    def __createPerezangleTrioObject(self, position1, position2, position3):
        idx = len(self.angles.keys()) + 1
        angle = position1.get_angle(position2, position3)
        angle_trio_object = angleTrio(idx,
                                      position1.index,
                                      position2.index,
                                      position3.index,
                                      angle,
                                      angleTrio.perez_epsilon)
        return angle_trio_object
    
    def __getPerezangleTrio(self, chain1, chain2, i, i_pair):
        s_b_b_angle = self.__createPerezangleTrioObject(self.positions[chain1][i],
                                                   self.positions[chain1][i+1],
                                                   self.positions[chain2][i_pair])
        self.__appendObject(s_b_b_angle, self.angles)
        
    def __getangleTrio(self, chain, i):
        p_s_b_angle = self.__createangleTrioObject(self.positions[chain][i],
                                                   self.positions[chain][i+1],
                                                   self.positions[chain][i+2])
        self.__appendObject(p_s_b_angle, self.angles)
        
        try:
            
            b_s_p_next_angle = self.__createangleTrioObject(
                                                    self.positions[chain][i+2],
                                                    self.positions[chain][i+1],
                                                    self.positions[chain][i+3])
            self.__appendObject(b_s_p_next_angle, self.angles)

            p_s_p_next_angle = self.__createangleTrioObject(
                                                    self.positions[chain][i],
                                                    self.positions[chain][i+1],
                                                    self.positions[chain][i+3])
            self.__appendObject(p_s_p_next_angle, self.angles)

            s_p_next_s_next_angle = self.__createangleTrioObject(
                                                    self.positions[chain][i+1],
                                                    self.positions[chain][i+3],
                                                    self.positions[chain][i+4])
            self.__appendObject(s_p_next_s_next_angle, self.angles)
        
        except:
            print(f'reached the end of the chain at {self.positions[chain][i+3]}')
    
    def __getdihedralQuart(self, chain, i):
        try:
            p_s_p_s_dihedral = self.__createdihedralQuart(
                                                 self.positions[chain][i],
                                                 self.positions[chain][i+1],
                                                 self.positions[chain][i+3],
                                                 self.positions[chain][i+4])
            b_s_p_s_dihedral = self.__createdihedralQuart(
                                                 self.positions[chain][i+2],
                                                 self.positions[chain][i+1],
                                                 self.positions[chain][i+3],
                                                 self.positions[chain][i+4])
        except:
            print(f'reached the end of the chain at {self.positions[chain][i+3]}')

    def __createdihedralQuart(self, pos1, pos2, pos3, pos4):
        idx = len(self.dihedrals.keys()) + 1
        dihedral = pos1.get_dihedral(pos2, pos3, pos4)
        dihedral_object = dihedralQuart(idx,
                                       pos1.index,
                                       pos2.index,
                                       pos3.index,
                                       pos4.index,
                                       dihedral,
                                       dihedralQuart.default_kd)
        self.__appendObject(dihedral_object, self.dihedrals)
        
    def __createPerezdihedralQuart(self, pos1, pos2, pos3, pos4):
        idx = len(self.dihedrals.keys()) + 1
        dihedral = pos1.get_dihedral(pos2, pos3, pos4)
        dihedral_object = dihedralQuart(idx,
                                       pos1.index,
                                       pos2.index,
                                       pos3.index,
                                       pos4.index,
                                       dihedral,
                                       dihedralQuart.perez_kd)
        self.__appendObject(dihedral_object, self.dihedrals)
    
    def __createDNAcontacts(self, dna_chains):
        self.__createDNAstacking(dna_chains)
        self.__createDNAbp(dna_chains)
        pass
        
    def __createDNAstacking(self, dna_chains):
        for chain in dna_chains:
            for key, value in self.positions[chain].items():
#                 if value.name == 'B':
                    i = key +1
                    while True:
                        if i >= len(self.positions[chain].keys()) or\
                            self.positions[chain][i].name != value.name:
                            break
                        else:
                            self.__createContactObject(value, self.positions[chain][i],
                                                        r_type = 'stacking')
                            i += 1
                            
    def __createDNAbp(self, dna_chains):
        assert len(dna_chains) <= 2
        for key1, value1 in self.positions[dna_chains[0]].items():
#             if value1.name == 'B':
            for key2, value2 in self.positions[dna_chains[1]].items():
#                     if value2.name == 'B':
                    self.__createContactObject(value1, value2)
                        
    def __createContactObject(self, position1, position2, r_type = 'bp'):
        idx = len(self.contacts.keys()) + 1
        distance = position1.get_distance(position2) ** 2
        if r_type != 'bp':
            threshold = contact.r_threshold
        else:
            threshold = contact.stacking_r_threshold
            
        if distance <= threshold:
            contact_object = contact(idx,
                                    position1.index,
                                    position2.index,
                                    distance,
                                    contact.default_epsilon)
            # append
#             print(f'FOUND CONTACT: {position1.index}-{position2.index}\n')
            if contact_object.key not in self.repulsive.keys():
                self.__appendObject(contact_object, self.contacts)
            
    def __appendObject(self, objct, dictionary):
        if objct.key not in dictionary.keys():
            dictionary[objct.key] = objct
#             print(f'{dictionary[objct.key]} was added\n{objct}')
        else:
            pass
#             print(f'{dictionary[objct.key]} is already in dictionary\n')
    
    def __createDNArepulsions(self, dna_chains):
        for chain in dna_chains:
            for key, value in self.positions[chain].items():
                for i_chain in dna_chains:
                    for i_value in self.positions[i_chain].values():
                        if (value.index, i_value.index) in [
                            self.bonds.keys(), self.contacts.keys()] or\
                            value.index in range(i_value.index - 5, i_value.index + 5):
                            pass
                        else:
                            repulsion_object = self.__createRepulsionsObject(value,
                            i_value,
                            ctype='dnadna')
                            self.__appendObject(repulsion_object, self.repulsive)
#                             print(f'found repuslion {(value.index, i_value.index)}')
                            
    def __createRepulsionsObject(self, position1, position2,ctype,distance=repulsivePair.default_distance):
        idx = len(self.repulsive.keys()) + 1
        if isinstance(position1,atomPosition):
            index1 = position1.index
            index2 = position2.index
        elif isinstance(position1,int):
            index1 = position1
            index2 = position2
        else:
            print('SOMETHING IS FUNNY ABOUT THE INPUT')
            return
        if ctype=='caca' or ctype=='dnadna':
            repulsion_distance = repulsivePair.default_distance
        elif ctype=='cadna':
            repulsion_distance = repulsivePair.protDNA_distance
        elif ctype=='custom':
            repulsion_distance = distance

        repulsion_object = repulsivePair(idx,
                                index1,
                                index2,
                                repulsion_distance,
                                repulsivePair.default_epsilon)
        return repulsion_object
    
    def updateSection(self, section, dictionary):
        section.update(dictionary)
        
        count=1
        for value in section.values():
            value.idx = count
            count+=1

    def removeContact(self, key,ctype='caca',distance=repulsivePair.default_distance):
        try:
            repulsion_object = self.__createRepulsionsObject(key[0],key[1],ctype,distance)
            self.repulsive[repulsion_object.idx] = repulsion_object
            del self.contacts[key]
        except:
            print('FAILED TO FIND THE KEY OR CREATE REPULSION OBJECT\nCHECK INPUT'.center(40))

    def writeDat(self, filename, writeEmpty=True):
        with open(filename, 'w') as outdat:
            for var in vars(self).keys():
                self.__writeSection(var, outdat, writeEmpty)
                
    
    def __writeSection(self, var, outdat, writeEmpty):
#         print('DEBUG: ',var)
        if var == 'n_positions_in_chains':
            return
            
        if var == 'positions':
#             self.__writeAtomPositions()
            total = sum(len(lst) for lst in vars(self)[var].values())
            chns = len(vars(self)[var].keys())
            outdat.write(f'{total:>8} atoms positions .\n')
#             print the number of chains
            outdat.write(f'{chns:>8} chains .\n')
    
            i = 1
            for n_chain in self.n_positions_in_chains:
                outdat.write(f'{n_chain:>8} is the size of chain{i:>12}.\n')
                i += 1
                
#             print the number of current chain positions
            for key in self.positions.keys():
                for value in self.positions[key].values():
#                     print(value)
                    outdat.write(value.__str__())
            return
            
        if isinstance(vars(self)[var],dict) and\
            len(vars(self)[var].keys()) != 0:
            outdat.write(f'{len(vars(self)[var].keys()):>8} {var}.\n')
            for value in vars(self)[var].values():
                outdat.write(value.__str__())
        elif isinstance(vars(self)[var],dict) and\
            len(vars(self)[var].keys()) == 0 and\
            writeEmpty:
                outdat.write(f'0 {var}.\n')
        elif isinstance(vars(self)[var],dict) and\
            len(vars(self)[var].keys()) == 0 and\
            not writeEmpty and\
            var in ['bonds','angles','dihedrals','contacts','repulsive','electrostatic']:
            outdat.write(f'0 {var}.\n')
                
    def returnDihedType(self, key):
        if len(key) == 2:
            dihedType = 'CACA'
        elif len(key) == 4:
            dihedType = ''
            for pos in key:
                for chain in self.positions.keys():
                    if pos in self.positions[chain].keys():
                        dihedType += self.positions[chain][pos].name
        return dihedType

    def returnBondType(self, key):
        bondType = ''
        for pos in key:
            for chain in self.positions.keys():
                if pos in self.positions[chain].keys():
                    bondType += self.positions[chain][pos].name
        if bondType == 'BB' and abs(key[0]-key[1]) == 3:
            bondType += '-WC'
        return bondType

    def returnAngleType(self, key):
        angleType = ''
        for pos in key:
            for chain in self.positions.keys():
                if pos in self.positions[chain].keys():
                    angleType += self.positions[chain][pos].name
        return angleType
    
    def adjustPerezForce(self):
        kds = {
            'BSPS': 0.3,
            'PSPS': 0.8,
            'PSBB': 0.1,
            'CACACACA': 1.0
        }
        eps_bonds = {
            'BB': 28.7,
            'BB-WC': 6.5,
            'CACA': 100,
            'SP': 50,
            'PS': 15,
            'SB': 40,
            'CACA': 100.0
        }
        eps_angles = {
            'PSB': 20,
            'BSP': 20,
            'PSP': 20,
            'SPS': 20,
            'SBB': 5,
            'CACACA': 20
        }
        for key, value in self.dihedrals.items():
            dihedType = self.returnDihedType(key)
            value.kd = kds[dihedType]
        for key, value in self.bonds.items():
            bondType = self.returnBondType(key)
            value.epsilon = eps_bonds[bondType]
        for key, value in self.angles.items():
            angleType = self.returnAngleType(key)
            value.epsilon = eps_angles[angleType]

# %%


# %%



