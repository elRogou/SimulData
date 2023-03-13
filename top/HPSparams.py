#%% Force-field parameters for HPS, taken from the following paper:
# =================================================================
# Dignon, Gregory L., et al. 
# "Sequence determinants of protein phase behavior from a coarse-grained model." 
# PLoS computational biology 14.1 (2018): e1005941.
# =================================================================
#AA     Mass    Charge  Sigma   Lambda  
hps = {
'ALA'  : {'sigma': 5.040,'lambda': 0.730, 'epsilon': 0.2},  
'ARG'  : {'sigma': 6.560,'lambda': 0.000, 'epsilon': 0.2},  
'ASN'  : {'sigma': 5.680,'lambda': 0.432, 'epsilon': 0.2},  
'ASP'  : {'sigma': 5.580,'lambda': 0.378, 'epsilon': 0.2},  
'CYS'  : {'sigma': 5.480,'lambda': 0.595, 'epsilon': 0.2},  
'GLN'  : {'sigma': 6.020,'lambda': 0.514, 'epsilon': 0.2},  
'GLU'  : {'sigma': 5.920,'lambda': 0.459, 'epsilon': 0.2},  
'GLY'  : {'sigma': 4.500,'lambda': 0.649, 'epsilon': 0.2},  
'HIS'  : {'sigma': 6.080,'lambda': 0.514, 'epsilon': 0.2},  
'ILE'  : {'sigma': 6.180,'lambda': 0.973, 'epsilon': 0.2},  
'LEU'  : {'sigma': 6.180,'lambda': 0.973, 'epsilon': 0.2},  
'LYS'  : {'sigma': 6.360,'lambda': 0.514, 'epsilon': 0.2},  
'MET'  : {'sigma': 6.180,'lambda': 0.838, 'epsilon': 0.2},  
'PHE'  : {'sigma': 6.360,'lambda': 1.000, 'epsilon': 0.2},  
'PRO'  : {'sigma': 5.560,'lambda': 1.000, 'epsilon': 0.2},  
'SER'  : {'sigma': 5.180,'lambda': 0.595, 'epsilon': 0.2},  
'THR'  : {'sigma': 5.620,'lambda': 0.676, 'epsilon': 0.2},  
'TRP'  : {'sigma': 6.780,'lambda': 0.946, 'epsilon': 0.2},  
'TYR'  : {'sigma': 6.460,'lambda': 0.865, 'epsilon': 0.2},  
'VAL'  : {'sigma': 5.860,'lambda': 0.892, 'epsilon': 0.2},
}

# %%
#%% The following table for the FB-HPS was taken from the following paper:
# ===================================================================
# Dannenhoffer-Lafage, Thomas, and Robert B. Best. 
# "A data-driven hydrophobicity scale for predicting 
# liquid–liquid phase separation of proteins." 
# The Journal of Physical Chemistry B 125.16 (2021): 4046-4056.
# ====================================================================
#AA     Mass    Charge  Sigma   Lambda
fb = {
'ALA'  : {'sigma': 5.040,'lambda': 0.51507, 'epsilon': 0.2}, 
'ARG'  : {'sigma': 6.560,'lambda': 0.24025, 'epsilon': 0.2}, 
'ASN'  : {'sigma': 5.680,'lambda': 0.78447, 'epsilon': 0.2}, 
'ASP'  : {'sigma': 5.580,'lambda': 0.30525, 'epsilon': 0.2}, 
'CYS'  : {'sigma': 5.480,'lambda': 0.46169, 'epsilon': 0.2}, 
'GLN'  : {'sigma': 6.020,'lambda': 0.29516, 'epsilon': 0.2}, 
'GLU'  : {'sigma': 5.920,'lambda': 0.42621, 'epsilon': 0.2}, 
'GLY'  : {'sigma': 4.500,'lambda': 1.24153, 'epsilon': 0.2}, 
'HIS'  : {'sigma': 6.080,'lambda': 0.55537, 'epsilon': 0.2}, 
'ILE'  : {'sigma': 6.180,'lambda': 0.83907, 'epsilon': 0.2}, 
'LEU'  : {'sigma': 6.180,'lambda': 0.51207, 'epsilon': 0.2}, 
'LYS'  : {'sigma': 6.360,'lambda': 0.47106, 'epsilon': 0.2}, 
'MET'  : {'sigma': 6.180,'lambda': 0.64648, 'epsilon': 0.2}, 
'PHE'  : {'sigma': 6.360,'lambda': 1.17854, 'epsilon': 0.2},
'PRO'  : {'sigma': 5.560,'lambda': 0.34128, 'epsilon': 0.2}, 
'SER'  : {'sigma': 5.180,'lambda': 0.11195, 'epsilon': 0.2}, 
'THR'  : {'sigma': 5.620,'lambda': 0.27538, 'epsilon': 0.2}, 
'TRP'  : {'sigma': 6.780,'lambda': 0.97588, 'epsilon': 0.2}, 
'TYR'  : {'sigma': 6.460,'lambda': 1.04266, 'epsilon': 0.2}, 
'VAL'  : {'sigma': 5.860,'lambda': 0.55645, 'epsilon': 0.2},
}

#%% The following table for the HPS-Urry was taken from the following paper:
# ===================================================================
##1)
# Urry, D. W.; Gowda, D. C.; Parker, T. M.; Luan, C. ?H; Reid, M. C. 
# Harris, C. M.; Pattanaik, A.; Harris, R. D. 
# Hydrophobicity Scale for Proteins Based Based on Inverse Temperature
# Transitions.
# Biopolymers 1992, 32 (9), 1243?1250.
# https://doi.org/10.1002/bip.360320913. 
##2)
#Regy, R. M., Thompson, J., Kim, Y. C., & Mittal, J. (2021).
# Improved coarse grained model for studying sequence dependent
# phase separation of disordered proteins.
# Protein Science, 30(7), 1371-1379.
# ====================================================================
"""
epsilon=[(lambda1*lambda2)/2]+delta
for fb-HPS and HPS delta=0 for Urry delta=0.08
"""

delta=-0.08
#AA   Mass    Charge  Sigma   Lambda
urry = {                                                      
'ALA'  : {'sigma': 5.040,'lambda': (delta+0.602942), 'epsilon': 0.2},
'ARG'  : {'sigma': 6.560,'lambda': (delta+0.558824), 'epsilon': 0.2},
'ASN'  : {'sigma': 5.680,'lambda': (delta+0.588236), 'epsilon': 0.2},
'ASP'  : {'sigma': 5.580,'lambda': (delta+0.294119), 'epsilon': 0.2},
'CYS'  : {'sigma': 5.480,'lambda': (delta+0.647060), 'epsilon': 0.2},
'GLN'  : {'sigma': 6.020,'lambda': (delta+0.558824), 'epsilon': 0.2},
'GLU'  : {'sigma': 5.920,'lambda': (delta+0.000000), 'epsilon': 0.2},
'GLY'  : {'sigma': 4.500,'lambda': (delta+0.573530), 'epsilon': 0.2},
'HIS'  : {'sigma': 6.080,'lambda': (delta+0.764707), 'epsilon': 0.2},
'ILE'  : {'sigma': 6.180,'lambda': (delta+0.705883), 'epsilon': 0.2},
'LEU'  : {'sigma': 6.180,'lambda': (delta+0.720589), 'epsilon': 0.2},
'LYS'  : {'sigma': 6.360,'lambda': (delta+0.382354), 'epsilon': 0.2},
'MET'  : {'sigma': 6.180,'lambda': (delta+0.676471), 'epsilon': 0.2},
'PHE'  : {'sigma': 6.360,'lambda': (delta+0.823530), 'epsilon': 0.2},
'PRO'  : {'sigma': 5.560,'lambda': (delta+0.758824), 'epsilon': 0.2},
'SER'  : {'sigma': 5.180,'lambda': (delta+0.588236), 'epsilon': 0.2},
'THR'  : {'sigma': 5.620,'lambda': (delta+0.588236), 'epsilon': 0.2},
'TRP'  : {'sigma': 6.780,'lambda': (delta+1.000000), 'epsilon': 0.2},
'TYR'  : {'sigma': 6.460,'lambda': (delta+0.897059), 'epsilon': 0.2},
'VAL'  : {'sigma': 5.860,'lambda': (delta+0.664707), 'epsilon': 0.2}
}

#%% The following table for the HPS-Urry-improved was taken from the following paper:
# ===================================================================
##1)
# Urry, D. W.; Gowda, D. C.; Parker, T. M.; Luan, C. ?H; Reid, M. C. 
# Harris, C. M.; Pattanaik, A.; Harris, R. D. 
# Hydrophobicity Scale for Proteins Based Based on Inverse Temperature
# Transitions.
# Biopolymers 1992, 32 (9), 1243?1250.
# https://doi.org/10.1002/bip.360320913. 
##2)
#Regy, R. M., Thompson, J., Kim, Y. C., & Mittal, J. (2021).
# Improved coarse grained model for studying sequence dependent
# phase separation of disordered proteins.
# Protein Science, 30(7), 1371-1379.
##3)perams were taken from here
#Joseph, J. A., Reinhardt, A., Aguirre, A., Chew, P. Y., Russell, K. O., Espinosa, J. R., ... & Collepardo-Guevara, R. (2021). 
#Physics-driven coarse-grained model for biomolecular phase separation with near-quantitative accuracy. 
#Nature Computational Science, 1(11), 732-743.
# ====================================================================
"""
epsilon=[(lambda1*lambda2)/2]+delta
for fb-HPS and HPS delta=0 for Urry delta=-0.08
"""

#AA   Mass    Charge  Sigma   Lambda
#
# urry_imp = {                                                      
# 'ALA': {'sigma': 5.04, 'lambda': 0.522942, 'epsilon': 0.2},
# 'ARG': {'sigma': 6.56, 'lambda': 0.478824, 'epsilon': 0.2},
# 'ASN': {'sigma': 5.68, 'lambda': 0.508236, 'epsilon': 0.2},
# 'ASP': {'sigma': 5.58, 'lambda': 0.214119, 'epsilon': 0.2},
# 'CYS': {'sigma': 5.48, 'lambda': 0.567060, 'epsilon': 0.2},
# 'GLN': {'sigma': 6.02, 'lambda': 0.478824, 'epsilon': 0.2},
# 'GLU': {'sigma': 5.92, 'lambda': 0.001000, 'epsilon': 0.2},
# 'GLY': {'sigma':  4.5, 'lambda': 0.493530, 'epsilon': 0.2},
# 'HIS': {'sigma': 6.08, 'lambda': 0.684707, 'epsilon': 0.2},
# 'ILE': {'sigma': 6.18, 'lambda': 0.625883, 'epsilon': 0.2},
# 'LEU': {'sigma': 6.18, 'lambda': 0.640589, 'epsilon': 0.2},
# 'LYS': {'sigma': 6.36, 'lambda': 0.302354, 'epsilon': 0.2},
# 'MET': {'sigma': 6.18, 'lambda': 0.596471, 'epsilon': 0.2},
# 'PHE': {'sigma': 6.36, 'lambda': 0.743530, 'epsilon': 0.2},
# 'PRO': {'sigma': 5.56, 'lambda': 0.678824, 'epsilon': 0.2},
# 'SER': {'sigma': 5.18, 'lambda': 0.508236, 'epsilon': 0.2},
# 'THR': {'sigma': 5.62, 'lambda': 0.508236, 'epsilon': 0.2},
# 'TRP': {'sigma': 6.78, 'lambda': 0.920000, 'epsilon': 0.2},
# 'TYR': {'sigma': 6.46, 'lambda': 0.817059, 'epsilon': 0.2},
# 'VAL': {'sigma': 5.86, 'lambda': 0.584707, 'epsilon': 0.2}
# }

#
#%% The following table for the M3 was taken from the following paper:
# ===================================================================
##1)
# Tesei, G., Schulze, T. K., Crehuet, R., & Lindorff-Larsen, K.
# (2021). Accurate model of liquid?liquid phase behavior of
# intrinsically disordered proteins from optimization
# of single-chain properties.
# Proceedings of the National Academy of Sciences, 118(44), e2111696118.?
# ====================================================================
"""
epsilon=[(lambda1*lambda2)/2]+delta
for fb-HPS and HPS delta=0 for Urry delta=0.08
"""
#AA   Mass    Charge  Sigma   Lambda
#
m3 = {                                                      
'ALA'  : {'sigma': 5.040,'lambda':	0.003075, 'epsilon': 0.2},
'ARG'  : {'sigma': 6.560,'lambda':	0.723334, 'epsilon': 0.2},
'ASN'  : {'sigma': 5.680,'lambda':	0.159612, 'epsilon': 0.2},
'ASP'  : {'sigma': 5.580,'lambda':	0.001706, 'epsilon': 0.2},
'CYS'  : {'sigma': 5.480,'lambda':	0.399824, 'epsilon': 0.2},
'GLN'  : {'sigma': 6.020,'lambda':	0.467838, 'epsilon': 0.2},
'GLU'  : {'sigma': 5.920,'lambda':	0.022450, 'epsilon': 0.2},
'GLY'  : {'sigma': 4.500,'lambda':	0.784127, 'epsilon': 0.2},
'HIS'  : {'sigma': 6.080,'lambda':	0.486967, 'epsilon': 0.2},
'ILE'  : {'sigma': 6.180,'lambda':	0.686738, 'epsilon': 0.2},
'LEU'  : {'sigma': 6.180,'lambda':	0.335171, 'epsilon': 0.2},
'LYS'  : {'sigma': 6.360,'lambda':	0.094808, 'epsilon': 0.2},
'MET'  : {'sigma': 6.180,'lambda':	0.992889, 'epsilon': 0.2},
'PHE'  : {'sigma': 6.360,'lambda':	0.870904, 'epsilon': 0.2},
'PRO'  : {'sigma': 5.560,'lambda':	0.470697, 'epsilon': 0.2},
'SER'  : {'sigma': 5.180,'lambda':	0.487225, 'epsilon': 0.2},
'THR'  : {'sigma': 5.620,'lambda':	0.273742, 'epsilon': 0.2},
'TRP'  : {'sigma': 6.780,'lambda':	0.752763, 'epsilon': 0.2},
'TYR'  : {'sigma': 6.460,'lambda':	0.984442, 'epsilon': 0.2},
'VAL'  : {'sigma': 5.860,'lambda':	0.427771, 'epsilon': 0.2},
}


# %%
HPSparams = {
'hps' : hps,
'fb' : fb,
'urry' : urry,
'm3' : m3
# 'urry_imp' : urry_imp
}
# %%
print(' THE FOLLOWING HPS-TYPE MODELS WERE LOADED: '.center(40,'='))
for key in HPSparams:
    print(f'\n{key}')
