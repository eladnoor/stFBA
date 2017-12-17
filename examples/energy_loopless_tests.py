# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:43:05 2016

@author: noore
"""

from cobra.flux_analysis import loopless
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.json import load_json_model
from cobra import Metabolite, Reaction
from find_energy_generating_cycle import construct_loopless_model
#from cobra.flux_analysis.loopless import construct_loopless_model

solver = 'cplex'


#model_fname = 'e_coli_core.xml.gz'
model_fname = 'iAF1260.xml.gz'
model = create_cobra_model_from_sbml_file(model_fname)

model_fname = 'iJO1366.json'
model = load_json_model(model_fname)

model.optimize(solver=solver)
print "original model and FBA: ", model.solution

# set the model to anaerobic conditions, in order to make ATP a limiting
# factor for the maximal yield
model.reactions.EX_o2_e.lower_bound = 0

model.optimize(solver=solver)
print "anaerobic model and FBA: ", model.solution


# add a futile cycle that generates ATP    
na1_e = Metabolite('M_na1_e', formula='Na', name='na1_e', compartment='C_e')
na1_c = Metabolite('M_na1_c', formula='Na', name='na1_c', compartment='C_c')
ser_L_e = Metabolite('M_ser_L_e', formula='C3H7NO3', name='ser_L_e', compartment='C_e')
ser_L_c = Metabolite('M_ser_L_c', formula='C3H7NO3', name='ser_L_c', compartment='C_c')
cys_L_e = Metabolite('M_cys_L_e', formula='C3H7NO2S', name='cys_L_e', compartment='C_e')
cys_L_c = Metabolite('M_cys_L_c', formula='C3H7NO2S', name='cys_L_c', compartment='C_c')
gly_e = Metabolite('M_gly_e', formula='C2H5NO2', name='gly_e', compartment='C_e')
gly_c = Metabolite('M_gly_c', formula='C2H5NO2', name='gly_c', compartment='C_c')

model.add_metabolites([na1_e, na1_c, ser_L_e, ser_L_c, cys_L_e, cys_L_c, gly_e, gly_c])
symporter1_reaction = Reaction('R_SER_CYS_SYMPORTER')
symporter1_reaction.add_metabolites({na1_c: -1, na1_e: 1, ser_L_e: -1, ser_L_c: 1, cys_L_e: 1, cys_L_c:-1})
symporter1_reaction.lower_bound = -1000
symporter2_reaction = Reaction('R_GLY_CYS_SYMPORTER')
symporter2_reaction.add_metabolites({gly_e: -1, gly_c: 1, cys_L_c:-1, cys_L_e: 1})
symporter2_reaction.lower_bound = -1000
symporter3_reaction = Reaction('R_GLY_SER_SYMPORTER')
symporter3_reaction.add_metabolites({gly_e: -1, gly_c: 1, ser_L_c:-1, ser_L_e: 1})
symporter3_reaction.lower_bound = -1000
atpase_reaction = Reaction('R_ATPASE_NA')
atpase_reaction.add_metabolites({na1_c: -1, na1_e:1,
                                 model.metabolites.adp_c:-1, model.metabolites.pi_c:-1,
                                 model.metabolites.atp_c:1, model.metabolites.h2o_c:1})
atpase_reaction.lower_bound = -1000
model.add_reactions([symporter1_reaction, symporter2_reaction,
                     symporter3_reaction, atpase_reaction])

#model.optimize(solver='cplex')
#print "model with futile cycle and FBA: ", model.solution

loopless_model = construct_loopless_model(model)
loopless_model.optimize(solver=solver)

print "model with futile cycle and ll-FBA: ", loopless_model.solution

R = 8.314e-3 # kJ/mol/K
T = 298 # K
loopless_model.reactions.thermo_var_atp_c.lower_bound = -2295.8/(R*T)
loopless_model.reactions.thermo_var_atp_c.upper_bound = -2295.8/(R*T)
loopless_model.reactions.thermo_var_adp_c.lower_bound = -1423.6/(R*T)
loopless_model.reactions.thermo_var_adp_c.upper_bound = -1423.6/(R*T)
loopless_model.reactions.thermo_var_pi_c.lower_bound = -1073.3/(R*T)
loopless_model.reactions.thermo_var_pi_c.upper_bound = -1073.3/(R*T)
loopless_model.reactions.thermo_var_h2o_c.lower_bound = -157.6/(R*T)
loopless_model.reactions.thermo_var_h2o_c.upper_bound = -157.6/(R*T)
loopless_model.optimize(solver=solver)
print "model with futile cycle and ll-FBA + ATP constraints: ", loopless_model.solution
