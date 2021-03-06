# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:43:05 2016

@author: noore
"""

from cobra.flux_analysis import loopless
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.json import load_json_model
from cobra import Metabolite, Reaction
from stfba.find_energy_generating_cycle import stFBA
from cobra.flux_analysis.loopless import construct_loopless_model
from stfba import settings

R = 8.314e-3 # kJ/mol/K
T = 298 # K

model_name_to_obj = {'e_coli_core': 'BIOMASS_Ecoli_core_w_GAM',
                     'iAF1260':     'BIOMASS_Ec_iAF1260_core_59p81M',
                     'iJO1366':     'BIOMASS_Ec_iJO1366_core_53p95M'}

model_name = 'iJO1366'

model = load_json_model('model/%s.json' % model_name)

sol = model.optimize(solver=settings.LP_SOLVER)
print("original model and FBA: %.3f" % sol.f)

# set the model to anaerobic conditions, in order to make ATP a limiting
# factor for the maximal yield
model.reactions.EX_o2_e.lower_bound = 0

sol = model.optimize(solver=settings.LP_SOLVER)
print("anaerobic model and FBA: %.3f" % sol.f)

loopless_model = construct_loopless_model(model)
sol = loopless_model.optimize(solver=settings.LP_SOLVER)

print("anaerobic model ll-FBA: %.3f" % sol.f)

loopless_model.reactions.thermo_var_atp_c.lower_bound = -2295.8/(R*T)
loopless_model.reactions.thermo_var_atp_c.upper_bound = -2295.8/(R*T)
loopless_model.reactions.thermo_var_adp_c.lower_bound = -1423.6/(R*T)
loopless_model.reactions.thermo_var_adp_c.upper_bound = -1423.6/(R*T)
loopless_model.reactions.thermo_var_pi_c.lower_bound = -1073.3/(R*T)
loopless_model.reactions.thermo_var_pi_c.upper_bound = -1073.3/(R*T)
loopless_model.reactions.thermo_var_h2o_c.lower_bound = -157.6/(R*T)
loopless_model.reactions.thermo_var_h2o_c.upper_bound = -157.6/(R*T)
sol = loopless_model.optimize(solver=settings.LP_SOLVER)
print("anaerobic model with ll-FBA + ATP constraints: %.3f" % 
      sol.f)

###############################################################################
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
                                 model.metabolites.adp_c:-1,
                                 model.metabolites.pi_c:-1,
                                 model.metabolites.atp_c:1,
                                 model.metabolites.h2o_c:1})
atpase_reaction.lower_bound = -1000
model.add_reactions([symporter1_reaction, symporter2_reaction,
                     symporter3_reaction, atpase_reaction])

sol = model.optimize(solver=settings.LP_SOLVER)
print("model with futile cycle and FBA: %.3f" % sol.f)

loopless_model = construct_loopless_model(model)
sol = loopless_model.optimize(solver=settings.LP_SOLVER)

print("model with futile cycle and ll-FBA: %.3f" % sol.f)

loopless_model.reactions.thermo_var_atp_c.lower_bound = -2295.8/(R*T)
loopless_model.reactions.thermo_var_atp_c.upper_bound = -2295.8/(R*T)
loopless_model.reactions.thermo_var_adp_c.lower_bound = -1423.6/(R*T)
loopless_model.reactions.thermo_var_adp_c.upper_bound = -1423.6/(R*T)
loopless_model.reactions.thermo_var_pi_c.lower_bound = -1073.3/(R*T)
loopless_model.reactions.thermo_var_pi_c.upper_bound = -1073.3/(R*T)
loopless_model.reactions.thermo_var_h2o_c.lower_bound = -157.6/(R*T)
loopless_model.reactions.thermo_var_h2o_c.upper_bound = -157.6/(R*T)
sol = loopless_model.optimize(solver=settings.LP_SOLVER)
print("model with futile cycle and ll-FBA + ATP constraints: %.3f" % 
      sol.f)
