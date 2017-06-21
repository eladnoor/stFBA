#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 12:07:00 2017

@author: noore
"""

from cobra.io.json import load_json_model
from cobra import Metabolite, Reaction
from find_energy_generating_cycle import find_egc, construct_stfba_model
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.flux_analysis.loopless import construct_loopless_model
from cobra.io.sbml import create_cobra_model_from_sbml_file

solver = 'cplex'

model_name_to_obj = {'e_coli_core': 'BIOMASS_Ecoli_core_w_GAM',
                     'iAF1260':     'BIOMASS_Ec_iAF1260_core_59p81M',
                     'iJO1366':     'BIOMASS_Ec_iJO1366_core_53p95M'}

model_name = 'iJO1366'

model = create_cobra_model_from_sbml_file('model/%s.xml.gz' % model_name)
model.objective = model_name_to_obj[model_name]
model.reactions.EX_glc__D_e.lower_bound = -10
model.reactions.EX_o2_e.lower_bound = 0
model.reactions.EX_o2_e.upper_bound = 0

if model_name == 'iJO1366':
    #model.reactions.SPODM.lower_bound = 0
    #model.reactions.SPODM.upper_bound = 0
    model.reactions.MOX.lower_bound = 0

print "Calculating anaerobic yield with original %s model" % model_name

print "FBA max yield: ",
fba_sol = model.optimize(solver=solver)
print fba_sol.f

#print "ll-FBA max yield: ",
#loopless_model = construct_loopless_model(model)
#llfba_sol = loopless_model.optimize(solver=solver)
#print llfba_sol.f
#
#print "st-FBA max yield: ",
#stfba_model = construct_stfba_model(model)
#stfba_sol = stfba_model.optimize(solver=solver)
#print stfba_sol.f

# Make sure there are no EGCs in the core model before we start
print "\nLooking for EGCs in the original model..."
sol = find_egc(model)

if sol is not None:
    print "There is already at least one EGC in the model!"
else:
    print "\nAdding an EGC and trying again..."
    na1_e = Metabolite('M_na1_e', formula='Na',
                       name='na1_e', compartment='C_e')
    na1_c = Metabolite('M_na1_c', formula='Na',
                       name='na1_c', compartment='C_c')
    ser_L_e = Metabolite('M_ser_L_e', formula='C3H7NO3',
                         name='ser_L_e', compartment='C_e')
    ser_L_c = Metabolite('M_ser_L_c', formula='C3H7NO3',
                         name='ser_L_c', compartment='C_c')
    cys_L_e = Metabolite('M_cys_L_e', formula='C3H7NO2S',
                         name='cys_L_e', compartment='C_e')
    cys_L_c = Metabolite('M_cys_L_c', formula='C3H7NO2S',
                         name='cys_L_c', compartment='C_c')
    gly_e = Metabolite('M_gly_e', formula='C2H5NO2',
                       name='gly_e', compartment='C_e')
    gly_c = Metabolite('M_gly_c', formula='C2H5NO2',
                       name='gly_c', compartment='C_c')

    model.add_metabolites([na1_e, na1_c, ser_L_e, ser_L_c, cys_L_e,
                           cys_L_c, gly_e, gly_c])
    symporter1_reaction = Reaction('SER_CYS_SYMPORTER')
    symporter1_reaction.add_metabolites({na1_c: -1, na1_e: 1,
                                         ser_L_e: -1, ser_L_c: 1,
                                         cys_L_e: 1, cys_L_c: -1})
    symporter2_reaction = Reaction('GLY_CYS_SYMPORTER')
    symporter2_reaction.add_metabolites({gly_e: -1, gly_c: 1,
                                         cys_L_c: -1, cys_L_e: 1})
    symporter3_reaction = Reaction('GLY_SER_SYMPORTER')
    symporter3_reaction.add_metabolites({gly_e: -1, gly_c: 1,
                                         ser_L_c: -1, ser_L_e: 1})
    atpase_reaction = Reaction('ATPASE_NA')
    atpase_reaction.add_metabolites({na1_c: -1, na1_e: 1,
                                     model.metabolites.adp_c: -1,
                                     model.metabolites.pi_c: -1,
                                     model.metabolites.atp_c: 1,
                                     model.metabolites.h2o_c: 1})
    symporter1_reaction.lower_bound = -1000
    symporter2_reaction.lower_bound = -1000
    symporter3_reaction.lower_bound = -1000
    atpase_reaction.lower_bound = -1000
    model.add_reactions([symporter1_reaction, symporter2_reaction,
                         symporter3_reaction, atpase_reaction])

    # Try to find the EGC that we just added to the model
    sol = find_egc(model)

print "\nAdding STFBA constraints..."
stfba_model = construct_stfba_model(model)
sol = find_egc(stfba_model, already_irreversible=True)

print "Calculating anaerobic yield with updated %s model" % model_name

print "FBA max yield: ",
fba_sol = model.optimize(solver=solver)
print fba_sol.f

print "ll-FBA max yield: ",
loopless_model = construct_loopless_model(model)
llfba_sol = loopless_model.optimize(solver=solver)
print llfba_sol.f

print "st-FBA max yield: ",
stfba_model = construct_stfba_model(model)
stfba_sol = stfba_model.optimize(solver=solver)
print stfba_sol.f
