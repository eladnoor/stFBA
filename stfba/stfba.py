#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 15:28:15 2017

@author: noore
"""
import pandas as pd
import numpy as np
from six import iteritems
from cobra import Reaction, Metabolite
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.manipulation.modify import convert_to_irreversible

def find_egc(cobra_model, solver='cplex'):
    """
        try to locate EGCs by blocking all transport reactions and
        maximizing the objective of creating ATP from ADP + Pi
    """
    model = cobra_model.copy()
    if "EGC_tester" not in model.reactions:
        obj = Reaction("EGC_tester", lower_bound=1e-3)
        obj.add_metabolites({model.metabolites.atp_c: -1,
                             model.metabolites.h2o_c: -1,
                             model.metabolites.adp_c: 1,
                             model.metabolites.pi_c: 1})
        model.add_reaction(obj)
    model.objective = "EGC_tester"

    # disable exchange reactions and ATM maintenance
    for r in model.reactions:
        if r.id[0:3] == 'EX_':
            r.lower_bound = 0
            r.upper_bound = 0
    model.reactions.ATPM.lower_bound = 0
    model.reactions.ATPM.upper_bound = 0

    # protons are sometimes not balanced well, so we ignore them
    model.reactions.EX_h_e.lower_bound = -1000
    model.reactions.EX_h_e.upper_bound = 1000

    FBA_sol = model.optimize(solver=solver)
    FBA_sol = optimize_minimal_flux(model)
    print_solution(FBA_sol)
    return FBA_sol


def print_solution(sol):
    if sol is None:
        print("Linear problem is not feasible")
    elif sol.status == 'optimal':
        print("maximal ATP futile production: %s" % sol.f)
        for rid, flux in sol.x_dict.iteritems():
            if abs(flux) > 1e-9:
                print('%20s: %6.2f' % (rid, flux))
    else:
        print("Linear problem is %s" % sol.status)


def construct_stfba_model(cobra_model, config_fname='stfba_config.tsv'):
    """construct a semi-thermodynamic model

    This adds MILP constraints to prevent internal flux cycling and
    energy generating cycles, as done in XXXXXX.

    This must be solved with an MILP capable solver.

    """
    config_df = pd.DataFrame.from_csv('stfba_config.csv')
    config_df['ln_min_conc'] = np.log(config_df['min_concentration']) - 1e-5
    config_df['ln_max_conc'] = np.log(config_df['max_concentration']) + 1e-5
    R = 8.314e-3  # kJ/mol/K
    T = 298  # K

    model = cobra_model.copy()
    convert_to_irreversible(model)
    max_ub = max(model.reactions.list_attr("upper_bound"))
    K = 10000  # used for the 'K' parameter

    # a dict for storing S^T
    thermo_stoic = {metabolite.id: {}
                    for metabolite in model.metabolites}
    # Slice operator is so that we don't get newly added metabolites
    original_metabolites = model.metabolites[:]

    for reaction in model.reactions[:]:
        # Boundary reactions are not subjected to these constraints
        if len(reaction._metabolites) == 1:
            continue
        # populate the S^T dict
        bound_id = "thermo_bound_" + reaction.id
        for met, stoic in iteritems(reaction._metabolites):
            thermo_stoic[met.id][bound_id] = stoic

        # encode the st-FBA constrainsts:
        #   1) v <= z*v_max
        #   2) S^T * dfG <= K - K*z
        # first, we create the boolean variable (z), in the form of a
        # special reactions that can have only two flux values (0 or 1)
        reaction_ind = Reaction(reaction.id + "_indicator")
        reaction_ind.variable_kind = "integer"
        reaction_ind.upper_bound = 1

        # For the first constraint we create a "fake" metabolite
        # in the indicator reaction it has a stoichiometric coeff = v_max
        # in the original reaction it has a stoichiometric coeff = -1
        # The mass-balance constraint is set to be >= 0
        # Mathematically, this translates to:
        #   -1 * v + v_max * z >= 0
        reaction_ub = Metabolite(reaction.id + "_ind_ub")
        reaction_ub._constraint_sense = "G"
        reaction.add_metabolites({reaction_ub: -1})
        reaction_ind.add_metabolites({reaction_ub: max_ub})

        # For the second constraint we create anotehr "fake" metabolite and
        # add it also to the indicator reaction with coefficint K
        # The mass-balance constraint is set to be <= K.
        # Mathematically, this translates to:
        #   K * z <= K
        # Later, we will add the reaction Gibbs energy to the left side.
        reaction_bound = Metabolite(bound_id)
        reaction_bound._constraint_sense = "L"
        reaction_bound._bound = K
        reaction_ind.add_metabolites({reaction_bound: K})
        model.add_reaction(reaction_ind)

    # Now we create the thermodynamic variables (dfG), each will be a new
    # reaction with the corresponding metabolite name.
    for metabolite in original_metabolites:
        metabolite_var = Reaction("thermo_var_" + metabolite.id)
        if metabolite.id in config_df.index:
            dfG0 = config_df.loc[metabolite.id, 'formation_energy']
            lb = dfG0 + R*T*config_df.loc[metabolite.id, 'min_concentration']
            ub = dfG0 + R*T*config_df.loc[metabolite.id, 'max_concentration']
            if lb > ub:
                raise ValueError('''lower concentration bound is higher
                                    than upper bound''')
            metabolite_var.bounds = (lb, ub)
        else:
            # if there are no specific bounds for this metabolite, use a very
            # wide range for its formation energy
            metabolite_var.bounds = (-K, K)
        model.add_reaction(metabolite_var)

        # Finally, we add the stoichiometric coefficients of each reaction that
        # creates or uses this metabolite (i.e. the transposed stoichiometric
        # matrix).
        for k, v in iteritems(thermo_stoic[metabolite.id]):
            m = model.metabolites.get_by_id(k)
            metabolite_var.add_metabolites({m: v})

    return model
