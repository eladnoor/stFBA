import pandas as pd
import pulp
import os
import csv
import numpy as np
from cobra.flux_analysis.parsimonious import pfba
from cobra.manipulation.modify import convert_to_irreversible
import sys

import settings
import json

DG_STD_MAX = 50 # kJ/mol

CYTOPLASMIC_PH = 7.4
CYTOPLASMIC_IONIC_STRENGTH = 0.25 # in M

LB_DICT = {'c': 1e-9, 'p': 1e-9, 'e': 1e-9} # in M
UB_DICT = {'c': 1e-1, 'p': 1e-1, 'e': 1} # in M

if False:
    # override constraints for testing purposes
    LB_DICT = {'c': 1e-99, 'p': 1e-99, 'e': 1e-99} # in M
    UB_DICT = {'c': 1e99, 'p': 1e99, 'e': 1e99} # in M


sys.path.append(os.path.expanduser('~/git/equilibrator-api/'))
from equilibrator_api import ComponentContribution, Reaction
equilibrator = ComponentContribution(pH=CYTOPLASMIC_PH,
                                     ionic_strength=CYTOPLASMIC_IONIC_STRENGTH)

USE_CORE = False

if USE_CORE:
    from cobra.io import read_sbml_model
    cobra_model = read_sbml_model(settings.CORE_SBML_FNAME)
    BM_RXN = 'BIOMASS_Ecoli_core_w_GAM'
else:
    from cobra.io import load_json_model
    cobra_model = load_json_model(settings.IJO1366_JSON_FNAME)
    BM_RXN = 'BIOMASS_Ec_iJO1366_core_53p95M'

###############################################################################

def get_metabolite_df():
    bigg2kegg = []
    
    # manually add nh4_c and nh4_e (they appear as nh3 in the text
    # file BIGG_METABOLITE_FNAME, but in some models nh4 is used)
    bigg2kegg += [('nh4_c', 'nh4', 'C00014')]
    bigg2kegg += [('nh4_p', 'nh4', 'C00014')]
    bigg2kegg += [('nh4_e', 'nh4', 'C00014')]
    
    # manually add q8_c q8h2_c (ubiquinone and ubiquinol)
    # note that although in BiGG, these have a specific length (8)
    # we map them to unspecific ones in KEGG (since they don't have a
    # a specific one for ubiquinol-8 in KEGG)
    bigg2kegg += [('q8_c', 'q8', 'C00399')]
    bigg2kegg += [('q8h2_c', 'q8h2', 'C00390')]
    
    with open(settings.BIGG_METABOLITE_FNAME, 'r') as fp:
        csv_reader = csv.reader(fp, delimiter='\t')
        next(csv_reader)
        for row in csv_reader:
            bigg_id = row[0]
            universal_bigg_id = row[1]
            database_links = json.loads(row[4])
            if 'KEGG Compound' in database_links:
                for d in database_links['KEGG Compound']:
                    bigg2kegg.append(
                        (bigg_id, universal_bigg_id, d['id']))

    df = pd.DataFrame(bigg2kegg, columns=['bigg.metabolite with suffix',
                                          'bigg.metabolite', 'KEGG ID'])
    
    return df.groupby('bigg.metabolite with suffix').first()


def get_kegg_dict(met_df, coeffs):
    _df = met_df.join(coeffs, on='bigg.metabolite')
    _df = _df.loc[_df[coeffs.name] != 0, ['KEGG ID', coeffs.name]]
    
    if pd.isnull(_df['KEGG ID']).any():
        # we don't have a KEGG mapping to one of the reactants
        # and therefore we cannot have a dG'0 value
        return None
    
    # convert the DataFrame to a dictionary
    return _df.groupby('KEGG ID').sum()[_df.columns[1]].to_dict()


def get_dGr0_prime(sparse):
    try:
        r = Reaction(sparse)
        if not r.check_full_reaction_balancing():
            return np.nan, np.nan, 'unbalanbed reaction'
        dG0_prime, dG0_std = equilibrator.dG0_prime(r)
        return dG0_prime, dG0_std, ''
    except ValueError as e:
        return np.nan, np.nan, 'value error: ' + str(e)
    except KeyError as e:
        return np.nan, np.nan, 'key error: ' + str(e)


def get_pfba_fluxes(cobra_model):
    convert_to_irreversible(cobra_model)
    pfba_sol = pfba(cobra_model, solver='cglpk', already_irreversible=True)
    x_dict = dict(pfba_sol.x_dict.items())
    print("FBA JSON: %s, %f" % (pfba_sol.status, x_dict[BM_RXN]))
    
    # revert the solution to a reversible one:
    reverse_reactions = [x for x in x_dict.keys() if x.endswith('_reverse')]
    
    for rxn_id in reverse_reactions:
        fwd_id = rxn_id.replace('_reverse', '')
        x_dict[fwd_id] -= x_dict[rxn_id]
        x_dict.pop(rxn_id)
    return x_dict

###############################################################################
bigg2kegg_df = get_metabolite_df()
bigg2kegg_dict = bigg2kegg_df['KEGG ID'].to_dict()

S = settings.get_stoichiometry_from_model(cobra_model)
flux_dict = get_pfba_fluxes(cobra_model)

#%%
# construct the LP for NEST (network-embedded semi-thermodynamic analysis)
pulp_solver = pulp.GLPK_CMD(msg=0, options=['--xcheck'])
ln_conc = pulp.LpVariable.dicts('lnC', indexs=S.index,
                                lowBound=-settings.M,
                                upBound=settings.M,
                                cat=pulp.LpContinuous)


# keep only reactions that carry flux and are internal
active_rxn = [x for x in flux_dict.keys() if x[0:3] != 'EX_']
active_rxn.remove(BM_RXN)
active_rxn = [x for x in active_rxn if np.abs(flux_dict[x]) > settings.eps]
S_int = S[active_rxn]

# keep only metabolites that are not involved and active internal reaction
S_int = S_int.loc[np.abs(S_int).sum(1) > 0, :]

#%%
dir_constraints = []
dgr = {}
rxn_df = pd.DataFrame(index=S_int.columns)
rxn_df['flux'] = rxn_df.index.map(flux_dict.get)
rxn_df['flux_direction'] = rxn_df['flux'].apply(lambda x: np.sign(x) * (np.abs(x) > 1e-4))
rxn_df['dG0_prime'] = np.nan
rxn_df['dG0_prime_std'] = np.nan
rxn_df['comment'] = ''
rxn_df['constrain_dg'] = False

#%%
coeff2lpsum = lambda c: pulp.lpSum(ln_conc[met] * coeff for met, coeff in c.iteritems())

for rxn in S_int.columns:
    if rxn == 'ATPS4rpp':
        # ATP synthase is complicated since it uses the proton motive force
        # to decrease the dG'0 of the ATP synthesis reaction. We just skip
        # it for now...
        rxn_df.at[rxn, 'comment'] = 'ATP synthase is currently ignored'
        continue
    
    coeffs = S_int.loc[S_int[rxn] != 0, rxn]
    kegg_ids = coeffs.index.map(bigg2kegg_dict.get)
    if None in kegg_ids:
        rxn_df.at[rxn, 'comment'] = 'could not find KEGG mappings for all compounds'
        continue

    kegg_coeffs = pd.DataFrame(list(zip(kegg_ids, coeffs)))
    kegg_coeffs = kegg_coeffs.groupby(0).sum()
    sparse_kegg = kegg_coeffs[1].to_dict()
    
    dgr0, dgr0_std, msg = get_dGr0_prime(sparse_kegg)
    rxn_df.at[rxn, 'dG0_prime'] = dgr0
    rxn_df.at[rxn, 'dG0_prime_std'] = dgr0_std
    rxn_df.at[rxn, 'comment'] = msg

    direction = rxn_df.at[rxn, 'flux_direction']
    if direction != 0 and np.isfinite(dgr0_std) and dgr0_std < DG_STD_MAX:
        rxn_df.at[rxn, 'constrain_dg'] = True
        
        dgr[rxn] = settings.RT * coeff2lpsum(coeffs)
        c = pulp.LpConstraint(direction * dgr[rxn],
                              pulp.LpConstraintLE,
                              name='dGr_%s' % rxn,
                              rhs=-settings.eps - direction*dgr0)
        dir_constraints.append(c)

#%%
# get the formation energies of all internal metabolites:
met_df = pd.DataFrame(data=S_int.index, columns=['bigg.metabolite'])
met_df['compartment'] = met_df['bigg.metabolite'].str.rsplit('_', 1).str[1]
met_df = met_df.join(bigg2kegg_df[['KEGG ID']], on='bigg.metabolite')
met_df['ln_lb'] = met_df['compartment'].apply(lambda s: np.log(LB_DICT.get(s)))
met_df['ln_ub'] = met_df['compartment'].apply(lambda s: np.log(UB_DICT.get(s)))
met_df = met_df.set_index('bigg.metabolite')

conc_df = pd.read_csv(os.path.join(settings.DATA_DIR, 'concentration_ranges.csv'),
                      index_col=0)

for met, row in conc_df.iterrows():
    if met in met_df.index:
        met_df.loc[met, 'ln_lb'] = np.log(row['conc_lb'])
        met_df.loc[met, 'ln_ub'] = np.log(row['conc_ub'])

#%%
conc_constraints = []
for met, row in met_df.iterrows():
    ln_lb = row['ln_lb']
    ln_ub = row['ln_ub']

    if ln_lb > ln_ub:
        raise Exception('lb > ub for metabolite %s' % met)
    if ln_ub - ln_lb < settings.eps:
        # add small margin to avoid numberical problems in LP
        ln_lb -= settings.eps/2.0
        ln_ub += settings.eps/2.0
    
    const_lb = pulp.LpConstraint(ln_conc[met], pulp.LpConstraintGE,
                                 name='conc_%s_lower' % met,
                                 rhs=ln_lb)
    const_ub = pulp.LpConstraint(ln_conc[met], pulp.LpConstraintLE,
                                 name='conc_%s_upper' % met,
                                 rhs=ln_ub)
    
    conc_constraints += [const_lb, const_ub]

#%%
lp_nodir = pulp.LpProblem('NEST_nodir')
for c in conc_constraints:
    lp_nodir.addConstraint(c)
lp_nodir.writeLP(os.path.join(settings.RESULT_DIR, 'nest_nodir.lp'))

lp_nest = pulp.LpProblem('NEST')
for c in conc_constraints + dir_constraints:
    lp_nest.addConstraint(c)
lp_nest.writeLP(os.path.join(settings.RESULT_DIR, 'nest.lp'))

#%% First, check that the NEST problem is feasible (i.e. that it has at least
#   one solution)
lp_nest.solve(pulp_solver)
if lp_nest.status != pulp.LpStatusOptimal:
    print('NEST problem cannot be solved: status = %d' % lp_nest.status)
    is_nest_feasible = False
else:
    print('NEST problem is feasible!')
    is_nest_feasible = True

#%%
def get_min_and_max(lp, lp_var, name=None):
    lp.setObjective(lp_var)
    lp.sense = pulp.LpMinimize
    lp.solve(pulp_solver)
    if lp.status != pulp.LpStatusOptimal:
        print('%s, LP status %s' % (name or '', lp.status))
        min_val = np.nan
    else:
        min_val = pulp.value(lp_var)
    
    lp.sense = pulp.LpMaximize
    lp.solve(pulp_solver)
    if lp.status != pulp.LpStatusOptimal:
        print('%s: LP status %s' % (name or '', lp.status))
        max_val = np.nan
    else:
        max_val = pulp.value(lp_var)
        
    return min_val, max_val

#%% calculate the min and max dGf for each metabolite:
if is_nest_feasible:
    met_res_df = []
    for met in S_int.index:
        lnc_nest_lb, lnc_nest_ub = get_min_and_max(lp_nest, ln_conc[met], met)
        met_res_df.append((met, lnc_nest_lb, lnc_nest_ub))
    
    met_res_df = pd.DataFrame(data=met_res_df,
                              columns=['bigg.metabolite', 'lnc_nest_lb', 'lnc_nest_up'])
    met_res_df = met_res_df.join(met_df, on='bigg.metabolite')
    met_res_df = met_res_df.set_index('bigg.metabolite')
    met_res_df.round(2).to_csv(os.path.join(settings.RESULT_DIR, 'nest_met.csv'))
else:
    met_df.round(2).to_csv(os.path.join(settings.RESULT_DIR, 'nest_met.csv'))
    
#%% calculate the min and max dGr for each reaction:
rxn_nest_results = []
for rxn in S_int.columns:
    # first calculate the dGr range when directionaliry constraints are not imposed
    # (i.e. as a isolated reaction)
    if rxn_df.at[rxn, 'constrain_dg']:
        dgr_nodir_lb, dgr_nodir_ub = get_min_and_max(lp_nodir, dgr[rxn], rxn)
        if is_nest_feasible:
            dgr_nest_lb, dgr_nest_ub = get_min_and_max(lp_nest, dgr[rxn], rxn)
        else:
            dgr_nest_lb, dgr_nest_ub = None, None
        
        rxn_nest_results.append((rxn, dgr_nodir_lb, dgr_nodir_ub,
                                 dgr_nest_lb, dgr_nest_ub))
        
rxn_res_df = pd.DataFrame(rxn_nest_results,
                          columns=['bigg.reaction', 'dgr_nodir_lb',
                                   'dgr_nodir_ub', 'dgr_nest_lb', 'dgr_nest_ub'])

rxn_res_df = rxn_res_df.join(rxn_df, on='bigg.reaction', how='right')
rxn_res_df.set_index('bigg.reaction')

# we now need to add the dGr'0 to the results, since the variable we use in the
# LP is actually the difference and not dGr' itself
rxn_res_df['dgr_nodir_lb'] += rxn_res_df['dG0_prime']
rxn_res_df['dgr_nodir_ub'] += rxn_res_df['dG0_prime']
rxn_res_df['dgr_nest_lb'] += rxn_res_df['dG0_prime']
rxn_res_df['dgr_nest_ub'] += rxn_res_df['dG0_prime']

forward = rxn_res_df['dgr_nodir_lb'] < 0
backward = rxn_res_df['dgr_nodir_ub'] > 0
rxn_res_df['allowed_direction'] = forward * 1 + backward * (-1)
rxn_res_df['infeasible'] = (rxn_res_df['allowed_direction'] * rxn_res_df['flux_direction'] == -1)

rxn_res_df.round(2).to_csv(os.path.join(settings.RESULT_DIR, 'nest_rxn.csv'))
