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

CYTOPLASMIC_PH = 7.4
CYTOPLASMIC_IONIC_STRENGTH = 0.25 # in M
CYTOPLASMIC_CONC_LB = 1e-12 # in M
CYTOPLASMIC_CONC_UB = 0.1 # in M

EXTERNAL_CONC_LB = 1e-20 # in M
EXTERNAL_CONC_UB = 1e20 # in M

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
    
    # manually add nh4_c and nh4_e
    bigg2kegg += [('nh4_c', 'nh4', 'C00014')]
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


def get_dGr0_prime(df):

    if pd.isnull(df['KEGG ID']).any():
        # we don't have a KEGG mapping to one of the reactants
        # and therefore we cannot have a dG'0 value
        return np.nan, np.nan
    
    # convert the DataFrame to a dictionary
    sparse = df.groupby('KEGG ID').sum()[df.columns[1]].to_dict()
    
    try:
        r = Reaction(sparse)
        dG0_prime, dG0_std = equilibrator.dG0_prime(r)
        return dG0_prime, dG0_std
    except ValueError as e:
        print('value error: ' + str(e))
        return np.nan, np.nan
    except KeyError as e:
        print('key error: ' + str(e))
        return np.nan, np.nan


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
# get the formation energies of all internal metabolites:
met_df = pd.DataFrame(data=S_int.index, columns=['bigg.metabolite'])
met_df['compartment'] = met_df['bigg.metabolite'].str.rsplit('_', 1).str[1]
met_df = met_df.join(bigg2kegg_df[['KEGG ID']], on='bigg.metabolite')

#%%
dir_constraints = []
dgr = {}
rxn_res_df = pd.DataFrame(index=S_int.columns)
rxn_res_df['flux'] = rxn_res_df.index.map(flux_dict.get)
rxn_res_df['flux_direction'] = rxn_res_df['flux'].apply(lambda x: np.sign(x) * (np.abs(x) > 1e-4))
rxn_res_df['dG0_prime'] = np.nan
rxn_res_df['dG0_prime_std'] = np.nan

for rxn, coeffs in S_int.iteritems():
    _df = met_df.join(coeffs, on='bigg.metabolite')
    _df = _df.loc[_df[coeffs.name] != 0, ['KEGG ID', coeffs.name]]
    dgr0, dgr0_std = get_dGr0_prime(_df)
    rxn_res_df.at[rxn, 'dG0_prime'] = dgr0
    rxn_res_df.at[rxn, 'dG0_prime_std'] = dgr0_std

rxn_res_df['constraint_dg'] = ((rxn_res_df['flux_direction'] != 0) &
                               (~pd.isnull(rxn_res_df['dG0_prime_std'])) &
                               (rxn_res_df['dG0_prime_std'] < 50))
#%
for rxn, coeffs in S_int.iteritems():
    direction = rxn_res_df.at[rxn, 'flux_direction']
    dgr0 = rxn_res_df.at[rxn, 'dG0_prime']
    dgr0_std = rxn_res_df.at[rxn, 'dG0_prime_std']
    
    if rxn_res_df.at[rxn, 'constraint_dg']:
        # convert the 'coeffs' series to a KEGG reaction format
        dgr[rxn] = dgr0 + settings.RT * direction * pulp.lpSum(ln_conc[met] * coeff
            for met, coeff in coeffs[coeffs != 0].iteritems())
        dir_constraints.append(pulp.LpConstraint(dgr[rxn], pulp.LpConstraintLE,
                                                 name='dGr_%s' % rxn,
                                                 rhs=-settings.eps))

# try to get formation energies from equilibrator-api
#met_df['dG0f_prime'] = met_df['KEGG ID'].apply(kegg_to_dg0f_prime)

conc_df = pd.read_csv(os.path.join(settings.DATA_DIR, 'concentration_ranges.csv'),
                      index_col=0)

met_df = met_df.join(conc_df, on='bigg.metabolite')

#%%
conc_constraints = []
for _, row in met_df.iterrows():
    met = row['bigg.metabolite']
    
    if pd.isnull(row['conc_lb']):
        if row['compartment'] == 'e':
            conc_lb = EXTERNAL_CONC_LB
        elif row['compartment'] == 'c':
            conc_lb = CYTOPLASMIC_CONC_LB
    else:
        conc_lb = row['conc_lb']
    
    if pd.isnull(row['conc_ub']):
        if row['compartment'] == 'e':
            conc_ub = EXTERNAL_CONC_UB
        elif row['compartment'] == 'c':
            conc_ub = CYTOPLASMIC_CONC_UB
    else:
        conc_ub = row['conc_ub']

    ln_lb = np.log(conc_lb)
    ln_ub = np.log(conc_ub)
    if ln_lb > ln_ub:
        raise Exception('lb > ub for metabolite %s' % met)
    if ln_ub - ln_lb < settings.eps:
        ln_lb -= settings.eps/2.0
        ln_ub += settings.eps/2.0
        
    conc_constraints.append(pulp.LpConstraint(ln_conc[met], pulp.LpConstraintGE,
                                              name='conc_%s_lower' % met,
                                              rhs=ln_lb))
    conc_constraints.append(pulp.LpConstraint(ln_conc[met], pulp.LpConstraintLE,
                                              name='conc_%s_upper' % met,
                                              rhs=ln_ub))

#%%
lp_nodir = pulp.LpProblem('NEST_nodir')
for c in conc_constraints:
    lp_nodir.addConstraint(c)
lp_nodir.writeLP(os.path.join(settings.RESULT_DIR, 'nest_nodir.lp'))

lp_nest = pulp.LpProblem('NEST')
for c in conc_constraints + dir_constraints:
    lp_nest.addConstraint(c)
lp_nest.writeLP(os.path.join(settings.RESULT_DIR, 'nest.lp'))

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
met_res_df = []
for met in S_int.index:
    lnc_nodir_lb, lnc_nodir_ub = get_min_and_max(lp_nodir, ln_conc[met], met)
    lnc_nest_lb, lnc_nest_ub = get_min_and_max(lp_nest, ln_conc[met], met)
    met_res_df.append((met, lnc_nodir_lb, lnc_nodir_ub, lnc_nest_lb, lnc_nest_ub))

met_res_df = pd.DataFrame(data=met_res_df,
                          columns=['bigg.metabolite', 'lnc_nodir_lb', 'lnc_nodir_ub',
                                   'lnc_nest_lb', 'lnc_nest_up'])
met_res_df = met_res_df.set_index('bigg.metabolite')
met_res_df.round(2).to_csv(os.path.join(settings.RESULT_DIR, 'nest_met.csv'))

#%% calculate the min and max dGr for each reaction:
rxn_nest_results = []
for rxn in S_int.columns:
    # first calculate the dGr range when directionaliry constraints are not imposed
    # (i.e. as a isolated reaction)
    if rxn_res_df.at[rxn, 'constraint_dg']:
        dgr_nodir_lb, dgr_nodir_ub = get_min_and_max(lp_nodir, dgr[rxn], rxn)
        dgr_nest_lb, dgr_nest_ub = get_min_and_max(lp_nest, dgr[rxn], rxn)
        rxn_nest_results.append((rxn, dgr_nodir_lb, dgr_nodir_ub, dgr_nest_lb, dgr_nest_ub))
        
rxn_nest_results = pd.DataFrame(rxn_nest_results,
                                columns=['bigg.reaction', 'dgr_nodir_lb',
                                         'dgr_nodir_ub', 'dgr_nest_lb', 'dgr_nest_ub'])
rxn_res_df = rxn_nest_results.join(rxn_res_df, on='bigg.reaction')
rxn_res_df.set_index('bigg.reaction')
    
forward = rxn_res_df['dgr_nodir_lb'] < 0
backward = rxn_res_df['dgr_nodir_ub'] > 0
rxn_res_df['allowed_direction'] = forward * 1 + backward * (-1)
rxn_res_df['infeasible'] = (rxn_res_df['allowed_direction'] * rxn_res_df['flux_direction'] == -1)

rxn_res_df.round(2).to_csv(os.path.join(settings.RESULT_DIR, 'nest_rxn.csv'))
