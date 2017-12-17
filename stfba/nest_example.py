import pandas as pd
import pulp
import os
import settings
import numpy as np
from cobra.flux_analysis.parsimonious import pfba
from cobra.manipulation.modify import convert_to_irreversible
from cobra.io import read_sbml_model, load_json_model

#cobra_model = load_json_model(settings.IJO1366_JSON_FNAME)
#BM_RXN = 'BIOMASS_Ec_iJO1366_core_53p95M'

cobra_model = read_sbml_model(settings.CORE_SBML_FNAME)
BM_RXN = 'BIOMASS_Ecoli_core_w_GAM'

S = settings.get_stoichiometry_from_model(cobra_model)

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

#%%
# construct the LP for NEST (network-embedded semi-thermodynamic analysis)
pulp_solver = pulp.GLPK_CMD(msg=0, options=['--xcheck'])
lp = pulp.LpProblem('NEST')
dgf = pulp.LpVariable.dicts('dGf', indexs=S.index,
                            lowBound=-settings.M,
                            upBound=settings.M,
                            cat=pulp.LpContinuous)


# keep only reactions that carry flux and are internal
active_rxn = [x for x in x_dict.keys() if x[0:3] != 'EX_']
active_rxn.remove(BM_RXN)
active_rxn = [x for x in active_rxn if np.abs(x_dict[x]) > settings.eps]
S_int = S[active_rxn]

# keep only metabolites that are not involved and active internal reaction
S_int = S_int.loc[np.abs(S_int).sum(1) > 0, :]

#%%
dgr = {}
for rxn, coeffs in S_int.iteritems():
    dgr[rxn] = \
        pulp.lpSum(dgf[met] * coeff
                   for met, coeff in coeffs[coeffs != 0].iteritems())
    
    if x_dict[rxn] > settings.eps:
        const = pulp.LpConstraint(dgr[rxn], pulp.LpConstraintLE,
                                  name='dGr_%s' % rxn,
                                  rhs=-settings.eps)
    elif x_dict[rxn] < -settings.eps:
        const = pulp.LpConstraint(dgr[rxn], pulp.LpConstraintGE,
                                  name='dGr_%s' % rxn,
                                  rhs=settings.eps)
    else:
        continue
    
    lp.addConstraint(const)

#%%
conc_df = pd.read_csv(os.path.join(settings.DATA_DIR, 'concentration_ranges.csv'),
                      index_col=0)

for met, row in conc_df.iterrows():
    if met in cobra_model.metabolites:
        ln_lb = np.log(row['lower_bound_mM']*1e-3) * settings.RT + row['dG0f'] - settings.eps
        ln_ub = np.log(row['upper_bound_mM']*1e-3) * settings.RT + row['dG0f'] + settings.eps
        
        if ln_ub - ln_lb < 1e-5:
            raise Exception('margin too small for %s' % met)
        
        lp.addConstraint(pulp.LpConstraint(dgf[met], pulp.LpConstraintGE,
                                           name='conc_%s_lower' % met,
                                           rhs=ln_lb))
        lp.addConstraint(pulp.LpConstraint(dgf[met], pulp.LpConstraintLE,
                                           name='conc_%s_upper' % met,
                                           rhs=ln_ub))

lp.writeLP(os.path.join(settings.RESULT_DIR, 'nest.lp'))

#%% calculate the min and max dGf for each metabolite:
results = []
for met in S_int.index:
    lp.setObjective(dgf[met])
    lp.sense = pulp.LpMaximize
    lp.solve(pulp_solver)
    if lp.status != pulp.LpStatusOptimal:
        print('%s, LP status %s' % (met, lp.status))
        continue
    
    dgf_ub = pulp.value(dgf[met])
    
    lp.sense = pulp.LpMinimize
    lp.solve(pulp_solver)
    if lp.status != pulp.LpStatusOptimal:
        print('%s, LP status %s' % (met, lp.status))
        continue
    dgf_lb = pulp.value(dgf[met])
    results.append((met, 'metabolite', dgf_lb, dgf_ub))
    print('%6.2f < %20s < %6.2f' % (dgf_lb, 'dGf(' + met + ')', dgf_ub))
    
#%% calculate the min and max dGr for each reaction:
for rxn in S_int.columns:
    if np.abs(x_dict[rxn]) > settings.eps:
        lp.setObjective(dgr[rxn])
        lp.sense = pulp.LpMaximize
        lp.solve(pulp_solver)
        if lp.status != pulp.LpStatusOptimal:
            print('%s, LP status %s' % (rxn, lp.status))
            continue
        
        dgr_ub = pulp.value(dgr[rxn])
        
        lp.sense = pulp.LpMinimize
        lp.solve(pulp_solver)
        if lp.status != pulp.LpStatusOptimal:
            print('%s, LP status %s' % (rxn, lp.status))
            continue
        dgr_lb = pulp.value(dgr[rxn])

        results.append((rxn, 'reaction', dgr_lb, dgr_ub))
        print('%6.2f < %20s < %6.2f' % (dgr_lb, 'dGr(' + rxn + ')', dgr_ub))
        
results_df = pd.DataFrame(data=results, columns=['bigg', 'type', 'lb', 'up'])
results_df.to_csv(os.path.join(settings.RESULT_DIR, 'nest.csv'))