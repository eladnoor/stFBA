from cobra.flux_analysis.loopless import construct_loopless_model
from cobra.io.json import load_json_model
from cobra.io.sbml import create_cobra_model_from_sbml_file

#model_name = 'model/e_coli_core'
#model_name = 'model/iAF1260'
model_name = 'model/iJO1366'

model_json = load_json_model(model_name + '.json')
fba_sol_json = model_json.optimize(solver='cplex')

#model_sbml = create_cobra_model_from_sbml_file(model_name + '.xml.gz')
#model_sbml.objective = model_sbml.reactions.BIOMASS_Ec_iJO1366_WT_53p95M
#model_sbml.reactions.EX_glc__D_e.lower_bound = -10
#model_json.reactions.ATPM.lower_bound = 0

#%%
#fba_sol_sbml = model_sbml.optimize(solver='cplex')
#print "FBA SBML", fba_sol_sbml.status, fba_sol_sbml.f

print "FBA JSON", fba_sol_json.status, fba_sol_json.f

#loopless_model_sbml = construct_loopless_model(model_sbml)
#ll_sol_sbml = loopless_model_sbml.optimize(solver='cplex')
#print "ll-FBA SBML", ll_sol_sbml.status, ll_sol_sbml.f
#
loopless_model_json = construct_loopless_model(model_json)
ll_sol_json = loopless_model_json.optimize(solver='cplex')
print "ll-FBA JSON", ll_sol_json.status, ll_sol_json.f


#%%
print "SBML solution"
for r, v in fba_sol_sbml.x_dict.iteritems():
    if v != 0 and r[0:3] == 'EX_':
        print "%20s : %.1f" % (r, v)
print "JSON solution"
for r, v in fba_sol_json.x_dict.iteritems():
    if v != 0 and r[0:3] == 'EX_':
        print "%20s : %.1f" % (r, v)
