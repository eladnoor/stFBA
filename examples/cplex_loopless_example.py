from cobra.flux_analysis.loopless import construct_loopless_model
from cobra.io.json import load_json_model

#model_name = 'model/e_coli_core'
#model_name = 'model/iAF1260'
model_name = 'model/iJO1366'

model_json = load_json_model(model_name + '.json')
fba_sol_json = model_json.optimize(solver='cglpk')

print("FBA JSON: %s, %f" % (fba_sol_json.status, fba_sol_json.f))

loopless_model_json = construct_loopless_model(model_json)
ll_sol_json = loopless_model_json.optimize(solver='cglpk')
print("ll-FBA JSON: %s, %f" % (ll_sol_json.status, ll_sol_json.f))
