import doranet.modules.enzymatic as enzymatic
import os
from equilibrator_api import ComponentContribution, Reaction
from equilibrator_assets.local_compound_cache import LocalCompoundCache


user_starters = [
    "CC(N)C(=O)O",
]


job_name = os.path.basename(__file__).removesuffix(".py")


cc = ComponentContribution()
lc = LocalCompoundCache()
lc.generate_local_cache_from_default_zenodo("compounds.sqlite")  ###
lc.load_cache("compounds.sqlite")


def rxn_dG(my_rxn):
    reas_pros = my_rxn["reactants"] + my_rxn["products"]
    cpd_results = lc.get_compounds(reas_pros)
    compound_dict = dict()
    for idx, i in enumerate(reas_pros):
        if cpd_results[idx].compound is None:
            return None
        if idx < len(my_rxn["reactants"]):
            compound_dict[cpd_results[idx].compound] = -1
        else:
            compound_dict[cpd_results[idx].compound] = 1
    reaction = Reaction(compound_dict)
    dG0_prime = cc.standard_dg_prime(reaction)
    return dG0_prime.value.magnitude / 4.184


forward_network = enzymatic.generate_network(
    job_name=job_name,
    starters=user_starters,
    gen=2,
    direction="forward",
    rxn_thermo_calculator=rxn_dG,
    max_rxn_thermo_change=0,
    allow_multiple_reactants=True,
)
