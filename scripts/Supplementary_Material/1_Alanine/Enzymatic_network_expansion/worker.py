#!/usr/bin/env python3
import json
import sys
from equilibrator_api import ComponentContribution, Reaction
from equilibrator_assets.local_compound_cache import LocalCompoundCache

def main():
    # Input: JSON on stdin with keys: reactants, products, sqlite_path
    payload = json.load(sys.stdin)
    reactants = payload["reactants"]
    products = payload["products"]
    sqlite_path = payload["sqlite_path"]

    cc = ComponentContribution()
    lc = LocalCompoundCache()
    lc.load_cache(sqlite_path)

    reas_pros = reactants + products
    cpd_results = lc.get_compounds(reas_pros)

    compound_dict = {}
    for idx, _ in enumerate(reas_pros):
        if cpd_results[idx].compound is None:
            print("None")
            return

        coeff = -1 if idx < len(reactants) else 1
        compound_dict[cpd_results[idx].compound] = coeff

    reaction = Reaction(compound_dict)
    dG0_prime = cc.standard_dg_prime(reaction)

    # equilibrator returns kJ/mol; convert to kcal/mol
    kcal_per_mol = dG0_prime.value.magnitude / 4.184
    print(kcal_per_mol)

if __name__ == "__main__":
    main()

