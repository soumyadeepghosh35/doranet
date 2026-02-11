import json
import collections.abc
import dataclasses
import re
import time
import typing
import argparse
import multiprocessing
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdqueries
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.rdmolops import GetFormalCharge
import doranet as dn
from doranet import interfaces, metadata


def clean_SMILES(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Chem.rdmolops.RemoveStereochemistry(mol)
    cpd_smiles = Chem.MolToSmiles(mol)
    return cpd_smiles


bio_rules_path = "JN3604IMT_rules.tsv"
cofactors_path = "all_cofactors.tsv"

bio_rules = pd.read_csv(bio_rules_path, sep="\t")
cofactors = pd.read_csv(cofactors_path, sep="\t")

excluded_cofactors = tuple()

cofactors_dict = dict()
cofactors_set = set()
for idx, x in enumerate(cofactors["SMILES"]):
    cofactors_dict[cofactors["#ID"][idx]] = Chem.MolToSmiles(Chem.MolFromSmiles(x))
    cofactors_set.add(Chem.MolToSmiles(Chem.MolFromSmiles(x)))

cofactors_clean_dict = dict()
cofactors_clean = set()
for idx, x in enumerate(cofactors["SMILES"]):
    if cofactors["#ID"][idx] not in excluded_cofactors:
        cofactors_clean_dict[cofactors["#ID"][idx]] = clean_SMILES(x)
        cofactors_clean.add(clean_SMILES(x))


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class SMILESCalculator(metadata.MolPropertyCalc[str]):
    # Calculate SMILES for molecules and save in network
    smiles_key: collections.abc.Hashable

    @property
    def key(self) -> collections.abc.Hashable:
        return self.smiles_key

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(molecule_keys={self.smiles_key})

    @property
    def resolver(self) -> metadata.MetaDataResolverFunc[str]:
        return metadata.TrivialMetaDataResolverFunc

    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.MolDatBase],
        prev_value: typing.Optional[str] = None,
    ) -> typing.Optional[str]:
        if prev_value is not None:
            return prev_value
        item = data.item
        if data.meta is not None and self.smiles_key in data.meta:
            return None
        if not isinstance(item, interfaces.MolDatRDKit):
            raise NotImplementedError(
                f"""Calculator only implemented for molecule type \
                    MolDatRDKit, not {type(item)}"""
            )
        return item.smiles


@typing.final
@dataclasses.dataclass(frozen=True)
class Reaction_Type_Filter(interfaces.RecipeFilter):
    # used for bio rxns, check reactants
    Allow_multi_reactants: bool

    def __call__(self, recipe: interfaces.RecipeExplicit) -> bool:
        if not self.Allow_multi_reactants:
            reas = set()
            for mol in recipe.reactants:
                if mol.meta is None:
                    raise RuntimeError("No molecule metadata found!")
                SMILES = mol.meta["SMILES"]
                if clean_SMILES(SMILES) not in cofactors_clean:
                    reas.add(clean_SMILES(SMILES))
            if len(reas) != 1:
                return False
        if recipe.operator.meta is None:
            raise RuntimeError("No operator metadata found!")
        reas_type = recipe.operator.meta["Reactants"].split(";")
        for idx, mol in enumerate(recipe.reactants):
            # check if reactant type is correct
            if mol.meta is None:
                raise RuntimeError("No molecule metadata found!")
            SMILES = mol.meta["SMILES"]
            if reas_type[idx] == "Any" and clean_SMILES(SMILES) in cofactors_clean:
                return False
            if (
                reas_type[idx] != "Any"
                and clean_SMILES(SMILES) != cofactors_clean_dict[reas_type[idx]]
            ):
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={"Reactants", "SMARTS"}, molecule_keys={"SMILES"}
        )


@typing.final
@dataclasses.dataclass(frozen=True)
class Product_Filter(metadata.ReactionFilterBase):
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if recipe.operator.meta is None:
            raise RuntimeError("No operator metadata found!")
        pros_type = recipe.operator.meta["Products"].split(";")
        for idx, mol in enumerate(recipe.products):
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"""Filter only implemented for molecule type \
                        MolDatRDKit, not {type(mol.item)}"""
                )
            if (
                pros_type[idx] != "Any"
                and clean_SMILES(mol.item.smiles)
                != cofactors_clean_dict[pros_type[idx]]
            ):
                return False
            if (
                pros_type[idx] == "Any"
                and clean_SMILES(mol.item.smiles) in cofactors_clean
            ):
                return False

            if (
                Descriptors.NumRadicalElectrons(Chem.MolFromSmiles(mol.item.smiles))
                != 0
            ):
                return False
            if "." in mol.item.smiles:
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(operator_keys={"Products", "Name"})


@typing.final
@dataclasses.dataclass(frozen=True)
class Check_balance_filter(metadata.ReactionFilterBase):
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if True:
            charge_diff = 0
            reactants_dict: dict[str, int | float] = dict()
            products_dict: dict[str, int | float] = dict()
            pattern = r"([A-Z][a-z]*)(\d*)"
            for mol in recipe.reactants:
                if not isinstance(mol.item, interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Calculator only implemented for molecule type \
                            MolDatRDKit, not {type(mol.item)}"""
                    )
                charge_diff += GetFormalCharge(mol.item.rdkitmol)
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches = re.findall(pattern, smiles)
                for match in matches:
                    element, count = match
                    count = int(count) if count else 1
                    reactants_dict[element] = reactants_dict.get(element, 0) + count
            for mol in recipe.products:
                if not isinstance(mol.item, interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Calculator only implemented for molecule type \
                            MolDatRDKit, not {type(mol.item)}"""
                    )
                charge_diff -= GetFormalCharge(mol.item.rdkitmol)
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches = re.findall(pattern, smiles)
                for match in matches:
                    element, count = match
                    count = int(count) if count else 1
                    products_dict[element] = products_dict.get(element, 0) + count
            if charge_diff and "H" in reactants_dict:
                reactants_dict["H"] = reactants_dict["H"] - charge_diff
            if reactants_dict != products_dict:
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={
                "reactants_stoi",
                "products_stoi",
                "ring_issue",
                "enthalpy_correction",
                "name",
            }
        )


def predictions(
    starters=False,
    allow_multiple_reactants=False,
    direction="forward",
):
    engine = dn.create_engine()
    network = engine.new_network()

    for key, value in cofactors_dict.items():
        if excluded_cofactors is None or key not in excluded_cofactors:
            network.add_mol(
                engine.mol.rdkit(value),
                meta={"SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(value))},
            )

    my_start_i = -1
    for smiles in starters:
        if my_start_i == -1:
            my_start_i = network.add_mol(
                engine.mol.rdkit(smiles),
                meta={"SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(smiles))},
            )
        else:
            network.add_mol(
                engine.mol.rdkit(smiles),
                meta={"SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(smiles))},
            )

    for idx, x in enumerate(bio_rules["SMARTS"]):
        reas_types = bio_rules["Reactants"][idx].split(";")
        pros_types = bio_rules["Products"][idx].split(";")

        if excluded_cofactors is None or (
            not set(excluded_cofactors) & set(reas_types)
            and not set(excluded_cofactors) & set(pros_types)
        ):
            network.add_op(
                engine.op.rdkit(x, kekulize=False, drop_errors=True),
                meta={
                    "name": bio_rules["Name"][idx],
                    "Reactants": bio_rules["Reactants"][idx],
                    "Products": bio_rules["Products"][idx],
                    "Comments": bio_rules["Comments"][idx],
                    "SMARTS": x,
                    "reactants_stoi": (1,) * len(x.split(">>")[0].split(".")),
                    "products_stoi": (1,) * len(x.split(">>")[1].split(".")),
                    "enthalpy_correction": 0,
                    "Reaction_type": "Enzymatic",
                    "Reaction_direction": direction,
                },
            )

    strat = engine.strat.cartesian(network)

    SMILES_Cal = SMILESCalculator("SMILES")
    Product_check = Product_Filter()
    coreactants_filter = engine.filter.bundle.coreactants(tuple(range(my_start_i)))

    reaction_plan = SMILES_Cal >> Product_check >> Check_balance_filter()

    Type_Filter = Reaction_Type_Filter(allow_multiple_reactants)

    strat.expand(
        num_iter=1,
        reaction_plan=reaction_plan,
        bundle_filter=coreactants_filter,
        recipe_filter=Type_Filter,
        save_unreactive=False,
    )

    # get reactions from network
    total_sets_of_products_smiles = []
    for rxn in network.rxns:
        products = rxn.products
        products_smiles = [network.mols[i].uid for i in products]
        total_sets_of_products_smiles.append(products_smiles)

    return_list = []
    for sets in total_sets_of_products_smiles:
        sets = list(set(sets))
        sets.sort()
        if sets not in return_list:
            return_list.append(sets)
    return return_list


def process_rhea_rxn(rxn):
    reas = rxn.split(">>")[0].split(".")
    pros = rxn.split(">>")[1].split(".")
    predicted_pros_list = predictions(starters=reas)
    for pro_set in predicted_pros_list:
        if set(pros).issubset(set(pro_set)):
            return True
    return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--procs",
        "-p",
        type=int,
        default=10,
        help="number of parallel processes to use (default: 1)",
    )
    args = parser.parse_args()
    num_procs = args.procs

    with open("Rhea_random_1000.json", encoding="utf-8") as f:
        rhea_rxns = json.load(f)
    start_time = time.time()

    if num_procs == 1:
        results = [process_rhea_rxn(rxn) for rxn in rhea_rxns]
    else:
        with multiprocessing.Pool(processes=num_procs) as pool:
            results = pool.map(process_rhea_rxn, rhea_rxns)

    checked_reactions = len(rhea_rxns)
    correct_num = sum(1 for r in results if r)
    total_calls = checked_reactions

    correction_rate = correct_num / checked_reactions if checked_reactions else 0
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print("success count:", correct_num)
    print("rhea reactions reproduction rate:", correction_rate)
    print("time used:", "{:.2f}".format(elapsed_time), " minutes")
    print("Total doranet calls:", total_calls)
