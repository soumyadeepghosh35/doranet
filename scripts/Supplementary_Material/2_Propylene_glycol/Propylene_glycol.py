import doranet.modules.synthetic as synthetic
import doranet.modules.post_processing as post_processing
from pathermo.properties import Hf


forward_helpers = {
    "O=O",
    "O",
    "S",
    "N",
    "ClCl",
    "[H][H]",
    "Cl",
}
forward_starters = {
    "CCCCCCCCC=CCCCCCCCC(=O)O",
    "CC(=O)O",
    "CCO",
    "CCCCCCCCCCCCCCCC(=O)O",
    "CSCCC(N)C(=O)O",
    "CCC",
    "CO",
    "C=C(C)C",
    "CCC(=O)O",
    "CC(=O)NC1C(O)OC(CO)C(O)C1O",
    "Cc1ccccc1C",
    "OCC1OC(O)C(O)C1O",
    "NC(N)=O",
    "NCCCCC(N)C(=O)O",
    "OC1COC(O)C(O)C1O",
    "C=CC",
    "OCCO",
    "COc1cc(C=CCO)ccc1O",
    "O=C(O)c1ccccc1",
    "C=CC=C",
    "CCCC(=O)O",
    "C=Cc1ccccc1",
    "O=C=O",
    "CCCCCC=CCC=CCCCCCCCC(=O)O",
    "O=CC(O)C(O)C(O)C(O)CO",
    "O=CC(O)C(O)C(O)CO",
    "O=[N+]([O-])O",
    "O=[N+][O-]",
    "C",
    "O=C(O)c1ccc(C(=O)O)cc1",
    "Cc1ccccc1",
    "C=CCC",
    "[N]=O",
    "OCC(O)C(O)C(O)C(O)CO",
    "CCCCCCCCCCCCCCCCCC(=O)O",
    "c1ccccc1",
    "[C-]#[O+]",
    "Cc1cccc(C)c1",
    "O=C(Cl)Cl",
    "Cc1ccc(C)cc1",
    "C=C",
    "OCC1OC(O)C(O)C(O)C1O",
    "CCCCC(=O)O",
    "CC=CC",
    "CCCCCC(=O)O",
    "O=CO",
    "CC(O)C(=O)O",
}

job_name = "Propylene_glycol"

retro_helpers = {
    "O",
    "O=C=O",
    "Cl",
    "[Cl][Cl]",
    "O=O",
    "N#N",
    "[C-]#[O+]",
    "N",
    "S",
    "CN(C)C",
}
user_target = {"CC(O)CO"}


if __name__ == "__main__":
    forward_network = synthetic.generate_network(
        job_name=job_name,
        starters=forward_starters,
        helpers=forward_helpers,
        gen=1,
        direction="forward",
        molecule_thermo_calculator=Hf,
        max_rxn_thermo_change=15,
        max_atoms={"C": 18},
    )

    retro_network = synthetic.generate_network(
        job_name=job_name,
        starters=user_target,
        helpers=retro_helpers,
        gen=3,
        direction="retro",
        molecule_thermo_calculator=Hf,
        max_rxn_thermo_change=15,
        max_atoms={"C": 10},
    )

    # Perform post-processing in a single step
    post_processing.one_step(
        networks={forward_network, retro_network},
        molecule_thermo_calculator=Hf,
        max_rxn_thermo_change=15,
        total_generations=4,
        search_depth=4,
        max_num_rxns=5,
        min_rxn_atom_economy=0.6,
        num_process=10,
        starters=forward_starters,
        helpers=forward_helpers,
        target=user_target,
        job_name=job_name,
        consider_name_difference=False,
    )

    # Alternatively, perform post-processing in multiple steps for more control (e.g., when using Reaxys results)
    # post_processing.pretreat_networks(
    #     networks={forward_network, retro_network},
    #     total_generations=4,
    #     starters=forward_starters,
    #     helpers=forward_helpers,
    #     job_name=job_name,
    #     molecule_thermo_calculator=Hf,
    # )

    # post_processing.pathway_finder(
    #     starters=forward_starters,
    #     helpers=forward_helpers,
    #     target=user_target,
    #     search_depth=4,
    #     max_num_rxns=5,
    #     min_rxn_atom_economy=0.6,
    #     job_name=job_name,
    #     consider_name_difference=False,
    # )

    # post_processing.pathway_ranking(
    #     starters=forward_starters,
    #     helpers=forward_helpers,
    #     target=user_target,
    #     num_process=10,
    #     job_name=job_name,
    #     molecule_thermo_calculator=Hf,
    #     max_rxn_thermo_change=15,
    # )

    # post_processing.pathway_visualization(
    #     starters=forward_starters,
    #     helpers=forward_helpers,
    #     num_process=10,
    #     job_name=job_name,
    # )
