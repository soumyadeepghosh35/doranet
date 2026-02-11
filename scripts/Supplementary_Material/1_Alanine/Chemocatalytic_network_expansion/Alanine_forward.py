import doranet.modules.synthetic as synthetic
import os
from pathermo.properties import Hf


user_starters = [
    "CC(N)C(=O)O",
]


user_helpers = (
    "O",
    "Cl",
    "ClCl",
    "O=O",
    "[C-]#[O+]",
    "[H][H]",
)

job_name = os.path.basename(__file__).removesuffix(".py")


forward_network = synthetic.generate_network(
    job_name=job_name,
    starters=user_starters,
    helpers=user_helpers,
    gen=3,
    direction="forward",
    molecule_thermo_calculator=Hf,
    max_rxn_thermo_change=40,
)
