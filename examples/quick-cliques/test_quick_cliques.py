import os

from evodesign.Chemistry.ChainFactory import ChainFactory
from evodesign.Chemistry.Geometry import compute_contacts_map
from evodesign.Metrics.RCC import RCC

home = os.environ["HOME"]
rcc_calc = RCC(
    path_to_quick_cliques=f"{home}/quick-cliques/bin/qc",
    distance_threshold=6.0,
    min_separation=3,
)
molecule = ChainFactory.create_from_pdb("../1y32.pdb")
contact_map = compute_contacts_map(
    molecule.ca_atoms,
    rcc_calc.distance_threshold,
    rcc_calc.min_separation,
)
rcc = rcc_calc.compute_rcc_vector(contact_map)
print(rcc)
