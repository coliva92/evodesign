import os

from evodesign.Evaluation.Misc.RCC import RCC
from evodesign.Protein.ChainFactory import ChainFactory
from evodesign.Protein.Geometry import compute_contacts_map

home = os.environ["HOME"]
example_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
target_pdb_path = os.path.abspath(os.path.join(example_folder, "1y32.pdb"))
output_folder = os.path.join(example_folder, "outputs")
rcc_calc = RCC(
    path_to_quick_cliques=os.path.join(home, "quick-cliques", "bin", "qc"),
    distance_threshold=6.0,
    min_separation=3,
)
graph_filename = os.path.join(output_folder, "1y32_graph.txt")
chain = ChainFactory.create_from_pdb(target_pdb_path)
contact_map = compute_contacts_map(
    chain.ca_atoms,
    rcc_calc.distance_threshold,
    rcc_calc.min_separation,
)
rcc = rcc_calc.compute_rcc_vector(contact_map, graph_filename)
print(rcc)  # [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
