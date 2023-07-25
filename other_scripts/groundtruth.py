from argparse import ArgumentParser
import evodesign.Chain as Chain
from evodesign.Metrics import Rmsd, Gdt





parser = ArgumentParser()
parser.add_argument('refPdbFilename')
parser.add_argument('modelPdbFilename')
args = parser.parse_args()
reference = Chain.load_structure_from_pdb_file(args.refPdbFilename)
reference = Chain.filter_backbone_atoms_from_chain(reference)
model = Chain.load_structure_from_pdb_file(args.modelPdbFilename)
model = Chain.filter_backbone_atoms_from_chain(model)
rmsd_calculator, gdt_calculator = Rmsd(), Gdt()
with open('groundtruth.txt', 'wt', encoding='utf-8') as output:
  rmsd = rmsd_calculator.compute(model, reference)
  gdt = gdt_calculator.compute(model, reference)
  output.write(f'model vs. ref\nrmsd: {rmsd}, gdt: {gdt}\n\n')
  rmsd = rmsd_calculator.compute(reference, reference)
  gdt = gdt_calculator.compute(reference, reference)
  output.write(f'ref vs. ref\nrmsd: {rmsd}, gdt: {gdt}\n')
