from evodesign.Chain import Chain
from Bio.PDB.Atom import Atom



structure = Chain.load_structure('../8hjc.pdb')
n = Chain.length(structure)
backbone = Chain.backbone_atoms(structure)

test_passed = type(n) == int and \
              n > 0 and \
              type(backbone) == list and \
              sum([ type(x) == Atom for x in backbone ]) == 4*n and \
              sum([ x.get_name() in Chain.BACKBONE_ATOMS for x in backbone ]) == 4*n
print('PASSED' if test_passed else 'FAILED')
