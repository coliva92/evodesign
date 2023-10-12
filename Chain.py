from typing import List
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom





BACKBONE_ATOM_NAMES = { 'N', 'CA', 'C', 'O' }



def count_chain_residues(structure: Structure, 
                         modelId: int = 0, 
                         chainId: str = 'A'
                         ) -> int:
  return len(list(structure[modelId][chainId].get_residues()))



def filter_backbone_atoms_in_chain(structure: Structure,
                                   modelId: int = 0, 
                                   chainId: str = 'A'
                                   ) -> List[Atom]:
  return list(filter(lambda atom: atom.get_name() in BACKBONE_ATOM_NAMES, 
                     structure[modelId][chainId].get_atoms()))



def filter_alpha_carbons_in_backbone(backbone: List[Atom]) -> List[Atom]:
  return filter(lambda atom: atom.get_name() == 'CA', backbone)
