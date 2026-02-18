from argparse import ArgumentParser
from evodesign.Metrics.RMSD import RMSD
from evodesign.Metrics.GDT import GDT
from evodesign.Metrics.TMScore import TMScore
from evodesign.Metrics.lDDT import lDDT
from evodesign.Chemistry.ChainFactory import ChainFactory


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("model_pdb_path", type=str)
    parser.add_argument("ref_pdb_path", type=str)
    args = parser.parse_args()
    rmsd_calc = RMSD()
    gdt_calc = GDT()
    tm_calc = TMScore()
    lddt_calc = lDDT()
    model = ChainFactory.create_from_pdb(args.model_pdb_path)
    reference = ChainFactory.create_from_pdb(args.ref_pdb_path)
    rmsd, _ = rmsd_calc.do(model.backbone_atoms, reference.backbone_atoms)
    gdt = gdt_calc.do(model.backbone_atoms, reference.backbone_atoms)
    tm_score = tm_calc.do(model.backbone_atoms, reference.backbone_atoms)
    lddt = lddt_calc.do(model.backbone_atoms, reference.backbone_atoms)
    print(f"RMSD: {rmsd}\nGDT: {gdt}\nTMScore: {tm_score}\nlDDT: {lddt}\n")
