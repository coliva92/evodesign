from argparse import ArgumentParser
import evodesign.Settings as Settings
from evodesign.Prediction.DirectoryManager import DirectoryManager
from evodesign.Utils.ChainFactory import ChainFactory
from evodesign.Metrics.TMScore import TMScore
import pandas as pd
import seaborn as sns
import os


def compute_sctm_scores(
    csv_path: str,
    predictor_json_path: str,
    output_dir: str,
    sequence_csv_col: str = "Sequence",
    ref_sequence_csv_col: str = "ReferenceSequence",
) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    predictor = Settings.load(predictor_json_path)
    manager = DirectoryManager(output_dir, output_dir, output_dir)
    calc = TMScore()
    for idx, row in df.iterrows():
        seq = row[sequence_csv_col]
        ref_seq = row[ref_sequence_csv_col]
        predictor.predict_single_pdb_file(seq, f"seq_{idx}", manager)
        pdb_path = os.path.join(
            manager.prediction_pdbs_dir, f"{manager.prefix}_seq_{idx}.pdb"
        )
        model_chain = ChainFactory.create_from_pdb(pdb_path)
        predictor.predict_single_pdb_file(ref_seq, f"ref_{idx}", manager)
        pdb_path = os.path.join(
            manager.prediction_pdbs_dir, f"{manager.prefix}_ref_{idx}.pdb"
        )
        ref_chain = ChainFactory.create_from_pdb(pdb_path)
        tm_score = calc.do(model_chain.backbone_atoms, ref_chain.backbone_atoms)
        df.loc[idx, "scTMScore"] = tm_score
    return df


if __name__ == "__main__":
    parser = ArgumentParser("sctm_from_csv")
    parser.add_argument("csv_path", type=str)
    parser.add_argument("predictor_jsons_dir", type=str)
    parser.add_argument("output_dir", type=str)
    parser.add_argument("--sequence_csv_col", "-s", type=str, default="Sequence")
    parser.add_argument(
        "--ref_sequence_csv_col", "-r", type=str, default="ReferenceSequence"
    )
    parser.add_argument("--color_palette_name", "-c", type=str, default="colorblind")
    parser.add_argument("--context", "-k", type=str, default="paper")
    parser.add_argument("--height", "-h", type=float, default=2)
    parser.add_argument("--aspect", "-a", type=float, default=1)
    parser.add_argument("--font_scale", "-f", type=float, default=1)
    parser.add_argument("--x_key", "--X", type=str, default="Predictor")
    parser.add_argument("--y_key", "--Y", type=str, default="scTMScore")
    parser.add_argument("--hue_key", "--Z", type=str, default="Designer")
    args = parser.parse_args()
    dfs = [
        compute_sctm_scores(
            args.csv_path,
            os.path.join(args.predictor_jsons_dir, filename),
            args.output_dir,
            args.sequence_csv_col,
            args.ref_sequence_csv_col,
        )
        for filename in os.listdir(args.predictor_jsons_dir)
        if filename[-5:] == ".json"
    ]
    out_df = pd.concat(dfs, ignore_index=True)
    out_df.to_csv(os.path.join(args.output_dir, "sctm_scores.csv"), index=False)
    g = sns.displot(
        data=out_df,
        x=args.x_key,
        y=args.y_key,
        hue=args.hue_key,
        height=args.height,
        aspect=args.aspect,
    )
    g.savefig(os.path.join(args.output_dir, "sctm_scores.svg", format="svg"))
