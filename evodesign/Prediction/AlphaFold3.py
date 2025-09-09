from .Predictor import Predictor
from ..Utils.Subprocess import run_subprocess
from ..Utils.mmCIF import convert_cif_to_pdb
from Bio.PDB import MMCIFParser, PDBIO
import os
import shutil
import json
from typing import List, Optional


class AlphaFold3(Predictor):

    _parser = MMCIFParser()
    _io = PDBIO()

    def __init__(
        self,
        path_to_run_alphafold_py: str,
        model_dir: str,
        output_dir: str,
        num_recycles: int = 1,
        num_diffusion_samples: int = 1,
        model_seeds: Optional[List[int]] = None,
        version: int = 2,
        run_data_pipeline: bool = False,
        max_template_date: str = "2021-09-30",
        gpu_device: int = 0,
        jackhmmer_binary_path: Optional[str] = None,
        nhmmer_binary_path: Optional[str] = None,
        hmmalign_binary_path: Optional[str] = None,
        hmmsearch_binary_path: Optional[str] = None,
        hmmbuild_binary_path: Optional[str] = None,
        db_dir: Optional[str] = None,
        small_bfd_database_path: Optional[str] = None,
        mgnify_database_path: Optional[str] = None,
        uniprot_cluster_annot_database_path: Optional[str] = None,
        uniref90_database_path: Optional[str] = None,
        ntrna_database_path: Optional[str] = None,
        rfam_database_path: Optional[str] = None,
        rna_central_database_path: Optional[str] = None,
        pdb_database_path: Optional[str] = None,
        seqres_database_path: Optional[str] = None,
        jax_compilation_cache_dir: Optional[str] = None,
    ):
        super().__init__()
        self.path_to_run_alphafold_py = os.path.abspath(path_to_run_alphafold_py)
        self.model_dir = os.path.abspath(model_dir)
        self.output_dir = os.path.abspath(output_dir)
        self.num_recycles = num_recycles
        self.num_diffusion_samples = num_diffusion_samples
        if model_seeds is None:
            import time

            model_seeds = [int(time.time())]
        self.model_seeds = model_seeds
        self.version = version
        self.run_data_pipeline = run_data_pipeline
        self.max_template_date = max_template_date
        self.gpu_device = gpu_device
        self.jackhmmer_binary_path = jackhmmer_binary_path
        self.nhmmer_binary_path = nhmmer_binary_path
        self.hmmalign_binary_path = hmmalign_binary_path
        self.hmmsearch_binary_path = hmmsearch_binary_path
        self.hmmbuild_binary_path = hmmbuild_binary_path
        self.db_dir = db_dir
        self.small_bfd_database_path = small_bfd_database_path
        self.mgnify_database_path = mgnify_database_path
        self.uniprot_cluster_annot_database_path = uniprot_cluster_annot_database_path
        self.uniref90_database_path = uniref90_database_path
        self.ntrna_database_path = ntrna_database_path
        self.rfam_database_path = rfam_database_path
        self.rna_central_database_path = rna_central_database_path
        self.pdb_database_path = pdb_database_path
        self.seqres_database_path = seqres_database_path
        self.jax_compilation_cache_dir = jax_compilation_cache_dir

    def _get_cmd_array(
        self, input_path: str, use_batch_inference: bool = False
    ) -> List[str]:
        cmd = [
            "python3",
            self.path_to_run_alphafold_py,
            (
                f"--input_dir={input_path}"
                if use_batch_inference
                else f"--json_path={input_path}"
            ),
            f"--model_dir={self.model_dir}",
            f"--output_dir={self.output_dir}",
            f"--run_data_pipeline={self.run_data_pipeline}",
            f"--num_recycles={self.num_recycles}",
            f"--num_diffusion_samples={self.num_diffusion_samples}",
            f"--max_template_date={self.max_template_date}",
            f"--gpu_device={self.gpu_device}",
        ]
        if self.jackhmmer_binary_path is not None:
            cmd.append(f"--jackhmmer_binary_path={self.jackhmmer_binary_path}")
        if self.nhmmer_binary_path is not None:
            cmd.append(f"--nhmmer_binary_path={self.nhmmer_binary_path}")
        if self.hmmalign_binary_path is not None:
            cmd.append(f"--hmmalign_binary_path={self.hmmalign_binary_path}")
        if self.hmmsearch_binary_path is not None:
            cmd.append(f"--hmmsearch_binary_path={self.hmmsearch_binary_path}")
        if self.hmmbuild_binary_path is not None:
            cmd.append(f"--hmmbuild_binary_path={self.hmmbuild_binary_path}")
        if self.db_dir is not None:
            cmd.append(f"--db_dir={self.db_dir}")
        if self.small_bfd_database_path is not None:
            cmd.append(f"--small_bfd_database_path={self.small_bfd_database_path}")
        if self.mgnify_database_path is not None:
            cmd.append(f"--mgnify_database_path={self.mgnify_database_path}")
        if self.uniprot_cluster_annot_database_path is not None:
            cmd.append(f"--uniprot_cluster_annot_database_path={self.uniprot_cluster_annot_database_path}")
        if self.uniref90_database_path is not None:
            cmd.append(f"--uniref90_database_path={self.uniref90_database_path}")
        if self.ntrna_database_path is not None:
            cmd.append(f"--ntrna_database_path={self.ntrna_database_path}")
        if self.rfam_database_path is not None:
            cmd.append(f"--rfam_database_path={self.rfam_database_path}")
        if self.rna_central_database_path is not None:
            cmd.append(f"--rna_central_database_path={self.rna_central_database_path}")
        if self.pdb_database_path is not None:
            cmd.append(f"--pdb_database_path={self.pdb_database_path}")
        if self.seqres_database_path is not None:
            cmd.append(f"--seqres_database_path={self.seqres_database_path}")
        if self.jax_compilation_cache_dir is not None:
            cmd.append(f"--jax_compilation_cache_dir={self.jax_compilation_cache_dir}")
        cmd.append("--force_output_dir=True")
        return cmd

    def create_input_json(
        self, sequence: str, prediction_jobname: str, structure_id: str = "A"
    ) -> dict:
        # Assume we're only predicting monomeric structures
        return {
            "name": prediction_jobname,
            "modelSeeds": self.model_seeds,
            "sequences": [
                {
                    "protein": {
                        "id": structure_id,
                        "sequence": sequence,
                        "unpairedMsa": "",
                        "pairedMsa": "",
                        "templates": [],
                    }
                }
            ],
            "dialect": "alphafold3",
            "version": self.version,
        }

    def run_alphafold(self, input_path: str, use_batch_inference: bool = False):
        run_subprocess(self._get_cmd_array(input_path, use_batch_inference))

    def predict_single_pdb_str(self, sequence: str) -> str:
        pdb_path = "tmp_prediction.pdb"
        self.predict_single_pdb_file(sequence, pdb_path)
        with open(pdb_path, "rt", encoding="utf-8") as pdb_file:
            prediction = pdb_file.read()
        os.remove(pdb_path)
        return prediction

    def predict_single_pdb_file(self, sequence: str, pdb_path: str) -> None:
        os.makedirs(self.output_dir, exist_ok=True)
        prediction_jobname = os.path.splitext(os.path.basename(pdb_path))[0]
        input_json = self.create_input_json(sequence, prediction_jobname)
        json_path = os.path.join(self.output_dir, f"{prediction_jobname}.json")
        with open(json_path, "wt", encoding="utf-8") as json_file:
            json.dump(input_json, json_file)
        self.run_alphafold(json_path)
        prediction_dir = os.path.join(self.output_dir, prediction_jobname)
        prediction_cif = os.path.join(prediction_dir, f"{prediction_jobname}_model.cif")
        prediction_pdb_path = convert_cif_to_pdb(
            prediction_cif, parser=self._parser, io=self._io
        )
        shutil.copyfile(prediction_pdb_path, pdb_path)
        os.remove(json_path)
        self.delete_folder(prediction_dir)

    def do(
        self, sequences: List[str], pdbs_dir: str, pdb_prefix: str = "tmp_prediction"
    ) -> None:
        os.makedirs(pdbs_dir, exist_ok=True)
        for i, sequence in enumerate(sequences):
            prediction_jobname = f"{pdb_prefix}_{i}"
            input_json = self.create_input_json(sequence, prediction_jobname)
            json_path = os.path.join(pdbs_dir, f"{prediction_jobname}.json")
            with open(json_path, "wt", encoding="utf-8") as json_file:
                json.dump(input_json, json_file)
        self.run_alphafold(pdbs_dir, use_batch_inference=True)
        for i in range(len(sequences)):
            prediction_jobname = f"{pdb_prefix}_{i}"
            prediction_dir = os.path.join(self.output_dir, prediction_jobname)
            prediction_cif = os.path.join(
                prediction_dir, f"{prediction_jobname}_model.cif"
            )
            prediction_pdb = convert_cif_to_pdb(
                prediction_cif, parser=self._parser, io=self._io
            )
            pdb_path = os.path.join(self.output_dir, os.path.basename(prediction_pdb))
            shutil.copyfile(prediction_pdb, pdb_path)
            self.delete_folder(prediction_dir)
        self.delete_folder(pdbs_dir)
