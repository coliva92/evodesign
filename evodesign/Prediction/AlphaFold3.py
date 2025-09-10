from .AlphaFoldInterface import AlphaFoldInterface
from .DirectoryManager import DirectoryManager
from ..Utils.mmCIF import convert_cif_to_pdb
from Bio.PDB import MMCIFParser, PDBIO
import os
import shutil
import json
from typing import List, Optional


class AlphaFold3(AlphaFoldInterface):

    _parser = MMCIFParser()
    _io = PDBIO()

    def __init__(
        self,
        path_to_run_alphafold_py: str,
        model_dir: str,
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

    def _create_model_input(
        self,
        sequence: str,
        protein_name: str,
        input_dir: str,
        output_dir: str,
    ) -> str:
        input_json = self.create_input_json(sequence, protein_name)
        json_path = os.path.join(input_dir, f"{protein_name}.json")
        with open(json_path, "wt", encoding="utf-8") as json_file:
            json.dump(input_json, json_file)
        return json_path

    def _prediction_pdb_path(
        self,
        protein_name: str,
        output_dir: str,
    ) -> str:
        prediction_dir = os.path.join(output_dir, protein_name)
        prediction_cif = os.path.join(prediction_dir, f"{protein_name}_model.cif")
        return convert_cif_to_pdb(prediction_cif, parser=self._parser, io=self._io)

    def _create_cmd_array(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ) -> List[str]:
        cmd = [
            "python3",
            self.path_to_run_alphafold_py,
            (
                f"--input_dir={input_path}"
                if do_batch_inference
                else f"--json_path={input_path}"
            ),
            f"--model_dir={self.model_dir}",
            f"--output_dir={output_dir}",
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
            cmd.append(
                f"--uniprot_cluster_annot_database_path={self.uniprot_cluster_annot_database_path}"
            )
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
        self,
        sequence: str,
        protein_name: str,
        structure_id: str = "A",
    ) -> dict:
        # Assume we're only predicting monomeric structures
        return {
            "name": protein_name,
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

    def do(
        self,
        sequences: List[str],
        directory: DirectoryManager,
    ) -> None:
        os.makedirs(directory.prediction_pdbs_dir, exist_ok=True)
        os.makedirs(directory.model_input_dir, exist_ok=True)
        directory.empty_folders_content()
        for i, sequence in enumerate(sequences):
            protein_name = f"{directory.prefix}_{i}"
            json_path = self._create_model_input(
                sequence,
                protein_name,
                directory.model_input_dir,
                directory.model_output_dir,
            )
        self.run_inference(json_path, directory.output_dir, do_batch_inference=True)
        for i in range(len(sequences)):
            protein_name = f"{directory.prefix}_{i}"
            prediction_pdb = self._prediction_pdb_path()
            pdb_path = os.path.join(
                directory.prediction_pdbs_dir, f"{protein_name}_{i}.pdb"
            )
            shutil.copyfile(prediction_pdb, pdb_path)
