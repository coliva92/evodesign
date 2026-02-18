from .AlphaFoldInterface import AlphaFoldInterface
from .DirectoryManager import DirectoryManager
from ..Files.mmCIF import convert_cif_to_pdb
from ..System.Subprocess import run_subprocess
from Bio.PDB import MMCIFParser, PDBIO
import os
import shutil
import json
from typing import List, Optional


class AlphaFold3Docker(AlphaFoldInterface):

    _parser = MMCIFParser()
    _io = PDBIO()

    def __init__(
        self,
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
        self.jackhmmer_binary_path = (
            os.path.abspath(jackhmmer_binary_path)
            if jackhmmer_binary_path is not None
            else jackhmmer_binary_path
        )
        self.nhmmer_binary_path = (
            os.path.abspath(nhmmer_binary_path)
            if nhmmer_binary_path is not None
            else nhmmer_binary_path
        )
        self.hmmalign_binary_path = (
            os.path.abspath(hmmalign_binary_path)
            if hmmalign_binary_path is not None
            else hmmalign_binary_path
        )
        self.hmmsearch_binary_path = (
            os.path.abspath(hmmsearch_binary_path)
            if hmmsearch_binary_path is not None
            else hmmsearch_binary_path
        )
        self.hmmbuild_binary_path = (
            os.path.abspath(hmmbuild_binary_path)
            if hmmbuild_binary_path is not None
            else hmmbuild_binary_path
        )
        self.db_dir = os.path.abspath(db_dir) if db_dir is not None else db_dir
        self.small_bfd_database_path = (
            os.path.abspath(small_bfd_database_path)
            if small_bfd_database_path is not None
            else small_bfd_database_path
        )
        self.mgnify_database_path = (
            os.path.abspath(mgnify_database_path)
            if mgnify_database_path is not None
            else mgnify_database_path
        )
        self.uniprot_cluster_annot_database_path = (
            os.path.abspath(uniprot_cluster_annot_database_path)
            if uniprot_cluster_annot_database_path is not None
            else uniprot_cluster_annot_database_path
        )
        self.uniref90_database_path = (
            os.path.abspath(uniref90_database_path)
            if uniref90_database_path is not None
            else uniref90_database_path
        )
        self.ntrna_database_path = (
            os.path.abspath(ntrna_database_path)
            if ntrna_database_path is not None
            else ntrna_database_path
        )
        self.rfam_database_path = (
            os.path.abspath(rfam_database_path)
            if rfam_database_path is not None
            else rfam_database_path
        )
        self.rna_central_database_path = (
            os.path.abspath(rna_central_database_path)
            if rna_central_database_path is not None
            else rna_central_database_path
        )
        self.pdb_database_path = (
            os.path.abspath(pdb_database_path)
            if pdb_database_path is not None
            else pdb_database_path
        )
        self.seqres_database_path = (
            os.path.abspath(seqres_database_path)
            if seqres_database_path is not None
            else seqres_database_path
        )
        self.jax_compilation_cache_dir = (
            os.path.abspath(jax_compilation_cache_dir)
            if jax_compilation_cache_dir is not None
            else jax_compilation_cache_dir
        )
        return

    def _create_model_input(
        self,
        sequence: str,
        protein_full_name: str,
        input_dir: str,
        output_dir: str,
    ) -> str:
        input_json = self.create_input_json(sequence, protein_full_name)
        json_path = os.path.join(input_dir, f"{protein_full_name}.json")
        with open(json_path, "wt", encoding="utf-8") as json_file:
            json.dump(input_json, json_file)
        return json_path

    def _prediction_pdb_path(
        self,
        protein_full_name: str,
        output_dir: str,
    ) -> str:
        prediction_dir = os.path.join(output_dir, protein_full_name)
        prediction_cif = os.path.join(prediction_dir, f"{protein_full_name}_model.cif")
        return convert_cif_to_pdb(prediction_cif, parser=self._parser, io=self._io)

    def _create_cmd_array(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ) -> List[str]:
        cmd = [
            "docker",
            "run",
            "-it",
            "--volume",
            (
                f"{input_path}:/root/af_input"
                if do_batch_inference
                else f"{input_path}:/root/af_input/input.json"
            ),
            "--volume",
            f"{output_dir}:/root/af_output",
            "--volume",
            f"{self.model_dir}:/root/models",
            "--gpus",
            "all",
            "alphafold3",
            "python",
            "run_alphafold.py",
            (
                f"--input_dir=/root/af_input"
                if do_batch_inference
                else f"--json_path=/root/af_input/input.json"
            ),
            "--model_dir=/root/models",
            "--output_dir=/root/af_output",
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
        protein_full_name: str,
        structure_id: str = "A",
    ) -> dict:
        # Assume we're only predicting monomeric structures
        return {
            "name": protein_full_name,
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
            self._create_model_input(
                sequence,
                protein_name,
                directory.model_input_dir,
                directory.model_output_dir,
            )
        self.run_inference(
            directory.model_input_dir,
            directory.model_output_dir,
            do_batch_inference=True,
        )
        for i in range(len(sequences)):
            protein_name = f"{directory.prefix}_{i}"
            prediction_pdb = self._prediction_pdb_path(
                protein_name, directory.model_output_dir
            )
            pdb_path = os.path.join(
                directory.prediction_pdbs_dir, f"{protein_name}.pdb"
            )
            shutil.copyfile(prediction_pdb, pdb_path)

    def run_inference(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ) -> None:
        run_subprocess(
            self._create_cmd_array(input_path, output_dir, do_batch_inference)
        )
        # AF3's docker image creates the output files as root;
        # we must change the owner
        cmd = [
            "docker",
            "run",
            "-it",
            "--volume",
            f"{output_dir}:/root/af_output",
            "alphafold3",
            "chown",
            "-R",
            "1000:1000",
            "/root/af_output",
        ]
        run_subprocess(cmd)
