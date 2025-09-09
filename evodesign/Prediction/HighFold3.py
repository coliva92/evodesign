from .AlphaFold3 import AlphaFold3
from typing import List, Optional


class HighFold3(AlphaFold3):

    def __init__(
        self,
        path_to_run_alphafold_py: str,
        model_dir: str,
        output_dir: str,
        head_to_tail: bool = True,
        disulfide_chain_res: Optional[List[List[int]]] = None,
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
        super().__init__(
            path_to_run_alphafold_py,
            model_dir,
            output_dir,
            num_recycles,
            num_diffusion_samples,
            model_seeds,
            version,
            run_data_pipeline,
            max_template_date,
            gpu_device,
            jackhmmer_binary_path,
            nhmmer_binary_path,
            hmmalign_binary_path,
            hmmsearch_binary_path,
            hmmbuild_binary_path,
            db_dir,
            small_bfd_database_path,
            mgnify_database_path,
            uniprot_cluster_annot_database_path,
            uniref90_database_path,
            ntrna_database_path,
            rfam_database_path,
            rna_central_database_path,
            pdb_database_path,
            seqres_database_path,
            jax_compilation_cache_dir,
        )
        self.head_to_tail = head_to_tail
        self.disulfide_chain_res = disulfide_chain_res

    def _get_cmd_array(
        self, input_path: str, use_batch_inference: bool = False
    ) -> List[str]:
        cmd = super()._get_cmd_array(input_path, use_batch_inference)
        # the following also removes --force_output_dirs:
        cmd[-1] = f"--head_to_tail={self.head_to_tail}"
        if self.disulfide_chain_res is not None:
            cmd.append(f"--disulfide_chain_res {self.disulfide_chain_res}")
        return cmd
