from abc import ABC, abstractmethod
from .Predictor import Predictor
from .DirectoryManager import DirectoryManager
from typing import List
from ..System.Subprocess import run_subprocess
import os
import shutil


class AlphaFoldInterface(Predictor, ABC):

    def predict_single_pdb_file(
        self,
        sequence: str,
        protein_name_suffix: str,
        directory: DirectoryManager,
    ) -> None:
        directory.create_folders()
        protein_full_name = directory.protein_full_name(protein_name_suffix)
        input_path = self._create_model_input(
            sequence,
            protein_full_name,
            directory.model_input_dir,
            directory.model_output_dir,
        )
        self.run_inference(
            input_path, directory.model_output_dir, do_batch_inference=False
        )
        prediction_pdb_path = self._prediction_pdb_path(
            protein_full_name, directory.model_output_dir
        )
        output_pdb_path = os.path.join(
            directory.prediction_pdbs_dir, f"{protein_full_name}.pdb"
        )
        shutil.copyfile(prediction_pdb_path, output_pdb_path)
        return

    @abstractmethod
    def _create_model_input(
        self,
        sequence: str,
        protein_full_name: str,
        input_dir: str,
        output_dir: str,
    ) -> str:
        raise NotImplementedError

    @abstractmethod
    def _prediction_pdb_path(
        self,
        protein_full_name: str,
        output_dir: str,
    ) -> str:
        raise NotImplementedError

    @abstractmethod
    def _create_cmd_array(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ) -> List[str]:
        raise NotImplemented

    def run_inference(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ) -> None:
        run_subprocess(
            self._create_cmd_array(input_path, output_dir, do_batch_inference)
        )
