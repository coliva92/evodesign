from abc import ABC, abstractmethod
from .Predictor import Predictor
from .DirectoryManager import DirectoryManager
from typing import List
from ..Utils.Subprocess import run_subprocess
import os
import shutil


class AlphaFoldInterface(Predictor, ABC):

    def predict_single_pdb_file(
        self,
        sequence: str,
        directory: DirectoryManager,
    ) -> None:
        os.makedirs(directory.model_input_dir, exist_ok=True)
        os.makedirs(directory.model_output_dir, exist_ok=True)
        input_path = self._create_model_input(
            sequence,
            directory.prefix,
            directory.model_input_dir,
            directory.model_output_dir,
        )
        self.run_inference(
            input_path, directory.model_output_dir, do_batch_inference=False
        )
        prediction_pdb = self._prediction_pdb_path(
            directory.model_output_dir, directory.prefix
        )
        shutil.copyfile(prediction_pdb, directory.prediction_pdbs_dir)

    @abstractmethod
    def _create_model_input(
        self,
        sequence: str,
        protein_name: str,
        input_dir: str,
        output_dir: str,
    ) -> str:
        raise NotImplementedError

    @abstractmethod
    def _prediction_pdb_path(
        self,
        protein_name: str,
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
