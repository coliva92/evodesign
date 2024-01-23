from .Predictor import Predictor
from ..Workspace import Workspace
import subprocess
import os





class AlphaFold(Predictor):

  @classmethod
  def name(cls) -> str:
    return 'Predictor_AlphaFold'
  


  def __init__(self,
               fakeMsaScript: str,
               alphafoldScript: str) -> None:
    """
    Interface for interacting with the AlphaFold 2 model for protein structure 
    prediction.

    Parameters
    ----------
    fakeMsaScript : str
        The path and name for the script provided from https://github.com/Zuricho/ParaFold_dev/blob/main/parafold/create_fakemsa.py.
        When doing "de novo" protein design with AlphaFold, we usually want
        to skip the MSA procedure. In order to do this, we provide AlphaFold
        with an empty MSA file. The aformentioned script allow us to create
        such MSA file.
    alphafoldScript : str
        The path and name for the script that allow us to run AlphaFold.
    """
    super().__init__()
    self.fake_msa_script = f'{fakeMsaScript}/create_fakemsa.py'
    self.alphafold_script = f'{alphafoldScript}/docker/run_docker.py'
    self.alphafold_outputs = f'{Workspace.instance().path}/alphafold_outputs'



  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    """
    Predicts the 3D structure of a given amino acid sequence using the 
    AlphaFold 2 model.

    Parameters
    ----------
    sequence : str
        The amino acid sequence which structure will be predicted. Each residue
        must be represented with a single letter corresponding to one of the
        20 essential amino acids.
    pdbFilename : str
        The path and name of the PDB file where the predicted structure will
        be stored.
    """
    os.makedirs(self.alphafold_outputs, exist_ok=True)
    protein_id = os.path.splitext(os.path.basename(pdbFilename))[0]
    workspace = Workspace.instance()
    fasta_filename = f'{workspace.path}/{protein_id}.fasta'
    with open(fasta_filename, 'wt', encoding='utf-8') as fasta_file:
      fasta_file.write(f'>{protein_id}\n{sequence}\n')
    # run the script for creating an empty MSA
    subprocess.call([
      'python3',
      self.empty_msa_filename,
      f'--fasta_paths={fasta_filename}',
      f'--output_dir={self.alphafold_outputs}'
    ])
    for line in self._run_alphafold_docker(fasta_filename):
      print(line, end='')
    prediction_pdb = f'{self.alphafold_outputs}/{protein_id}/ranked_0.pdb'
    os.symlink(prediction_pdb, pdbFilename)
    os.remove(fasta_filename)



  def _run_alphafold_docker(self, fastaFilename: str) -> None:
    #TODO let these options be input from `settings.json`
    cmd = [
      'python3',
      f'{self.alphafold_script}',
      '--use_precomputed_msas=True',
      f'--fasta_paths={fastaFilename}',
      '--max_template_date=2020-05-14',
      '--model_preset=monomer_casp14',
      '--db_preset=reduced_dbs',
      f'--output_dir={self.alphafold_outputs}',
      '--mgnify_database_path=/media/biocomp/My Passport/mgnify/mgy_clusters_2018_12.fa',
      '--data_dir=/media/biocomp/My Passport/reduced_dbs'
    ]
    popen = subprocess.Popen(cmd, 
                             stdout=subprocess.PIPE, 
                             universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
      yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
      raise subprocess.CalledProcessError(return_code, cmd)
    