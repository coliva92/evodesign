from .Predictor import Predictor
import subprocess
import os





class AlphaFold(Predictor):

  @classmethod
  def name(cls) -> str:
    return 'Predictor_AlphaFold'
  


  def __init__(self, 
               workspaceRoot: str,
               emptyMsaScriptFolder: str,
               alphafoldFolder: str) -> None:
    super().__init__()
    # para ejecutar AF2 sin MSA, hay que ejecutar el script disponible en:
    # https://github.com/Zuricho/ParaFold_dev/blob/main/parafold/create_fakemsa.py
    # para crear un MSA vacío y alimentarlo a AF2; suponemos que este
    # script ya está disponible
    self.empty_msa_filename = os.path.join(emptyMsaScriptFolder, 
                                           'create_fakemsa.py')
    self.run_docker_filename = os.path.join(alphafoldFolder, 
                                            'docker', 
                                            'run_docker.py')
    self.workspace_root = workspaceRoot
    self.pdbs_folder = os.path.join(workspaceRoot, 'pdbs')



  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    name = f'prot_{sequence}'
    fasta_filename = os.path.join(self.workspace_root, f'{name}.fasta')
    with open(fasta_filename, 'wt', encoding='utf-8') as fasta_file:
      fasta_file.write(f'>prot_{sequence}\n{sequence}\n')
    # ejecutamos el script para crear el MSA vacío
    subprocess.call([
      'python3',
      self.empty_msa_filename,
      f'--fasta_paths={fasta_filename}',
      f'--output_dir={self.pdbs_folder}'
    ])
    for line in self.run_docker(fasta_filename):
      print(line, end='')
    os.remove(fasta_filename)
    prediction_pdb = os.path.join(self.pdbs_folder, name, 'ranked_0.pdb')
    os.symlink(prediction_pdb, pdbFilename)



  def run_docker(self, fastaFilename: str) -> None:
    # ejecutamos AF2 a través de Docker
    #TODO todas las opciones deben obtenerse desde `settings.json`
    cmd = [
      'python3',
      f'{self.run_docker_filename}',
      '--use_precomputed_msas=True',
      f'--fasta_paths={fastaFilename}',
      '--max_template_date=2020-05-14',
      '--model_preset=monomer_casp14',
      '--db_preset=reduced_dbs',
      f'--output_dir={self.pdbs_folder}',
      '--mgnify_database_path=/media/biocomp/"My Passport"/mgnify/mgy_clusters_2018_12.fa',
      '--data_dir=/media/biocomp/"My Passport"/reduced_dbs'
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
    