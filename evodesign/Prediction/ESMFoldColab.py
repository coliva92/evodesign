from .Predictor import Predictor





class ESMFoldColab(Predictor):

  _model = None
  


  @classmethod
  def name(cls) -> str:
    return 'Predictor_ESMFold_Colab'
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbPath: str
                        ) -> None:
    """
    Predicts the 3D structure of a given amino acid sequence using the 
    ESMFold v1 model. This predictor can only be used within a Google Colab 
    notebook.

    Parameters
    ----------
    sequence : str
        The amino acid sequence which structure will be predicted. Each residue
        must be represented with a single letter corresponding to one of the
        20 essential amino acids.
    pdbPath : str
        The path to the PDB file where the predicted structure will
        be stored.
    """
    import torch
    if not ESMFoldColab._model:
      ESMFoldColab._model = torch.load("esmfold.model")
      ESMFoldColab._model.eval().cuda().requires_grad_(False)
      ESMFoldColab._model.set_chunk_size(128)
    prediction = ESMFoldColab._model.infer_pdb(sequence)
    torch.cuda.empty_cache()
    with open(pdbPath, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(prediction)
