from .Predictor import Predictor





class ESMFold(Predictor):

  _model = None
  


  @classmethod
  def name(cls) -> str:
    return 'Predictor_ESMFold'
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    """
    Predicts the 3D structure of a given amino acid sequence using the 
    ESMFold v1 model.

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
    import torch
    import esm
    if not ESMFold._model:
      ESMFold._model = esm.pretrained.esmfold_v1()
      ESMFold._model.eval().cuda().requires_grad_(False)
      ESMFold._model.set_chunk_size(128)
    prediction = ESMFold._model.infer_pdb(sequence)
    torch.cuda.empty_cache()
    with open(pdbFilename, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(prediction)
