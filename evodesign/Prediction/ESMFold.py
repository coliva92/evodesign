from .Predictor import Predictor





class ESMFold(Predictor):

  _model = None
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbPath: str
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
    pdbPath : str
        The path to the PDB file where the predicted structure will
        be stored.
    """
    import torch
    import esm
    if not self.model:
      self.model = esm.pretrained.esmfold_v1()
      self.model.eval().cuda()
      self.model.set_chunk_size(128)
    with torch.no_grad():
      prediction = self.model.infer_pdb(sequence)
    torch.cuda.empty_cache()
    with open(pdbPath, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(prediction)
