from .Predictor import Predictor





class ESMFoldColab(Predictor):

  _model = None
  


  @classmethod
  def name(cls) -> str:
    return 'Predictor_ESMFold_Colab'
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    import torch
    if not ESMFoldColab._model:
      ESMFoldColab._model = torch.load("esmfold.model")
      ESMFoldColab._model.eval().cuda().requires_grad_(False)
      ESMFoldColab._model.set_chunk_size(128)
    prediction = ESMFoldColab._model.infer_pdb(sequence)
    torch.cuda.empty_cache()
    with open(pdbFilename, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(prediction)
