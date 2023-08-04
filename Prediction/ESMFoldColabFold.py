from .Predictor import Predictor





class ESMFoldColabFold(Predictor):

  _model = None
  


  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_ESMFold_ColabFold'
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    import torch
    if not ESMFoldColabFold._model:
      ESMFoldColabFold._model = torch.load("esmfold.model")
      ESMFoldColabFold._model.eval().cuda().requires_grad_(False)
      ESMFoldColabFold._model.set_chunk_size(128)
    torch.cuda.empty_cache()
    prediction = ESMFoldColabFold._model.infer_pdb(sequence)
    with open(pdbFilename, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(prediction)
    # from numba import cuda
    # cuda.select_device(0) 
    # cuda.close()
