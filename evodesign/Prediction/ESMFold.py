from . import Predictor





class ESMFold(Predictor):

  _model = None
  


  @classmethod
  def name(cls) -> str:
    return 'Predictor_ESMFold'
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
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
