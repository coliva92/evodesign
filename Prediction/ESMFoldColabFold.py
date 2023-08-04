from .Predictor import Predictor





class ESMFoldColabFold(Predictor):
  
  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_ESMFold_ColabFold'
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    import torch
    model = torch.load("esmfold.model")
    model = model.eval().cuda().requires_grad_(False)
    model.set_chunk_size(128)
    torch.cuda.empty_cache()
    prediction = model.infer_pdb(sequence, 
                                 num_recycles=3,
                                 residue_index_offset=512)
    with open(pdbFilename, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(model.output_to_pdb(prediction)[0])
