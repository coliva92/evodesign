from .Predictor import Predictor





class ESMFold(Predictor):

    model = None



    def _params(self) -> dict:
        return { 'use_gpu': self._use_gpu }
  


    def __init__(self, use_gpu: bool = True) -> None:
        """
        Interface for interacting with the ESMFold v1 model for protein structure 
        prediction.

        Parameters
        ----------
        use_gpu : bool
            If true, the model will run in the first GPU available. Otherwise, the
            model will run entirely on the CPU.
        """
        super().__init__()
        self._use_gpu = use_gpu
  


    def predict_raw_pdb(self, sequence: str) -> str:
        """
        Predicts the 3D structure of a given amino acid sequence using the 
        ESMFold v1 model.

        Parameters
        ----------
        sequence : str
            The amino acid sequence which structure will be predicted. Each residue
            must be represented with a single letter corresponding to one of the
            20 essential amino acids.
        
        Returns
        -------
        str
            The predicted 3D structure codified in a string as a PDB file.
        """
        import torch
        if self.model is None:
            import esm
            self.model = esm.pretrained.esmfold_v1()
            self.model.eval()
        if self._use_gpu and torch.cuda.is_available():
            self.model = self.model.cuda()
            self.model.set_chunk_size(128)
        with torch.no_grad():
            prediction = self.model.infer_pdb(sequence)
        return prediction
