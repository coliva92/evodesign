from Descriptor import Descriptor





class CTDC(Descriptor):

  @classmethod
  def column_name(cls) -> str:
    return 'CTDC'
  


  def __init__(self, ilearnDir: str) -> None:
    super().__init__(ilearnDir, 'CTDC')
