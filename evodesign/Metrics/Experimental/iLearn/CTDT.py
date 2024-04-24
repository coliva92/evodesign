from Descriptor import Descriptor





class CTDT(Descriptor):

  @classmethod
  def column_name(cls) -> str:
    return 'CTDC'
  


  def __init__(self, ilearnDir: str) -> None:
    super().__init__(ilearnDir, 'CTDT')
