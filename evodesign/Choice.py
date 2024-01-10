from typing import Tuple
import random





def flip_coin(weights: Tuple[float] = (0.5, 0.5)) -> bool:
  return random.choices([ True, False ], weights)[0]
