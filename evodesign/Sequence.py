import numpy as np



AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')





def create_random(rng: np.random.Generator, length: int) -> str:
    """
    Generates a random sequence of amino acids with uniform probability for
    all amino acids at each residue position. 

    Parameters
    ----------
    rng : numpy.random.Generator
        The RNG to be used to produce the result of the coin toss.
    length : int
        The number of residues in the produced sequence.

    Returns
    -------
    str
        The generated sequence. Each residue is represented as a
        single letter corresponding to one of the 20 standard amino acids.
    """
    return ''.join(rng.choice(AMINO_ACIDS, size=length).tolist())



def swap_letter(rng: np.random.Generator, letter: str) -> str:
    """
    Produce a single random amino acid letter that's different from the one
    specified by `letter`.

    Parameters
    ----------
    rng : numpy.random.Generator
        The RNG to be used to produce the result of the coin toss.
    letter : str
        A single letter representing one of the 20 standard amino acids.

    Returns
    -------
    str
        A single letter that is different from the one specified by `letter`. 
        The produced letter also corresponds to one of the 20 standard amino 
        acids.
    """
    new_letter = letter
    while new_letter == letter:
        new_letter = AMINO_ACIDS[rng.choice(len(AMINO_ACIDS))]
    return new_letter
