import evodesign.Random as r



AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')





def create_random(length: int) -> str:
    """
    Generates a random sequence of amino acids with uniform probability for
    all amino acids at each residue position. 

    Parameters
    ----------
    length : int
        The number of residues in the produced sequence.

    Returns
    -------
    str
        The generated sequence. Each residue is represented as a
        single letter corresponding to one of the 20 standard amino acids.
    """
    rng = r.generator()
    return ''.join(rng.choice(AMINO_ACIDS, size=length).tolist())



def switch_letter(letter: str) -> str:
    """
    Produce a single random amino acid letter that's different from the one
    specified by `letter`.

    Parameters
    ----------
    letter : str
        A single letter representing one of the 20 standard amino acids.

    Returns
    -------
    str
        A single letter that is different from the one specified by `letter`. 
        The produced letter also corresponds to one of the 20 standard amino 
        acids.
    """
    rng = r.generator()
    new_letter = letter
    while new_letter == letter:
        new_letter = AMINO_ACIDS[rng.choice(len(AMINO_ACIDS))]
    return new_letter
