import evodesign.Sequence as Sequence
from evodesign import Random



sequence = Sequence.random_sequence(10)
test1_passed = type(sequence) == str and \
               len(sequence) == 10 and \
               sum([ i in Sequence.AMINO_ACIDS for i in sequence ]) == len(sequence)

rng = Random.get_rng()
old_letter = rng.choice(len(Sequence.AMINO_ACIDS))
new_letter = Sequence.switch_residue(old_letter)
test2_passed = type(new_letter) == str and \
               len(new_letter) == 1 and \
               new_letter != old_letter

print('PASSED' if test1_passed and test2_passed else 'FAILED')
