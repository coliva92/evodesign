from argparse import ArgumentParser
from evodesign.Prediction import Predictor_ESMFold_RemoteApi





parser = ArgumentParser()
parser.add_argument('fastaFilename')
parser.add_argument('pdbFilename')
args = parser.parse_args()
skip_line = True
for line in open(args.fastaFilename, 'rt', encoding='utf-8'):
  if skip_line:
    skip_line = False
    continue
  sequence = line.strip()
  break
predictor = Predictor_ESMFold_RemoteApi()
predictor.predict_structure(sequence, args.pdbFilename)
