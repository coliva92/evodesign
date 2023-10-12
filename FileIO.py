from typing import Tuple
from .Population import Population 
from .Individual import Individual
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
import json
import csv
import os





def load_structure_from_pdb(filename: str) -> Structure:
  structure_id = os.path.splitext(os.path.basename(filename))[0]
  parser = PDBParser()
  return parser.get_structure(structure_id, filename)



def save_population_csv(population: Population, filename: str) -> bool:
  def serialize_metrics(jsonData: dict) -> dict:
    a = { key: value for key, value in jsonData.items() if key != 'metrics' }
    b = { key: value for key, value in jsonData['metrics'].items() }
    return { **a, **b }
  
  is_new_file = not os.path.isfile(filename)
  with open(filename, 'wt', encoding='utf-8') as csv_file:
    rows = list(map(serialize_metrics, population.as_json()))
    a, b = rows[-1].keys(), rows[0].keys()
    fieldnames = a if len(a) > len(b) else b
    writer = csv.DictWriter(csv_file, 
                            fieldnames, 
                            dialect='unix', 
                            quoting=csv.QUOTE_NONE)
    writer.writeheader()
    for row in rows:
      writer.writerow(row)
  return is_new_file



def load_population_csv(filename: str, iterationId: int = 0) -> Population:
  def deserialize_metrics(csv_row: dict) -> Individual:
    a = {
      'sequence': csv_row['sequence'], 
      'fitness': float(csv_row['fitness']) if csv_row['fitness'] else None
    }
    b = {
      'metrics': { 
        key: float(value) if len(value) > 0 else None
        for key, value in csv_row.items()
        if key != 'sequence' and key != 'fitness' 
      }
    }
    return Individual(**{ **a, **b })

  with open(filename, 'rt', encoding='utf-8') as csv_file:
    rows = csv.DictReader(csv_file, dialect='unix')
    individuals = list(map(deserialize_metrics, rows))
  return Population(individuals, iterationId)



def save_population_json(population: Population, filename: str) -> bool:
  is_new_file = not os.path.isfile(filename)
  with open(filename, 'wt', encoding='utf-8') as json_file:
    json_file.write(json.dumps(population.as_json(), indent=2) + '\n')
  return is_new_file



def load_population_json(filename: str) -> Population:
  with open(filename, 'rt', encoding='utf-8') as json_file:
    pop_data = json.load(json_file)
  return Population([ Individual(**params) for params in pop_data ])



def save_population_fasta(population: Population, filename: str) -> bool:
  def serialize_metrics(item: Tuple[int, dict]) -> str:
    idx, jsonData = item
    metrics = [ 
      f'{value}' 
      for _, value in jsonData['metrics'].items() 
    ]
    return f'>{idx}|{population.iteration_id}|' + \
      f'{jsonData["fitness"]}|' + \
      '|'.join(metrics) + \
      f'\n{jsonData["sequence"]}'
  
  is_new_file = not os.path.isfile(filename)
  json_data = population.as_json()
  content = '\n'.join(list(map(serialize_metrics, enumerate(json_data))))
  with open(filename, 'wt', encoding='utf-8') as fasta_file:
    fasta_file.write(content + '\n')
  return is_new_file
