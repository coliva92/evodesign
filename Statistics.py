"""Colección de funciones auxiliares para obtener estadísticas poblacionales 
del algoritmo evolutivo, y guardarlas en un archivo.
"""
from typing import List, Dict, Optional
from .Population import Individual
import math
import os
import ast
import matplotlib.pyplot as plt





def compute_statistics(population: List[Individual]) -> Dict[str, float]:
  """Calcula varias estadísticas de la población especificada por `population` 
  y las retorna en un diccionario. Dicho diccionario contiene como mínimo los 
  siguientes datos:
  - `min_fitness`: la aptitud mínima.
  - `fitness_mean`: la aptitud promedio.
  - `max_fitness`: la aptitud máxima.
  - `best_sequence_id`: el identificador único universal (UUID) de aquél 
    individuo contenido en `population` que tiene la aptitud más alta.
  - `best_sequence_fitness`: el valor de aptitud más alto encontrado.
  - `best_sequence`: la secuencia de aminoácidos correspondiente al individuo 
    de `population` cuya aptitud es la más alta. 
  
  El diccionario también puede contener campos adicionales que contienen otras métricas de calidad o atributos correspondientes al individuo con la aptitud más alta.
  """
  stats = {
    'min_fitness': math.inf,
    'fitness_mean': None,
    'max_fitness': -math.inf,
    'best_sequence_id': None,
    'best_sequence_fitness': None,
    'best_sequence': None
  }
  s = 0.0
  for individual in population:
    if individual.fitness < stats['min_fitness']:
      stats['min_fitness'] = individual.fitness
    if individual.fitness > stats['max_fitness']:
      stats['max_fitness'] = individual.fitness
      stats['best_sequence_id'] = individual.id
      stats['best_sequence_fitness'] = individual.fitness
      stats['best_sequence'] = f'"{individual.sequence}"'
    s += individual.fitness
  stats['fitness_mean'] = s / len(population)
  return stats



def save_statistics_to_csv_file(stats: Dict[str, float],
                                iterationId: int,
                                filename: str
                                ) -> None:
  """Calcula varias estadísticas de la población especificada por `population` 
  y las guarda en un archivo CSV, en la locación especificada por `filename`.
  - `iterationId`: el identificador de la iteración que corresponde a la 
    población a respaldar.
  - `population`: la colección de individuos que van a respaldarse. 
  - `filename`: el nombre del archivo donde se guardarán las estadísticas 
    calculadas.
  """
  file_exists = os.path.isfile(filename)
  with open(filename, 'at', encoding='utf-8') as file:
    if not file_exists:
      temp = [ name for name in stats.keys() ]
      file.write('iteration_id,' + ','.join(temp) + '\n')
    temp = [ f'{data}' for data in stats.values() ]
    file.write(f'{iterationId},' + ','.join(temp) + '\n')



def load_statistics_from_csv_file(filename: str) -> List[Dict[str, float]]:
  """Carga las estadísticas poblacionales almacenadas en el archivo CSV 
  especificado por `filename`.
  """
  if not os.path.isfile(filename): return []
  stats = []
  is_heading = True
  for line in open(filename, 'rt', encoding='utf-8'):
    if is_heading:
      column_names = line.strip().split(',')
      is_heading = False
      continue
    values = line.strip().split(',')
    temp = {}
    for i in range(len(column_names)):
      temp[column_names[i]] = ast.literal_eval(values[i])
    stats.append(temp)
  return stats



def plot_fitness_over_iterations(csvFilename: str, 
                                 graphFilename: str,
                                 title: str
                                 ) -> None:
  """Crea un archivo PNG que contiene una gráfica de línea donde se muestra la 
  aptitud mínima, la promedio y la máxima con respecto a cada iteración 
  ejecutada por el algoritmo evolutivo.
  - `csvFilename`: la locación del archivo CSV de donde se cargarán las 
    estadísticas de aptitud.
  - `graphFilename`: la locación donde se guardará la gráfica generada.
  - `title`: el título que se desplegará en la gráfica.
  """
  iterations = []
  minimums = []
  averages = []
  maximums = []
  best_solutions = []
  is_heading = True
  for line in open(csvFilename, 'rt', encoding='utf-8'):
    values = line.strip().split(',')
    if is_heading:
      for i, header in enumerate(values):
        if header == 'iteration_id':
          iteration_id = i
          continue
        if header == 'min_fitness':
          min_fitness = i
          continue
        if header == 'fitness_mean':
          fitness_mean = i
          continue
        if header == 'max_fitness':
          max_fitness = i
          continue
        if header == 'best_sequence_fitness':
          best_solution = i
          continue
      is_heading = False
      continue
    iterations.append(ast.literal_eval(values[iteration_id]))
    minimums.append(ast.literal_eval(values[min_fitness]))
    averages.append(ast.literal_eval(values[fitness_mean]))
    maximums.append(ast.literal_eval(values[max_fitness]))
    best_solutions.append(ast.literal_eval(values[best_solution]))
  fig, ax = plt.subplots(1)
  ax.plot(iterations, averages, label='Fitness mean')
  ax.fill_between(iterations, minimums, maximums, alpha=0.3)
  ax.plot(iterations, best_solutions, label='Best solution found')
  ax.set_title(title)
  ax.set_xlabel('Iterations')
  ax.set_ylabel('Fitness')
  ax.legend(loc='lower right')
  ax.grid()
  fig.savefig(graphFilename)
