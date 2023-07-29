"""Colección de funciones auxiliares para obtener estadísticas poblacionales 
del algoritmo evolutivo, y guardarlas en un archivo.
"""
from typing import List, Dict
from .Population import Individual
import math
import os
import ast
import matplotlib.pyplot as plt
from dataclasses import dataclass





@dataclass()
class PopulationStats:

  min_fitness: float = math.inf
  fitness_mean: float = -math.inf
  max_fitness: float = -math.inf
  best_sequence_id: int = -1
  best_sequence_fitness: float = -math.inf
  best_sequence: str = ''





def compute_statistics(population: List[Individual]) -> PopulationStats:
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
  stats = PopulationStats()
  s = 0.0
  for individual in population:
    if individual.fitness < stats.min_fitness:
      stats.min_fitness = individual.fitness
    if individual.fitness > stats.max_fitness:
      stats.max_fitness = individual.fitness
      stats.best_sequence_id = individual.id
      stats.best_sequence_fitness = individual.fitness
      stats.best_sequence = f'"{individual.sequence}"'
    s += individual.fitness
  stats.fitness_mean = s / len(population)
  return stats



def save_statistics(iterationId: int,
                    stats: PopulationStats,
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
  with open(filename, 'at', encoding='utf-8') as the_file:
    if not file_exists:
      temp = [ name for name, _ in vars(stats).items() ]
      the_file.write('iteration_id,' + ','.join(temp) + '\n')
    temp = [ f'{data}' for _, data in vars(stats).items() ]
    the_file.write(f'{iterationId},' + ','.join(temp) + '\n')



def load_statistics(filename: str) -> List[PopulationStats]:
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
    data = {}
    for i in range(len(column_names)):
      if column_names[i] != 'iteration_id':
        data[column_names[i]] = ast.literal_eval(values[i])
    stats.append(PopulationStats(**data))
  return stats



def plot_fitness(csvFilename: str, 
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
  ax.legend(loc='best')
  ax.grid()
  fig.savefig(graphFilename)
