from .SimpleGeneticAlgorithm import SimpleGeneticAlgorithm
from typing import Callable
from ..Fitness import FitnessFunction
from ..Prediction import Predictor, Predictor_Null
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from ..Population import Population





class SurrogatedGeneticAlgorithm(SimpleGeneticAlgorithm):

  _null_predictor = Predictor_Null()



  @classmethod
  def get_name(cls) -> str:
    return 'GA_Simple_Surrogated'
  

  
  def __init__(self, 
               workspaceRoot: str, 
               targetPdbFilename: str, 
               predictor: Predictor, 
               fitnessFunction: FitnessFunction, 
               surrogateFitnessFunction: FitnessFunction,
               populationSize: int, 
               numIterations: int, 
               numPredictions: int,
               selection: Selection, 
               recombination: Recombination, 
               mutation: Mutation
               ) -> None:
    super().__init__(workspaceRoot, 
                     targetPdbFilename, 
                     predictor, 
                     fitnessFunction, 
                     populationSize, 
                     numIterations, 
                     selection, 
                     recombination, 
                     mutation)
    self._num_predictions = numPredictions
    self._surrogate_fitness_fn = surrogateFitnessFunction
    # TODO: buscar una manera más limpia de realizar esta inicializacion
    if 'energy_score' in surrogateFitnessFunction._metric_calculators:
      temp = surrogateFitnessFunction._metric_calculators['energy_score']
      temp.initialize(self.workspace.reference_filename,
                      self.workspace.root_folder,
                      self.workspace.pdbs_folder)
  


  def _update_fitness(self, 
                      population: Population, 
                      saveFunction: Callable[[Population], None]
                      ) -> bool:
    population.update_fitness(self._surrogate_fitness_fn, 
                              SurrogatedGeneticAlgorithm._null_predictor, 
                              self._reference_backbone)
    try:
      # deberíamos definir aquí un criterio para seleccionar los
      # individuos a quienes se les predecirá la estructura, en base al fitness
      # surrogado; sin embargo, por motivos de tiempo y porque esta no
      # es una característica inmediatamente necesaria, se omitirá por ahora
      top = Population(list(filter(lambda ind: 'rmsd' not in ind.metrics or \
                                               'lddt' not in ind.metrics, 
                                   population[-self._num_predictions:])))
      for individual in top:
        if 'surrogate_fitness' not in individual.metrics:
          individual.metrics['surrogate_fitness'] = individual.fitness
          individual.fitness = None
      saveFunction(population)
      result = top.update_fitness(self._fitness_fn, 
                                  self._predictor, 
                                  self._reference_backbone, 
                                  self.workspace.pdbs_folder)
    except BaseException as e:
      saveFunction(population)
      self.workspace.save_rng_settings()
      raise e
    population.individuals = population[:-self._num_predictions] + \
                             top.individuals
    saveFunction(population)
    return result
  