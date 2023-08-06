# EvoDesign

Este es un programa y un _framework_ rudimentario, escrito en Python 3.8, para implementar diferentes [algoritmos evolutivos](https://en.wikipedia.org/wiki/Evolutionary_algorithm) para el _diseño de proteínas_. 
En el diseño de proteínas, el objetivo es encontrar la secuencia de aminoácidos que se pliega en una estructura objetivo determinada. 
Para verificar que una secuencia particular se pliega en la estructura objetivo, los algoritmos evolutivos implementados con EvoDesign utilizan un algoritmo de _predicción de la estructura de una proteína_, como por ejemplo, [ESMFold](https://github.com/facebookresearch/esm).

## Contenido

1. [Instalación](#instalacion)
2. [Instrucciones de uso desde la consola](#instrucciones-consola)
  1. [Ejemplo](#ejemplo)
  2. [Catálogo completo de las diferentes opciones de configuración](#opciones-catalogo)
  3. [Condición de finalización del algoritmo evolutivo](#condicion-finalizacion)
  4. [Reanudar la ejecución del programa a partir de una interrupción anterior](#reanudar-ejecucion)
  5. [Extender la ejecución del programa](#extender-ejecucion)
3. [Instrucciones de uso desde Python](#instrucciones-python)

<a name="instalacion"></a>
## Instalación

Para utilizar EvoDesign, se recomienda instalar [Anaconda](https://www.anaconda.com/).
Alternativamente, las dependencias pueden instalarse de manera individual:

- [Instrucciones de instalación de BioPython](https://biopython.org/wiki/Download).
- [Instrucciones de instalación de Matplotlib](https://matplotlib.org/stable/users/getting_started/index.html#installation-quick-start).

<a name="instrucciones-consola"></a>
## Instrucciones de uso desde la consola

Para correr un algoritmo evolutivo utilizando EvoDesign, se deben seguir los siguientes pasos:

1. Crear la carpeta donde se guardarán los datos producidos durante la ejecución del algoritmo. En esta documentación, nos referiremos a esta carpeta como el _workspace_ del algoritmo.
2. Dentro de la carpeta _workspace_, crear un archivo de texto llamado `settings.json`.
3. Llenar el archivo `settings.json` con la configuración del algoritmo que se desea ejecutar. En el ejemplo siguiente se describe el formato en que debe escribirse dicha configuración.
4. Comenzar la ejecución del programa usando el comando: `python3 -m evodesign <workspace>/settings.json`. 

<a name="ejemplo"></a>
### Ejemplo

Supóngase que se desea ejecutar un algoritmo genético que cumpla con las siguientes características:

- **Estructura objetivo**: aquella identificada por el PDBID 1Y32.
- **Representación de los individuos**: secuencias de aminoácidos, en donde cada cadena lateral se representa por medio de una letra.
- **Función de aptitud**: el GDT de la superposición contra la estructura objetivo. 
- **Tamaño de la población**: 10 individuos.
- **Selección**: por [torneo](https://en.wikipedia.org/wiki/Tournament_selection) entre 3 individuos.
- **Número de padres e hijos**: 2.
- **Recombinación**: por [_crossover_ de dos puntos](https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)).
- **Mutación**: cambiar un aminoácido por otro de manera aleatoria.
- **Probabilidad de mutación**: 10%.
- **Reemplazo**: se preservan los 10 individuos de mayor aptitud.
- **Máximo número de iteraciones**: 20.

Todas estas características pueden alimentarse a EvoDesign por medio de un archivo JSON, denominado `settings.json`, que debe contener lo siguiente:

```json
{
  "algorithmType": "GA_Steady",
  "algorithmParams": {
    "workspaceName": "example",
    "targetPdbFilename": "example/1Y32.pdb",
    "predictor": "Predictor_ESMFold_RemoteApi",
    "fitnessFunction": "Fitness_GDT",
    "populationSize": 10,
    "numIterations": 20,
    "selection": "GA_Selection_Tournament",
    "selectionParams": {
      "selectionSize": 2,
      "tournamentSize": 3
    },
    "recombination": "GA_Recombination_TwoPointsCrossover",
    "mutation": "GA_Mutation_SingleSwitch",
    "mutationParams": {
      "probability": 0.1
    },
    "replacement": "GA_Replacement_WorseOut"
  }
}
```

Una vez escrita este archivo de configuración, el programa puede ejecutarse ingresando el siguiente comando en la consola: 

```
python -m evodesign example/settings.json
```

Al finalizar la ejecución del programa, la carpeta _workspace_ debería contener la siguiente estructura:

```
.
├─ example/
|  ├─ pdbs/
|  |  ├─ prot_ABC.pdb
|  |  ├─ ...
|  ├─ populations/
|  |  ├─ pop_0.json
|  |  ├─ ...
|  ├─ ~children.tmp
|  ├─ fitness.png
|  ├─ settings.json
|  ├─ statistics.csv
```

- La carpeta `pdbs` contiene los archivos PDB de las estructuras que fueron predichas para todas las secuencias producidas durante la ejecución del programa (y que cuya aptitud fue evaluada).
- La carpeta `populations` contiene las poblaciones (esto es, el conjunto de secuencias) producidas en cada iteración del algoritmo evolutivo. En cada uno de estos archivos se registra la aptitud y demás métricas de cada secuencia.
- El archivo `~children.tmp` es un archivo temporal que se usa durante la ejecución del programa. **Este archivo puede ignorarse**.
- El archivo `fitness.png` contiene la gráfica que muestra la aptitud de la población en cada iteración del algoritmo.
- El archivo `settings.json` es el mismo archivo de configuración que se mencionó anteriormente. No se modifica durante la ejecución del algoritmo. 
- El archivo `statistics.csv` contiene datos relevantes de cada iteración del algoritmo. Estos datos son los siguientes:
  - La aptitud mínima, promedio y máxima de la población.
  - El valor de aptitud y la secuencia correspondiente al individuo de mayor aptitud encontrado por el algoritmo.

<a name="reanudar-ejecucion"></a>
### Reanudar la ejecución del programa a partir de una interrupción anterior

Es posible que en algún punto se vea interrumpida la ejecución del programa (ya sea porque ocurrió un error, porque el API remoto de ESMFold falló en responder a una petición, o porque el usuario interrumpió el programa deliveradamente). 
En estos casos, es posible reanudar la ejecución del programa, continuando desde la iteración que fue interrumpida. Para ello, simplemente hay que volver a ingresar en la consola el mismo comando mostrado anteriormente: 

```
python -m evodesign example/settings.json
```

<a name="extender-ejecucion"></a>
### Extender la ejecución del programa

Supóngase que el algoritmo utilizado en el ejemplo anterior terminó de ejecutar las 20 iteraciones que fueron especificadas en `example/settings.json`, y supóngase que se desea correr este mismo algoritmo por 10 iteraciones adicionales. Para lograr esto, simplemente debe modificarse el archivo `example/settings.json` para cambiar el campo que originalmente decía `"numIterations": 20` por `"numIterations": 30`. Luego, se vuelve a ejecutar el comando `python -m evodesign example/settings.json`.

<a name="instrucciones-consola"></a>
## Instrucciones de uso desde Python

El archivo `settings.json` mostrado en el [ejemplo de la sección anterior](#ejemplo), es equivalente a escribir el siguiente programa en Python:

```python
from evodesign.GA import GA_Steady
from evodesign.Prediction import Predictor_ESMFold_RemoteApi
from evodesign.Fitness import Fitness_GDT
from evodesign.GA.Selection import GA_Selection_Tournament
from evodesign.GA.Recombination import GA_Recombination_TwoPointsCrossover
from evodesign.GA.Mutation import GA_Mutation_SingleSwitch
from evodesign.GA.Replacement import GA_Replacement_WorseOut

algorithm = GA_Steady(workspaceName='example',
                      targetPdbFilename='example/1Y32.pdb',
                      predictor=Predictor_ESMFold_RemoteApi(),
                      fitnessFunction=Fitness_GDT(),
                      populationSize=10,
                      numIterations=20,
                      selection=GA_Selection_Tournament(selectionSize=2,
                                                        tournamentSize=3),
                      recombination=GA_Recombination_TwoPointsCrossover(),
                      mutation=GA_Mutation_SingleSwitch(probability=0.1),
                      replacement=GA_Replacement_WorseOut())
algorithm() # iniciamos la ejecución del algoritmo
```

Alternativamente, es posible replicar la ejecución del algoritmo anterior utilizando las clases individuales, como se muestra en el siguiente ejemplo:

```python
from evodesign.Prediction import Predictor_ESMFold_RemoteApi
from evodesign.Fitness import Fitness_GDT
from evodesign.GA.Selection import GA_Selection_Tournament
from evodesign.GA.Recombination import GA_Recombination_TwoPointsCrossover
from evodesign.GA.Mutation import GA_Mutation_SingleSwitch
from evodesign.GA.Replacement import GA_Replacement_WorseOut
import evodesign.Chain as Chain
from evodesign import Population, Statistics

fitnessFn = Fitness_GDT()
predictor = Predictor_ESMFold_RemoteApi()
selection = GA_Selection_Tournament(selectionSize=2, tournamentSize=3)
recombination = GA_Recombination_TwoPointsCrossover()
mutation = GA_Mutation_SingleSwitch(probability=0.1)
replacement = GA_Replacement_WorseOut()

target_structure = Chain.load_structure_from_pdb('example/1Y32.pdb')
sequence_length = Chain.count_chain_residues(target_structure)
target_backbone = Chain.filter_backbone_atoms_in_chain(target_structure)

population = Population.new_random(size=10, 
                                   sequenceLength=sequence_length,
                                   iterationId=1)
population.update_fitness(fitnessFn, predictor, target_backbone)
# después de calcular la aptitud, la población se ordena de manera ascendente
print(Statistics.new_from_population(population))
while True:
  if population.iterationId == 20: # max. 20 iteraciones
    break
  if population[-1].fitness >= fitnessFn.upper_bound():
    break
  parents = selection(population)
  children = recombination(parents)
  mutation(children)
  children.update_fitness(fitnessFn, predictor, target_backbone)
  population = replacement(population, children)
  print(Statistics.new_from_population(population))
print(population[-1]) # mejor solución
```
