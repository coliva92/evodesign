# EvoDesign

Esta es un programa y un _framework_ rudimentario, escrito en Python 3.8, para ejecutar diferentes algoritmos evolutivos para el diseño de proteínas. En el diseño de proteínas, el objetivo es encontrar la secuencia de aminoácidos que se pliega en una estructura objetivo determinada.

## 1. Instalación

Para utilizar EvoDesign, primero se debe instalar BioPython y Matplotlib en sus versiones más recientes. **Estos dos paquetes de software ya vienen incluídos en Anaconda Python**.

- [Instrucciones de instalación de BioPython](https://biopython.org/wiki/Download).
- [Instrucciones de instalación de Matplotlib](https://matplotlib.org/stable/users/getting_started/index.html#installation-quick-start).

## 2. Instrucciones de uso desde la línea de comandos

Para correr un algoritmo evolutivo utilizando EvoDesign, se deben seguir los siguientes pasos:

1. Crear la carpeta donde se guardarán los datos producidos durante la ejecución del algoritmo. En esta documentación, nos referiremos a esta carpeta como el _workspace_ del algoritmo.
2. Dentro de la carpeta _workspace_, crear un archivo de texto llamado `setup.json`.
3. Llenar el archivo `setup.json` con la configuración del algoritmo que se desea ejecutar. En el ejemplo siguiente se describe detalladamente el formato en que debe escribirse dicha configuración.
4. Comenzar la ejecución del programa usando el comando: `python3 -m evodesign <workspace>/setup.json`. 

### Ejemplo

Supóngase que se desea ejecutar un algoritmo genético que cumpla con las siguientes características:

- **Representación de los individuos**: secuencias de aminoácidos (cada cadena lateral se representa por medio de una letra).
- **Función de aptitud**: el RMSD de la superposición contra la estructura objetivo. 
- **Tamaño de la población**: 10 individuos.
- **Selección**: por torneo entre 3 individuos.
- **Número de padres e hijos**: 2.
- **Recombinación**: por _crossover_ de dos puntos.
- **Mutación**: cambiar un aminoácido por otro de manera aleatoria.
- **Número de aminoácidos a mutar**: 1.
- **Probabilidad de mutación**: 10%.
- **Reemplazo**: se preservan los 10 individuos de mayor aptitud.
- **Máximo número de iteraciones**: 20.

Para lograr esto, el archivo de configuración `setup.json` debe contener lo siguiente:

```json
{
  "algorithmType": "GA",
  "algorithmParams": {
    "taskName": "example",
    "referencePdbFilename": "example/1Y32.pdb",
    "predictor": "Predictor_ESMFold",
    "fitnessFunction": "Fitness_NegativeRmsd",
    "populationSize": 10,
    "numIterations": 20,
    "selection": "GA_Selection_Tournament",
    "selectionParams": {
      "selectionSize": 2,
      "tournamentSize": 3
    },
    "recombination": "GA_Recombination_TwoPointCrossover",
    "recombinationParams": {},
    "mutation": "GA_Mutation_MultipleSwitches",
    "mutationParams": {
      "probability": 0.1,
      "numSwitches": 1
    },
    "replacement": "GA_Replacement_WorseOut",
    "replacementParams": {}
  }
}
```

La configuración anterior puede interpretarse de la siguiente manera:

- El campo `"algorithmType": "GA"` indica que se desea usar un algoritmo genético tradicional.
- El campo `"workspaceName": "example"` indica el nombre de la carpeta _workspace_ donde se guardarán los datos producidos por el algoritmo.
- El campo `"referencePdbFilename": "example/1Y32.pdb"` indica la ubicación del archivo PDB que contiene la estructura de entrada. **No es necesario que dicho archivo esté localizado en la carpeta _workspace_**.
- El campo `"fitnessFunction": "Fitness_NegativeRmsd"` indica que la función de aptitud a utilizar es el RMSD.
- El campo `"populationSize": 10` indica que el tamaño de la población es 10.
- El campo `"numIterations: 20"` indica que el algoritmo correrá por 20 iteraciones.
- El campo `"selection": "GA_Selection_Tournament"` indica que se usará selección por torneo. Así mismo, los campos `"selectionSize": 2` y `"tournamentSize": 3`, en la sección `"selectionParams"`, indican respectivamente que se seleccionarán dos padres y que el torneo se realizará entre 3 individuos.
- El campo `"recombination": "GA_Recombination_TwoPointCrossover"` indica que se usará la [recombinación _crossover_ de dos puntos](https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)).
- El campo `"mutation": "GA_Recombation_MultipleSwitches"` indica que la mutación consistirá en cambiar un aminoácido de la secuencia por otro aminoácido elegido de manera aleatoria. Así mismo, `"probability"`
- El campo `""`
