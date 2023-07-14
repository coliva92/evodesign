# EvoDesign

Esta es un programa y un _framework_ rudimentario, escrito en Python 3.8, para ejecutar diferentes algoritmos evolutivos para el diseño de proteínas. En el diseño de proteínas, el objetivo es encontrar la secuencia de aminoácidos que se pliega en una estructura objetivo determinada.

## Contenido

1. [Instalación](#instalacion)
2. [Instrucciones de uso desde la consola](#instrucciones-consola)
  1. 1. [Ejemplo](#ejemplo)
  2. [Condición de finalización del algoritmo evolutivo](#condicion-finalizacion)
  3. [Reanudar la ejecución del programa a partir de una interrupción anterior](#reanudar-ejecucion)
  4. [Extender la ejecución del programa](#extender-ejecucion)

<a name="instalacion"></a>
## Instalación

Para utilizar EvoDesign, primero se debe instalar BioPython y Matplotlib en sus versiones más recientes. **Si usted usa Anaconda Python, no es necesario que instale estos paquetes de software**.

- [Instrucciones de instalación de BioPython](https://biopython.org/wiki/Download).
- [Instrucciones de instalación de Matplotlib](https://matplotlib.org/stable/users/getting_started/index.html#installation-quick-start).

<a name="instrucciones-consola"></a>
## Instrucciones de uso desde la consola

Para correr un algoritmo evolutivo utilizando EvoDesign, se deben seguir los siguientes pasos:

1. Crear la carpeta donde se guardarán los datos producidos durante la ejecución del algoritmo. En esta documentación, nos referiremos a esta carpeta como el _workspace_ del algoritmo.
2. Dentro de la carpeta _workspace_, crear un archivo de texto llamado `setup.json`.
3. Llenar el archivo `setup.json` con la configuración del algoritmo que se desea ejecutar. En el ejemplo siguiente se describe el formato en que debe escribirse dicha configuración.
4. Comenzar la ejecución del programa usando el comando: `python3 -m evodesign <workspace>/setup.json`. 

<a name="ejemplo"></a>
### Ejemplo

Supóngase que se desea ejecutar un algoritmo genético que cumpla con las siguientes características:

- **Representación de los individuos**: secuencias de aminoácidos (cada cadena lateral se representa por medio de una letra).
- **Función de aptitud**: el RMSD de la superposición contra la estructura objetivo. 
- **Tamaño de la población**: 10 individuos.
- **Selección**: por [torneo](https://en.wikipedia.org/wiki/Tournament_selection) entre 3 individuos.
- **Número de padres e hijos**: 2.
- **Recombinación**: por [_crossover_ de dos puntos](https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)).
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
    "workspaceName": "example",
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
- El campo `"recombination": "GA_Recombination_TwoPointCrossover"` indica que se usará la recombinación _crossover_ de dos puntos.
- El campo `"mutation": "GA_Recombation_MultipleSwitches"` indica que la mutación consistirá en cambiar un aminoácido de la secuencia por otro aminoácido elegido de manera aleatoria. Así mismo, `"probability"` y `"numSwitches"`, en la sección `"mutationParams"`, indican respectivamente que la mutación se aplicará con 10% de probabilidad y que solo se mutará 1 posición de la secuencia.
- El campo `"replacement": "GA_Replacement_WorseOut"` indica que, después de producir los hijos, se eliminarán aquellos individuos del conjunto total cuya aptitud es la menor, de tal manera que el tamaño de la población se mantenga constante en cada iteración. 

Una vez escrita la configuración en el archivo `setup.json`, el programa se puede ejecutar ingresando el siguiente comando en la consola: 

```
python -m evodesign example/setup.json
```

**Durante la ejecución del programa, no se despliega ningún mensaje en la consola a menos que haya ocurrido un error.** Sin embargo, uno puede verificar que el programa está corriendo correctamente al inspeccionar el contenido de la carpeta _workspace_—nuevos archivos deberían crearse dentro de dicha carpeta periódicamente.

Al finalizar la ejecución del programa, la carpeta _workspace_ debería contener la siguiente estructura:

```
.
├─ example/
|  ├─ pdbs/
|  |  ├─ prot_0.pdb
|  |  ├─ ...
|  ├─ populations/
|  |  ├─ pop_0.json
|  |  ├─ ...
|  ├─ ~children.bkp
|  ├─ fitness.png
|  ├─ setup.json
|  ├─ statistics.csv
```

- La carpeta `pdbs` contiene los archivos PDB de las estructuras predichas para todas las secuencias producidas durante la ejecución del programa (y que cuya aptitud fue evaluada).
- La carpeta `populations` contiene las poblaciones (esto es, el conjunto de secuencias) producidas en cada iteración del algoritmo evolutivo. En cada uno de estos archivos se registra la aptitud y demás métricas de cada secuencia.
- El archivo `~children.bkp` es un archivo temporal que se usa durante la ejecución del programa. **Este archivo puede ignorarse**.
- El archivo `fitness.png` contiene la gráfica que muestra la aptitud de la población en cada iteración del algoritmo.
- El archivo `setup.json` es el mismo archivo de configuración que se mencionó anteriormente. No se modifica durante la ejecución del algoritmo. 
- El archivo `statistics.csv` contiene datos relevantes de cada iteración del algoritmo. Estos datos son los siguientes:
  - La aptitud mínima, promedio y máxima de la población.
  - El identificador, el valor de aptitudo, y los valores individuales de las demás métricas correspondientes al individuo de mayor aptitud.
  - La secuencia de aminoácidos correspondiente al individuo de mayor aptitud.

<a name="condicion-finalizacion"></a>
### Condición de finalización del algoritmo evolutivo

La condición de finalización está preestablecida y es particular según el algoritmo evolutivo seleccionado en la configuración del programa. Sin embargo, en general el programa finalizará su ejecución al momento que se cumpla al menos una de las siguientes condiciones: 

- se encuentra la solución óptima (esto es, cuando se encuentra un individuo cuya aptitud es el valor máximo posible);
- se alcanzó el número máximo de iteraciones a ejecutar (cantidad que debió ser especificada por el usuario en el archivo de configuración);
- hubo cierta cantidad de iteraciones consecutivas donde la aptitud poblacional no experimentó cambios relevantes (es decir, cuando se ha llegado a un óptimo local). 

<a name="reanudar-ejecucion"></a>
### Reanudar la ejecución del programa a partir de una interrupción anterior

Es posible que en algún punto se vea interrumpida la ejecución del programa (ya sea porque ocurrió un error, porque el API remoto de ESMFold falló en responder a una petición, o porque el usuario interrumpió el programa deliveradamente). En estos casos, es posible reanudar la ejecución del programa, continuando desde la iteración que fue interrumpida. Para ello, simplemente hay que volver a ingresar en la consola el mismo comando mostrado anteriormente: 

```
python -m evodesign example/setup.json
```

<a name="extender-ejecucion"></a>
### Extender la ejecución del programa

Supóngase que el algoritmo utilizado en el ejemplo anterior terminó de ejecutar las 20 iteraciones que fueron especificadas en `example/setup.json`, y supóngase que se desea correr este mismo algoritmo por 10 iteraciones adicionales. Para lograr esto, simplemente debe modificarse el archivo `example/setup.json` para cambiar el campo que originalmente decía `"numIterations": 20` por `"numIterations": 30`. Luego, se vuelve a ejecutar el comando `python -m evodesign example/setup.json`.
