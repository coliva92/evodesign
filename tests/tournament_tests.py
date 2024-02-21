from evodesign.Population import Population
from evodesign.GA.Selection.Tournament import Tournament
from evodesign.Random import Random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# The purpose of this program is to verify if the order of the individuals
# in a population has any impacto on the selection process in a GA through
# binary tournament.





num_experiments = 10000
pop_size = 250
pop = Population.create_random(pop_size, 10)
pop['survivor'] = True
rng = Random.generator()
fitness = pd.DataFrame({ 
  'fitness': rng.random(pop_size),
  'idx': list(range(pop_size)),
  'count': pop_size * [ 0 ]
})
unsorted_pop = pd.concat([ pop, fitness ], axis=1)
sorted_pop = unsorted_pop.sort_values(by='fitness', 
                                      ascending=False, 
                                      inplace=False,
                                      ignore_index=True)
sorted_pop['idx'] = list(range(pop_size))
selection = Tournament(tournamentSize=2,
                       fitnessColumns=[ 'fitness' ],
                       ascendingSort=[ False ],
                       elitism=False)
init_state = rng.bit_generator.state

pops = [ unsorted_pop, sorted_pop ]
for P in pops:
  for k in range(num_experiments):
    pop = selection(P)
    for i in pop['idx']:
      P.at[i, 'count'] += 1
    print(k)

unsorted_pop.sort_values(by='fitness', 
                         ascending=True, 
                         inplace=True, 
                         ignore_index=True)
unsorted_pop.to_csv('tests/unsorted_pop.csv')
sorted_pop.to_csv('tests/sorted_pop.csv')
print((unsorted_pop['sequence_id'] == sorted_pop['sequence_id']).all())

x = np.arange(pop_size)
plt.ioff()
fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(20, 6))
ax[0].bar(x, unsorted_pop['count'], align='center', color='C0', label='Unsorted')
ax[1].bar(x, sorted_pop['count'], align='center', color='C1', label='Sorted')
plt.savefig('tests/tournament_histograms.png')
