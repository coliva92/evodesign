from evodesign.Fitness.Rmsd import Rmsd
from evodesign.Fitness.Gdt import Gdt
from evodesign.Fitness.Rastrigin import Rastrigin
from evodesign.Fitness.Experimental.Cyclization import Cyclization
from evodesign.Fitness.Experimental.SideChainPacking import SideChainPacking
from evodesign.Workspace import Workspace
import math

Workspace('foo', 'bar.pdb') # SideChainPacking requires a valid Workspace instance
rmsd = Rmsd()
gdt = Gdt()
rastrigin = Rastrigin()
cyclization = Cyclization()
sc_packing = SideChainPacking()

params = rmsd.settings()
test1_passed = type(params) == dict and \
               'Fitness.Rmsd' in params and \
               type(params['Fitness.Rmsd']) == dict and \
               'upperBound' in params['Fitness.Rmsd'] and \
               type(params['Fitness.Rmsd']['upperBound']) == float and \
               params['Fitness.Rmsd']['upperBound'] == -0.5

params = gdt.settings()
test2_passed = type(params) == dict and \
               'Fitness.Gdt' in params and \
               type(params['Fitness.Gdt']) == dict and \
               'upperBound' in params['Fitness.Gdt'] and \
               type(params['Fitness.Gdt']['upperBound']) == float and \
               params['Fitness.Gdt']['upperBound'] == 0.95 and \
               'cutoffs' in params['Fitness.Gdt'] and \
               type(params['Fitness.Gdt']['cutoffs']) == list and \
               type(params['Fitness.Gdt']['cutoffs'][0]) == float and \
               params['Fitness.Gdt']['cutoffs'][0] == 1.0 and \
               params['Fitness.Gdt']['cutoffs'][1] == 2.0 and \
               params['Fitness.Gdt']['cutoffs'][2] == 4.0 and \
               params['Fitness.Gdt']['cutoffs'][3] == 8.0

params = rastrigin.settings()
test3_passed = type(params) == dict and \
               'Fitness.Rastrigin' in params and \
               type(params['Fitness.Rastrigin']) == dict and \
               'upperBound' in params['Fitness.Rastrigin'] and \
               type(params['Fitness.Rastrigin']['upperBound']) == float and \
               params['Fitness.Rastrigin']['upperBound'] == 0.0 and \
               'windowWidth' in params['Fitness.Rastrigin'] and \
               type(params['Fitness.Rastrigin']['windowWidth']) == int and \
               params['Fitness.Rastrigin']['windowWidth'] == 3

params = cyclization.settings()
test4_passed = type(params) == dict and \
               'Fitness.Experimental.Cyclization' in params and \
               type(params['Fitness.Experimental.Cyclization']) == dict and \
               'upperBound' in params['Fitness.Experimental.Cyclization'] and \
               type(params['Fitness.Experimental.Cyclization']['upperBound']) == float and \
               params['Fitness.Experimental.Cyclization']['upperBound'] == -1.4

params = sc_packing.settings()
test5_passed = type(params) == dict and \
               'Fitness.Experimental.SideChainPacking' in params and \
               type(params['Fitness.Experimental.SideChainPacking']) == dict and \
               'upperBound' in params['Fitness.Experimental.SideChainPacking'] and \
               type(params['Fitness.Experimental.SideChainPacking']['upperBound']) == float and \
               params['Fitness.Experimental.SideChainPacking']['upperBound'] == -math.inf and \
               'scwrlExecutablePath' in params['Fitness.Experimental.SideChainPacking'] and \
               type(params['Fitness.Experimental.SideChainPacking']['scwrlExecutablePath']) == str and \
               params['Fitness.Experimental.SideChainPacking']['scwrlExecutablePath'] == './scwrl4/Scwrl4'

all_passed = test1_passed and \
             test2_passed and \
             test3_passed and \
             test4_passed and \
             test5_passed
print('PASSED' if all_passed else 'FAILED')
