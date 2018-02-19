'''
Created on 13 dec. 2016

@author: fremorti
'''

from AdaptNGB import Metapopulation as Metapopulation
import numpy as np
import os
import sys
import time
start = time.clock()
default_path = os.getcwd()


#runs model with fixed dispersal trait 
def LH_dispersal(disp, cost = 0):
    runall(disp, None, 0, 1, sys._getframe().f_code.co_name)

#run model with fixed niche width trait    
def LH_varT(disp, cost = 0):
    runall(None, disp, 1, 0, sys._getframe().f_code.co_name)

#run model without fixed traits, both traits evolve
def LH_both(cost = 0):
    runall(None, None, 1, 1, sys._getframe().f_code.co_name)


def runall(initialthreshold, initialvarT, mutable_threshold, mutable_variability, mode):
    meta = run(MAXTIME, dim,R_res,K_res, maxd, initialvarT, initialthreshold, mutable_threshold, mutable_variability, departure, settlement, cost)
    
    #create an array where we calculate 10 metrics of this simulation
    data = np.zeros(10)
    diversity = [ind.muT for ind in meta.population]
    thresholds = [ind.threshold for ind in meta.population]
    nichebr = [ind.varT for ind in meta.population]
    habitatmatch = [abs(ind.muT-meta.environment[ind.x][ind.y]) for ind in meta.population]
    
    #mean optimal environment trait    
    data[0] = sum(diversity)/len(diversity)
    #mean dispersal trait
    data[1] = sum(thresholds)/len(thresholds)
    #mean niche width
    data[2] = sum(nichebr)/len(nichebr)
    #mean habitat mismatch
    data[3] = abs(sum(habitatmatch)/len(habitatmatch))
    #metapopulation size
    data[4] = len(meta.population)
    #dispersal propensity
    data[5] = meta.disp_prop
    #prospection propensity
    data[6] = meta.pros_prop

    
    #temporal standard deviation of local population sizes during the last 5 generations at each location
    localstddev = [[pow(np.var([size[x,y] for size in meta.localsizes[-5:]]), 0.5) for y in range(dim)] for x in range(dim)]
    #temporal mean local population sizes during the last 5 generations at each location
    localmean = [[np.mean([size[x,y] for size in meta.localsizes[-5:]]) for y in range(dim)] for x in range(dim)]
    #temporal variance of total metapopulation size during the last 5 generations
    globalvar = np.var([np.sum(size) for size in meta.localsizes[-5:]])
    #temporal mean of total metapopulation size during the last 5 generations
    globalmean = np.mean([np.sum(size) for size in meta.localsizes[-5:]])
    
    #local population variability (alpha variability)
    data[7] = pow(np.sum(localstddev)/np.sum(localmean), 2)
    #metapopulation variability (gamma variability)
    data[8] = globalvar/pow(globalmean, 2)
    #metapopulation synchrony (beta variability)
    data[9] = pow(np.sum(localstddev), 2)/globalvar
    
    #save metrics
    if not os.path.exists('{}/data/{}/departure{}/settlement{}/{}/{}'.format(default_path, mode, str(departure), str(settlement), str(cost), str(trait))):
        os.makedirs('{}/data/{}/departure{}/settlement{}/{}/{}'.format(default_path, mode, str(departure), str(settlement), str(cost), str(trait)))
    np.save('{}/data/{}/departure{}/settlement{}/{}/{}/rep{}'.format(default_path, mode, str(departure), str(settlement), str(cost), str(trait), str(rep)), data)
    #save final local population sizes
    localsizes = [meta.localsizes[-1][x][y] for x in range(dim) for y in range(dim)]
    np.save('{}/data/{}/departure{}/settlement{}/{}/{}/locrep{}'.format(default_path, mode, str(departure), str(settlement), str(cost), str(trait), str(rep)), localsizes)
    
    
def run(MAXTIME, dim,R_res,K_res, maxd, initialvarT, initialthreshold, mutable_dispersal, mutable_variability, departure, directed, cost):
    '''
    helper function that initiates a metapopulation and let it evolve for a given amount of generations
    dispersal and niche width are either fixed or evolvable, with initial values passed for fixed traits
    all other parameters are passed to the metapopulation object
    '''
    #initialize metapopulation
    meta = Metapopulation(dim,dim,R_res,K_res, maxd, initialvarT, initialthreshold, mutable_dispersal, mutable_variability, departure, directed, cost)
    
    #initialize landscape
    meta.loadlandscape()
    
    #simulate MAXTIME generations (print generation time and metapopulation size for quickly checking during runs)
    for timer in range(MAXTIME): 
        meta.lifecycle()
    print('generation ',timer)
    print("popsize: {}\n".format(len(meta.population)))
    return(meta)


    
'''
PARAMETERS
'''


MAXTIME=500                     #generation time, default:50
dim = 32                        #grid size (dim X dim gridcells), default:32
R_res = 0.25                    #optimal growth resources, default:0.25
K_res = 1                       #carrying capacity resources, default:1
maxd = 2                        #maximum dispersal length, default:2
func = LH_dispersal            #fixed trait
#departure =int(sys.argv[4])     #departure decision, int(sys.argv[5])
#settlement =int(sys.argv[3])    #settlement desicion, int(sys.argv[4]) 
cost = float(sys.argv[1])       #cost of directed dispersal, float(sys.argv[1])  
rep  = int(sys.argv[3])                  #replicate, sys.argv[3]

for departure in [0, 1]:
    for settlement in [0, 1]:
        trait = (1 if func == LH_varT else 5 if departure else 1)*float(sys.argv[2]) 
        func(trait)
print(str(time.clock()))