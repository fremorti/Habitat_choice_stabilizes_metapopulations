'''
Created on Aug 9, 2016

@author: frederik
added incomplete habitat choice
'''

import random as rnd
import numpy as np
import math as math

class Individual:
    '''Class that regulates individuals and their properties'''
    def __init__(self,
                 x,
                 y,
                 muT,
                 varT,
                 maxd,
                 d,
                 settlement,
                 cost_of_disp):
        '''Initialization'''
        self.x = x                  #location (x, y)
        self.y = y
        self.muT=muT                #optimal environment trait
        self.varT=varT              #Niche width,
        self.maxd = maxd            #max dispersal distance
        self.d = d                  #dispersal trait
        self.amax=0.05              #maximum encounter rate, default:0.05
        self.ct=1                   #strength of the generalism trade-off, default:1
        self.h=0.2                  #handling time, default: 0.2
        self.settlement = settlement                    #absence or presence of settlement choice
        self.sigma=300 - settlement*cost_of_disp        #conversion factor, default: 300
        
        
    def move(self, max_x, max_y,  env):
        '''an individual moves within the boundaries of the landscape
        where settlement is or is not dependent on the environment env
        toroidal landscape
        '''
        
        #all possible new x coordinates (g) and all possible y coordinates (h) to disperse to    
        g = [x%len(env) for x in range(self.x-self.maxd, self.x+self.maxd+1)]
        h = [y%len(env[0]) for y in range(self.y-self.maxd, self.y+self.maxd+1)]
        rnd.shuffle(g)
        rnd.shuffle(h)
        
        if self.settlement:
            #with settlement decision
            diff = 1
            #loop through every combination of xs and ys   
            for x in g:
                for y in h:
                    diff_ = abs(self.muT - env[x][y])
                    if diff_ < diff and (self.x != x or self.y != y):
                        '''if the mismatch of the individuals optimal environment trait 
                        with the environmental value of this location is smaller than
                        with the environmental value of the previously proposed location
                        (and it is not the individual's current location)
                        propose a new settlement location'''
                        diff = diff_
                        x_ = x
                        y_ = y
            self.x, self.y = x_, y_
        
        else:
            #without settlement decision
            while self.x == g[0] and self.y == h[0]:
                '''take the combination of the first coordinates as the new settlement location 
                unless it is the individual's current location. Since the coordinates are shuffled, 
                this is a random coordinate'''
                rnd.shuffle(g)
                rnd.shuffle(h)
            self.x, self.y = g[0], h[0]
        

    def mutation(self,rate, md, mv):
        #an individual mutates according to a probability (rate)
        
        #dispersal mutates only if it is mutable (md)
        if rnd.random()<rate and md :
            self.threshold=abs(np.random.normal(self.threshold,0.1))
            self.threshold = 2-self.threshold if self.threshold > 1 else self.threshold
            
        if rnd.random()<rate:
            self.muT=np.random.normal(self.muT,0.1)
        if rnd.random()<rate and mv :
            self.varT=abs(np.random.normal(self.varT,0.1))
            self.varT = 2-self.varT if self.varT > 1 else self.varT
        
                                                                                                       
    def resource_use(self,localhabitat,R):
        #calculates the expected resources an individual will use in a certain location, derived from chaianunporn and hoverstadt 2012
        Gamma=math.exp(-self.ct*self.varT)    #factor due to trade off function 
        Alpha_ij=self.amax*Gamma*math.exp(-(math.pow((self.muT-localhabitat),2)/math.pow((self.varT),2)))  #habitat match and niche width determine how performant the individual is in this environment
        Ri=(Alpha_ij*R/(1+(self.h+self.h*Alpha_ij*R)))  #resource use is determined by how performant the individual is and how many resources there are locally, a
        
        return Ri
        
                  
    def fitness(self,localhabitat,R):        
        #draw a random number of offspring with an average proportionate with the resources used
        return np.random.poisson(self.sigma*self.resource_use(localhabitat, R))
        
   
class Metapopulation:
    '''Contains the whole population, regulates daily affairs'''
    def __init__(self,
                 max_x,
                 max_y,
                 res_R,
                 res_K,
                 initialmaxd,
                 fixedvarT,
                 fixedd,
                 mutable_dispersal,
                 mutable_varT,
                 departure, 
                 settlement,
                 cod):
        '''Initialization'''           
        self.max_x = max_x                                      #number of grid cells along the first dimension of the landscape
        self.max_y = max_y                                      #number of grid cells along the second dimension of the landscape
        self.res_R = res_R                                      #Optimal growth rate of the resources
        self.res_K = res_K                                      #Carrying capacity of the resources
        self.initialmaxd = initialmaxd                          #maximum dispersal distance: for each dimension, maximum difference between the departure and settlement location
        self.fixedvarT = fixedvarT                              #value of the fixed niche width, None in scenarios of fixed dispersal trait
        self.fixedd = fixedd                                    #value of the fixed dispersal trait, None in scenarios of fixed varT
        self.mv = mutable_varT                                  #0 if varT is fixed, 1 if varT evolves
        self.md = mutable_dispersal                             #0 if d is fixed, 1 if d evolves
        self.departure = departure                              #indicates whether all individuals have a departure choice
        self.settlement = settlement                            #indicates whether all individuals have a settlement choice
        self.environment = np.zeros((self.max_x,self.max_y))    #the environmental-value map
        self.resources = np.zeros((self.max_x,self.max_y))      #the local resources map
        self.population = []
        self.localsizes = []                                    #list of population sizes at each location for each generation
        self.cod = cod                                          #cost of settlement choice
        self.initialize_pop()
        
    def initialize_pop(self):
        '''Initialize individuals'''
        startpop = 70000  #initial metapopulation size
        
        for _ in range(startpop):
            x, y, muT = rnd.randint(0,(self.max_x-1)), rnd.randint(0,(self.max_y-1)), rnd.random()
            self.population.append(Individual(x,
                                              y,
                                              muT,
                                              (self.fixedvarT if self.fixedvarT else 0.5*rnd.random()) , 
                                              self.initialmaxd,
                                              (self.fixedd if self.fixedd else (5 if self.departure else 1)*rnd.random()),
                                              self.settlement,
                                              self.cod))

                                             
    def lifecycle(self):   
        '''all actions during one generation for the metapopulation'''
        
        #resources grow
        self.resources += self.res_R*(1-self.resources/self.res_K)

        #replace generation with new one
        oldpop = self.population[:]
        del self.population[:]
        
        #randomize the order in which individuals will perfom their actions
        rnd.shuffle(oldpop)
        
        movenumber, prospectnumber = 0, 0
        oldpopsize = len(oldpop)        #old metapopulation size
        newlocalsizes= np.zeros((self.max_x,self.max_y))
        
        for ind in oldpop:
            
            #mutate
            ind.mutation(0.01, self.md, self.mv)
                     
            #move
            #calculate how much resources the individual needs to reproduce
            necessary_resources=ind.resource_use(self.environment[ind.x,ind.y],self.resources[ind.x, ind.y])
            #decide to move: according to available resources when there is a departure decision, random when there is not
            if (ind.sigma*necessary_resources if self.departure else rnd.random()) < ind.d:
                x_, y_ = ind.x, ind.y 
                ind.move(self.max_x, self.max_y, self.environment)
                prospectnumber += 1
                if not(ind.x == x_ and ind.y == y_):
                    movenumber += 1
                
            #reproduce
            necessary_resources=ind.resource_use(self.environment[ind.x,ind.y],self.resources[ind.x, ind.y])
            #if there are enough resources present locally...
            if necessary_resources<self.resources[ind.x,ind.y]:
                #...deplete resources
                self.resources[ind.x,ind.y]-=necessary_resources   
                #reproduce according to the fitness value
                Fitness=ind.fitness(self.environment[ind.x,ind.y],self.resources[ind.x, ind.y])
                newlocalsizes[ind.x, ind.y] += Fitness
                for _ in range(Fitness):
                    #add a new individual with the same traits as its parent to the new population
                    self.population.append(Individual(ind.x,
                                                      ind.y,
                                                      ind.muT,
                                                      ind.varT,
                                                      ind.maxd,
                                                      ind.d,
                                                      self.settlement,
                                                      self.cod))
    
            else:
                #deplete resources, but no reproduction (fitness dependent on environmnet only, not resources)
                self.resources[ind.x, ind.y] = 0
        
        #resp. the dispersaland prospecting propensity
        self.disp_prop = movenumber/oldpopsize
        self.pros_prop = prospectnumber/oldpopsize
        #calculate local population sizes of this generation
        self.localsizes.append(newlocalsizes)
         
    def loadlandscape(self):
        'Initialize envionmental values in the landscape'
        for x in range(self.max_x):
            for y in range(self.max_y):
                self.environment[x,y]=rnd.random()