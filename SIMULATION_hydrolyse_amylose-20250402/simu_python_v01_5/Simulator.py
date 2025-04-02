# Framework for simulation of molecular diffusion - INSA
# Catherine Pothier and Christophe Rigotti

import numpy as np

from Molecule import Molecule


class Simulator:
    def __init__(self, nbEnzymes, nbPolymers,
                nbIterations, xMax, yMax):
        
        # creation of the attributes
        
        self.nbEnzymes = nbEnzymes
        self.nbPolymers = nbPolymers
        self.xMax = xMax
        self.yMax = yMax
        self.nbIterations = nbIterations #?
        
        self.nbInternalObstacleTrials = 200 #nombre d'obstacle 
        
        self.rng = np.random.default_rng(seed=42) #la racine de test aléatoire est initialisée à 42
        
        self.iteration = 0 #?
        self.enzymes = [] #?
        self.polymers = [] #?
        
        self.grid = np.zeros((xMax+1, yMax+1), dtype=object) #on peut gérer la taille de grille ici, elle est remplie de zéro type objet au début 
        self.gridBackground = np.zeros((xMax+1, yMax+1), dtype=int) #on peut gérer la taille de grille ici, elle est remplie de zéro type entier au début, celle dont on se sert en grille de contexte  
        self.universeViews = np.zeros((nbIterations, xMax+1, yMax+1), dtype=int)
        
        self.intialization()
        
        
    def intialization(self):
                
        # initializing both grids        
        for x in range(self.xMax+1):
            for y in range(self.yMax+1):
                self.grid[x, y] = None;
                
                if((x<2) or (x>self.xMax-2) or
                   (y<2) or (y>self.yMax-2)):
                    self.gridBackground[x, y] = -1; # an obstacle, déclaration de la frontière de la grille pour que les molécules soient contenues 
                else:
                    self.gridBackground[x, y] = -2; # solvent, zone de travaille 

        # random trials to add obstacles
        for i in range(self.nbInternalObstacleTrials):
            x = self.rng.integers(low=0, high=self.xMax-4) + 2
            y = self.rng.integers(low=0, high=self.yMax-4) + 2
            self.gridBackground[x, y] = -1; # an obstacle
                
        # enzyme creation
        for i in range(self.nbEnzymes):
            freeCell = False;
            while not(freeCell):
                # Draw a location in a region in the center.
                # Retry if the location already contains
                # a molecule or an obstacle.
                x = self.rng.integers(low=self.xMax//4,
                                       high=(self.xMax*3)//4)
                y = self.rng.integers(low=self.yMax//4,
                                       high=(self.yMax*3)//4)

                if ((self.grid[x, y] == None) and
                    (self.gridBackground[x, y] == -2)):
                    freeCell = True
                
            molId = Molecule(x, y, species="E", size=1)
            self.enzymes.append(molId)
            self.grid[x, y] = molId;
        
        # polymer creation
        for i in range(self.nbPolymers):
            freeCell = False;
            while not(freeCell):
                # Draw a location in a region in the center.
                # Retry if the location already contains
                # a molecule or an obstacle.
                x = self.rng.integers(low=self.xMax//4,
                                       high=(self.xMax*3)//4)
                y = self.rng.integers(low=self.yMax//4,
                                       high=(self.yMax*3)//4)

                if ((self.grid[x, y] == None) and
                    (self.gridBackground[x, y] == -2)):
                    freeCell = True
                
            molId = Molecule(x, y, species="P", size=1)
            self.polymers.append(molId)
            self.grid[x, y] = molId;
            
            
    def moveOneMolecule(self, molId):
        dx = self.rng.integers(low=-1, high=2)
        dy = self.rng.integers(low=-1, high=2)
        x = molId.x;
        y = molId.y;
        if ((self.grid[x+dx, y+dy] == None)
            and (self.gridBackground[x+dx, y+dy] == -2)):
            self.grid[x, y] = None;
            self.grid[x+dx, y+dy] = molId;
            molId.x = x+dx;
            molId.y = y+dy;
        
        
    def execOneIteration(self):
        
        # moving each enzyme
        for molId in self.enzymes:
            self.moveOneMolecule(molId)
        
        # moving each polymer
        for molId in self.polymers:
            self.moveOneMolecule(molId)
            
          
    def computeUniverseView(self):
        universeView = self.gridBackground.copy()
        universeView[0, 0] = -3 # to indicate the (0,0) corner
        for molId in self.enzymes:
            universeView[molId.x, molId.y] = 1
            
        for molId in self.polymers:
            universeView[molId.x, molId.y] = 2
            
        return(universeView)
    
    
    def run(self):
        iteration = 1
        while iteration <= self.nbIterations:
            self.execOneIteration()
            
            self.universeViews[iteration-1, :, :]  = self.computeUniverseView()
            
            print(f"Iteration {iteration} done")
            iteration += 1
            
        return(self.universeViews)

		
