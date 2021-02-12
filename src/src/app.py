import pandas as pd
import numpy as np
import itertools, yaml, sys, pathlib
import logging, logging.config
import networkx as nx

class triornot():
    def __init__(self) -> None:
        p = pathlib.Path('./src/src/logging.yaml')
        with open(p,'r') as f:
            logconf = yaml.safe_load(f)
        logging.config.dictConfig(logconf)
        self.logger = logging.getLogger('TriOrNot')

        p = pathlib.Path('./src/src/config.yaml')
        with open(p,'r') as f:
            self.cfg = yaml.safe_load(f)
            
        self.G = nx.Graph()
        self.logger.debug(self.cfg)

    def readResFile(self) -> None:
        p = pathlib.Path(self.cfg['resfile'])
        self.df = pd.read_csv(p, sep='\s+', skiprows=8, engine='python',
                 names=['atom','b','x','y','z','f','g'], index_col='atom')

        self.df = self.df.drop(axis=1, columns=['b','f','g'])
        self.df = self.df.drop(axis=0, index=['HKLF','END'])
        self.logger.debug(self.df)

    def extractCell(self) -> None:
        '''
        This is a mess. Would rather keep everything as dataframes
        or series. Meh. Works for now. 
        '''
        p = pathlib.Path(self.cfg['resfile'])
        with open(p,'r') as f:
            string = f.readlines()

        cell = {}
        _cell = string[3].split()
        cell['a'] = float(_cell[2])
        cell['b'] = float(_cell[3])
        cell['c'] = float(_cell[4])
        cell['alpha'] = float(_cell[5])
        cell['beta'] = float(_cell[6])
        cell['gamma'] = float(_cell[7])
        self.cell = cell
        self.logger.debug(cell)

    def makeCombinations(self) -> None:
        '''
        Make list of every possible pair of atoms
        '''
        self.listToCheck = itertools.combinations(self.df.index,2)

    def calc(self,cell:dict,pair:list) -> float:
        '''
        Turn absolute difference of pair of coordinates (x,y,z), and turn it into
        a distance between the two atoms. Works directly in unit cell space.
        '''
        delta = abs(self.df.loc[pair[0]] - self.df.loc[pair[1]])

        d2 = (delta['x'] * cell['a'])**2 + \
             (delta['y'] * cell['b'])**2 + \
             (delta['z'] * cell['c'])**2 + \
             (2 * cell['b']*cell['c']*delta['y']*delta['z']*np.cos(np.radians(cell['alpha']))) + \
             (2 * cell['a']*cell['c']*delta['x']*delta['z']*np.cos(np.radians(cell['beta']))) +  \
             (2 * cell['a']*cell['b']*delta['x']*delta['y']*np.cos(np.radians(cell['gamma'])))

        return d2**(1/2)
    
    def getDistBetweenNodes(self) -> None:
        '''
        Calls calc() on every pair of atoms, and creates a new dataframe. 
        '''
        dist = []
        for i in self.listToCheck:
            dist.append([i[0], i[1], self.calc(self.cell,i)])
        self.dfdist = pd.DataFrame(dist,columns=['atom1','atom2','dist'])
        self.logger.debug(self.dfdist)

    def rejectOutliers(self) -> None:
        '''
        First round of rejections based purely on distance between atoms. 
        '''
        triangle = self.cfg['triangle']
        min, max = triangle['side']-triangle['error'],\
                   triangle['side']+triangle['error']
        _data = self.dfdist

        hits = _data[pd.Series(_data['dist']).between(min,max)]
        self.hits = hits.sort_values(by='dist').reset_index(drop=True)
        self.logger.debug(self.hits)

    def buildGraph(self) -> None:
        '''
        Takes remanining pairs of atoms and builds a graph from their relationships. 

        In a graphical evironment, this one could be plotted to check for triangles
        using matplotlib by calling somethinglike nx.plot(self.G)
        '''
        for i,_ in enumerate(self.hits.index):
            a,b,_ = self.hits.iloc[i].values
            self.G.add_edges_from([(a,b)])
        self.logger.debug(f'added {list(self.G).__len__()} nodes to graph.')

    def pruneGraph(self) -> None:
        '''
        Looks for cyclical (i.e. closed, undirectional parts of the graph) elements,
        and checks if they are triangles (i.e. hat drei Ecken). 

        Note that nx.find_cycle only returns the FIRST cycle it finds, not ALL.
        So, look, report, delete node, then check again for another cycle. 
        '''
        while len(list(self.G)) > 0:
            try:
                cycle = nx.find_cycle(self.G)
            except:
                self.logger.info(f'No more cyclical graphs found; exiting.')
                break
                
            nodes = len(list(cycle))
            if nodes == 3:
                self.logger.info(f'Triangle found: {cycle[0][0],cycle[1][0],cycle[2][0]}!!!' )
            else:
                self.logger.debug(f'Found a polygon with more edges than three. Cool, but do not want.')
            
            self.logger.debug(f'Removing node {cycle[0][0]} to break current cycle.')
            self.G.remove_node(cycle[0][0])


    def report(self):
        '''
        Could be used to clean up results. Not in use.
        '''
        pass

    

tri = triornot()
tri.readResFile()
tri.extractCell()
tri.makeCombinations()
tri.getDistBetweenNodes()
tri.rejectOutliers()
tri.buildGraph()
tri.pruneGraph()
#tri.report()