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

    def extractSG(self) -> None:
        p = pathlib.Path(self.cfg['resfile'])
        with open(p,'r') as f:
            string = f.readlines()

        sg = {}
        _sg = string[3].split()
        sg['a'] = float(_sg[2])
        sg['b'] = float(_sg[3])
        sg['c'] = float(_sg[4])
        sg['alpha'] = float(_sg[5])
        sg['beta'] = float(_sg[6])
        sg['gamma'] = float(_sg[7])
        self.sg = sg
        self.logger.debug(sg)

    def makeCombinations(self) -> None:
        self.listToCheck = itertools.combinations(self.df.index,2)

    def calc(self,sg:dict,pair:list) -> float:
        delta = self.df.loc[pair[0]] - self.df.loc[pair[1]]

        d2 = (delta['x'] * sg['a'])**2 + \
             (delta['y'] * sg['b'])**2 + \
             (delta['z'] * sg['c'])**2 + \
             (2 * sg['b']*sg['c']*delta['y']*delta['z']*np.cos(np.radians(sg['alpha']))) + \
             (2 * sg['a']*sg['c']*delta['x']*delta['z']*np.cos(np.radians(sg['beta']))) +  \
             (2 * sg['a']*sg['b']*delta['x']*delta['y']*np.cos(np.radians(sg['gamma'])))

        return d2**(1/2)
    
    def getDistBetweenNodes(self) -> None:
        dist = []
        for i in self.listToCheck:
            dist.append([i[0], i[1], self.calc(self.sg,i)])
        self.dfdist = pd.DataFrame(dist,columns=['atom1','atom2','dist'])
        self.logger.debug(self.dfdist)

    def rejectOutliers(self) -> None:
        triangle = self.cfg['triangle']
        min, max = triangle['side']-triangle['error'],\
                   triangle['side']+triangle['error']
        _data = self.dfdist

        hits = _data[pd.Series(_data['dist']).between(min,max)]
        self.hits = hits.sort_values(by='dist').reset_index(drop=True)
        self.logger.debug(self.hits)

    def buildGraph(self) -> None:
        for i,_ in enumerate(self.hits.index):
            a,b,_ = self.hits.iloc[i].values
            self.G.add_edges_from([(a,b)])
        self.logger.debug(f'added {list(self.G).__len__()} nodes to graph.')

    def pruneGraph(self) -> None:
        while len(list(self.G)) > 0:
            try:
                cycle = nx.find_cycle(self.G)
            except:
                self.logger.info(f'No more cyclical graphs found; exiting.')
                break
                
            nodes = len(list(cycle))
            if nodes == 3:
                self.logger.info(f'Triangle found: {cycle[0][0],cycle[1][0],cycle[2][0]} !!!' )
            else:
                self.logger.debug(f'Found a polygon with more edges than three. Cool, but do not want.')
            
            self.logger.debug(f'Removing node {cycle[0][0]} to break current cycle.')
            self.G.remove_node(cycle[0][0])


    def report(self):
        pass

    

tri = triornot()
tri.readResFile()
tri.extractSG()
tri.makeCombinations()
tri.getDistBetweenNodes()
tri.rejectOutliers()
tri.buildGraph()
tri.pruneGraph()
#tri.report()