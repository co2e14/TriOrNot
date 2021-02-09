import pandas as pd
import numpy as np
import itertools, yaml, sys, pathlib
import logging, logging.config

class triornot():
    def __init__(self) -> None:
        p = pathlib.Path('./src/src/logging.yaml')
        with open(p,'r') as f:
            logconf = yaml.safe_load(f)
        logging.config.dictConfig(logconf)
        self.logger = logging.getLogger('simpleExample')

        p = pathlib.Path('./src/src/config.yaml')
        with open(p,'r') as f:
            self.cfg = yaml.safe_load(f)

        self.logger.debug(self.cfg)

    def readResFile(self):
        p = pathlib.Path(self.cfg['resfile'])
        self.df = pd.read_csv(p, sep='\s+', skiprows=8, engine='python',
                 names=['atom','b','x','y','z','f','g'], index_col='atom')

        self.df = self.df.drop(axis=1, columns=['b','f','g'])
        self.df = self.df.drop(axis=0, index=['HKLF','END'])
        self.logger.debug(self.df)

    def extractSG(self):
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

    def makeCombinations(self):
        self.listToCheck = itertools.combinations(self.df.index,2)

    def calc(self,sg,pair):
        delta = self.df.loc[pair[0]] - self.df.loc[pair[1]]

        d2 = (delta['x'] * sg['a'])**2 + \
             (delta['y'] * sg['b'])**2 + \
             (delta['z'] * sg['c'])**2 + \
             (2 * sg['b']*sg['c']*delta['y']*delta['z']*np.cos(np.radians(sg['alpha']))) + \
             (2 * sg['a']*sg['c']*delta['x']*delta['z']*np.cos(np.radians(sg['beta']))) +  \
             (2 * sg['a']*sg['b']*delta['x']*delta['y']*np.cos(np.radians(sg['gamma'])))

        return d2**(1/2)
    
    def findTriangles(self):
        dist = []
        for i in self.listToCheck:
            dist.append([i[0], i[1], self.calc(self.sg,i)])
        self.dfdist = pd.DataFrame(dist,columns=['atom1','atom2','dist'])
        self.logger.debug(self.dfdist)

    def rejectNonTriangles(self):
        triangle = self.cfg['triangle']
        min, max = triangle['side']-triangle['error'],\
                   triangle['side']+triangle['error']
        _data = self.dfdist

        hits = _data[pd.Series(_data['dist']).between(min,max)]
        self.hits = hits.sort_values(by='dist').reset_index(drop=True)
        self.logger.debug(self.hits)

    def report(self):
        pass

    

tri = triornot()
tri.readResFile()
tri.extractSG()
tri.makeCombinations()
tri.findTriangles()
tri.rejectNonTriangles()
#tri.report()