import pandas as pd
import numpy as np
import itertools, yaml, sys, pathlib
import logging, logging.config

class triornot():
    def __init__(self) -> None:
        p1 = pathlib.Path('./src/src/logging.yaml')
        with open(p1,'r') as f:
            logconf = yaml.safe_load(f)
        logging.config.dictConfig(logconf)
        self.logger = logging.getLogger('simpleExample')

        p2 = pathlib.Path('./src/src/config.yaml')
        with open(p2,'r') as f:
            self.cfg = yaml.safe_load(f)

        self.logger.debug(self.cfg)

tri = triornot()