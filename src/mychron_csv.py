#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 12:42:47 2025

@author: bflynt
"""

import numpy as np
import pathlib


class MychronCSV:
    
    def __init__(self, 
                 filename: pathlib.PosixPath):
        self.filename  = filename
        self.header    = {}
        self.units     = {}
        self.data      = {}
        self.lap_times = []
        self.load_file_()
        
        
    def load_file_(self):
        """
        Load the CSV data from a MyChron *.csv file

        Returns
        -------
        None.

        """
        # Check for errors
        if not self.filename.exists():
            print(f"ERROR: File does not exist: {self.filename}")
        
        # Use text reader for header & units
        with open(self.filename, "r") as fin:      
            
            # Load Header
            for i, line in enumerate(fin):
                key, value = line.replace('"','').replace('\n','').split(',', 1)
                self.header[key] = value
                if 'Segment Times' in key:
                    header_length = i
                    break
        
            # Load Units
            line = fin.readline()  # Blank line
            line = fin.readline()
            keys = line.replace('"','').replace('\n','').split(',')
            
            line = fin.readline()
            unit = line.replace('"','').replace('\n','').split(',')
    
            for key, value in zip(keys, unit):
                self.units[key] = value
                                
        # Use Numpy CSV reader for raw data
        header_length = header_length + 4
        rawdata = np.loadtxt(self.filename, delimiter=",", skiprows=header_length, dtype=str)        
        for i, key in enumerate(keys):
            self.data[key] = np.array([np.float64(t.replace('"','')) for t in rawdata[:, i]], dtype=np.float64)
            
        return        
   