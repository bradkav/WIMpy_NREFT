# Experiment.py
# Last modified: 07/06/2017

#Most recently modified for compatibility with Xenon1T data

#Summary of changes:
#   Adapted to accept single-entry lists of runs
#   Now reads in list of isotopes

#

# Bradley Kavanagh

import numpy as np
import sys
from scipy.interpolate import interp1d
from scipy.ndimage.filters import median_filter

base_dir = '../data/'

def GetFormFactors(root, A, Z):
    print " Loading Form Factor for (A, Z) = (" + str(int(A)) + ", " + str(int(Z)) + ")..."
    return np.loadtxt(root +'/FormFactors/FormFactors_Z='\
                 + str(int(Z)) + '_A=' + str(int(A)) +'.dat')

class Experiment:
    
    #Initialisation
    def __init__(self, expt_name, run_label):
        
        self.expt_name = expt_name
        self.run_label = run_label
        
        #Read in data
        runs_data = np.atleast_1d(np.genfromtxt(base_dir + expt_name +'/' + expt_name +'_runs.dat',
            names=True, comments='-', missing_values='-1', \
            dtype=None))
        
        #List column headings
        names = runs_data.dtype.names
            
        #Determine which row we want
        current_run = []

        row_num = -1
        Nrows = len(runs_data)
        for i in range(Nrows):
            if (runs_data[i]['Label'] == run_label):
                row_num = i
                break
        #except TypeError:
        #    print "Gets here..."
        #    Nrows = 1
        #    if (runs_data['Label'] == run_label):
        #        current_run = runs_data
 
        if (row_num == -1):
            print " Experiment.py: Run label '" + run_label+ "' not found..."
            sys.exit()

        #NOTE - I might not have to do this, it might be done automatically 
        #by genfromtxt
 
        #Create a dictionary for the different column headings
        self.run_info = dict((x, y) for x, y in zip(names, runs_data[row_num]))
        #You can now access 'exposure', 'E_min', 'E_max', 'N_obs' as you like
        #using run_info['exposure'], etc.

        #Isotopes
        self.Z, self.A, self.frac = np.loadtxt(base_dir + expt_name + '/' + expt_name + "_isotopes.dat", unpack=True)
        self.N_iso = len(self.Z)
        
        #Efficiency
        eff_data = np.loadtxt(base_dir + expt_name + '/' + self.run_info["Efficiency_file"])
        #eff_data[:,1] = median_filter(eff_data[:,1],4)
        self.efficiency = interp1d(eff_data[:,0], np.clip(eff_data[:,1], 0, 1.0), bounds_error=False, fill_value=0.0)
    
        #Load formfactors
        self.FF = [GetFormFactors(base_dir + expt_name, self.A[i], self.Z[i]) for i in range(self.N_iso)]
        #print self.FF[0]
            

    