#!/usr/bin/env python
import numpy as np
import sys
import os
CURRENT_PATH = os.getcwd()

#   ----    initialize scripts     ----    #
sys.path.append(CURRENT_PATH+'/src/md/build_file/')
from buildMD import *

lmp_bin_path 	= '~/Program/bin/'	# need to be provided
workDir			= CURRENT_PATH + '/sim/'; create_folder(workDir)
dataFolder		= CURRENT_PATH + '/src/md/init_structure/'
paramsFolder	= CURRENT_PATH + '/src/md/energy_file/'
simFolder		= CURRENT_PATH + '/src/md/sim_file/'

runId = 0
prepare_simulation(lmp_bin_path,workDir,runId,dataFolder,paramsFolder,simFolder)
run_simulation(workDir,runId)