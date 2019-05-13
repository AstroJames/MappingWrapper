import numpy as np
import os

# Density Experiments
os.system("python AutoMapping.py -numOfModels 150 -initDens 1e0  -finalDens 1e15 -initTemp 5e3")
os.system("mv slab* ./T5e3")
os.system("python AutoMapping.py -numOfModels 150 -initDens 1e0  -finalDens 1e15 -initTemp 1e4")
os.system("mv slab* ./T1e4")
os.system("python AutoMapping.py -numOfModels 150 -initDens 1e0  -finalDens 1e15 -initTemp 2e4")
os.system("mv slab* ./T20e3")

# Temperature Experiments
#os.system("python AutoMapping.py -numOfModels 150 -initDens 1e2 -initTemp 5e3 -finalTemp 25e3")
#os.system("mv slab* ./N1e2")
#os.system("python AutoMapping.py -numOfModels 150 -initDens 1e4 -initTemp 5e3 -finalTemp 25e3")
#os.system("mv slab* ./N1e4")
