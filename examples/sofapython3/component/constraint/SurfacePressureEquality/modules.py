import Sofa
import Sofa.Core
import Sofa.Simulation
import SofaRuntime

# 1. Importar las definiciones base de los Solvers (¡El paso que faltaba!)
from Sofa import SofaConstraintSolver

# 2. Cargar los plugins C++ al motor
SofaRuntime.importPlugin("SoftRobots")
SofaRuntime.importPlugin("SoftRobots.Inverse")

# 3. Ahora la inyección no fallará y podremos importar el módulo
import Sofa.SoftRobotsInverse

import os
import numpy as np