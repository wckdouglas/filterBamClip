import numpy as np
import pyximport
import pysam
pyximport.install(setup_args = {'include_dirs':[np.get_include() ]+pysam.get_include()})

from .bam_splitter import *
