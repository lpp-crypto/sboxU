# Time-stamp: <2021-08-12 16:45:37 lperrin>

import os
HOME = os.path.expanduser('~')

if os.path.exists(HOME + '/.fftw'):
    from .cpp_diff_lin import *
else:
    from .cpp_diff_lin_no_fp_lat import *

from .cpp_utils import *
from .cpp_equiv import *
from .cpp_equiv_approx import *
from .cpp_ccz import *
