#!/usr/bin/env sage

"""sboxU is a package providing both all the algorithms relevant to the study of vectorial Boolean functions, in particular those intended to be used as S-boxes within symmetric cryptographic primitives.

A particular focus is put on efficiency: large parts of this library are written in C++, and multi-threading is used as much as possible.

"""



from sboxU.core import *

from sboxU.display import *

from sboxU.algorithms import *

from sboxU.biblio import *
    
from sboxU.statistics import *

from sboxU.ccz import *

from sboxU.apn import *

from sboxU.databases import *

from sboxU.random_objects import *

