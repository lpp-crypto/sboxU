#!/usr/bin/env sage

"""sboxUv2 is a package providing both all the algorithms relevant to the study of vectorial Boolean functions, in particular those intended to be used as S-boxes within symmetric cryptographic primitives.

A particular focus is put on efficiency: large parts of this library are written in C++, and multi-threading is used as much as possible.

"""



from .core import *

from .display import *

from .algorithms import *

from .biblio import *
    
from .statistics import *

from .ccz import *

from .apn import *

from .databases import *

from .random_objects import *
