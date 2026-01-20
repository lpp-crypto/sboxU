#!/usr/bin/env sage

"""sboxUv2 is a package providing both all the algorithms relevant to the study of vectorial Boolean functions, in particular those intended to be used as S-boxes within symmetric cryptographic primitives.

A particular focus is put on efficiency: large parts of this library are written in C++, and multi-threading is used as much as possible.

"""



from sboxUv2.core import *

from sboxUv2.display import *

from sboxUv2.algorithms import *

from sboxUv2.biblio import *
    
from sboxUv2.statistics import *

from sboxUv2.ccz import *

from sboxUv2.apn import *

from sboxUv2.databases import *

from sboxUv2.random_objects import *

from sboxUv2.exploration import *
