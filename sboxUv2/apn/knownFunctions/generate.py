from sboxUv2.apn.database import *
from sboxUv2.apn.cython_functions import *


DEFAULT_APN_DB_NAME = "apnDB_{}.db"

def generate_apn_ea_classes_database(
        n,
        db_path=None
):
    """!TODO! write docstring
    """
    # handling inputs
    if n == 6:
        from .reprs6 import ccz_class_representatives
    elif n == 7:
        from .reprs7 import ccz_class_representatives
    if db_path == None:
        db_path = DEFAULT_APN_DB_NAME.format(n)
    # generating the DB
    result = []
    with APNFunctions(db_path) as db:
        for s in ccz_class_representatives:
            result.append(s)
    return result
