from sage.all import *
from sboxUv2 import *
from time import *
import hashlib
from tqdm import *
from ast import *



if __name__ == "__main__":
    # Use the db with the previous functions generated in /scripts/apnDB in 'ccz_compact' mode
    # Run : sboxU_apn_db_generation -n 8 -p apn8.db -mode ccz_compact 
    n = 8
    with APNQuadraticFunctions_ccz_only('apn8.db') as db:
        for k in range(38):
            print("{}th Batch".format(k))
            path = "/home/pg/code/APN/databases/8_bit_3m_apn/8_bit_3m_apn_"+str(k)
            file = open(path,'r')
            lines = file.readlines()
            print('Read {} functions'.format(len(lines)))
            file.close()
            functions = [ literal_eval(l) for l in lines]
            functions_entries = []
            print('Imported {} functions'.format(len(lines)))
            for i in trange(len(functions)):
                temp = {}
                temp["qcr"] = bytearray(quadratic_compact_representation(functions[i]))
                # We hash the mugshot for memory concerns
                mug = apn_ea_mugshot(functions[i])
                h = hashlib.sha256()
                h.update(mug)
                mug = h.digest()
                temp["mugshot"] = mug
                functions_entries.append(temp)
            db.insert_many_functions(functions_entries)
            print('Added {} functions'.format(len(functions_entries)))
        print("DB size : {}".format(db.number_of_functions))

        