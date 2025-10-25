from sage.all import Integer as sage_Integer
from sboxUv2.core import Sb

import sqlite3


class FunctionsDB:
    """This idea of this class is to factor away the interaction with
    any database of S-boxes.

    In particular, it handles the generation of the SELECT queries,
    and the (admitedly small) boilerplate needed by the `with` syntax.

    """



    # !SECTION! Initialization 

    def __init__(self, db_file, row_structure):
        """Initializes a database of S-boxes. It contains a table called "functions" storing S-boxes, and a table called "bibliography" which contains bibliography entries. The exact structure of the rows of the "functions" table is specified by the `row_structure` argument, which must be a dictionary where keys are the identifiers of the rows of the database, and the values are the tinysql type of the corresponding row.

        This class should be thought of as a virtual class since it doesn't provide *all* that is needed. In particular, the specifics of how to parse the content of a row and how to generate one are not provided: this is a job too specific for a general purpose class, and is left to its children.
        
        Args:
            db_file (str): the path to the file that will/already does store the database.
            row_structure (dict): a description of the structure of the rows of the "functions" table.
        """
        self.db_file = db_file
        self.row_structure = row_structure
        if "id" not in row_structure:
            self.row_structure["id"] = "INTEGER"
        self.functions_table = "functions"
        self.bibliography_table = "bibliography"
        # preparing queries
        self.function_insertion_query = "INSERT INTO {} VALUES ({} ?)".format(
            self.functions_table,
            "?, " * (len(row_structure) - 1) # the last question mark
                                             # is already in the
                                             # string above
        )
        self.connection = sqlite3.connect(self.db_file)
        self.cursor = self.connection.cursor()
        # if the file exists, we initialize the length
        try:
            self.cursor.execute("SELECT COUNT(id) FROM {}".format(self.functions_table))
            self.number_of_functions = self.cursor.fetchall()[0][0]
        # otherwise we create the DB file
        except:
            self.create()

        
    def create(self):
        creation_query = "CREATE TABLE IF NOT EXISTS {} (".format(self.functions_table)
        for column in sorted(self.row_structure.keys()):
            creation_query += "{} {},".format(column, self.row_structure[column])
        creation_query = creation_query[:-1] + ")"
        self.cursor.execute(creation_query)
        self.number_of_functions = 0
        
        

    # !SECTION! Handling queries 

    
    def parse_function_from_row(self, row):
        raise Exception("virtual method that shouldn't be called")

    
    def query_functions(self, query_description):
        where_clause = "SELECT * from functions WHERE "
        values = []
        for constraint in query_description.keys():
            if constraint not in self.row_structure.keys():
                raise Exception("unrecognized parameter in query : {} (={})".format(
                    constraint,
                    query_description[constraint]
                    ))
            else:
                where_clause += constraint + " = ? AND "
                # c_type = self.row_structure[constraint]
                # if c_type  == "BLOB":
                #     values.append("x'{}'".format(query_description[constraint].hex()))
                # elif c_type == "INTEGER":
                #     values.append(int(query_description[constraint]))
                # else:
                values.append(query_description[constraint])
        values = tuple(values)
        where_clause = where_clause[:-4]
        self.cursor.execute(where_clause, values)
        result = self.cursor.fetchall()
        if len(result) == 0:
            return []
        else:
            return (
                self.parse_function_from_row(row)
                for row in result
            )
    
    

    # !SECTION! Handling insertions 

    def insert_function(self, entry):
        entry["id"] = self.number_of_functions
        inserted_list = [entry[k] for k in sorted(self.row_structure.keys())]
        try:
            self.cursor.execute(self.function_insertion_query, tuple(inserted_list))
            self.number_of_functions += 1
            return entry["id"]
        except:
            raise Exception("Insertion failed for \n {}\n".format(entry))


    def __len__(self):
        return self.number_of_functions
    
    
    # !SECTION!  handling the `with` syntax
        
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.connection.commit()
        self.connection.close()
