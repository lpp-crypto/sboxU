from sage.all import Integer as sage_Integer
import sqlite3
import json


from sboxU.core import get_sbox



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
        self.journal_table = "journal"
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
        self.new_db = False
        try:
            self.cursor.execute("SELECT COUNT(id) FROM {}".format(self.functions_table))
            self.number_of_functions = self.cursor.fetchall()[0][0]
        # otherwise we create the DB file
        except Exception:
            self.create()

        # CREATE TABLE IF NOT EXISTS for compatibility with db without the journal
        self.cursor.execute(
            "CREATE TABLE IF NOT EXISTS {} (id INTEGER, timestamp TEXT, operation TEXT, comment TEXT)".format(
                self.journal_table
            )
        )
        self.cursor.execute("SELECT COUNT(id) FROM {}".format(self.journal_table))
        self.number_of_journal_entries = self.cursor.fetchone()[0]


    def create(self):
        creation_query = "CREATE TABLE IF NOT EXISTS {} (".format(self.functions_table)
        for column in sorted(self.row_structure.keys()):
            creation_query += "{} {},".format(column, self.row_structure[column])
        creation_query = creation_query[:-1] + ")"
        self.cursor.execute(creation_query)
        self.number_of_functions = 0
        self.new_db = True



    # !SECTION! Journal

    def log_journal(self, operation, comment=None):
        """Appends an entry to this database's journal, unless `comment` is None.

        Args:
            operation: short string naming the operation performed (e.g. "update_database",
                       "insert_quadratic", "add_quadratic_ccz_class").
            comment:   free-text note chosen by the caller describing the change, or None.

        Returns:
            The row id of the inserted journal entry, or None if `comment` is None.
        """
        if comment is None:
            return None
        entry_id = self.number_of_journal_entries
        self.cursor.execute(
            "INSERT INTO {} VALUES (?, datetime('now'), ?, ?)".format(self.journal_table),
            (entry_id, operation, comment)
        )
        self.number_of_journal_entries += 1
        return entry_id


    def get_journal(self):
        """Returns the full journal as a list of dicts, ordered by id."""
        self.cursor.execute(
            "SELECT id, timestamp, operation, comment FROM {} ORDER BY id".format(self.journal_table)
        )
        return [
            {"id": row[0], "timestamp": row[1], "operation": row[2], "comment": row[3]}
            for row in self.cursor.fetchall()
        ]


    # !SECTION! Handling queries

    
    def parse_function_from_row(self, row):
        raise Exception("virtual method that shouldn't be called")

    
    def query_functions(self, query_description):
        """Queries the database using a dictionary of selectors and returns the result as a list.

        Be careful with the queries depending on the database: **.query_functions({}) returns the entire database and can flood your RAM    
        !TODO! finish the docstring
        
        """
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
                # the following is needed to ensure an integer is indeed a python-style `int` and not a SAGE-style `Integer`
                if self.row_structure[constraint] == "INTEGER":
                    values.append(int(query_description[constraint]))
                else:
                    values.append(query_description[constraint])
        values = tuple(values)
        where_clause = where_clause[:-4]
        self.cursor.execute(where_clause, values)
        result = self.cursor.fetchall()
        if len(result) == 0:
            return []
        else:
            return [
                self.parse_function_from_row(row)
                for row in result
            ]

        
    def __getitem__(self, index):
        if not isinstance(index, (int, sage_Integer)):
            raise Exception("db[ index ] can only work if `index` is an integer")
        else:
            counter = 0
            for entry in self.query_functions({"id" : index}):
                if counter > 0:
                    raise Exception("ERROR: the database contains multiple entries with the same identifier ({})".format(index))
                else:
                    return entry
    

    # !SECTION! Handling insertions 

    def insert_function(self, entry):
        entry["id"] = self.number_of_functions
        inserted_list = [entry[k] for k in sorted(self.row_structure.keys())]
        try:
            self.cursor.execute(self.function_insertion_query, tuple(inserted_list))
            self.number_of_functions += 1
            return entry["id"]
        except Exception as e:
            raise Exception("Insertion failed for \n {}\n".format(entry)) from e


    def batch_insert_function(self, entries):
        """Batch analog of `insert_function`: inserts a list of entries with a
        single `executemany` call instead of one `execute` per entry.

        Args:
            entries: a list of dicts, one per row, using the same column keys
                     as `insert_function` ("id" is added/overwritten automatically).

        Returns:
            A list of row ids, one per element of `entries`, in the same order.
        """
        start_id = self.number_of_functions
        columns = sorted(self.row_structure.keys())
        ids = list(range(start_id, start_id + len(entries)))
        rows = []
        for entry, row_id in zip(entries, ids):
            entry["id"] = row_id
            rows.append(tuple(entry[k] for k in columns))
        try:
            self.cursor.executemany(self.function_insertion_query, rows)
            self.number_of_functions = start_id + len(entries)
            return ids
        except Exception as e:
            raise Exception("Batch insertion failed for {} entries\n".format(len(entries))) from e




    def __len__(self):
        return self.number_of_functions
    
    
    # !SECTION!  handling the `with` syntax
        
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.connection.commit()
        self.connection.close()

