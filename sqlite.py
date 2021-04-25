import sqlite3
db = sqlite3.connect('results.db')
cursor = db.cursor()
cursor.execute('''
    CREATE TABLE algorithms(id INTEGER PRIMARY KEY, name TEXT NOT NULL)
''')
cursor.execute('''
    CREATE TABLE workloads(id INTEGER PRIMARY KEY, name TEXT NOT NULL)
''')
cursor.execute('''
    CREATE TABLE column_distributions(id INTEGER PRIMARY KEY, name TEXT NOT NULL)
''')
cursor.execute('''
    CREATE TABLE experiments(id INTEGER PRIMARY KEY, algorithm_id INTEGER NOT NULL, workload_id INTEGER NOT NULL, column_size INTEGER NOT NULL, query_selectivity REAL NOT NULL, start_updates INTEGER NOT NULL,update_frequency INTEGER NOT NULL,update_size INTEGER NOT NULL,
    FOREIGN KEY (workload_id) REFERENCES workload(id), FOREIGN KEY (algorithm_id) REFERENCES algorithm(id))
''')

cursor.execute('''
    CREATE TABLE queries(experiment_id INTEGER , query_number INTEGER NOT NULL , swap_time REAL NOT NULL,index_insert_time REAL NOT NULL,scan_time REAL NOT NULL,lookup_time REAL NOT NULL, prune_time REAL NOT NULL, update_time REAL NOT NULL, FOREIGN KEY (experiment_id) REFERENCES experiment(id),  PRIMARY KEY(experiment_id,query_number))
''')

algorithms = [('AVL',),('SPST',)]
cursor.executemany(''' INSERT INTO algorithms(name) VALUES(?)''', algorithms)

workload = [('Random',),('Sequential',),('Skew',)]
cursor.executemany(''' INSERT INTO workloads(name) VALUES(?)''', workload)


db.commit()
db.close()