import os
import inspect
import sqlite3
from pathlib import Path

SCRIPT_PATH =  os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
os.chdir(SCRIPT_PATH)

# Setting Values For Workload Patterns
RANDOM = "1"
SEQUENTIAL = "2"
SKEWED = "3"

#Cracker Index Types
AVL= "1"
SPST="2"


os.system("python sqlite.py")
db = sqlite3.connect('results.db')
cursor = db.cursor()

def compile():
    print("Compiling")
    os.environ['OPT'] = 'true'
    if os.system('cmake -DCMAKE_BUILD_TYPE=Release && make') != 0:
        print("Make Failed")
        exit()

def run_experiment(column_size,query,selectivity,algorithm,num_queries,update_start,update_frequency,update_size):
    codestr ="./main --num-queries=" + str(num_queries) + " --column-size=" + str(column_size) + \
             " --workload-pattern="+str(query)+ " --update-frequency=" + str(update_frequency) + \
             " --update-size="+str(update_size)+ " --query-selectivity=" + str(selectivity) + \
             " --algorithm="+str(algorithm) + " --start-updates="+str(update_start)

    print(codestr)
    result = os.popen(codestr).read()
    cursor.execute('''INSERT INTO experiments(algorithm_id, workload_id, column_size, query_selectivity,start_updates,update_frequency,update_size)
                  VALUES(:algorithm_id,:workload_id, :column_size, :query_selectivity,:start_updates,:update_frequency,:update_size)''',
                  {'algorithm_id':algorithm, 'workload_id':query, 'column_size':column_size, 'query_selectivity':selectivity, 'start_updates':update_start, 'update_frequency':update_frequency, 'update_size':update_size })
    experiment_id =  cursor.execute('''
               SELECT id FROM experiments where algorithm_id = (?) and workload_id=(?) and column_size=(?) and query_selectivity=(?) and start_updates=(?) and update_frequency = (?) and update_size = (?)
            ''', (algorithm,query,column_size,selectivity,update_start,update_frequency,update_size))
    experiment_id = cursor.fetchone()
    experiment_id = experiment_id[0]
    result = result.split("\n")
    for query_number in range(0, len(result)-1):
        query_result = result[query_number].split(";")
        cursor.execute('''INSERT INTO queries(experiment_id, query_number, swap_time,index_insert_time,scan_time,lookup_time,prune_time,update_time)
              VALUES(:experiment_id,:query_number, :swap_time, :index_insert_time, :scan_time,:lookup_time,:prune_time,:update_time)''',
              {'experiment_id':experiment_id, 'query_number':query_number, 'swap_time':query_result[0], 'index_insert_time':query_result[1], 'scan_time':query_result[2], 'lookup_time': query_result[3], 'prune_time': query_result[4], 'update_time': query_result[5]})

def template_run(update_frequency,update_size):
    compile()
    ALGORITHM_LIST=[AVL,SPST]
    WORKLOAD_LIST = [RANDOM,SEQUENTIAL,SKEWED]
    column_size = 100000000
    selectivity = 1
    num_queries = 10000
    update_start = 0
    for query in WORKLOAD_LIST:
        for algorithm in ALGORITHM_LIST:
            cursor.execute('''
               SELECT id FROM experiments where algorithm_id = (?) and workload_id=(?) and column_size=(?) and query_selectivity=(?) and start_updates=(?) and update_frequency = (?) and update_size = (?)
            ''', (algorithm,query,column_size,selectivity,update_start,update_frequency,update_size))
            experiment_exists = cursor.fetchone()
            if experiment_exists is None:
                run_experiment(column_size,query,selectivity,algorithm,num_queries,update_start,update_frequency,update_size)
            db.commit()

template_run(0,0)
# template_run(1000,1000)
# template_run(10,10)
