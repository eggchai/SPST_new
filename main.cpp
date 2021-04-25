#include <vector>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <ctime>
#include <ios>
#include <iomanip>
#include <climits>
#include <cassert>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <queue>
#include <math.h>
#include <map>
#include "include/structs.h"
#include "include/cracking/standard_cracking.h"
#include "include/cracking/updates.h"
#include "index/skiplist_with_rank.h"


#pragma clang diagnostic ignored "-Wformat"

int64_t COLUMN_SIZE, NUM_QUERIES, UPDATES_SIZE, UPDATE_FREQUENCY, UPDATE_START_AFTER;
int ALGORITHM, WORKLOAD_PATTERN;
float QUERY_SELECTIVITY, ZIPF_ALPHA = 2.0;
TotalTime query_times;
size_t current_query;

using namespace std;

int64_t zipf(double alpha, int64_t n) {
    static int first = true;      // Static first time flag
    static double c = 0;          // Normalization constant
    double z;                     // Uniform random number (0 < z < 1)
    double sum_prob;              // Sum of probabilities
    double zipf_value = 0.0;      // Computed exponential value to be returned
    int64_t i;                     // Loop counter

    // Compute normalization constant on first call only
    if (first == true) {
        for (i = 1; i <= n; i++)
            c = c + (1.0 / pow((double) i, alpha));
        c = 1.0 / c;
        first = false;
    }

    // Pull a uniform random number (0 < z < 1)
    do {
        z = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    } while ((z == 0) || (z == 1));

    // Map z to the value
    sum_prob = 0;
    for (i = 1; i <= n; i++) {
        sum_prob = sum_prob + c / pow((double) i, alpha);
        if (sum_prob >= z) {
            zipf_value = i;
            break;
        }
    }

    // Assert that zipf_value is between 1 and N
    assert((zipf_value >= 1) && (zipf_value <= n));

    return zipf_value;
}

void generate_skewed_data(vector<int64_t> *data, int64_t maxLeftQueryVal) {
    // the focus should be in the center of the dataset
    int64_t hotspot = COLUMN_SIZE / 2;

    // compute zipf distribution
    typedef map<int64_t, int64_t> result_t;
    typedef result_t::iterator result_iterator_t;

    result_t result;
    for (size_t i = 0; i < NUM_QUERIES; ++i) {
        int64_t nextValue = zipf(ZIPF_ALPHA, COLUMN_SIZE);
        result_iterator_t it = result.find(nextValue);
        if (it != result.end()) {
            ++it->second;
        } else {
            result.insert(make_pair(nextValue, 1));
        }
    }

    int64_t zoneSize = hotspot / result.size();

    int64_t zone = 0;
    for (result_iterator_t it = result.begin(); it != result.end(); ++it) {
        for (int i = 0; i < it->second; ++i) {
            int64_t direction = rand() % 2 == 0 ? 1 : -1;
            int64_t zoneBegin = hotspot + (zone * zoneSize * direction);
            int64_t zoneEnd = zoneBegin + (zoneSize * direction);
            if (direction == -1) {
                int64_t tmp = zoneBegin;
                zoneBegin = zoneEnd;
                zoneEnd = tmp;
            }
            int64_t predicate = rand() % (zoneEnd - zoneBegin + 1) + zoneBegin;
            while (predicate > maxLeftQueryVal) {
                direction = rand() % 2 == 0 ? 1 : -1;
                zoneBegin = hotspot + (zone * zoneSize * direction);
                zoneEnd = zoneBegin + (zoneSize * direction);
                if (direction == -1) {
                    int64_t tmp = zoneBegin;
                    zoneBegin = zoneEnd;
                    zoneEnd = tmp;
                }
                predicate = rand() % (zoneEnd - zoneBegin + 1) + zoneBegin;
            }
            data->push_back(predicate);
        }
        ++zone;
    }
    random_shuffle(data->begin(), data->end());

}

int64_t scanQuery(IndexEntry *c, int64_t from, int64_t to) {
    int64_t sum = 0;
    for (int64_t i = from; i <= to; i++) {
        sum += c[i].m_key;
    }

    return sum;
}

int64_t fullScan(IndexEntry *c, int64_t low, int64_t high) {
    int64_t sum = 0;
    for (int64_t i = 0; i < COLUMN_SIZE; i++) {
        if (c[i].m_key >= low && c[i].m_key < high)
            sum += c[i].m_key;
    }
    return sum;
}

void generateColumn(vector<int64_t> &column) {
    for (int i = 0; i < COLUMN_SIZE; i++) {
        column.push_back(i);
    }
    random_shuffle(column.begin(), column.end());
}

void generateWorkload(vector<int64_t> &leftQuery, vector<int64_t> &rightQuery) {
    int64_t s = COLUMN_SIZE / 100 * QUERY_SELECTIVITY;
    int64_t j = COLUMN_SIZE / NUM_QUERIES;
    int64_t maxLeftQueryVal = COLUMN_SIZE - s;
    if (WORKLOAD_PATTERN == 1) {
        for (int i = 0; i < NUM_QUERIES; i++) {
            leftQuery.push_back(rand() % (COLUMN_SIZE - s));
            rightQuery.push_back(leftQuery[i] + s);
        }
    } else if (WORKLOAD_PATTERN == 2) {
        for (int i = 0; i < NUM_QUERIES; i++) {
            leftQuery.push_back(i * j +1);
            rightQuery.push_back(leftQuery[i] + s);
        }
    } else if (WORKLOAD_PATTERN == 3) {
        vector<int64_t> aux;
        generate_skewed_data(&aux, maxLeftQueryVal);
        for (int i = 0; i < aux.size(); i++) {
            leftQuery.push_back(aux[i]);
            rightQuery.push_back(leftQuery[i] + s);
        }
    } else
        fprintf(stderr, "This workload does not exist");
}

void generateUpdates(vector<int64_t> &column) {
    for (int i = 0; i < NUM_QUERIES/UPDATE_FREQUENCY * UPDATES_SIZE; i++)
        column.push_back(rand() % (COLUMN_SIZE));
}

int64_t binary_search(int64_t *c, int64_t key, int64_t lower, int64_t upper, bool *foundKey) {

    *foundKey = false;
    upper--;
    while (lower <= upper) {
        int middle = (lower + upper) / 2;
        auto middleElement = c[middle];

        if (middleElement < key) {
            lower = middle + 1;
        } else if (middleElement > key) {
            upper = middle - 1;
        } else {
            *foundKey = true;
            return middle;
        }
    }
    return upper;
}


int64_t binary_search_gte(int64_t *c, int64_t key, int64_t start, int64_t end) {
    bool found = false;
    int pos = binary_search(c, key, start, end, &found);
    if (found) {
        while (--pos >= start && c[pos] == key);
    }
    ++pos;
    return pos;
}


void standard_cracking_avl(vector<int64_t> &column, vector<int64_t> &leftQuery, vector<int64_t> &rightQuery,
                       vector<int64_t> &to_do_updates) {
    chrono::time_point<chrono::system_clock> start, end;


    start = chrono::system_clock::now();
    size_t capacity = COLUMN_SIZE;
    IndexEntry *crackercolumn = (IndexEntry *) malloc(COLUMN_SIZE * 2 * sizeof(int64_t));
    //Creating Cracker Column
    for (size_t i = 0; i < COLUMN_SIZE; i++) {
        crackercolumn[i].m_key = column[i];
        crackercolumn[i].m_rowId = i;
    }
    end = chrono::system_clock::now();
//    query_times.swap_time[0] += chrono::duration<double>(end - start).count();
    //Initialitizing Cracker Index
    AvlTree T = NULL;
    vector<int64_t> updates;
    size_t update_index = 0;
    size_t update_count;
    if (UPDATE_FREQUENCY)
        update_count = to_do_updates.size() / UPDATE_FREQUENCY;
    for (current_query = 0; current_query < NUM_QUERIES; current_query++) {
//        if (UPDATE_FREQUENCY){
//            start = chrono::system_clock::now();
//            // Time to fill append list
//            if (((current_query + 1) % UPDATE_FREQUENCY) == 0 && current_query > UPDATE_START_AFTER) {
//                for (size_t j = 0; j < update_count && update_index + j < to_do_updates.size(); j++) {
//                    updates.push_back(to_do_updates[update_index + j]);
//                }
//                update_index = std::min(update_index + update_count, to_do_updates.size());
//            }
//            // Perform MRI
//            if (updates.size()) {
//                int64_t initial_offset = 0;
//                int64_t final_offset = 0;
//                sort(begin(updates), std::end(updates));
//                initial_offset = binary_search_gte(&updates[0], leftQuery[current_query], 0, updates.size() - 1);
//                final_offset = binary_search_gte(&updates[0], rightQuery[current_query], 0, updates.size() - 1);
//                if (initial_offset < final_offset) {
//                    while (final_offset < updates.size() && updates[final_offset] <= rightQuery[current_query]) {
//                        final_offset++;
//                    }
//                    merge_ripple(crackercolumn, capacity, T, updates, initial_offset, final_offset - 1,
//                                 leftQuery[current_query], rightQuery[current_query]);
//                } else if (final_offset == updates.size() - 1 && initial_offset == updates.size() - 1
//                           && updates[final_offset] >= leftQuery[current_query] &&
//                           updates[final_offset] <= rightQuery[current_query])
//                    merge_ripple(crackercolumn, capacity, T, updates, initial_offset, final_offset,
//                                 leftQuery[current_query], rightQuery[current_query]);
//            }
//            end = chrono::system_clock::now();
//            query_times.update_time[current_query] += chrono::duration<double>(end - start).count();
//        }
        T = standardCracking(crackercolumn, COLUMN_SIZE, T, leftQuery[current_query], rightQuery[current_query]);

        //Querying
        start = chrono::system_clock::now();
        IntPair p1 = FindNeighborsGTE(leftQuery[current_query], (AvlTree) T, COLUMN_SIZE - 1);
        IntPair p2 = FindNeighborsLT(rightQuery[current_query], (AvlTree) T, COLUMN_SIZE - 1);
        int offset1 = p1->first;
        int offset2 = p2->second;
        free(p1);
        free(p2);
        end = chrono::system_clock::now();
        query_times.lookup_time[current_query] += chrono::duration<double>(end - start).count();
        start = chrono::system_clock::now();
        int64_t sum = scanQuery(crackercolumn, offset1, offset2);
        end = chrono::system_clock::now();
        query_times.scan_time[current_query] += chrono::duration<double>(end - start).count();
        int64_t scan_sum = 0;

//        int64_t scan_sum = fullScan(crackercolumn, leftQuery[current_query], rightQuery[current_query]);
        if (sum != scan_sum) {
            fprintf(stderr, "", current_query, scan_sum, sum);
        }
    }
    free(crackercolumn);
}

void standard_cracking_spst(vector<int64_t> &column, vector<int64_t> &leftQuery, vector<int64_t> &rightQuery,
                           vector<int64_t> &to_do_updates) {
    chrono::time_point<chrono::system_clock> start, end,start_prune, end_prune;


    start = chrono::system_clock::now();
    size_t capacity = COLUMN_SIZE;
    IndexEntry *crackercolumn = (IndexEntry *) malloc(COLUMN_SIZE * 2 * sizeof(int64_t));
    //Creating Cracker Column
    for (size_t i = 0; i < COLUMN_SIZE; i++) {
        crackercolumn[i].m_key = column[i];
        crackercolumn[i].m_rowId = i;
    }
    end = chrono::system_clock::now();
//    query_times.swap_time[0] += chrono::duration<double>(end - start).count();
    //Initialitizing Cracker Index
    SPSTTree T = NULL;
    vector<int64_t> updates;
    size_t update_index = 0;
    size_t update_count;
    if (UPDATE_FREQUENCY)
        update_count = to_do_updates.size() / UPDATE_FREQUENCY;
    for (current_query = 0; current_query < NUM_QUERIES; current_query++) {
//        if (UPDATE_FREQUENCY){
//            start = chrono::system_clock::now();
//            // Time to fill append list
//            if (((current_query + 1) % UPDATE_FREQUENCY) == 0 && current_query > UPDATE_START_AFTER) {
//                for (size_t j = 0; j < update_count && update_index + j < to_do_updates.size(); j++) {
//                    updates.push_back(to_do_updates[update_index + j]);
//                }
//                update_index = std::min(update_index + update_count, to_do_updates.size());
//            }
//            // Perform MRI
//            if (updates.size()) {
//                int64_t initial_offset = 0;
//                int64_t final_offset = 0;
//                sort(begin(updates), std::end(updates));
//                initial_offset = binary_search_gte(&updates[0], leftQuery[current_query], 0, updates.size() - 1);
//                final_offset = binary_search_gte(&updates[0], rightQuery[current_query], 0, updates.size() - 1);
//                if (initial_offset < final_offset) {
//                    while (final_offset < updates.size() && updates[final_offset] <= rightQuery[current_query]) {
//                        final_offset++;
//                    }
//                    start_prune = chrono::system_clock::now();
//                    T = PruneLeaves(T);
//                    end_prune = chrono::system_clock::now();
//                    query_times.prune_time[current_query] += chrono::duration<double>(end_prune - start_prune).count();
//                    merge_ripple(crackercolumn, capacity, T, updates, initial_offset, final_offset - 1,
//                                 leftQuery[current_query], rightQuery[current_query]);
//                } else if (final_offset == updates.size() - 1 && initial_offset == updates.size() - 1
//                           && updates[final_offset] >= leftQuery[current_query] &&
//                           updates[final_offset] <= rightQuery[current_query]){
//                    start_prune = chrono::system_clock::now();
//                    T = PruneLeaves(T);
//                    end_prune = chrono::system_clock::now();
//                    query_times.prune_time[current_query] += chrono::duration<double>(end_prune - start_prune).count();
//                    merge_ripple(crackercolumn, capacity, T, updates, initial_offset, final_offset,
//                                 leftQuery[current_query], rightQuery[current_query]);
//                }
//
//            }
//            end = chrono::system_clock::now();
//            query_times.update_time[current_query] += chrono::duration<double>(end - start).count() - query_times.prune_time[current_query];
//        }
        T = standardCracking(crackercolumn, COLUMN_SIZE, T, leftQuery[current_query], rightQuery[current_query]);

        //Querying
        start = chrono::system_clock::now();
        T = IntervalSplay(leftQuery[current_query],rightQuery[current_query],T);
        IntPair p1 = FindNeighborsGTE(leftQuery[current_query], (SPSTTree) T, COLUMN_SIZE - 1);
        IntPair p2 = FindNeighborsLT(rightQuery[current_query], (SPSTTree) T, COLUMN_SIZE - 1);
        int offset1 = p1->first;
        int offset2 = p2->second;
        free(p1);
        free(p2);
        end = chrono::system_clock::now();
        query_times.lookup_time[current_query] += chrono::duration<double>(end - start).count();
        start = chrono::system_clock::now();
        int64_t sum = scanQuery(crackercolumn, offset1, offset2);
        end = chrono::system_clock::now();
        query_times.scan_time[current_query] += chrono::duration<double>(end - start).count();
        int64_t scan_sum = 0;
        fprintf(stderr, "", current_query, scan_sum, sum);
//        cout << current_query << "\n";
//        int64_t scan_sum = fullScan(crackercolumn, leftQuery[current_query], rightQuery[current_query]);
        if (sum != scan_sum) {
            fprintf(stderr, "", current_query, scan_sum, sum);
        }
    }
    free(crackercolumn);
}

void standard_cracking_skiplist(vector<int64_t> &column, vector<int64_t> &leftQuery, vector<int64_t> &rightQuery,
                            vector<int64_t> &to_do_updates) {
    chrono::time_point<chrono::system_clock> start, end,start_prune, end_prune;


    start = chrono::system_clock::now();
    size_t capacity = COLUMN_SIZE;
    IndexEntry *crackercolumn = (IndexEntry *) malloc(COLUMN_SIZE * 2 * sizeof(int64_t));
    //Creating Cracker Column
    for (size_t i = 0; i < COLUMN_SIZE; i++) {
        crackercolumn[i].m_key = column[i];
        crackercolumn[i].m_rowId = i;
    }
    end = chrono::system_clock::now();
//    query_times.swap_time[0] += chrono::duration<double>(end - start).count();
    //Initialitizing Cracker Index
    skiplist* T = skiplist_new();
    for (current_query = 0; current_query < NUM_QUERIES; current_query++) {
        T = standardCracking(crackercolumn, COLUMN_SIZE, T, leftQuery[current_query], rightQuery[current_query]);

        //Querying
        start = chrono::system_clock::now();
        IntPair p1 = slFindNeighborsGTE(leftQuery[current_query], (skiplist*) T, COLUMN_SIZE - 1);
        IntPair p2 = slFindNeighborsLT(rightQuery[current_query], (skiplist*) T, COLUMN_SIZE - 1);
        int offset1 = p1->first;
        int offset2 = p2->second;
        free(p1);
        free(p2);
        end = chrono::system_clock::now();
        query_times.lookup_time[current_query] += chrono::duration<double>(end - start).count();
        start = chrono::system_clock::now();
        int64_t sum = scanQuery(crackercolumn, offset1, offset2);
        end = chrono::system_clock::now();
        query_times.scan_time[current_query] += chrono::duration<double>(end - start).count();
        int64_t scan_sum = 0;
        fprintf(stderr, "", current_query, scan_sum, sum);
//        cout << current_query << "\n";
//        int64_t scan_sum = fullScan(crackercolumn, leftQuery[current_query], rightQuery[current_query]);
        if (sum != scan_sum) {
            fprintf(stderr, "", current_query, scan_sum, sum);
        }
    }
    free(crackercolumn);
}

void print_help(int argc, char **argv) {
    fprintf(stderr, "Unrecognized command line option.\n");
    fprintf(stderr, "Usage: %s [args]\n", argv[0]);
    fprintf(stderr, "   --num-queries\n");
    fprintf(stderr, "   --workload-pattern\n");
    fprintf(stderr, "   --column-size\n");
    fprintf(stderr, "   --algorithm\n");
    fprintf(stderr, "   --update-frequency\n");
    fprintf(stderr, "   --update-size\n");
    fprintf(stderr, "   --start-updates\n");
    fprintf(stderr, "   --query-selectivity\n");

}

pair<string, string> split_once(string delimited, char delimiter) {
    auto pos = delimited.find_first_of(delimiter);
    return {delimited.substr(0, pos), delimited.substr(pos + 1)};
}

int main(int argc, char **argv) {
    NUM_QUERIES = 1000;
    COLUMN_SIZE = 10000000;
    ALGORITHM = 1;
    UPDATE_FREQUENCY = 0;
    UPDATES_SIZE = 0;
    UPDATE_START_AFTER = 0;
    QUERY_SELECTIVITY = 1.0;

    int repetition = 1;
    for (int i = 1; i < argc; i++) {
        auto arg = string(argv[i]);
        if (arg.substr(0, 2) != "--") {
            print_help(argc, argv);
            exit(EXIT_FAILURE);
        }
        arg = arg.substr(2);
        auto p = split_once(arg, '=');
        auto &arg_name = p.first;
        auto &arg_value = p.second;
        if (arg_name == "num-queries") {
            NUM_QUERIES = atoi(arg_value.c_str());
        } else if (arg_name == "workload-pattern") {
            WORKLOAD_PATTERN = atoi(arg_value.c_str());
        } else if (arg_name == "column-size") {
            COLUMN_SIZE = atoi(arg_value.c_str());
        } else if (arg_name == "algorithm") {
            ALGORITHM = atoi(arg_value.c_str());
        } else if (arg_name == "update-frequency") {
            UPDATE_FREQUENCY = atoi(arg_value.c_str());
        } else if (arg_name == "update-size") {
            UPDATES_SIZE = atoi(arg_value.c_str());
        } else if (arg_name == "start-updates") {
            UPDATE_START_AFTER = atoi(arg_value.c_str());
        } else if (arg_name == "query-selectivity") {
            QUERY_SELECTIVITY = atof(arg_value.c_str());
        } else {
            print_help(argc, argv);
            exit(EXIT_FAILURE);
        }
    }

    vector<int64_t> column;
    vector<int64_t> leftQueries;
    vector<int64_t> rightQueries;
    vector<int64_t> updates;
    generateColumn(column);
    generateWorkload(leftQueries, rightQueries);
//    generateUpdates(updates);

    query_times.Initialize(NUM_QUERIES);
    vector<double> times(NUM_QUERIES);
    for (size_t i = 0; i < repetition; i++) {
        current_query = 0;
        if(ALGORITHM == 1)
            standard_cracking_avl(column, leftQueries, rightQueries, updates);
        else if (ALGORITHM == 2)
            standard_cracking_spst(column, leftQueries, rightQueries, updates);
        else if(ALGORITHM == 3){
            //显示cracking时间120, insert 0.068, scan103.736, lookup0.071
            //与SPST比较分别是4.15, 0.008, 1.39, 0.035
            //与AVL比较分别是3.18, 0.019, 1.56, 0.019
            standard_cracking_skiplist(column, leftQueries, rightQueries, updates);
        }
        else
            fprintf(stderr, "Algorithm does not exist");

    }
    double swap=0;
    double insert=0;
    double scan=0;
    double lookup=0;
    double prune=0;
    double update=0;
    skiplist *list = skiplist_new();
    for (size_t i = 0; i < NUM_QUERIES; i++) {
        swap += query_times.swap_time[i];
        insert += query_times.index_insert_time[i];
        scan += query_times.scan_time[i];
        lookup += query_times.lookup_time[i];
        prune += query_times.prune_time[i];
        update += query_times.update_time[i];
        //swap_time;index_insert_time;scan_time;lookup_time;prune_time;update_time
//        cout << query_times.swap_time[i] / repetition << ";" << query_times.index_insert_time[i] / repetition << ";" <<
//             query_times.scan_time[i] / repetition << ";" << query_times.lookup_time[i] / repetition << ";"
//             << query_times.prune_time[i] / repetition << ";" << query_times.update_time[i] / repetition << "\n";
    }
    cout <<swap <<"," <<insert <<","<<scan<<","<<lookup<<","<<prune<<","<<update<<","<<endl;
}
