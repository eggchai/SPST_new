#include "../include/cracking/standard_cracking.h"
#include <chrono>

extern TotalTime query_times;
extern size_t current_query;

using namespace std;

IntPair cracking(IndexEntry *&c, IntPair p1, IntPair p2, int lowKey, int highKey){
    IntPair pivot_pair = (IntPair)malloc(sizeof(struct int_pair));
    if (p1->first == p2->first && p1->second == p2->second) {
        pivot_pair = crackInThreeItemWise(c, p1->first, p1->second, lowKey, highKey);
    } else {
        // crack in two
        pivot_pair = (IntPair) malloc(sizeof(struct int_pair));
        pivot_pair->first = crackInTwoItemWise(c, p1->first, p1->second, lowKey);
        pivot_pair->second = crackInTwoItemWise(c, pivot_pair->first, p2->second, highKey);
    }
    return pivot_pair;
}

AvlTree standardCracking(IndexEntry *&c, int64_t dataSize, AvlTree T, int64_t lowKey, int64_t highKey) {
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    IntPair p1, p2;
    p1 = FindNeighborsLT(lowKey, T, dataSize - 1);
    p2 = FindNeighborsLT(highKey, T, dataSize - 1);
    IntPair pivot_pair = NULL;
    end = chrono::system_clock::now();
    query_times.lookup_time[current_query]+= chrono::duration<double>(end - start).count();
    start = chrono::system_clock::now();
//    if (p1->first == p2->first && p1->second == p2->second) {
//        pivot_pair = crackInThreeItemWise(c, p1->first, p1->second, lowKey, highKey);
//    } else {
//        // crack in two
//        pivot_pair = (IntPair) malloc(sizeof(struct int_pair));
//        pivot_pair->first = crackInTwoItemWise(c, p1->first, p1->second, lowKey);
//        pivot_pair->second = crackInTwoItemWise(c, pivot_pair->first, p2->second, highKey);
//    }
    pivot_pair = cracking(c, p1, p2, lowKey, highKey);
    end = chrono::system_clock::now();
    query_times.swap_time[current_query]+= chrono::duration<double>(end - start).count();
    start = chrono::system_clock::now();
    T = Insert(pivot_pair->first, lowKey, T);
    T = Insert(pivot_pair->second, highKey, T);
    end = chrono::system_clock::now();
    query_times.index_insert_time[current_query]+= chrono::duration<double>(end - start).count();
    free(p1);
    free(p2);
    if (pivot_pair) {
        free(pivot_pair);
        pivot_pair = NULL;
    }
    return T;
}

SPSTTree standardCracking(IndexEntry *&c, int64_t dataSize, SPSTTree T, int64_t lowKey, int64_t highKey) {
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    IntPair p1, p2;
    p1 = FindNeighborsLT(lowKey, T, dataSize - 1);
    p2 = FindNeighborsLT(highKey, T, dataSize - 1);
    IntPair pivot_pair = NULL;
    end = chrono::system_clock::now();
    query_times.lookup_time[current_query]+= chrono::duration<double>(end - start).count();
    start = chrono::system_clock::now();
//    if (p1->first == p2->first && p1->second == p2->second) {
//        pivot_pair = crackInThreeItemWise(c, p1->first, p1->second, lowKey, highKey);
//    } else {
//        // crack in two
//        pivot_pair = (IntPair) malloc(sizeof(struct int_pair));
//        pivot_pair->first = crackInTwoItemWise(c, p1->first, p1->second, lowKey);
//        pivot_pair->second = crackInTwoItemWise(c, pivot_pair->first, p2->second, highKey);
//    }
    pivot_pair = cracking(c,p1, p2, lowKey, highKey);
    end = chrono::system_clock::now();
    query_times.swap_time[current_query]+= chrono::duration<double>(end - start).count();
    start = chrono::system_clock::now();
    T = Insert(pivot_pair->first, lowKey, T);
    T = Insert(pivot_pair->second, highKey, T);
    end = chrono::system_clock::now();
    query_times.index_insert_time[current_query]+= chrono::duration<double>(end - start).count();
    free(p1);
    free(p2);
    if (pivot_pair) {
        free(pivot_pair);
        pivot_pair = NULL;
    }
    return T;
}

skiplist* standardCracking(IndexEntry *&c, int64_t dataSize, skiplist* T, int64_t lowKey, int64_t highKey) {
    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();
    IntPair p1, p2;
    //p1，p2可以正确获取
    p1 = slFindNeighborsLT(lowKey, T, dataSize - 1);
    p2 = slFindNeighborsLT(highKey, T, dataSize - 1);
    IntPair pivot_pair = NULL;
    end = chrono::system_clock::now();
    query_times.lookup_time[current_query]+= chrono::duration<double>(end - start).count();
    start = chrono::system_clock::now();
    if (p1->first == p2->first && p1->second == p2->second) {
        //pivot_pair可以正确获取
        pivot_pair = crackInThreeItemWise(c, p1->first, p1->second, lowKey, highKey);
    } else {
        // crack in two
        pivot_pair = (IntPair) malloc(sizeof(struct int_pair));
        pivot_pair->first = crackInTwoItemWise(c, p1->first, p1->second, lowKey);
        pivot_pair->second = crackInTwoItemWise(c, pivot_pair->first, p2->second, highKey);
    }
    end = chrono::system_clock::now();
    query_times.swap_time[current_query]+= chrono::duration<double>(end - start).count();
    start = chrono::system_clock::now();
    //前面已经正确获取了pivot_pair,下面的参数依次为skiplist*, key, value
    T = skiplist_insert(T, lowKey, pivot_pair->first);
    T = skiplist_insert(T, highKey, pivot_pair->second);
    end = chrono::system_clock::now();
    query_times.index_insert_time[current_query]+= chrono::duration<double>(end - start).count();
    free(p1);
    free(p2);
    if (pivot_pair) {
        free(pivot_pair);
        pivot_pair = NULL;
    }
    return T;
}