#ifndef PROGRESSIVEINDEXING_CRACKING_UPDATES_H
#define PROGRESSIVEINDEXING_CRACKING_UPDATES_H
#include "stdio.h"
#include <cstdlib>
#include "../index/avl_tree.h"
#include "../index/spst.h"
#include "../structs.h"
#include "cracking_util.h"
#include <vector>

using namespace std;


void merge(IndexEntry *&column, size_t &capacity, AvlTree T, vector<int64_t>  &updates, int64_t posL, int64_t posH, int64_t _next = -1);
void merge_ripple(IndexEntry *&column,size_t &capacity, AvlTree T, vector<int64_t>  &updates, int64_t posL, int64_t posH, int64_t low, int64_t high);
void merge(IndexEntry *&column, size_t &capacity, SPSTTree T, vector<int64_t>  &updates, int64_t posL, int64_t posH, int64_t _next = -1);
void merge_ripple(IndexEntry *&column,size_t &capacity, SPSTTree T, vector<int64_t>  &updates, int64_t posL, int64_t posH, int64_t low, int64_t high);
#endif //PROGRESSIVEINDEXING_CRACKING_UPDATES_H
