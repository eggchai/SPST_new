#ifndef PROGRESSIVEINDEXING_STANDARD_CRACKING_H
#define PROGRESSIVEINDEXING_STANDARD_CRACKING_H

#include "stdio.h"
#include <cstdlib>
#include "../index/avl_tree.h"
#include "../index/spst.h"
#include "../../index/skiplist_with_rank.h"
#include "../structs.h"
#include "cracking_util.h"

AvlTree standardCracking(IndexEntry *&c, int64_t dataSize, AvlTree T, int64_t lowKey, int64_t highKey);

SPSTTree standardCracking(IndexEntry *&c, int64_t dataSize, SPSTTree T, int64_t lowKey, int64_t highKey);

skiplist* standardCracking(IndexEntry *&c, int64_t dataSize, skiplist* T, int64_t lowkey, int64_t highkey);
#endif //PROGRESSIVEINDEXING_STANDARD_CRACKING_H
