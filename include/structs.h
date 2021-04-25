#ifndef PROGRESSIVEINDEXING_STRUCTS_H_
#define PROGRESSIVEINDEXING_STRUCTS_H_

#pragma once

#include <vector>
#include <cstdlib>
#include <iostream>

#include <iomanip>

#include <cmath>
#include <queue>
#include <climits>
#include <algorithm>
#include <cstring>
#include "basestructs.h"

struct Column {
    std::vector<int64_t> data;
};

class IndexEntry {
public:


    int64_t m_key;
    int64_t m_rowId;

    IndexEntry(int64_t i)
            : m_key(i), m_rowId(-1) {
    }

    IndexEntry()
            : m_key(-1), m_rowId(-1) {}

    IndexEntry(int64_t key, int64_t rowId)
            : m_key(key), m_rowId(rowId) {}

//Query comparisons
    bool operator>(int64_t &other) const { return m_key > other; }

    bool operator>=(int64_t &other) const { return m_key >= other; }

    bool operator<(int64_t &other) const { return m_key < other; }

    bool operator<=(int64_t &other) const { return m_key <= other; }

    bool operator!=(int64_t &other) const { return m_key != other; }

    bool operator==(int64_t &other) const { return m_key == other; }

    bool operator>(const IndexEntry &other) const { return m_key > other.m_key; }

    bool operator>=(const IndexEntry &other) const { return m_key >= other.m_key; }

    bool operator<(const IndexEntry &other) const { return m_key < other.m_key; }

    bool operator<=(const IndexEntry &other) const { return m_key <= other.m_key; }

    bool operator!=(const IndexEntry &other) const { return m_key != other.m_key; }

};


struct AvlNode;

typedef struct AvlNode *PositionAVL;
typedef struct AvlNode *AvlTree;

struct AvlNode {
    ElementType Element;
    int64_t offset;

    AvlTree Left;
    AvlTree Right;
    int64_t Height;
};

struct SPSTNode;

typedef struct SPSTNode *PositionSPST;
typedef struct SPSTNode *SPSTTree;
struct SPSTNode{
    ElementType Element;
    int64_t offset;
     SPSTTree lchild;
     SPSTTree rchild;
};

struct TotalTime{

    std::vector<double> swap_time;
    std::vector<double> index_insert_time;
    std::vector<double> scan_time;
    std::vector<double> lookup_time;
    std::vector<double> prune_time;
    std::vector<double> update_time;

    void Initialize(size_t query_number) {
        swap_time= std::vector<double>(query_number);
        index_insert_time= std::vector<double>(query_number);
        scan_time= std::vector<double>(query_number);
        lookup_time= std::vector<double>(query_number);
        prune_time= std::vector<double>(query_number);
        update_time= std::vector<double>(query_number);
    };
};


#endif