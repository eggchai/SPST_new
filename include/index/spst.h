#ifndef SPST
#define SPST
#include "../structs.h"

SPSTTree New_Node(int64_t offset,ElementType key);
inline SPSTTree RR_Rotate(SPSTTree k2);
inline SPSTTree LL_Rotate(SPSTTree k2);
SPSTTree Splay(int key, SPSTTree root);
SPSTTree IntervalSplay(int keyL,int keyR, SPSTTree root);
SPSTTree PruneLeaves(SPSTTree node);
SPSTTree Insert(int64_t offset,ElementType key, SPSTTree root);
SPSTTree Delete(ElementType key, SPSTTree root);
SPSTTree Search(ElementType key, SPSTTree root);
void InOrder(SPSTTree root);
void PostOrder(SPSTTree p, int indent);
SPSTTree FindLastPiece(SPSTTree T);
SPSTTree FindNeighborsLTFinal(ElementType X, SPSTTree T);
std::vector<SPSTTree> GetNodesInOrder(SPSTTree T);
SPSTTree FindNeighborsGTFinal(ElementType X, SPSTTree T);
SPSTTree  FindFirstPiece(SPSTTree T);
IntPair FindNeighborsLT(ElementType X, SPSTTree T, ElementType limit);
IntPair  FindNeighborsGTE(ElementType X, SPSTTree T, ElementType limit);
#endif