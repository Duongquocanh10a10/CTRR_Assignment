#ifndef PATHNODE_H
#define PATHNODE_H

#include <string>

// Common linked-list node required by all five tasks (exact fields, no extras).
struct PathNode {
    std::string name;
    double f;
    double g;
    double h;
    PathNode* next;
};

// Allocates and initializes a single PathNode.
PathNode* createPathNode(const std::string& name, double f, double g, double h, PathNode* next = nullptr);

// Releases every node in a linked list returned by any pathfinding function.
void freePath(PathNode* head);

#endif // PATHNODE_H
