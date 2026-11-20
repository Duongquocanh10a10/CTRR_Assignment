#include "PathNode.h"

PathNode* createPathNode(const std::string& name, double f, double g, double h, PathNode* next) {
    PathNode* node = new PathNode;
    node->name = name;
    node->f = f;
    node->g = g;
    node->h = h;
    node->next = next;
    return node;
}

void freePath(PathNode* head) {
    while (head != nullptr) {
        PathNode* nextNode = head->next;
        delete head;
        head = nextNode;
    }
}
