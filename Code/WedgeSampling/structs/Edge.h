#ifndef EDGE_H
#define EDGE_H

#include "../util/Base.h"

// forward declare
class Node;

class Edge {

public:
    Node *from;
    Node *to;
    double weight;

    Edge();

    ~Edge();
};

#endif //EDGE_H
