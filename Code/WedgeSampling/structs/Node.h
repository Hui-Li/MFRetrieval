#ifndef NODE_H
#define NODE_H

#include "../util/Base.h"
#include "../util/AliasSamplingGenerator.h"

// forward declare
class Edge;

class Node {

public:

    double W;
    vector<Edge *> to_next_level;
    int next_level_size;

    Node();
    void init(const int next_level_size);
    ~Node();

};

#endif //NODE_H
