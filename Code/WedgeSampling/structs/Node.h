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
    vector<double> transistion;
    int next_level_size;
    bool set_random_generator;
    AliasSamplingGenerator *random_generator = nullptr;

    Node();
    void init(const int next_level_size);
    void clear_generator();
    ~Node();

};

#endif //NODE_H
