#include "Node.h"
#include "Edge.h"

Node::Node(){
    this->W = 1;
    this->next_level_size = 0;
}

void Node::init(const int next_level_size) {
    this->W = 1;
    this->next_level_size = next_level_size;

    if (this->next_level_size > 0) {
        this->to_next_level.resize(next_level_size);
    }
}

Node::~Node(){
    for (Edge* e:to_next_level) {
        delete e;
    }

}