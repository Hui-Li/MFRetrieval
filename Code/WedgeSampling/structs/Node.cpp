#include "Node.h"
#include "Edge.h"

Node::Node(){
    this->W = 1;
    this->next_level_size = 0;
    set_random_generator = false;
    random_generator = nullptr;
}

void Node::init(const int next_level_size) {
    this->next_level_size = next_level_size;
    this->set_random_generator = false;
    this->W = 1;
    this->random_generator = nullptr;

    if (this->next_level_size > 0) {
        this->to_next_level.resize(next_level_size);
        transistion.resize(next_level_size);

        for (int j = 0; j < next_level_size; j++) {
            to_next_level[j] = new Edge();
        }
    }
}

void Node::clear_generator(){
    delete random_generator;
    random_generator = nullptr;
    set_random_generator = false;
}

Node::~Node(){
    for (Edge* e:to_next_level) {
        delete e;
    }

    delete random_generator;
}