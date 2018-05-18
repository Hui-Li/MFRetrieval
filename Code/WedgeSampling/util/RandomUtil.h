#ifndef RANDOMUTIL_H
#define RANDOMUTIL_H

#include <random>
#include "../alg/Node.h"

namespace DiscreteRandomGenerator{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1.0);

    int DiscreteRandom(Node *node){

        if(!node->set_random_generator) {
            node->random_generator = new AliasSamplingGenerator(node->transistion.size(), node->transistion.data());
            node->set_random_generator = true;
        }

        return node->random_generator->sample(dis(gen), dis(gen));
    }

    int DiscreteRandom(vector<double> &weights){

        AliasSamplingGenerator generator(weights.size(), weights.data());
        return generator.sample(dis(gen), dis(gen));
    }
};
#endif //RANDOMUTIL_H
