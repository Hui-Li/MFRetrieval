#ifndef TOPK_H
#define TOPK_H

#include "../util/Base.h"
#include "MaxHeap.h"

/**
 * Top-K heap.
 *
 * At this point topk should contain the nearest k query keys and distances.
 */
class Topk {
private:
    unsigned K;
    MaxHeap <std::pair<double, unsigned> > heap;
    std::vector<std::pair<double, unsigned>> tops;
public:
    Topk() : K(0) {}

    /**
     * reset K value.
     * @param _K the K value in TopK.
     */
    void reset(int _K) {
        K = _K;
        tops.resize(0);
    }

    /**
     * push a value into the maxHeap.
     * @param key  the key.
     * @param dist the distance.
     */
    void push(unsigned key, double dist) {
        std::pair<double, unsigned> item(dist, key);
        unsigned S = heap.size();
        if (S < K) {
            heap.insert(item);
        } else if (item < heap.findMax()) {
            heap.deleteMax();
            heap.insert(item);
        }
    }

    /**
     * generate TopK.
     */
    void genTopk() {
        unsigned total = heap.size();
        for (unsigned i = 0; i != total; ++i) {
            std::pair<double, unsigned> top;
            heap.deleteMax(top);
            tops.push_back(top);
        }
        std::reverse(tops.begin(), tops.end());
    }

    /**
     * Get the std::vector<std::pair<double, unsigned> > instance which contains the nearest keys and distances.
     */
    const std::vector<std::pair<double, unsigned>> &getTopk() const {
        return tops;
    }

    /**
     * Get the std::vector<std::pair<double, unsigned> > instance which contains the nearest keys and distances.
     */
    std::vector<std::pair<double, unsigned>> &getTopk() {
        return tops;
    }

    /**
     * Calculate the recall vale with another heap.
     * @param  topk another TopK.
     */
    const double recall(const Topk &topk) const {
        std::vector<std::pair<double, unsigned>> tops = getTopk();
        std::vector<std::pair<double, unsigned>> benchTops = topk.getTopk();
        unsigned matched = 0;
        for (std::vector<std::pair<double, unsigned> >::iterator i = tops.begin(); i != tops.end();
             ++i) {
            for (std::vector<std::pair<double, unsigned> >::iterator j = benchTops.begin(); j != benchTops.end();
                 ++j) {
                if (i->second == j->second) {
                    ++matched;
                    break;
                }
            }
        }
        return double(matched + 1) / double(benchTops.size() + 1);
    }

    /**
     * Calculate the precision vale with another heap.
     * @param  topk another TopK.
     */
    const double precision(const Topk &topk) const {
        std::vector<std::pair<double, unsigned>> tops = getTopk();
        std::vector<std::pair<double, unsigned>> benchTops = topk.getTopk();
        unsigned matched = 0;
        for (std::vector<std::pair<double, unsigned> >::iterator i = tops.begin(); i != tops.end();
             ++i) {
            for (std::vector<std::pair<double, unsigned> >::iterator j = benchTops.begin(); j != benchTops.end();
                 ++j) {
                if (i->second == j->second) {
                    ++matched;
                    break;
                }
            }
        }
        return double(matched + 1) / double(tops.size() + 1);
    }
};



#endif //TOPK_H
