#ifndef SCANNER_H
#define SCANNER_H

#include "../util/Base.h"
#include "Metric.h"
#include "Topk.h"

/**
 * Top-K scanner.
 *
 * Scans keys for top-K query, this is the object passed into the LSH query interface.
 */
template<typename ACCESSOR>
class Scanner {
public:
    typedef typename ACCESSOR::Value Value;

    /**
     * Constructor for this class.
     *
     * @param accessor The scanner use accessor to retrieva values from keys.
     * @param metric The distance metric.
     * @param K Value used to reset internal Topk class.
     * @param R Value used to reset internal Topk class.
     */
    Scanner(
            const ACCESSOR &accessor,
            const Metric &metric,
            unsigned K
    ) : accessor_(accessor), metric_(metric), K_(K), cnt_(0) {}

    /**
      * Reset the query, this function should be invoked before each query.
      */
    void reset(Value query) {
        query_ = query;
        accessor_.reset();
        topk_.reset(K_);
        cnt_ = 0;
    }

    /**
     * Number of points scanned for the current query.
     */
    unsigned cnt() const {
        return cnt_;
    }

    /**
     * TopK results.
     */
    const Topk &topk() const {
        return topk_;
    }

    /**
     * TopK results.
     */
    Topk &topk() {
        return topk_;
    }

    /**
     * Update the current query by scanning key, this is normally invoked by the LSH
     * index structure.
     */
    void operator()(unsigned key) {
        if (accessor_.mark(key)) {
            ++cnt_;
            topk_.push(key, metric_.dist(query_, accessor_(key)));
        }
    }

private:
    ACCESSOR accessor_;
    Metric metric_;
    Topk topk_;
    Value query_;
    unsigned K_;
    unsigned cnt_;
};

#endif //SCANNER_H
