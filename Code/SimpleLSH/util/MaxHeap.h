#ifndef MAXHEAP_H
#define MAXHEAP_H

#include "../util/Base.h"

/**
 * Max Heap.
 *
 * This is a max heap for TopK.
 */
template<typename Comparable>
class MaxHeap {
public:
    explicit MaxHeap(int capacity = 100) : array(capacity + 1), currentSize(0) {}

    explicit MaxHeap(const std::vector<Comparable> &items) : array(items.size() + 10), currentSize(items.size()) {
        for (int i = 0; i < items.size(); ++i) {
            array[i + 1] = items[i];
        }
        buildHeap();
    }

    /**
     * check if the heap is empty.
     */
    bool isEmpty() const {
        return currentSize == 0;
    }

    /**
     * get the max value.
     */
    const Comparable &findMax() const {
        if (isEmpty()) {
            std::cout << "UnderflowException() ..." << std::endl;
        }
        return array[1];
    }

    /**
     * inser a value.
     * @param x the insert value.
     */
    void insert(const Comparable &x) {
        if (currentSize == array.size() - 1) {
            array.resize(array.size() * 2);
        }
        int hole = ++currentSize;
        for (; hole > 1 && x > array[hole / 2]; hole /= 2) {
            array[hole] = array[hole / 2];
        }
        array[hole] = x;
    }

    /**
     * delete the max value.
     */
    void deleteMax() {
        if (isEmpty()) {
            std::cout << "UnderflowException() ..." << std::endl;
        }
        array[1] = array[currentSize--];
        percolateDown(1);
    }

    /**
     * delete the max value.
     * @param minItem the max value.
     */
    void deleteMax(Comparable &maxItem) {
        if (isEmpty()) {
            std::cout << "UnderflowException() ..." << std::endl;
        }
        maxItem = array[1];
        array[1] = array[currentSize--];
        percolateDown(1);
    }

    /**
     * make the heap empty.
     */
    void makeEmpty() {
        currentSize = 0;
    }

    /**
     * the size of max heap.
     */
    int size() {
        return currentSize;
    }

private:
    int currentSize;
    std::vector<Comparable> array;

    /**
     * build heap.
     */
    void buildHeap() {
        for (int i = currentSize / 2; i > 0; --i) {
            percolateDown(i);
        }
    }

    /**
     * percolate down the binary heap.
     * @param hole the index.
     */
    void percolateDown(int hole) {
        int child;
        Comparable tmp = array[hole];
        for (; hole * 2 <= currentSize; hole = child) {
            child = hole * 2;
            if (child != currentSize && array[child + 1] > array[child]) {
                child++;
            }
            if (array[child] > tmp) {
                array[hole] = array[child];
            } else {
                break;
            }
        }
        array[hole] = tmp;
    }
};

#endif //MAXHEAP_H
