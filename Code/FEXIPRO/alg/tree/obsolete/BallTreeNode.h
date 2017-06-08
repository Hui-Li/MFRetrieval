#pragma once

#include "../../util/Matrix1D.h"

class BallTreeNode {
protected:
    double *center;
    int dimension;
    vector<Matrix1D> points;
    double radius;
    BallTreeNode *leftNode;
    BallTreeNode *rightNode;
    bool isLeaf;

public:
    BallTreeNode(){

    }

    BallTreeNode(vector<Matrix1D> &points, const int dimension) {
        this->points = points;
        this->dimension = dimension;
        this->isLeaf = false;
        this->center = NULL;
        this->leftNode = NULL;
        this->rightNode = NULL;

        center = new double[dimension];
        for (int dimIndex = 0; dimIndex < dimension; dimIndex++) {
            center[dimIndex] = 0;
        }

        for (auto point:points) {
            for (int dimIndex = 0; dimIndex < dimension; dimIndex++) {
                center[dimIndex] += point[dimIndex];
            }
        }

        for (int dimIndex = 0; dimIndex < dimension; dimIndex++) {
            center[dimIndex] /= points.size();
        }

        this->radius = DBL_MIN;
        for (auto point:points) {
            double distance = DataUtil::euclideanDisArray(point, center);
            if (distance > this->radius) {
                this->radius = distance;
            }
        }
    }

    ~BallTreeNode() {
        if (this->center)
            delete[]center;

        if (this->leftNode)
            delete this->leftNode;

        if (this->rightNode)
            delete this->rightNode;

    }

    int getDim() {
        return dimension;
    }

    int getSize() {
        return points.size();
    }
//
//    Matrix1D &getPoint(const int index) {
//        return points[index];
//    }

    const Matrix1D &operator[](const int &index) {
        return points[index];
    }

    virtual void splitNode(vector<Matrix1D> &leftPoints, vector<Matrix1D> &rightPoints) {
        leftNode = new BallTreeNode(leftPoints, this->dimension);
        rightNode = new BallTreeNode(rightPoints, this->dimension);
    }

    virtual BallTreeNode &getLeftNode() {
        return *leftNode;
    }

    virtual BallTreeNode &getRightNode() {
        return *rightNode;
    }

    double *getMean() {
        return center;
    }

    void setLeafFlag() {
        this->isLeaf = true;
    }

    bool isLeafNode() {
        return isLeaf;
    }

    vector<Matrix1D> &getPoints() {
        return this->points;
    }

    double getRadius() {
        return this->radius;
    }
};



