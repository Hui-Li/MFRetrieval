#pragma once

#include "BallTreeNode.h"
#include "../../util/DataUtil.h"

class AngleBallTreeNode : public BallTreeNode {

private:
	AngleBallTreeNode *leftNode;
	AngleBallTreeNode *rightNode;
public:
	AngleBallTreeNode(vector<Matrix1D> &points, const int dimension) {
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
		for (int i = 0; i < points.size(); i++) {
			double cosine = DataUtil::cosineArray(points[i], center);

			if (cosine + 1 > this->radius) {
				this->radius = cosine + 1;
			}
		}

	}

	virtual void splitNode(vector<Matrix1D> &leftPoints, vector<Matrix1D> &rightPoints) {
		leftNode = new AngleBallTreeNode(leftPoints, this->dimension);
		rightNode = new AngleBallTreeNode(rightPoints, this->dimension);
	}

	virtual AngleBallTreeNode &getLeftNode() {
		return *leftNode;
	}

	virtual AngleBallTreeNode &getRightNode() {
		return *rightNode;
	}
};