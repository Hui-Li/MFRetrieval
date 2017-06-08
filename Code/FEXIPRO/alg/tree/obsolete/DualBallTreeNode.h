#pragma once

#include "BallTreeNode.h"

class DualBallTreeNode : public BallTreeNode {
private:
	double minIP;
	DualBallTreeNode *leftNode;
	DualBallTreeNode *rightNode;
protected:
	double centerNorm;

public:
	DualBallTreeNode(vector<Matrix1D> &points, const int dimension) : BallTreeNode(points, dimension) {
		this->minIP = DBL_MIN;
		this->centerNorm = 0;
		this->leftNode = NULL;
		this->rightNode = NULL;
		for (int dimIndex = 0; dimIndex < dimension; dimIndex++) {
			this->centerNorm += center[dimIndex] * center[dimIndex];
		}
		this->centerNorm = sqrt(this->centerNorm);
	}

	~DualBallTreeNode() {

		if (this->leftNode)
			delete this->leftNode;

		if (this->rightNode)
			delete this->rightNode;

	}

	double getMinIP() {
		return this->minIP;
	}

	double setMinIP(double minIP) {
		this->minIP = minIP;
	}

	double getCenterNorm() {
		return this->centerNorm;
	}

	virtual DualBallTreeNode &getLeftNode() {
		return *(this->leftNode);
	}

	virtual DualBallTreeNode &getRightNode() {
		return *(this->rightNode);
	}

	virtual void splitNode(vector<Matrix1D> &leftPoints, vector<Matrix1D> &rightPoints) {
		this->leftNode = new DualBallTreeNode(leftPoints, dimension);
		this->rightNode = new DualBallTreeNode(rightPoints, dimension);
	}

};