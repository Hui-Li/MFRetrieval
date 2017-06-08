#pragma once

#include "DualBallTreeNode.h"
#include "../../util/DataUtil.h"

class ConeTreeNode : public DualBallTreeNode {
private:
	double cosineW_q;
	double sineW_q;
	ConeTreeNode *leftNode;
	ConeTreeNode *rightNode;


public:
	ConeTreeNode(vector<Matrix1D> &points, const int dimension) : DualBallTreeNode(points, dimension) {
		this->leftNode = NULL;
		this->rightNode = NULL;
		this->cosineW_q = DBL_MAX;

		for (int i = 0; i < points.size(); i++) {
			double cosine = DataUtil::cosineArray(points[i], center);

			if (cosine < this->cosineW_q) {
				this->cosineW_q = cosine;
			}
		}

		this->sineW_q = DataUtil::calSine(this->cosineW_q);
	}

	~ConeTreeNode() {
		if (this->leftNode)
			delete this->leftNode;

		if (this->rightNode)
			delete this->rightNode;
	}

	double getCosineW_q() {
		return this->cosineW_q;
	}

	double getSineW_q() {
		return this->sineW_q;
	}

	virtual ConeTreeNode &getLeftNode() {
		return *(this->leftNode);
	}

	virtual ConeTreeNode &getRightNode() {
		return *(this->rightNode);
	}

	virtual void splitNode(vector<Matrix1D> &leftPoints, vector<Matrix1D> &rightPoints) {
		this->leftNode = new ConeTreeNode(leftPoints, dimension);
		this->rightNode = new ConeTreeNode(rightPoints, dimension);
	}
};