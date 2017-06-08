#ifndef BALLTREENODE_H
#define BALLTREENODE_H

#include "BasicBallTreeNode.h"

class BallTreeNode : public BasicBallTreeNode {

protected:
	BallTreeNode *leftNode;
	BallTreeNode *rightNode;

public:
	inline BallTreeNode(){
		this->leftNode = NULL;
		this->rightNode = NULL;
	}

	inline BallTreeNode(vector<int> &pointIDs, vector<const double*> &pointPtrs, int
	dimension) : BasicBallTreeNode(pointIDs, pointPtrs, dimension) {
		this->leftNode = NULL;
		this->rightNode = NULL;
		this->constrain = -1;
		for (int i = 0; i < pointPtrs.size(); i++) {
			double value = 0;

			for (int colIndex = 0; colIndex < Conf::dimension; colIndex++) {
				value += pow(pointPtrs[i][colIndex] - center[colIndex], 2);
			}

			value = sqrt(value);

			if (value > this->constrain) {
				this->constrain = value;
			}
		}
	}

	inline ~BallTreeNode() {

		if (this->leftNode)
			delete this->leftNode;

		if (this->rightNode)
			delete this->rightNode;

	}

	inline virtual void splitNode(vector<int> &leftPointIDs, vector<const double *> &leftPointPtrs,
	                              vector<int> &rightPointIDs, vector<const double *> &rightPointPtrs) {
		leftNode = new BallTreeNode(leftPointIDs, leftPointPtrs, dimension);
		rightNode = new BallTreeNode(rightPointIDs, rightPointPtrs, dimension);
	}

	inline virtual BallTreeNode *getLeftNode() {
		return this->leftNode;
	}

	inline virtual BallTreeNode *getRightNode() {
		return this->rightNode;
	}

};

#endif //BEATLEMP_BALLTREENODE_H