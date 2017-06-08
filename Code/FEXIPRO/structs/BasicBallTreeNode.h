#ifndef BASICBALLTREENODE_H
#define BASICBALLTREENODE_H

#include "../util/Base.h"
#include "../util/Conf.h"

class BasicBallTreeNode {
protected:
	double *center;
	vector<int> pointIDs;
	vector<const double*> pointPtrs;
	double constrain;
	bool isLeaf;
	int dimension;
public:

	inline BasicBallTreeNode(){}

	inline BasicBallTreeNode(vector<int> &pointIDs, vector<const double*> &pointPtrs, int dimension) {

		this->pointIDs = pointIDs;
		this->pointPtrs = pointPtrs;
		this->isLeaf = false;
		this->constrain = 0;
		this->dimension = dimension;
		this->center = new double[dimension];

		for (int dimIndex = 0; dimIndex < dimension; dimIndex++) {
			center[dimIndex] = 0;
		}

		for (int i = 0; i < pointPtrs.size(); i++) {
			for (int dimIndex = 0; dimIndex < dimension; dimIndex++) {
				center[dimIndex] += pointPtrs[i][dimIndex];
			}
		}

		for (int dimIndex = 0; dimIndex < dimension; dimIndex++) {
			center[dimIndex] /= pointPtrs.size();
		}
	}

	inline ~BasicBallTreeNode() {
		if (this->center)
			delete[]center;
	}

	virtual void splitNode(vector<int> &leftPointIDs, vector<const double *> &leftPointPtrs,
	vector<int> &rightPointIDs, vector<const double *> &rightPointPtrs) = 0;

	inline virtual BasicBallTreeNode *getLeftNode() = 0;

	inline virtual BasicBallTreeNode *getRightNode() = 0;

	inline int getSize() const {
		return this->pointIDs.size();
	}

	inline int getID(const int index) const {
		return pointIDs[index];
	}

	inline const double *getPtr(const int index) const {
		return pointPtrs[index];
	}

	inline double *getMean() const {
		return center;
	}

	inline void setLeafFlag() {
		this->isLeaf = true;
	}

	inline bool isLeafNode() const {
		return isLeaf;
	}

	inline vector<const double *> &getPointPtrs() {
		return this->pointPtrs;
	}

	inline double getConstrain() const {
		return this->constrain;
	}
};
#endif //BASICBALLTREENODE_H