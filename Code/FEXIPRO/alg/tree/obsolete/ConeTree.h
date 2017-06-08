#pragma once

#include "ConeTreeNode.h"
#include "../../util/DataUtil.h"

void MakeConeTree(ConeTreeNode &node, const int nodeSize);
void VisitConeTree(ConeTreeNode &node);

void MakeConeTree(ConeTreeNode &node, const int nodeSize) {

	if (node.getSize() > nodeSize) {

		srand(time(NULL));
		int randomIndex = rand() % node.getSize();

		int firstIDIndex = -1;
		double f_minCosine = DBL_MAX;

		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == randomIndex) {
				continue;
			}
			double cosine = DataUtil::cosine(node[randomIndex], node[pointIndex]);

			if (cosine < f_minCosine) {
				firstIDIndex = pointIndex;
				f_minCosine = cosine;
			}
		}

		double s_minCosine = DBL_MAX;
		int secondIDIndex = -1;

		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == firstIDIndex) {
				continue;
			}
			double cosine = DataUtil::cosine(node[firstIDIndex], node[pointIndex]);

			if (cosine < s_minCosine) {
				secondIDIndex = pointIndex;
				s_minCosine = cosine;
			}
		}

		vector<Matrix1D> left;
		vector<Matrix1D> right;
		left.push_back(node[firstIDIndex]);
		right.push_back(node[secondIDIndex]);
		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == firstIDIndex || pointIndex == secondIDIndex) {
				continue;
			}

			double firstCosine = DataUtil::cosine(node[firstIDIndex], node[pointIndex]);
			double secondCosine = DataUtil::cosine(node[secondIDIndex], node[pointIndex]);

			if (firstCosine > secondCosine) {
				left.push_back(node[pointIndex]);
			} else {
				right.push_back(node[pointIndex]);
			}
		}

		node.splitNode(left, right);
		MakeConeTree(node.getLeftNode(), nodeSize);
		MakeConeTree(node.getRightNode(), nodeSize);
	} else {
		node.setLeafFlag();
	}
}

void MakeConeTree2(ConeTreeNode &node, const double cosineThreshold) {

	if (node.getCosineW_q() < cosineThreshold) {

		srand(time(NULL));
		int randomIndex = rand() % node.getSize();

		int firstIDIndex = -1;
		double f_minCosine = DBL_MAX;

		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == randomIndex) {
				continue;
			}
			double cosine = DataUtil::cosine(node[randomIndex], node[pointIndex]);

			if (cosine < f_minCosine) {
				firstIDIndex = pointIndex;
				f_minCosine = cosine;
			}
		}

		double s_minCosine = DBL_MAX;
		int secondIDIndex = -1;

		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == firstIDIndex) {
				continue;
			}
			double cosine = DataUtil::cosine(node[firstIDIndex], node[pointIndex]);

			if (cosine < s_minCosine) {
				secondIDIndex = pointIndex;
				s_minCosine = cosine;
			}
		}

		vector<Matrix1D> left;
		vector<Matrix1D> right;
		left.push_back(node[firstIDIndex]);
		right.push_back(node[secondIDIndex]);
		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == firstIDIndex || pointIndex == secondIDIndex) {
				continue;
			}

			double firstCosine = DataUtil::cosine(node[firstIDIndex], node[pointIndex]);
			double secondCosine = DataUtil::cosine(node[secondIDIndex], node[pointIndex]);

			if (firstCosine > secondCosine) {
				left.push_back(node[pointIndex]);
			} else {
				right.push_back(node[pointIndex]);
			}
		}

		node.splitNode(left, right);
		MakeConeTree2(node.getLeftNode(), cosineThreshold);
		MakeConeTree2(node.getRightNode(), cosineThreshold);
	} else {
		node.setLeafFlag();
	}
}

void VisitConeTree(ConeTreeNode &node) {
	if (!node.isLeafNode()) {
		VisitConeTree(node.getLeftNode());
		VisitConeTree(node.getRightNode());
	} else {
		cout << node.getSize() << ":" << node.isLeafNode() << endl;
	}
}