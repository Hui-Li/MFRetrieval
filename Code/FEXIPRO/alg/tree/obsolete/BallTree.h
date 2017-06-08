#pragma once

#include "BallTreeNode.h"
#include "../../util/DataUtil.h"

void MakeBallTree(BallTreeNode &node, const int nodeSize);
void VisitBallTree(BallTreeNode &node);

void MakeBallTree(BallTreeNode &node, const int nodeSize) {

	if (node.getSize() > nodeSize) {

		srand(time(NULL));
		int randomIndex = rand() % node.getSize();

		int firstIDIndex = -1;
		double f_maxDis = DBL_MIN;

		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == randomIndex) {
				continue;
			}
			double distance = DataUtil::euclideanDis(node[randomIndex], node[pointIndex]);

			if (distance > f_maxDis) {
				firstIDIndex = pointIndex;
				f_maxDis = distance;
			}
		}

		double s_maxDis = DBL_MIN;
		int secondIDIndex = -1;

		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == firstIDIndex) {
				continue;
			}
			double distance = DataUtil::euclideanDis(node[firstIDIndex], node[pointIndex]);

			if (distance > s_maxDis) {
				secondIDIndex = pointIndex;
				s_maxDis = distance;
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

			double firstDis = DataUtil::euclideanDis(node[firstIDIndex], node[pointIndex], false);
			double secondDis = DataUtil::euclideanDis(node[secondIDIndex], node[pointIndex], false);

			if (firstDis <= secondDis) {
				left.push_back(node[pointIndex]);
			} else {
				right.push_back(node[pointIndex]);
			}
		}

		node.splitNode(left, right);
		MakeBallTree(node.getLeftNode(), nodeSize);
		MakeBallTree(node.getRightNode(), nodeSize);
	} else {
		node.setLeafFlag();
	}
}

void MakeBallTreeSort(BallTreeNode &node, const int nodeSize) {

	if (node.getSize() > nodeSize) {

		srand(time(NULL));
		int randomIndex = rand() % node.getSize();

		int firstIDIndex = -1;
		double f_maxDis = DBL_MIN;

		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == randomIndex) {
				continue;
			}
			double distance = DataUtil::euclideanDis(node[randomIndex], node[pointIndex]);

			if (distance > f_maxDis) {
				firstIDIndex = pointIndex;
				f_maxDis = distance;
			}
		}

		double s_maxDis = DBL_MIN;
		int secondIDIndex = -1;

		for (int pointIndex = 0; pointIndex < node.getSize(); pointIndex++) {
			if (pointIndex == firstIDIndex) {
				continue;
			}
			double distance = DataUtil::euclideanDis(node[firstIDIndex], node[pointIndex]);

			if (distance > s_maxDis) {
				secondIDIndex = pointIndex;
				s_maxDis = distance;
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

			double firstDis = DataUtil::euclideanDis(node[firstIDIndex], node[pointIndex], false);
			double secondDis = DataUtil::euclideanDis(node[secondIDIndex], node[pointIndex], false);

			if (firstDis <= secondDis) {
				left.push_back(node[pointIndex]);
			} else {
				right.push_back(node[pointIndex]);
			}
		}

		sort(left.begin(), left.end(), greater<Matrix1D>());
		sort(right.begin(), right.end(), greater<Matrix1D>());
		node.splitNode(left, right);
		MakeBallTree(node.getLeftNode(), nodeSize);
		MakeBallTree(node.getRightNode(), nodeSize);
	} else {
		node.setLeafFlag();
	}
}

void VisitBallTree(BallTreeNode &node) {
	if (!node.isLeafNode()) {
		VisitBallTree(node.getLeftNode());
		VisitBallTree(node.getRightNode());
	} else {
		cout << node.getSize() << ":" << node.isLeafNode() << endl;
	}
}