#pragma once

#include "../../util/Matrix2D.h"
#include "../../util/Monitor.h"
#include "../../util/DataUtil.h"
#include "BallTree.h"
#include "DualBallTreeNode.h"

void IPDualBallTreeTopK(const int k, Matrix2D &userData, Matrix2D &itemData);
void IPDualBallTreeSearch(const int k, DualBallTreeNode &queryNode, DualBallTreeNode &refNode,
                          vector<priority_queue<pair<int, double>, vector<pair<int, double> >,
		                          Comparator::PairSmallCMP> > &topkLists,
                          const int dimension);
double CalBound(DualBallTreeNode &queryNode, DualBallTreeNode &refNode, const int dimension);

void IPDualBallTreeTopK(const int k, Matrix2D &userData, Matrix2D &itemData) {

	double step1;
	double step2;
	const int nodeSize = 20;
	Monitor tt;

	// Step1: build two ball trees
	tt.start();

	vector<Matrix1D> userPoints;
	for (int userIndex = 0; userIndex < userData.num_row; userIndex++) {
		userPoints.push_back(userData[userIndex]);
	}

	vector<Matrix1D> itemPoints;
	for (int itemIndex = 0; itemIndex < itemData.num_row; itemIndex++) {
		itemPoints.push_back(itemData[itemIndex]);
	}

	DualBallTreeNode userBallTreeRoot(userPoints, userData.num_col);
	MakeBallTree(userBallTreeRoot, nodeSize);

	DualBallTreeNode itemBallTreeRoot(itemPoints, itemData.num_col);
	MakeBallTree(itemBallTreeRoot, nodeSize);

	vector<priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> > topkLists
			(userData.num_row, priority_queue<pair<int, double>, vector<pair<int, double> >,
					Comparator::PairSmallCMP>());

	tt.stop();
	step1 = tt.getElapsedTime();

	DataUtil::logger.Log("step 1 build two ball tree: " + to_string(step1) + " secs", true);

	// Step2: search the tree
	tt.start();
	IPDualBallTreeSearch(k, userBallTreeRoot, itemBallTreeRoot, topkLists, itemData.num_col);
	tt.stop();
	step2 = tt.getElapsedTime();

	DataUtil::globalTimer.stop();

	DataUtil::logger.Log("step 2 search the tree: " + to_string(step2) + " secs", true);

	if (Conf::outputResult) {
		FileUtil::outputResult(topkLists);
	}
}

void IPDualBallTreeSearch(const int k, DualBallTreeNode &queryNode, DualBallTreeNode &refNode,
                          vector<priority_queue<pair<int, double>, vector<pair<int, double> >,
		                          Comparator::PairSmallCMP> > &topkLists,
                          const int dimension) {

	double minIP = CalBound(queryNode, refNode, dimension);

	if (queryNode.getMinIP() < minIP) {
		if (queryNode.isLeafNode() && refNode.isLeafNode()) {
			double minIP = DBL_MAX;
			for (auto queryPoint:queryNode.getPoints()) {

				for (auto refPoint:refNode.getPoints()) {
					double innerProductValue = DataUtil::innerProduct(queryPoint, refPoint);
					topkLists[queryPoint.id].push(make_pair(refPoint.id, innerProductValue));
				}

				while (topkLists[queryPoint.id].size() > k) {
					topkLists[queryPoint.id].pop();
				}

				if (topkLists[queryPoint.id].top().second < minIP) {
					minIP = topkLists[queryPoint.id].top().second;
				}
			}

			queryNode.setMinIP(minIP);

		} else if (refNode.isLeafNode()) {

			IPDualBallTreeSearch(k, queryNode.getLeftNode(), refNode, topkLists, dimension);
			IPDualBallTreeSearch(k, queryNode.getRightNode(), refNode, topkLists, dimension);
			double minIP = queryNode.getLeftNode().getMinIP() < queryNode.getRightNode().getMinIP() ? queryNode
					.getLeftNode().getMinIP() : queryNode.getRightNode().getMinIP();
			queryNode.setMinIP(minIP);

		} else if (queryNode.isLeafNode()) {

			double leftBound = CalBound(queryNode, refNode.getLeftNode(), dimension);
			double rightBound = CalBound(queryNode, refNode.getRightNode(), dimension);

			if (leftBound <= rightBound) {
				IPDualBallTreeSearch(k, queryNode, refNode.getRightNode(), topkLists, dimension);
				IPDualBallTreeSearch(k, queryNode, refNode.getLeftNode(), topkLists, dimension);
			} else {
				IPDualBallTreeSearch(k, queryNode, refNode.getLeftNode(), topkLists, dimension);
				IPDualBallTreeSearch(k, queryNode, refNode.getRightNode(), topkLists, dimension);
			}

		} else {
			double leftBound = CalBound(queryNode.getLeftNode(), refNode.getLeftNode(), dimension);
			double rightBound = CalBound(queryNode.getLeftNode(), refNode.getRightNode(), dimension);

			if (leftBound <= rightBound) {
				IPDualBallTreeSearch(k, queryNode.getLeftNode(), refNode.getRightNode(), topkLists, dimension);
				IPDualBallTreeSearch(k, queryNode.getLeftNode(), refNode.getLeftNode(), topkLists, dimension);
			} else {
				IPDualBallTreeSearch(k, queryNode.getLeftNode(), refNode.getLeftNode(), topkLists, dimension);
				IPDualBallTreeSearch(k, queryNode.getLeftNode(), refNode.getRightNode(), topkLists, dimension);
			}

			leftBound = CalBound(queryNode.getRightNode(), refNode.getLeftNode(), dimension);
			rightBound = CalBound(queryNode.getRightNode(), refNode.getRightNode(), dimension);

			if (leftBound <= rightBound) {
				IPDualBallTreeSearch(k, queryNode.getRightNode(), refNode.getRightNode(), topkLists, dimension);
				IPDualBallTreeSearch(k, queryNode.getRightNode(), refNode.getLeftNode(), topkLists, dimension);
			} else {
				IPDualBallTreeSearch(k, queryNode.getRightNode(), refNode.getLeftNode(), topkLists, dimension);
				IPDualBallTreeSearch(k, queryNode.getRightNode(), refNode.getRightNode(), topkLists, dimension);
			}

			double minIP = queryNode.getLeftNode().getMinIP() < queryNode.getRightNode().getMinIP() ? queryNode
					.getLeftNode().getMinIP() : queryNode.getRightNode().getMinIP();
			queryNode.setMinIP(minIP);

		}
	}
}

double CalBound(DualBallTreeNode &queryNode, DualBallTreeNode &refNode, const int dimension) {
	return DataUtil::innerProductDualArray(queryNode.getMean(), refNode.getMean(), dimension) +
	       queryNode.getRadius() * refNode.getRadius() + queryNode.getCenterNorm() * refNode.getRadius() +
	       refNode.getCenterNorm() * queryNode.getRadius();
}