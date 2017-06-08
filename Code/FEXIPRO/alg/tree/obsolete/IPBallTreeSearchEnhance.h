#pragma once

#include "../../util/Base.h"
#include "../../util/Matrix2D.h"
#include "../../util/DataUtil.h"
#include "../../util/Monitor.h"
#include "BallTree.h"

void IPBallTreeSearchEnhanceTopK(const int k, const int nodeSize, Matrix2D &userData, Matrix2D &itemData);

void IPBallTreeSearchEnhance(const int k, const Matrix1D &query, BallTreeNode &node,
                priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> &topk,
                const int dimension);


void IPBallTreeSearchEnhanceTopK(const int k, const int nodeSize, Matrix2D &userData, Matrix2D &itemData) {

	double step1;
	double step2;
	double step3;
	Monitor tt;

	// Step1: calculate norm
	tt.start();
	userData.CalNorms();
	itemData.CalNorms();
	tt.stop();
	step1 = tt.getElapsedTime();

	DataUtil::logger.Log("step 1 calculate norm: " + to_string(step1) + " secs", true);

	// Step2: build ball tree
	tt.start();
	vector<Matrix1D> points;
	for (int itemIndex = 0; itemIndex < itemData.num_row; itemIndex++) {
		points.push_back(itemData[itemIndex]);
	}

	BallTreeNode root(points, itemData.num_col);
	MakeBallTreeSort(root, nodeSize);
	tt.stop();
	step2 = tt.getElapsedTime();

	DataUtil::logger.Log("step 2 build ball tree: " + to_string(step2) + " secs", true);

	// Step3: search the tree
	tt.start();
	vector<priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> > topkLists;
	for (int userID = 0; userID < userData.num_row; userID++) {
		Matrix1D user = userData[userID];
		priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> topkSetInTheLoop;
		IPBallTreeSearchEnhance(k, user, root, topkSetInTheLoop, userData.num_col);

		topkLists.push_back(topkSetInTheLoop);
	}
	tt.stop();
	step3 = tt.getElapsedTime();
	DataUtil::globalTimer.stop();

	DataUtil::logger.Log("step 3 search the tree: " + to_string(step3) + " secs", true);

	if (Conf::outputResult) {
		FileUtil::outputResult(topkLists);
	}
}

void IPBallTreeSearchEnhance(const int k, const Matrix1D &query, BallTreeNode &node,
                priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> &topk,
                const int dimension) {

	if (node.isLeafNode()) {
		for (auto point : node.getPoints()) {
			if (topk.size() < k) {
				double innerProductValue = DataUtil::innerProduct(query, point);
				topk.push(make_pair(point.id, innerProductValue));
			} else {
				if (point.norm * query.norm >= topk.top().second) {
					double innerProductValue = DataUtil::innerProduct(query, point);
					topk.push(make_pair(point.id, innerProductValue));
					if (topk.size() > k) {
						topk.pop();
					}
				} else {
					break;
				}
			}
		}
	} else {

		double leftIP = DataUtil::innerProductArray(query, node.getLeftNode().getMean()) +
		                node.getLeftNode().getRadius() * query.norm;
		double rightIP = DataUtil::innerProductArray(query, node.getRightNode().getMean()) +
		                 node.getRightNode().getRadius() * query.norm;

		BallTreeNode *firstNode;
		BallTreeNode *secondNode;
		bool continueSearch = false;

		if (topk.size() < k) {
			firstNode = &node.getLeftNode();
			secondNode = &node.getRightNode();
			continueSearch = true;
		} else {
			if (leftIP >= rightIP && leftIP >= topk.top().second) {
				firstNode = &node.getLeftNode();
				secondNode = &node.getRightNode();
				continueSearch = true;
			} else if (rightIP >= topk.top().second) {
				firstNode = &node.getRightNode();
				secondNode = &node.getLeftNode();
				continueSearch = true;
			}
		}

		if (!continueSearch) {
			return;
		}

		IPBallTreeSearch(k, query, *firstNode, topk, dimension);

		if (topk.size() < k || DataUtil::innerProductArray(query, secondNode->getMean()) +
		                       secondNode->getRadius() * query.norm >= topk.top().second) {
			IPBallTreeSearch(k, query, *secondNode, topk, dimension);
		}
	}
}

