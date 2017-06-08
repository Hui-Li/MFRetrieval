#pragma once

#include "../../util/Base.h"
#include "../../util/Matrix2D.h"
#include "../../util/DataUtil.h"
#include "../../util/Monitor.h"
#include "BallTree.h"

void IPBallTreeTopK2(const int k, Matrix2D &userData, Matrix2D &itemData);
void IPBallTreeSearch2(const int k, const Matrix1D &query, BallTreeNode &node,
                       priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> &topk,
                       const int dimension);

// mapping and then do inner product search
void IPBallTreeTopK2(const int k, Matrix2D &userData, Matrix2D &itemData) {

	double step1;
	double step2;
	double step3;
	double step4;
	double step5;
	const int nodeSize = 20;
	Monitor tt;

	// Step1: map data
	tt.start();
	double maxUserNorm = DBL_MIN;
	double maxItemNorm = DBL_MIN;
	userData.CalNorms();
	itemData.CalNorms();

	for (int userID = 0; userID < userData.num_row; userID++) {
		if (userData[userID].norm > maxUserNorm) {
			maxUserNorm = userData[userID].norm;
		}
	}

	for (int itemID = 0; itemID < itemData.num_row; itemID++) {
		if (itemData[itemID].norm > maxItemNorm) {
			maxItemNorm = itemData[itemID].norm;
		}
	}

	Matrix2D transferedUserData;
	transferedUserData.Init(userData.num_row, userData.num_col + 2);
	for (int userID = 0; userID < userData.num_row; userID++) {
		for (int dimIndex = 0; dimIndex < userData.num_col; dimIndex++) {
			transferedUserData[userID][dimIndex] = userData[userID][dimIndex] / maxUserNorm;
		}
		transferedUserData[userID][userData.num_col] = 0;
		transferedUserData[userID][userData.num_col + 1] = sqrt(1 - (userData[userID].norm * userData[userID].norm) /
		                                                            (maxUserNorm * maxUserNorm));
	}

	vector<Matrix1D> points(itemData.num_row);
	Matrix2D transferedItemData;
	transferedItemData.Init(itemData.num_row, itemData.num_col + 2);
	for (int itemID = 0; itemID < itemData.num_row; itemID++) {
		for (int dimIndex = 0; dimIndex < itemData.num_col; dimIndex++) {
			transferedItemData[itemID][dimIndex] = itemData[itemID][dimIndex] / maxItemNorm;
		}
		transferedItemData[itemID][itemData.num_col] = sqrt(1 - (itemData[itemID].norm * itemData[itemID].norm) /
		                                                        (maxItemNorm * maxItemNorm));
		transferedItemData[itemID][itemData.num_col + 1] = 0;

		points[itemID] = transferedItemData[itemID];
	}

	tt.stop();
	step1 = tt.getElapsedTime();

	DataUtil::logger.Log("step 1 map data: " + to_string(step1) + " secs", true);

	// Step2: build ball tree
	tt.start();

	BallTreeNode root(points, transferedItemData.num_col);
	MakeBallTree(root, nodeSize);
	tt.stop();
	step2 = tt.getElapsedTime();

	DataUtil::logger.Log("step 2 build ball tree: " + to_string(step3) + " secs", true);

	// Step3: search the tree
	tt.start();
	vector<priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> > topkLists
			(transferedUserData.num_row, priority_queue<pair<int, double>, vector<pair<int, double> >,
					Comparator::PairSmallCMP>());
	for (int userID = 0; userID < transferedUserData.num_row; userID++) {
		Matrix1D user = transferedUserData[userID];
		IPBallTreeSearch2(k, user, root, topkLists[userID], transferedUserData.num_col);
	}
	tt.stop();
	step3 = tt.getElapsedTime();

	DataUtil::logger.Log("step 3 search the tree: " + to_string(step3) + " secs", true);

	// Step4: recover to original rating
	tt.start();
	vector<vector<pair<int, double> > > originalTopkLists(topkLists.size(), vector<pair<int, double> >());
	for (int userIndex = 0; userIndex < topkLists.size(); userIndex++) {
		while (!topkLists[userIndex].empty()) {
//			originalTopkLists[userIndex].push(make_pair(topkLists[userIndex].top().first, topkLists[userIndex].top()
//					                                                                              .second *
//					maxItemNorm * maxUserNorm));
			originalTopkLists[userIndex].push_back(make_pair(topkLists[userIndex].top().first, DataUtil::innerProduct
					(itemData[topkLists[userIndex].top().first], userData[userIndex])));
			topkLists[userIndex].pop();
		}
	}
	tt.stop();
	step4 = tt.getElapsedTime();

	DataUtil::logger.Log("step 4 recover to original rating: " + to_string(step3) + " secs", true);

	DataUtil::globalTimer.stop();

	if (Conf::outputResult) {
		FileUtil::outputResult(originalTopkLists);
	}
}

void IPBallTreeSearch2(const int k, const Matrix1D &query, BallTreeNode &node,
                      priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> &topk,
                      const int dimension) {

	if (node.isLeafNode()) {
		for (auto point : node.getPoints()) {
			double innerProductValue = DataUtil::innerProduct(query, point);
			topk.push(make_pair(point.id, innerProductValue));
			if (topk.size() > k) {
				topk.pop();
			}
		}
	} else {

		double leftIP = DataUtil::innerProductArray(query, node.getLeftNode().getMean()) +
		                node.getLeftNode().getRadius();
		double rightIP = DataUtil::innerProductArray(query, node.getRightNode().getMean()) +
		                 node.getRightNode().getRadius();

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

		IPBallTreeSearch2(k, query, *firstNode, topk, dimension);

		if (topk.size() < k || DataUtil::innerProductArray(query, secondNode->getMean()) +
		                       secondNode->getRadius() >= topk.top().second) {
			IPBallTreeSearch2(k, query, *secondNode, topk, dimension);
		}
	}
}

