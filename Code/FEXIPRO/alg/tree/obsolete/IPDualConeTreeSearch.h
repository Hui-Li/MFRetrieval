#pragma once

#include "../../util/Matrix2D.h"
#include "../../util/Monitor.h"
#include "../../util/DataUtil.h"
#include "BallTree.h"
#include "ConeTree.h"
#include "DualBallTreeNode.h"
#include "ConeTreeNode.h"

void IPDualConeTreeSearch(const int k, ConeTreeNode &queryNode, DualBallTreeNode &refNode,
                              vector<priority_queue<pair<int, double>, vector<pair<int, double> >,
		                              Comparator::PairSmallCMP> > &topkLists,
                              const int dimension);
double CalBound(ConeTreeNode &queryNode, ConeTreeNode &refNode, const int dimension);

void IPDualConeTreeTopK(const int k, Matrix2D &userData, Matrix2D &itemData) {

	double step1;
	double step2;
	double step3;
	double step4;

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

	// Step2: build dual tree
	tt.start();

	vector<Matrix1D> userPoints;
	for (int userIndex = 0; userIndex < userData.num_row; userIndex++) {
		userPoints.push_back(userData[userIndex]);
	}

	vector<Matrix1D> itemPoints;
	for (int itemIndex = 0; itemIndex < itemData.num_row; itemIndex++) {
		itemPoints.push_back(itemData[itemIndex]);
	}

	ConeTreeNode userConeTreeRoot(userPoints, userData.num_col);
	MakeConeTree(userConeTreeRoot, nodeSize);

	ConeTreeNode itemConeTreeRoot(itemPoints, itemData.num_col);
	MakeConeTree(itemConeTreeRoot, nodeSize);

	vector<priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> > topkLists
			(userData.num_row, priority_queue<pair<int, double>, vector<pair<int, double> >,
					Comparator::PairSmallCMP>());

	tt.stop();
	step2 = tt.getElapsedTime();

	DataUtil::logger.Log("step 2 build dual cone tree: " + to_string(step2) + " secs", true);

	// Step3: search the tree
	tt.start();
	IPDualConeTreeSearch(k, userConeTreeRoot, itemConeTreeRoot, topkLists, itemData.num_col);
	tt.stop();
	step3 = tt.getElapsedTime();

	DataUtil::logger.Log("step 3 search the tree: " + to_string(step3) + " secs", true);

	// Step4: recover to original rating
	tt.start();
	vector<priority_queue<pair<int, double>, vector<pair<int, double> >, Comparator::PairSmallCMP> > originalTopkLists
			(userData.num_row, priority_queue<pair<int, double>, vector<pair<int, double> >,
					Comparator::PairSmallCMP>());
	for (int userIndex = 0; userIndex < topkLists.size(); userIndex++) {
//		if(topkLists[userIndex].size()!=10){
//			cout << userIndex << endl;
//		}
		while (!topkLists[userIndex].empty()) {
			originalTopkLists[userIndex].push(make_pair(topkLists[userIndex].top().first, topkLists[userIndex].top()
					.second * userData[userIndex].norm));
			topkLists[userIndex].pop();
		}
	}
	tt.stop();
	step4 = tt.getElapsedTime();

	DataUtil::globalTimer.stop();

	DataUtil::logger.Log("step 4 recover to original rating: " + to_string(step4) + " secs", true);

	if (Conf::outputResult) {
		FileUtil::outputResult(originalTopkLists);
	}
}

void IPDualConeTreeSearch(const int k, ConeTreeNode &queryNode, ConeTreeNode &refNode,
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

			IPDualConeTreeSearch(k, queryNode.getLeftNode(), refNode, topkLists, dimension);
			IPDualConeTreeSearch(k, queryNode.getRightNode(), refNode, topkLists, dimension);
			double minIP = queryNode.getLeftNode().getMinIP() < queryNode.getRightNode().getMinIP() ? queryNode
					.getLeftNode().getMinIP() : queryNode.getRightNode().getMinIP();
			queryNode.setMinIP(minIP);

		} else if (queryNode.isLeafNode()) {

			double leftBound = CalBound(queryNode, refNode.getLeftNode(), dimension);
			double rightBound = CalBound(queryNode, refNode.getRightNode(), dimension);

			if (leftBound <= rightBound) {
				IPDualConeTreeSearch(k, queryNode, refNode.getRightNode(), topkLists, dimension);
				IPDualConeTreeSearch(k, queryNode, refNode.getLeftNode(), topkLists, dimension);
			} else {
				IPDualConeTreeSearch(k, queryNode, refNode.getLeftNode(), topkLists, dimension);
				IPDualConeTreeSearch(k, queryNode, refNode.getRightNode(), topkLists, dimension);
			}

		} else {
			double leftBound = CalBound(queryNode.getLeftNode(), refNode.getLeftNode(), dimension);
			double rightBound = CalBound(queryNode.getLeftNode(), refNode.getRightNode(), dimension);

			if (leftBound <= rightBound) {
				IPDualConeTreeSearch(k, queryNode.getLeftNode(), refNode.getRightNode(), topkLists, dimension);
				IPDualConeTreeSearch(k, queryNode.getLeftNode(), refNode.getLeftNode(), topkLists, dimension);
			} else {
				IPDualConeTreeSearch(k, queryNode.getLeftNode(), refNode.getLeftNode(), topkLists, dimension);
				IPDualConeTreeSearch(k, queryNode.getLeftNode(), refNode.getRightNode(), topkLists, dimension);
			}

			leftBound = CalBound(queryNode.getRightNode(), refNode.getLeftNode(), dimension);
			rightBound = CalBound(queryNode.getRightNode(), refNode.getRightNode(), dimension);

			if (leftBound <= rightBound) {
				IPDualConeTreeSearch(k, queryNode.getRightNode(), refNode.getRightNode(), topkLists, dimension);
				IPDualConeTreeSearch(k, queryNode.getRightNode(), refNode.getLeftNode(), topkLists, dimension);
			} else {
				IPDualConeTreeSearch(k, queryNode.getRightNode(), refNode.getLeftNode(), topkLists, dimension);
				IPDualConeTreeSearch(k, queryNode.getRightNode(), refNode.getRightNode(), topkLists, dimension);
			}

			double minIP = queryNode.getLeftNode().getMinIP() < queryNode.getRightNode().getMinIP() ? queryNode
					.getLeftNode().getMinIP() : queryNode.getRightNode().getMinIP();
			queryNode.setMinIP(minIP);

		}
	}
}

double CalBound(ConeTreeNode &queryNode, ConeTreeNode &refNode, const int dimension) {
	double cosineBetweenPerimeters = DataUtil::cosineDualArray(queryNode.getMean(), refNode.getMean(), dimension);
	double sumOfDelta = queryNode.getCosineW_q() * refNode.getCosineW_q() - queryNode.getSineW_q() * refNode
			                                                                                                 .getSineW_q();
	if (cosineBetweenPerimeters < queryNode.getCosineW_q()) {
		double cosine = cosineBetweenPerimeters * queryNode.getCosineW_q() + DataUtil::calSine(cosineBetweenPerimeters)
		                                                                     * queryNode.getSineW_q();
		return refNode.getCenterNorm() * cosine + refNode.getRadius();

	} else {
		return refNode.getCenterNorm() + refNode.getRadius();
	}

}
