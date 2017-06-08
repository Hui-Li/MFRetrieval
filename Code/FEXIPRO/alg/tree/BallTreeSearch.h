#ifndef BALLTREESEARCH_H
#define BALLTREESEARCH_H

#include "../../util/Base.h"
#include "../../util/Monitor.h"
#include "../../structs/VectorElement.h"
#include "../../structs/Matrix.h"
#include "../../util/Calculator.h"
#include "../../util/Logger.h"
#include "../../structs/BallTreeNode.h"

void ballTreeTopK(const int k, Matrix &q, Matrix &p);
void ballTreeSearch(const int k, const double *qRow, const double qNorm, const int dim, BallTreeNode *node, vector<VectorElement> &heap);
void makeBallTree(BallTreeNode *node, Matrix &p, const int nodeSize);

#ifdef TIME_IT
uint64_t counter = 0;
#endif

void ballTreeTopK(const int k, Matrix &q, Matrix &p) {

	const int nodeSize = 20;
	Monitor tt;
	double offlineTime;


	tt.start();

	// Step1 (offline): build ball tree

	vector<int> pointIDs(p.rowNum, 0);
	vector<const double *> pointPtrs(p.rowNum, NULL);
	for (int id = 0; id < p.rowNum; id++) {
		pointIDs[id] = id;
		pointPtrs[id] = p.getRowPtr(id);
	}

	BallTreeNode root(pointIDs, pointPtrs, p.colNum);
	makeBallTree(&root, p, nodeSize);

	tt.stop();
	offlineTime = tt.getElapsedTime();


	// Step2 (online): search the tree
	tt.start();

	vector<vector<VectorElement> > results(q.rowNum, vector<VectorElement>(k, VectorElement(0, -DBL_MAX)));
	double qNorm;

	for (int qID = 0; qID < q.rowNum; qID++) {
		vector<VectorElement> &currentHeap = results[qID];
		make_heap(currentHeap.begin(), currentHeap.end(), greater<VectorElement>());
		const double *qPtr = q.getRowPtr(qID);
		Calculator::calSingleNorm(qPtr, q.colNum, qNorm);
		ballTreeSearch(k, qPtr, qNorm, q.colNum, &root, results[qID]);
	}

	tt.stop();

#ifdef TIME_IT
	Logger::Log("Avg Num of Inner Product Calculation Before Stop: " + to_string((counter + 0.0) / q.rowNum));
#endif
	Logger::Log("offline time: " + to_string(offlineTime) + " secs");
	Logger::Log("online time: " + to_string(tt.getElapsedTime()) + " secs");

	if (Conf::outputResult) {
		string resultFileName = Conf::resultPathPrefix + "-" + Conf::dataset + "-" + Conf::algName + "-" + to_string(Conf::k) + ".txt";
		FileUtil::outputResult(k, results, resultFileName);
	}
}


void ballTreeSearch(const int k, const double *qRow, const double qNorm, const int dim, BallTreeNode *node,
					vector<VectorElement> &heap) {

	if (node->isLeafNode()) {
		for (int index = 0; index < node->getSize(); index++) {
#ifdef TIME_IT
			counter++;
#endif
			double innerProductValue = Calculator::innerProduct(node->getPtr(index), qRow, dim);
			if (innerProductValue > heap.front().data) {
				pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
				heap.pop_back();
				heap.push_back(VectorElement(node->getID(index), innerProductValue));
				push_heap(heap.begin(), heap.end(), greater<VectorElement>());
			}
		}
	} else {

		double leftIP = Calculator::innerProduct(node->getLeftNode()->getMean(), qRow, dim) +
		                node->getLeftNode()->getConstrain() * qNorm;
		double rightIP = Calculator::innerProduct(node->getRightNode()->getMean(), qRow, dim) +
		                 node->getRightNode()->getConstrain() * qNorm;

		BallTreeNode *firstNode;
		BallTreeNode *secondNode;
		bool continueSearch = false;

		if (leftIP >= rightIP && leftIP >= heap.front().data) {
			firstNode = node->getLeftNode();
			secondNode = node->getRightNode();
			continueSearch = true;
		} else if (rightIP >= heap.front().data) {
			firstNode = node->getRightNode();
			secondNode = node->getLeftNode();
			continueSearch = true;
		}

		if (!continueSearch) {
			return;
		}

		ballTreeSearch(k, qRow, qNorm, dim, firstNode, heap);

		if (Calculator::innerProduct(secondNode->getMean(), qRow, dim) + secondNode->getConstrain() *
		                                                                             qNorm >= heap.front().data) {
			ballTreeSearch(k, qRow, qNorm, dim, secondNode, heap);
		}
	}
}

inline void makeBallTree(BallTreeNode *node, Matrix &p, const int nodeSize) {

	if (node->getSize() > nodeSize) {
		srand(time(NULL));
		int randomIndex = rand() % node->getSize();
		const double* randomPointPtr = p.getRowPtr(node->getID(randomIndex));

		int firstIDIndex = -1;

		double f_maxDis = -1;
		double distance;
		for (int pointIndex = 0; pointIndex < node->getSize(); pointIndex++) {
			if (pointIndex == randomIndex) {
				continue;
			}
			const double* pointPtr = p.getRowPtr(node->getID(pointIndex));
			 distance = Calculator::l2Distance(pointPtr, randomPointPtr, p.colNum);

			if (distance > f_maxDis) {
				firstIDIndex = pointIndex;
				f_maxDis = distance;
			}
		}

		int firstPointID = node->getID(firstIDIndex);
		const double* firstPointPtr = p.getRowPtr(firstPointID);

		double s_maxDis = -1;
		int secondIDIndex = -1;

		for (int pointIndex = 0; pointIndex < node->getSize(); pointIndex++) {
			if (pointIndex == firstIDIndex) {
				continue;
			}
			const double* pointPtr = p.getRowPtr(node->getID(pointIndex));
			double distance = Calculator::l2Distance(pointPtr, firstPointPtr, p.colNum);

			if (distance > s_maxDis) {
				secondIDIndex = pointIndex;
				s_maxDis = distance;
			}
		}

		vector<int> leftPointIDs;
		vector<const double *> leftPointPtrs;
		vector<int> rightPointIDs;
		vector<const double *> rightPointPtrs;

		int secondPointID = node->getID(secondIDIndex);
		const double* secondPointPtr = p.getRowPtr(secondPointID);
		leftPointIDs.push_back(firstPointID);
		leftPointPtrs.push_back(firstPointPtr);
		rightPointIDs.push_back(secondPointID);
		rightPointPtrs.push_back(secondPointPtr);

		for (int pointIndex = 0; pointIndex < node->getSize(); pointIndex++) {
			if (pointIndex == firstIDIndex || pointIndex == secondIDIndex) {
				continue;
			}

			int pointID = node->getID(pointIndex);
			const double *pointPtr = p.getRowPtr(pointID);
			double firstDis = Calculator::l2Distance(firstPointPtr, pointPtr, p.colNum);
			double secondDis = Calculator::l2Distance(secondPointPtr, pointPtr, p.colNum);

			if (firstDis <= secondDis) {
				leftPointIDs.push_back(pointID);
				leftPointPtrs.push_back(pointPtr);
			} else {
				rightPointIDs.push_back(pointID);
				rightPointPtrs.push_back(pointPtr);
			}
		}

		node->splitNode(leftPointIDs, leftPointPtrs, rightPointIDs, rightPointPtrs);
		makeBallTree(node->getLeftNode(), p, nodeSize);
		makeBallTree(node->getRightNode(), p, nodeSize);
	} else {
		node->setLeafFlag();
	}
}

#endif //BALLTREESEARCH_H