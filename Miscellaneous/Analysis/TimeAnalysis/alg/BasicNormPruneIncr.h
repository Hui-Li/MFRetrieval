#ifndef BASICNORMPRUNEINCR_H
#define BASICNORMPRUNEINCR_H

#include "../util/Base.h"
#include "../util/FileUtil.h"
#include "../util/Calculator.h"
#include "../util/Monitor.h"
#include "../structs/VectorElement.h"
#include "../structs/Matrix.h"
#include "../structs/ExtendMatrix.h"
#include "../structs/ExtendMatrixRow.h"

// original order
void basicNormPruneIncr(const int k, const int checkDim, Matrix &q, Matrix &p) {

	Monitor tt;
	double offlineTime = 0;

	tt.start();
	// Step 1 (offline): Vectors with larger data are placed in the front.

	vector<VectorElement> pNorms(p.rowNum);
	vector<double> subNorms(p.rowNum);

	Calculator::calNormAndSubNorm(p, checkDim, pNorms, subNorms);
	sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());

	ExtendMatrix<ExtendMatrixRow> preprocessedP;
	preprocessedP.initExtendMatrix(p, pNorms, subNorms);

	tt.stop();
	offlineTime = tt.getElapsedTime();

	tt.start();

	vector<vector<VectorElement> > results(q.rowNum, vector<VectorElement>());
	double qNorm;
	double subQNorm;

	for (int qID = 0; qID < q.rowNum; qID++) {

		const double *qRow = q.getRowPtr(qID);
		Calculator::calSingleNorm(qRow, q.colNum, checkDim, qNorm, subQNorm);

		vector<VectorElement> &heap = results[qID];
		heap.resize(k);

		for (int pIndex = 0; pIndex < k; pIndex++) {
			const ExtendMatrixRow * pRow = preprocessedP.getRowPtr(pIndex);
			double innerProduct = Calculator::innerProduct(qRow, pRow->rawData, q.colNum);
			heap[pIndex] = VectorElement(pRow->gRowID, innerProduct);
		}

		make_heap(heap.begin(), heap.end(), greater<VectorElement>());
		double lowerBound = heap.front().data;

		Calculator::calSingleNorm(q.getRowPtr(qID), q.colNum, qNorm);

		for (int pIndex = k; pIndex < preprocessedP.rowNum; pIndex++) {

			const ExtendMatrixRow * pRow = preprocessedP.getRowPtr(pIndex);

			if (pRow->norm * qNorm <= lowerBound) {
				break;
			}
			else {

				const double *pPtr = pRow->rawData;

				double innerProduct = 0;
				for (int colIndex = 0; colIndex < checkDim; colIndex++) {
					innerProduct += qRow[colIndex] * pPtr[colIndex];
				}

				if (innerProduct + pRow->subNorm * subQNorm <=
					lowerBound) {
					continue;
				} else {
					for (int colIndex = checkDim; colIndex < p.colNum; colIndex++) {
						innerProduct += qRow[colIndex] * pPtr[colIndex];
					}

					if (innerProduct > lowerBound) {
						pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
						heap.pop_back();
						heap.push_back(VectorElement(pRow->gRowID, innerProduct));
						push_heap(heap.begin(), heap.end(), greater<VectorElement>());
						lowerBound = heap.front().data;
					}
				}

			}
		}
	}

	tt.stop();

	cout << "preprocess time: " << to_string(offlineTime) << " secs" << endl;
	cout << "online time: " << to_string(tt.getElapsedTime()) << " secs" << endl;

	string resultFileName = "result.txt";

	FileUtil::outputResult(results, resultFileName);

}

// original order
void basicNormPruneIncrRandomOrder(const int k, const int checkDim, Matrix &q, Matrix &p) {

	Monitor tt;
	double offlineTime = 0;

	tt.start();
	// Step 1 (offline): Vectors with larger data are placed in the front.

	vector<VectorElement> pNorms(p.rowNum);
	vector<double> subNorms(p.rowNum);

	Calculator::calNormAndSubNorm(p, checkDim, pNorms, subNorms);
	sort(pNorms.begin(), pNorms.end(), greater<VectorElement>());

	ExtendMatrix<ExtendMatrixRow> preprocessedP;
	preprocessedP.initExtendMatrix(p, pNorms, subNorms);

	// use for randomizing the sequence
	vector<int> rr = vector<int>(p.colNum);
	for (int i = 0; i < p.colNum; i++) {
		rr[i] = i;
	}
	random_shuffle(rr.begin(), rr.end());

	tt.stop();
	offlineTime = tt.getElapsedTime();

	tt.start();

	vector<vector<VectorElement> > results(q.rowNum, vector<VectorElement>());
	double qNorm;
	double subQNorm;
	double *reorderQ = new double[q.colNum]();

	for (int qID = 0; qID < q.rowNum; qID++) {

		const double *qRow = q.getRowPtr(qID);
		Calculator::calSingleNorm(qRow, q.colNum, checkDim, qNorm, subQNorm);

		random_shuffle(rr.begin(), rr.end());

		for (int i = 0; i < q.colNum; i++) {
			reorderQ[i] = qRow[rr[i]];
		}

		vector<VectorElement> &heap = results[qID];
		heap.resize(k);

		for (int pIndex = 0; pIndex < k; pIndex++) {
			const ExtendMatrixRow * pRow = preprocessedP.getRowPtr(pIndex);
			double innerProduct = Calculator::innerProduct(qRow, pRow->rawData, q.colNum);
			heap[pIndex] = VectorElement(pRow->gRowID, innerProduct);
		}

		make_heap(heap.begin(), heap.end(), greater<VectorElement>());
		double lowerBound = heap.front().data;

		Calculator::calSingleNorm(q.getRowPtr(qID), q.colNum, qNorm);

		for (int pIndex = k; pIndex < preprocessedP.rowNum; pIndex++) {

			const ExtendMatrixRow * pRow = preprocessedP.getRowPtr(pIndex);

			if (pRow->norm * qNorm <= lowerBound) {
				break;
			}
			else {

				const double *pPtr = pRow->rawData;

				double innerProduct = 0;
				for (int colIndex = 0; colIndex < checkDim; colIndex++) {
					innerProduct += reorderQ[colIndex] * pPtr[rr[colIndex]];
				}

				if (innerProduct + pRow->subNorm * subQNorm <=
					lowerBound) {
					continue;
				} else {
					for (int colIndex = checkDim; colIndex < p.colNum; colIndex++) {
						innerProduct += reorderQ[colIndex] * pPtr[rr[colIndex]];
					}

					if (innerProduct > lowerBound) {
						pop_heap(heap.begin(), heap.end(), greater<VectorElement>());
						heap.pop_back();
						heap.push_back(VectorElement(pRow->gRowID, innerProduct));
						push_heap(heap.begin(), heap.end(), greater<VectorElement>());
						lowerBound = heap.front().data;
					}
				}

			}
		}
	}

	tt.stop();

	delete [] reorderQ;

	cout << "preprocess time: " << to_string(offlineTime) << " secs" << endl;
	cout << "online time: " << to_string(tt.getElapsedTime()) << " secs" << endl;

	string resultFileName = "result.txt";

	FileUtil::outputResult(results, resultFileName);

}

#endif //BASICNORMPRUNEINCR_H