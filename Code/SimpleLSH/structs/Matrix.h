#ifndef MATRIX_H
#define MATRIX_H

#include "../util/Base.h"
#include "../util/Calculator.h"

class Matrix {

public:
    double *rawData;
    int rowNum;
    int colNum;

    inline Matrix() {
        this->rawData = NULL;
        this->rowNum = 0;
        this->colNum = 0;
    }

    inline ~Matrix() {
        delete[] rawData;
    }

    inline double *getRowPtr(const int rowIndex) const {
        return &rawData[rowIndex * colNum];
    }

    const double *operator[](int rowIndex) const {
        return rawData + rowIndex * colNum;
    }

    double *operator[](int rowIndex) {
        return rawData + rowIndex * colNum;
    }

    inline void init(double *rawData, const int rowNum, const int colNum) {
        this->rowNum = rowNum;
        this->colNum = colNum;
        this->rawData = rawData;
    }

    inline void init(const int rowNum, const int colNum) {
        this->rowNum = rowNum;
        this->colNum = colNum;
        this->rawData = new double[rowNum * colNum];
    }

    void transformPMatrix(Matrix &matrix) {

        double maxLen = 0;
        for (int i = 0; i < matrix.rowNum; i++) {
            double len = cal_len(matrix.getRowPtr(i), matrix.colNum);
            if (len > maxLen) {
                maxLen = len;
            }
        }
        double invMaxLen = 1 / maxLen;
        this->rowNum = matrix.rowNum;
        this->colNum = matrix.colNum + 1;
        this->rawData = new double[this->rowNum * this->colNum];

        for (int row_index = 0; row_index < rowNum; row_index++) {
            double *row_ptr = getRowPtr(row_index);
            double *ori_row_ptr = matrix.getRowPtr(row_index);
            for (int col_index = 0; col_index < matrix.colNum; col_index++) {
                row_ptr[col_index] = ori_row_ptr[col_index] * invMaxLen;
            }
            double len = cal_len(ori_row_ptr, matrix.colNum);
            row_ptr[matrix.colNum] = sqrt(1 - len * len / (maxLen * maxLen));
        }

    }

    void transformQMatrix(Matrix &matrix) {

        this->rowNum = matrix.rowNum;
        this->colNum = matrix.colNum + 1;
        this->rawData = new double[this->rowNum * this->colNum];

        for (int row_index = 0; row_index < matrix.rowNum; row_index++) {
            double *row_ptr = getRowPtr(row_index);
            double *ori_row_ptr = matrix.getRowPtr(row_index);
            double inv_ori_len = 1 / cal_len(ori_row_ptr, matrix.colNum);
            for (int col_index = 0; col_index < matrix.colNum; col_index++) {
                row_ptr[col_index] = ori_row_ptr[col_index] * inv_ori_len;
            }

            row_ptr[matrix.colNum] = 0;
        }

    }

    void readData(string dataFilePath) {
        vector <string> lines;
        string line;

        ifstream fin(dataFilePath.c_str());

        int rowNum = 0;
        int colNum = 0;

        while (getline(fin, line)) {
            if (line.length() == 0) {
                continue;
            }
            lines.push_back(line);
        }

        if (lines.size() == 0) {
            return;
        }

        fin.close();

        stringstream test(lines[0]);
        double tempValue;
        while (test >> tempValue) {
            if (test.peek() == ',') {
                test.ignore();
                colNum++;
            }
        }
        colNum++;
        rowNum = lines.size();

        this->rowNum = rowNum;
        this->colNum = colNum;
        this->rawData = new double[rowNum * colNum];
        int colIndex = 0;

        for (int rowIndex = 0; rowIndex < lines.size(); rowIndex++) {
            stringstream ss(lines[rowIndex]);
            colIndex = 0;
            while (ss >> this->rawData[rowIndex * colNum + colIndex]) {
                if (ss.peek() == ',') {
                    ss.ignore();
                    colIndex++;
                }
            }
        }
    }

    class Accessor {
        const Matrix &matrix_;
        std::vector<bool> flags_;
    public:
        typedef unsigned Key;
        typedef const double *Value;
        typedef double DATATYPE;

        Accessor(const Matrix &matrix) : matrix_(matrix) {
            flags_.resize(matrix_.rowNum);
        }

        void reset() {
            flags_.clear();
            flags_.resize(matrix_.rowNum);
        }

        bool mark(unsigned key) {
            if (flags_[key]) {
                return false;
            }
            flags_[key] = true;
            return true;
        }

        const double *operator()(unsigned key) {
            return matrix_[key];
        }
    };
};

#endif //MATRIX_H