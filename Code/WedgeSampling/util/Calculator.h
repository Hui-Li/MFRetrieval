#ifndef CALCULATOR_H
#define CALCULATOR_H

namespace Calculator{

    /**
     * Matrix Transportation
     * @param origin
     * @param target
     * @param num_of_row
     * @param num_of_col
     */
    inline void transpose(double *origin, double *target, const int num_of_row, const int num_of_col) {
        for (int row = 0; row < num_of_row; row++) {

            double *origin_row_ptr = origin + row * num_of_col;

            for (int col = 0; col < num_of_col; col++) {
                target[col * num_of_row + row] = origin_row_ptr[col];
            }
        }
    }

}
#endif //CALCULATOR_H
