#ifndef ALIASSAMPLINGGENERATOR_H
#define ALIASSAMPLINGGENERATOR_H

#include <cstdlib>
#include <cstdio>
#include <cstring>

class AliasSamplingGenerator {

private:
    int *alias;
    double *prob;
    int n;

public:

    AliasSamplingGenerator() {
        alias = nullptr;
        prob = nullptr;
    }

    AliasSamplingGenerator(const int n, double *weight) : n(n) {
        alias = (int *) malloc(n * sizeof(int));
        prob = (double *) malloc(n * sizeof(double));
        if (alias == NULL || prob == NULL) {
            printf("Error: memory allocation failed!\n");
            exit(1);
        }

        double *norm_prob = (double *) malloc(n * sizeof(double));
        int *large_block = (int *) malloc(n * sizeof(int));
        int *small_block = (int *) malloc(n * sizeof(int));
        if (norm_prob == NULL || large_block == NULL || small_block == NULL) {
            printf("Error: memory allocation failed!\n");
            exit(1);
        }

        double sum = 0;
        int cur_small_block, cur_large_block;
        int num_small_block = 0, num_large_block = 0;

        for (int k = 0; k != n; k++) {
            sum += weight[k];
        }

        for (int k = 0; k != n; k++) {
            norm_prob[k] = weight[k] * n / sum;
        }

        for (int k = n - 1; k >= 0; k--) {
            if (norm_prob[k] < 1) {
                small_block[num_small_block] = k;
                num_small_block++;
            } else {
                large_block[num_large_block] = k;
                num_large_block++;
            }
        }

        while (num_small_block > 0 && num_large_block > 0) {
            num_small_block--;
            cur_small_block = small_block[num_small_block];
            num_large_block--;
            cur_large_block = large_block[num_large_block];
            prob[cur_small_block] = norm_prob[cur_small_block];
            alias[cur_small_block] = cur_large_block;
            norm_prob[cur_large_block] = norm_prob[cur_large_block] + norm_prob[cur_small_block] - 1;
            if (norm_prob[cur_large_block] < 1) {
                small_block[num_small_block] = cur_large_block;
                num_small_block++;
            } else {
                large_block[num_large_block] = cur_large_block;
                num_large_block++;
            }
        }

        while (num_large_block) {
            num_large_block--;
            prob[large_block[num_large_block]] = 1;
        }

        while (num_small_block) {
            num_small_block--;
            prob[small_block[num_small_block]] = 1;
        }

        free(norm_prob);
        free(small_block);
        free(large_block);
    }

    int sample(double rand_value1, double rand_value2) {
        int k = (int) n * rand_value1;
        return rand_value2 < prob[k] ? k : alias[k];
    }

    ~AliasSamplingGenerator() {
        free(alias);
        free(prob);
    }
};

#endif //ALIASSAMPLINGGENERATOR_H
