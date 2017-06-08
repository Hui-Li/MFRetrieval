#ifndef SIMDUTIL_H
#define SIMDUTIL_H

#include <immintrin.h>
#include "Base.h"
#include "../structs/Matrix.h"

namespace SIMDUtil {


    inline int int32ForIntDotProduct(const int N, const int lowerN, const int16_t *a, const int16_t *b) {
#ifdef __AVX2__
        // 256 bits = 8 * int_32 = 16 * int_16 = 32 * int_8
        __m256i temp_sum = _mm256_setzero_si256();
        __m256i temp_1;
        __m256i temp_2;

        int index = 0;
        for (index = 0; index < lowerN; index += 16) {

            temp_1 = _mm256_loadu_si256((__m256i *)(a + index));
            temp_2 = _mm256_loadu_si256((__m256i *)(b + index));

            temp_sum = _mm256_add_epi32(temp_sum, _mm256_madd_epi16(temp_1, temp_2));
        }

        temp_sum = _mm256_hadd_epi32(temp_sum, temp_sum);
        temp_sum = _mm256_hadd_epi32(temp_sum, temp_sum);

        int sum = _mm256_extract_epi32(temp_sum, 0) + _mm256_extract_epi32(temp_sum, 4);

        for (; index < N; index++) {
            sum += a[index] * b[index];
        }

        return sum;
#else
        exit(-1);
#endif
    }

    inline double doubleDotProduct(const int N, const int lowerN, const double *a, const double *b) {
#ifdef __AVX2__
        __m256d sum = _mm256_setzero_pd();
        __m256d xmm1;
        __m256d xmm2;

        int i;
        for (i = 0; i < lowerN; i += 4) {
            xmm1 = _mm256_loadu_pd(a + i);
            xmm2 = _mm256_loadu_pd(b + i);
            sum = _mm256_add_pd(sum, _mm256_mul_pd(xmm1, xmm2));
        }

        sum = _mm256_hadd_pd(sum, sum);

        double result = ((double*)&sum)[0] + ((double*)&sum)[2];

        for(;i<N;i++) {
            result+=a[i] * b[i];
        }

        return result;
#else
        exit(-1);
#endif
    }

    inline double doubleSingleNorm(const int N, const int lowerN, const double *a) {
#ifdef __AVX2__
        __m256d sum = _mm256_setzero_pd();
        __m256d xmm1;
        __m256d xmm2;

        int i;
        for (i = 0; i < lowerN; i += 4) {
            xmm1 = _mm256_loadu_pd(a + i);
            xmm2 = _mm256_loadu_pd(a + i);
            sum = _mm256_add_pd(sum, _mm256_mul_pd(xmm1, xmm2));
        }

        sum = _mm256_hadd_pd(sum, sum);

        double result = ((double*)&sum)[0] + ((double*)&sum)[2];

        for(;i<N;i++) {
            result+=a[i] * a[i];
        }

        return sqrt(result);
#else
        exit(-1);
#endif
    }

    inline void doubleNorms(const int N, const int lowerN, const Matrix &m, vector<VectorElement> &norms) {

        norms.resize(m.rowNum);

        for (int rowID = 0; rowID < m.rowNum; rowID++) {
            const double *row = m.getRowPtr(rowID);
            norms[rowID] = VectorElement(rowID, doubleSingleNorm(N, lowerN, row));
        }
    }

    // http://stackoverflow.com/questions/15238978/sse3-intrinsics-how-to-find-the-maximum-of-a-large-array-of-floats
    inline double maxValue(const int N, const int lowerN, const double *a) {
#ifdef __AVX2__
        __m256d maxVector = _mm256_loadu_pd(a);
        __m256d cur;
        int i;
        for (i = 4; i < lowerN; i += 4) {
            cur = _mm256_loadu_pd(a + i);
            maxVector = _mm256_max_pd(maxVector, cur);
        }
        double maxValue = max(((double*)&maxVector)[3] , max(((double*)&maxVector)[2], max(((double*)&maxVector)[0], ((double*)&maxVector)[1])));

        for(;i<N;i++){
            maxValue = max(maxValue, a[i]);
        }

        return maxValue;
#else
        exit(-1);
#endif
    }

    inline double minValue(const int N, const int lowerN, const double *a) {
#ifdef __AVX2__
        __m256d minVector = _mm256_loadu_pd(a);
        __m256d cur;
        int i;
        for (i = 4; i < lowerN; i += 4) {
            cur = _mm256_loadu_pd(a + i);
            minVector = _mm256_min_pd(minVector, cur);
        }
        double minValue = min(((double*)&minVector)[3] , min(((double*)&minVector)[2], min(((double*)&minVector)[0], ((double*)&minVector)[1])));

        for(;i<N;i++){
            minValue = min(minValue, a[i]);
        }

        return minValue;
#else
        exit(-1);
#endif
    }

    inline void minMaxValue(const int N, const int lowerN, const double *a, double &maxValue, double &minValue) {
#ifdef __AVX2__
        __m256d maxVector = _mm256_loadu_pd(a);
        __m256d minVector = _mm256_loadu_pd(a);
        __m256d cur;
        int i;
        for (i = 4; i < lowerN; i += 4) {
            cur = _mm256_loadu_pd(a + i);
            minVector = _mm256_min_pd(minVector, cur);
            maxVector = _mm256_max_pd(maxVector, cur);
        }

        maxValue = max(((double*)&maxVector)[3] , max(((double*)&maxVector)[2], max(((double*)&maxVector)[0], ((double*)&maxVector)[1])));
        minValue = min(((double*)&minVector)[3] , min(((double*)&minVector)[2], min(((double*)&minVector)[0], ((double*)&minVector)[1])));

        for(;i<N;i++){
            maxValue = max(maxValue, a[i]);
            minValue = min(minValue, a[i]);
        }
#else
        exit(-1);
#endif

    }
}
#endif //SIMDUTIL_H