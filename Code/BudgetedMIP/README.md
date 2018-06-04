# BudgetedMIP

C++ implementation for papers: 

*Hsiang-Fu Yu, Cho-Jui Hsieh, Qi Lei, and Inderjit S. Dhillon. 2017. A Greedy Approach for Budgeted Maximum Inner Product Search. In NIPS. 5459–5468.*

See [https://github.com/rofuyu/exp-gmips-nips17](https://github.com/rofuyu/exp-gmips-nips17) for the details.

### Parameters
- `k`: top k (Default 10).
- `b`: budget for top-k search (Default 1024).
- `q_file`: file path to Q data file.
- `p_file`: file path to P data file.
- `outputFilePath`: file path to result file.
- `groundTruth`: file path to ground truth file.


You can check scripts/runBudgetedMIP.sh for a running example.