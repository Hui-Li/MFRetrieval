# WedgeSampling

C++ implementation for papers: 


+ *Edith Cohen and David D. Lewis. Approximating matrix multiplication for pat- tern recognition tasks. J. Algorithms, 30(2):211–252, 1999.*

+ *Edith Cohen and David D. Lewis. Approximating matrix multiplication for pat- tern recognition tasks. In SODA, pages 682–691, 1997.*

### Parameters
- `compareWithNaive`: the results file will include results from naive method (Default true).
- `k`: top k (Default 10).
- `s`: number of samples (Default 1000).
- `q_file`: file path to Q data file
- `p_file`: file path to P data file
- `outputFilePath`: file path to result file


You can check scripts/runWedgeSampling.sh for a running example.