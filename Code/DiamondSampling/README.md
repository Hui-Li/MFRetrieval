# DiamondSampling

C++ implementation for papers: 

*G. Ballard, T. G. Kolda, A. Pinar, and C. Seshadhri. Diamond sampling for approximate maximum all-pairs dot-product (MAD) search. In ICDM, pages 11–20, 2015.*

We modified DiamondSampling for maximum inner product problem, which is different from its original purpose, i.e., maximum all-pairs inner product problem. See *Section VI - E* in the paper.
### Parameters
- `verify`: postprocess the approximate results (Default true).
- `k`: top k (Default 10).
- `search_k`: budget for top-k search in postprocessing (Default 100).
- `s`: number of samples (Default 1000).
- `q_file`: file path to Q data file
- `p_file`: file path to P data file
- `outputFilePath`: file path to result file
- `groundTruth`: file path to ground truth file.


You can check scripts/runDiamondSampling.sh for a running example.