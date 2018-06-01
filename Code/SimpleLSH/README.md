# SimpleLSH

Use the transformation in the following paper to transform inner product search problem to KNN search problem:

*Yoram Bachrach, Yehuda Finkelstein, Ran Gilad-Bachrach, Liran Katzir, Noam Koenigstein, Nir Nice, and Ulrich Paquet. Speeding up the xbox recommender system using a euclidean transformation for inner-product spaces. In RecSys, pages 257–264, 2014.*

Then use LSH to solve the equivalent KNN problem. The lsh implementation is from [LSHBOX](https://github.com/RSIA-LIESMARS-WHU/LSHBOX/blob/master/include/lshbox/lsh/psdlsh.h), which use the algorithm in the following paper:

*MayurDatar,NicoleImmorlica,PiotrIndyk,andVahabS.Mirrokni.2004.Locality- sensitive hashing scheme based on p-stable distributions. In Symposium on Com- putational Geometry. 253–262.*