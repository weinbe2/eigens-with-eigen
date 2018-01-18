# eigens-with-eigen
A set of reference code showing how to calculate eigenvalues and eigenvectors with Eigen. Dense Eigen matrices are formed by taking matrix elements of arbitrary linear operators.

This code, unsurprisingly, depends on the Eigen library. If you already have it downloaded, great! If you don't, run the script:

```
./download_eigen.sh [directory to put it in]
```

And it will use wget to download it for you. If you do not pass a command line argument, it will simply download Eigen to the current directory. The hard-coded URL in that shell script may become out of date if/when the developers release a new version. Right now it downloads Eigen 3.3.4. If they do release a new version, let me know (see my e-mail below) and I'll update this repo.

By default, the Makefile in each directory will assume the downloaded Eigen library lives in the root directory of this repository, i.e., the same directory as this README. If you already have Eigen downloaded, or you chose to download it to a different directory, you'll have to update the path after the "-I" flag in the Makefile to the appropriate absolute or relative directory. 

I've included reference code for real symmetric, real indefinite, complex Hermitian, and complex indefinite matrices. Eigen can perform different optimizations in each case, so be sure to refer to the most specific reference code for your problem domain! Of course, it's possible to include more than one type of eigensolver in your code... just be warned that the executables get huge in that case.

Only remark of note: Eigen depends on compiling against at least C++11.

If you have any questions, feel free to kick me an e-mail at weinbe2@bu.edu or at evansweinberg@gmail.com


