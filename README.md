# eigens-with-eigen
A set of reference code showing how to calculate eigenvalues and eigenvectors with Eigen. Dense Eigen matrices are formed by taking matrix elements of arbitrary linear operator.

This code, unsurprisingly, depends on the Eigen library. If you already have it downloaded, great! If you don't, run the script:

./download_eigen.sh [directory to put it in]

And it'll use wget to download it for you. The hard-coded URL in that shell script may become out of date if/when the developers release a new version. Right now it downloads Eigen 3.3.4. If they do release a new version, let me know (see my e-mail below) and I'll update this repo.

Anyway, once you're done with that, pop the absolute or relative directory into a Makefile in the various directories as appropriate. Hopefully the directory names are self explanatory. 

Only remark of note: Eigen depends on compiling against at least C++11.

If you have any questions, feel free to kick me an e-mail at weinbe2@bu.edu or at evansweinberg@gmail.com


