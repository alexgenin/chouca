//
// This file defines some constants that are useful for all code
//

// Some typedefs for better legibility. ushort typically contains an integer
// for the state number, which is often way lower than its maximum value (2^16). We
// cannot use a shorter type here (e.g. unsigned char) because matrices of uchars are not
// supported by armadillo. We use c-style array in the 'compiled' backend for better
// performance which do use those, but not here because it adds complexity to the code.
//
// Note that in function signatures we can't use this typedef because they are copied
// by Rcpp in the RcppExports.cpp file, but without typedef and this makes the
// compilation fail on Windows.
typedef unsigned short ushort;

// We need to define the types of probability normalization functions as
// integers here
constexpr ushort NORMF_IDENTITY = 0;
constexpr ushort NORMF_SOFTMAX  = 1;

