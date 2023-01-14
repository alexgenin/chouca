// 
// This is a header file that defines a few things needed for compilation to make it on 
// all platforms. 
//

// We need this on windows for some reason as we cannot declare arma::Mat<ushort> on that
// platform 
#ifdef _WIN32
typedef short ushort; 
#endif
