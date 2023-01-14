// 
// This is a header file that defines a few things needed for compilation to make it on 
// all platforms. 
//

// We need this on windows for some reason to be able to declare arma::Mat with unsigned 
// shorts, on linux it works well. 
typedef unsigned short ushort; 
