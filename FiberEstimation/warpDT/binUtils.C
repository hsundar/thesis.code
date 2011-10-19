#include "binUtils.h"


//namespace binOp {

bool isPowerOfTwo(unsigned int n) {
  return !(n & (n - 1)) && n;
}

// compute the next highest power of 2 of 32-bit v
int getNextHighestPowerOfTwo(unsigned int n) {
  unsigned int v = n;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

// compute the prev highest power of 2 of 32-bit v
int getPrevHighestPowerOfTwo(unsigned int n) {
  unsigned int v = n;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return  v >> 1;
}

//}//end namespace

