#include "Util.h"

using namespace std;
namespace util {

// About mulAdd below (see Util.h): mulAdd(a,b,c) is a*b+c, as a fused multiply-add (a single rounding)
// on 64-bit ARM and as a multiplication followed by an addition on other processors. These are exactly
// the places where the earlier versions of Chrom3D, built with Clang on Apple Silicon (e.g. with the
// Xcode compiler), used fused multiply-adds, because Clang fuses a*b+c written in one expression when
// the processor has the instruction. Builds for x86-64 (e.g. with g++ on Linux) used none. Compilers
// fuse in different places on their own (GCC also across statements), so the build turns that off
// (-ffp-contract=off, see the Makefile) and the fused operations are written out here instead. That way
// the results are the same as before on each type of processor, with every compiler.

// Replaces boost::numeric::ublas::norm_2 (the Euclidean length of a vector)
static double norm_2(const CoordinateVector &vec) {
  double t = 0;
  for(uint i=0; i!=vec.size(); i++) {
    t = mulAdd(vec(i), vec(i), t); // t += vec(i) * vec(i)
  }
  return sqrt(t);
}

double euclideanDistance(CoordinateVector vec1, CoordinateVector vec2) {
  return norm_2(vec1-vec2);
}


CoordinateMatrix repeatVector(CoordinateVector vec, uint n) {  
  CoordinateMatrix mat(n, vec.size());
  for (uint i = 0; i < n; i++) {
    mat(i,0) = vec(0);
    mat(i,1) = vec(1);
    mat(i,2) = vec(2);
  }
  return mat;
}

CoordinateVector normalizeVector(CoordinateVector vec) {
  return vec / norm_2(vec);
}


CoordinateMatrix rotatePointsAroundArbitraryAxis(double theta, CoordinateVector startPoint, CoordinateVector endPoint, CoordinateMatrix mat) {
  CoordinateVector rU = normalizeVector(startPoint-endPoint); // Normalized rotation vector

  double S;
  double C;
  sinCos(theta, S, C); // S = sin(theta), C = cos(theta)
  double t = 1 - C;
  double ux = rU(0);
  double uy = rU(1);
  double uz = rU(2); 

  CoordinateMatrix rotMat(3,3);
  
  // In each entry the first product and the addition are done with mulAdd, see the top of this file:
  rotMat(0,0) = mulAdd(t*ux, ux, C);       // t*ux*ux + C
  rotMat(1,0) = mulAdd(t*ux, uy, S*uz);    // t*ux*uy + S*uz
  rotMat(2,0) = mulAdd(t*ux, uz, -(S*uy)); // t*ux*uz - S*uy
  rotMat(0,1) = mulAdd(t*ux, uy, -(S*uz)); // t*ux*uy - S*uz
  rotMat(1,1) = mulAdd(t*uy, uy, C);       // t*uy*uy + C
  rotMat(2,1) = mulAdd(t*uy, uz, S*ux);    // t*uy*uz + S*ux
  rotMat(0,2) = mulAdd(t*ux, uz, S*uy);    // t*ux*uz + S*uy
  rotMat(1,2) = mulAdd(t*uy, uz, -(S*ux)); // t*uy*uz - S*ux
  rotMat(2,2) = mulAdd(t*uz, uz, C);       // t*uz*uz + C


  // Same calculation as the earlier ublas version, trans(prod(rotMat, trans(mat - startPoint))) + startPoint,
  // with the sums done in the same order:
  CoordinateMatrix res(mat.size1(),3);
  for(uint j=0; j!=mat.size1(); j++) {
    for(uint i=0; i!=3; i++) {
      double rotated = 0; // row j of the rotated matrix at origo
      for(uint k=0; k!=3; k++) {
        rotated = mulAdd(rotMat(i,k), mat(j,k) - startPoint(k), rotated); // rotated += rotMat(i,k) * (mat(j,k) - startPoint(k))
      }
      res(j,i) = rotated + startPoint(i); // rotated matrix transformed back to original coordinate system
    }
  }
  return res;
}


// unif_real and unif_int use the same algorithms as boost::random::uniform_real_distribution and
// boost::random::uniform_int_distribution (used in earlier versions), so a given seed gives the
// same random numbers as before. (The std:: distributions are not guaranteed to do that, and give
// different numbers with different compilers.)

double unif_real(double min, double max, ENG &eng) {
  for(;;) {
    double numerator = static_cast<double>(eng() - ENG::min());
    double divisor = static_cast<double>(ENG::max() - ENG::min()) + 1;
    double result = mulAdd(numerator / divisor, max - min, min); // numerator / divisor * (max - min) + min
    if(result < max) return result;
  }
}

int unif_int(uint min, uint max, ENG &eng) {
  int minValue = min;
  int maxValue = max;
  assert(minValue <= maxValue);
  const uint32_t range = uint32_t(maxValue) - uint32_t(minValue);
  const uint32_t engineRange = ENG::max() - ENG::min();
  assert(range < engineRange);
  if(range == 0) {
    return minValue;
  }
  // Split the generator's output into range+1 equally sized buckets; reject numbers beyond the last bucket.
  uint32_t bucketSize = engineRange / (range+1);
  if(engineRange % (range+1) == range) {
    bucketSize++;
  }
  for(;;) {
    uint32_t result = (eng() - ENG::min()) / bucketSize;
    if(result <= range) return minValue + int(result);
  }
}

bool coordClash(double x1, double y1, double z1, double x2, double y2, double z2, double radius1, double radius2) {
  /*cout.precision(16);
  cout << "Is " <<  (pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2)) << " < " << pow(radius1+radius2,2) << endl;
  cout << ((pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2)) < pow(radius1+radius2,2))  << endl;*/
  return (square(x1-x2) + square(y1-y2) + square(z1-z2)) < square(radius1+radius2)-PRECISIONCONST;
}

bool coordClash(CoordinateVector vec1, CoordinateVector vec2, double radius1, double radius2) {
  return coordClash(vec1(0), vec1(1), vec1(2), vec2(0), vec2(1), vec2(2), radius1,radius2); 
}
  
  
// Does the bead at vec clash with any of the first nRows beads in mat? (Replaces passing a
// ublas::project of the first rows.)
bool coordClash(CoordinateVector vec, const CoordinateMatrix &mat, uint nRows, double radius, const vector<double> &radii) {
  for(uint i=0; i!=nRows; i++) {

    if(coordClash(vec(0), vec(1), vec(2), mat(i,0), mat(i,1), mat(i,2), radius, radii[i])) {
      return true;
    }
  }
  return false;
}


#if defined(__APPLE__)
// sin() and cos() in functions of their own, so that the compiler cannot merge the two calls into
// __sincos_stret (which sometimes differs in the last bit). See createSelfAvoidingWalk.
__attribute__((noinline)) static double sinNotMerged(double x) { return sin(x); }
__attribute__((noinline)) static double cosNotMerged(double x) { return cos(x); }
#endif

CoordinateMatrix createSelfAvoidingWalk(vector<double> radii, CoordinateVector origin, ENG& eng, uint maxIter) {
  assert(radii.size() > 0);
  double prevRadius=radii[0];
  CoordinateVector prevCoords = origin;
  CoordinateMatrix newCoords(radii.size(), 3); newCoords(0,0) = origin(0);   newCoords(0,1) = origin(1); newCoords(0,2) = origin(2);
  uint j;
  
  for(uint i=1; i!=radii.size(); i++) {
    double R = radii[i] + prevRadius;
    assert(R != 0);
    double x;
    double y;
    double z;
    for(j=0; j != maxIter; j++) {
      double theta = unif_real(0, 2*PI, eng);
      z = unif_real(-R,R, eng);
      // sin and cos as computed here by the earlier versions: on macOS, Clang used separate sin() and cos()
      // calls (not __sincos_stret as in rotatePointsAroundArbitraryAxis). A difference in the last bit can
      // change the results, because each new bead is placed exactly touching the previous one, so the
      // clash test below is decided by rounding. Elsewhere compilers merge the calls as before (GCC: sincos).
#if defined(__APPLE__)
      double sinTheta = sinNotMerged(theta);
      double cosTheta = cosNotMerged(theta);
#else
      double sinTheta = sin(theta);
      double cosTheta = cos(theta);
#endif
      x = sqrt(square(R) - square(z)) * cosTheta; // OBS: Speed up by doing sqrt only once??
      y = sqrt(square(R) - square(z)) * sinTheta; // Archimedes method (see http://tinyurl.com/nquq2d9)

      x += prevCoords(0);
      y += prevCoords(1);
      z += prevCoords(2);

      CoordinateVector proposedCoords(3); proposedCoords(0) = x; proposedCoords(1) = y; proposedCoords(2) = z;

      cout.precision(15);

      if(not coordClash(proposedCoords, newCoords, i, radii[i], radii)) { // compare with the i beads placed so far
	       break;
      }
    }
    
    if (j==maxIter) {
      // reached the end of maxIter without finding a structure that does not clash
      CoordinateMatrix newCoordsEmpty(0,0);
      return newCoordsEmpty;
    }
    
    prevCoords(0) = x; prevCoords(1) = y;  prevCoords(2) = z;
    prevRadius = radii[i];
    newCoords(i,0) = x;
    newCoords(i,1) = y;
    newCoords(i,2) = z;    
  }
  return newCoords;
}


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
  
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::vector<double> &splitDbl(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss(s);
    std::string item;

    while (std::getline(ss, item, delim)) {
        elems.push_back(atof(item.c_str()));
    }
    return elems;
}

std::vector<double> splitDbl(const std::string &s, char delim) {
    std::vector<double> elems;
    splitDbl(s, delim, elems);
    return elems;
}

std::map<std::string, std::string> makeMap(std::vector<std::string> vec1, std::vector<std::string> vec2) { // Make a template function here!
  //assert(vec1.size() == vec2.size());
  if(vec1.size() != vec2.size())
  {
    throw std::runtime_error("Header is missing in the GTrack file");
  }
  std::map<string,string> res;
  for(uint i=0; i != vec1.size(); i++) {
    res[vec1[i]] = vec2[i];
  }
  return res;
}

std::string errorLine(uint line_counter, std::string filename)
{
  std::ostringstream oss;
  oss << "Error in the line #" << line_counter << " of the GTrack file " << filename << endl;
  std::string error_line = oss.str();
  return error_line;	
}
  

  
MoveException::MoveException() { };


  
}

