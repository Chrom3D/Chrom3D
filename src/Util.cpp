#include "Util.h"

using namespace std;
namespace util {

// Replaces boost::numeric::ublas::norm_2 (the Euclidean length of a vector)
static double norm_2(const CoordinateVector &vec) {
  double t = 0;
  for(uint i=0; i!=vec.size(); i++) {
    t += vec(i) * vec(i);
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

  double C = cos(theta);
  double S = sin(theta);
  double t = 1 - C;
  double ux = rU(0);
  double uy = rU(1);
  double uz = rU(2); 

  CoordinateMatrix rotMat(3,3);
  
  rotMat(0,0) = t*ux*ux + C;
  rotMat(1,0) = t*ux*uy + S*uz;
  rotMat(2,0) = t*ux*uz - S*uy;
  rotMat(0,1) = t*ux*uy - S*uz;
  rotMat(1,1) = t*uy*uy + C;
  rotMat(2,1) =  t*uy*uz + S*ux;
  rotMat(0,2) = t*ux*uz + S*uy;
  rotMat(1,2) = t*uy*uz - S*ux;
  rotMat(2,2) =  t*uz*uz + C;


  // Same calculation as the earlier ublas version, trans(prod(rotMat, trans(mat - startPoint))) + startPoint,
  // with the sums done in the same order:
  CoordinateMatrix res(mat.size1(),3);
  for(uint j=0; j!=mat.size1(); j++) {
    for(uint i=0; i!=3; i++) {
      double rotated = 0; // row j of the rotated matrix at origo
      for(uint k=0; k!=3; k++) {
        rotated += rotMat(i,k) * (mat(j,k) - startPoint(k));
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
    double result = numerator / divisor * (max - min) + min;
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
  return (pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2)) < pow(radius1+radius2,2)-PRECISIONCONST;
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
      x = sqrt(pow(R,2) - pow(z,2)) * cos(theta); // OBS: Speed up by doing sqrt only once??
      y = sqrt(pow(R,2) - pow(z,2)) * sin(theta); // Archimedes method (see http://tinyurl.com/nquq2d9)

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

