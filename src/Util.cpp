#include "Util.h"

using namespace std;
namespace util {

double euclideanDistance(CoordinateVector vec1, CoordinateVector vec2) {
  return ublas::norm_2(vec1-vec2);
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
  return vec / boost::numeric::ublas::norm_2(vec);
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


  CoordinateMatrix res(mat.size1(),3);
  
  res = ublas::trans(ublas::prod(rotMat, ublas::trans(mat-repeatVector(startPoint,mat.size1())))); // rotated matrix at origo

  return res + repeatVector(startPoint,mat.size1()); // rotated matrix transformed back to original coordinate system
}


double unif_real(double min, double max, ENG &eng) {
  boost::random::uniform_real_distribution<double> urd = boost::random::uniform_real_distribution<double>(min,max);
  return urd(eng);
}

int unif_int(uint min, uint max, ENG &eng) {
  boost::random::uniform_int_distribution<int> urd = boost::random::uniform_int_distribution<int>(min,max);
  return urd(eng);
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
  
  
bool coordClash(CoordinateVector vec, CoordinateMatrix mat, double radius, vector<double> radii) {
  for(uint i=0; i!=mat.size1(); i++) {
    
    if(coordClash(vec,ublas::row(mat,i), radius, radii[i])) {
      return true;
    }
  }
  return false;
}



CoordinateMatrix createSelfAvoidingWalk(vector<double> radii, CoordinateVector origin, ENG& eng, uint maxIter) {
  assert(radii.size() > 0);
  CoordinateVector tmpRadii;
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
      vector<double> tmpRadii(radii.begin(), radii.begin()+i+1);
     
      cout.precision(15);
      CoordinateVector newc = ublas::row(newCoords,i);

      if(not coordClash(proposedCoords, project(newCoords, boost::numeric::ublas::range(0, i), boost::numeric::ublas::range(0,3)), radii[i], tmpRadii)) {
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
  oss << "Error in the line #" << boost::lexical_cast<std::string>(line_counter) << " of the GTrack file " << filename << endl;
  std::string error_line = oss.str();
  return error_line;	
}
  

  
MoveException::MoveException() { };


  
}

