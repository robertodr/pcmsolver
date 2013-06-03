#ifndef WAVELETCAVITY
#define WAVELETCAVITY

#include <iostream>
#include <string>

#include "Config.h"

#include <Eigen/Dense>

#include "vector3.h"
#include "Getkw.h"
#include "Cavity.h"

class WaveletCavity : public Cavity {
 public:
    WaveletCavity(){}
    //    WaveletCavity(string &filename);
    WaveletCavity(const Getkw & input, const string path = "Cavity");
    WaveletCavity(const Section & cavity);
    WaveletCavity(const std::vector<Sphere> & _spheres, double _probeRadius, int _patchLevel = 2, double _coarsity = 0.5) :
	   Cavity(_spheres), probeRadius(_probeRadius), patchLevel(_patchLevel), coarsity(_coarsity) 
           {
		uploadedDyadic = false;
		makeCavity();
           }
    ~WaveletCavity(){};
    static Cavity* Create(const std::vector<Sphere> & _spheres, double _area, double _probeRadius = 0.0, 
		    bool _addSpheres = false, int _patchLevel = 2, double _coarsity = 0.5)
    {
	    return new WaveletCavity(_spheres, _probeRadius, _patchLevel, _coarsity);
    }
    static bool registered;
    static bool Register();
    void makeCavity();
    void readCavity(const string & filename);
    void uploadPoints(int quadLevel, vector3 **** T_, bool isPWL);
//    VectorXd & getTessRadius(){return tessRadius;};
//    VectorXd & getSphereRadius(){return sphereRadius;};
//    int getNSpheres(){return nSpheres;};
//    Matrix3Xd & getSphereCenter(){return sphereCenter;};
//    Matrix3Xd & getTessSphereCenter(){return tessSphereCenter;};
//    double getTessRadius(int i){return tessRadius(i);};
    unsigned int getNPatches(){return nPatches;}
    unsigned int getNLevels(){return nLevels;}
    unsigned int getNPoints(){return nPoints;}
    Vector3d getNodePoint(int i){return nodePoint[i];}
    Vector3i getNodeIndex(int i){return nodeIndex[i];}
    friend std::ostream& operator<<(std::ostream &o, const WaveletCavity &c);
    void compFakePotential();
 private:
    void uploadPointsPWC(int quadLevel, vector3 **** T_);
    void uploadPointsPWL(int quadLevel, vector3 **** T_);
    std::vector<Vector3d> nodePoint;
    std::vector<Vector3i> nodeIndex;
    unsigned int nPatches;
    unsigned int nLevels;
    unsigned int nPoints;
    bool uploadedDyadic;
    void writeInput(string &fileName);
 //   int nSpheres;
 //   Matrix3Xd sphereCenter;
 //   Matrix3Xd tessSphereCenter;
 //   VectorXd sphereRadius;
 //   VectorXd tessRadius;
    int patchLevel;
    double probeRadius;
    double coarsity;
};

#ifndef REGISTER
#define REGISTER static const bool anonW = WaveletCavity::Register();
REGISTER
#endif
#undef REGISTER

#endif
