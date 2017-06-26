#ifndef UtilClasses_h
#define UtilClasses_h

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include "TH1.h"
#include "TCanvas.h"

/* #include "RecoNtuples/HGCalAnalysis/interface/AEvent.h" */
/* #include "RecoNtuples/HGCalAnalysis/interface/AObData.h" */

//#include "RecoNtuples/HGCalAnalysis/interface/NumberToString.h"

using namespace std;

/*
float deltaPhi( float phi1, float phi2) {
        
    float dPhi = phi1 - phi2;
    float pi = 3.14159265;
    if     ( dPhi <=-pi) dPhi += 2.0*pi;
    else if( dPhi >  pi) dPhi -= 2.0*pi;
        
    return dPhi;
}    
*/

/* namespace edm { */
/*   class Event; */
/*   class EventSetup; */
/* } */

class UtilClasses {    
 public:
  UtilClasses(float genEta, float genPhi);
  ~UtilClasses();
  
  
  float deltaEta(float eta1, float eta2);
  float deltaR(float eta1, float eta2, float phi1, float phi2);
  float deltaX(float x1, float x2, float y1, float y2);
  
  // need correct sign for z
  float etaPhiZtoX( float eta, float phi, float z );
  // need correct sign for z
  float etaPhiZtoY( float eta, float phi, float z );
  // go from HGC layer to z in cm
  float layerToZ( int layer, float eta );
  int Ztolayer(float z, float eta );
  
  float dsGenRecHit(float genEta, float genPhi, int recHitLayer, float recHitX, float recHitY );
  float dsGenRecoObj(float genEta, float genPhi, float recoObjZ, float recoObjX, float recoObjY );

  void layerIntersection(std::array<double,3> &to, const std::array<double,3> &from, const std::array<double,3> &fromB) const;
  
  vector<TH1*> initHists( string name );
  
  void writeHists( vector<TH1*> hists );
  void drawAndSaveHists( vector<TH1*> hists, TCanvas *c, bool doLogs = false );


 private:
  float _genEta;
  float _genPhi;

};

#endif
