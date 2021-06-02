#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"

// user include files
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>
#include <vector>
#include <string>

class SupportScintillator {
 public:

  SupportScintillator();

  void setValues(int type, int layer, float posX, float posY, float posZ, float posEta, float posPhi);
  int get_Ring();
  int get_iPhi();
  void get_ScintID(HGCScintillatorDetId& detId);

 private:
  int type_m;
  int layer_m;
  int ring_m;
  int iphi_m;
  bool trigger_m;
  int sipm_m;
  HGCScintillatorDetId myDet_m;
  const float pigreco = 3.1415926536;
  const float radiiMax[42] = {1060.98, 1084.12, 1108.30,
			      1132.48, 1157.73, 1182.99, 1209.37,
			      1235.75, 1263.31, 1290.87, 1319.66, 
			      1348.45, 1378.52, 1408.59, 1440.00, 
			      1471.42, 1504.23, 1537.05, 1571.33, 
			      1605.60, 1641.41, 1677.22, 1714.62, 
			      1752.03, 1791.10, 1830.17, 1870.99, 
			      1911.80, 1954.44, 1997.08, 2041.61, 
			      2086.15, 2132.67, 2179.20, 2227.80, 
			      2276.40, 2327.16, 2377.93, 2430.96, 
			      2483.99, 2539.39, 2594.79};



  const float radiiMin[42] = {1037.83, 1060.98, 1084.12, 
			      1108.30, 1132.48, 1157.73, 1182.99, 
			      1209.37, 1235.75, 1263.31, 1290.87, 
			      1319.66, 1348.45, 1378.52, 1408.59, 
			      1440.00, 1471.42, 1504.23, 1537.05, 
			      1571.33, 1605.60, 1641.41, 1677.22, 
			      1714.62, 1752.03, 1791.10, 1830.17,
			      1870.99, 1911.80, 1954.44, 1997.08,
			      2041.61, 2086.15, 2132.67, 2179.20,
			      2227.80, 2276.40, 2327.16, 2377.93, 
			      2430.96, 2483.99, 2539.39};
};
