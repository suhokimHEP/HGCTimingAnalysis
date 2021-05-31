#include "HGCTimingAnalysis/HGCTiming/plugins/SupportScintillator.h"


SupportScintillator::SupportScintillator():
  type_m(0), layer_m(0), ring_m(0), iphi_m(0)
{};

void SupportScintillator::setValues(int type, int layer, float posX, float posY, float posZ, float posEta, float posPhi){
  //  if(type == DetId::HGCalHSc) 
  type_m = type;
  layer_m = layer - (28);

  float cellSide = (2 * pigreco / 288); 

  float posR0 = 103.783;
  float etaZero = asinh(posZ/posR0);
  float Ring = std::abs((posEta - etaZero) / cellSide) + 1;

  float RingX = 43;
  float posR = sqrt(posX*posX + posY*posY);
  for(int ij=0; ij<43; ++ij){
    if(radiiLimits[ij] > posR * 10.){
      RingX = ij;
      //      std::cout << " radiiLimits[ij] = " << radiiLimits[ij] << " posR = " << posR << " Ring = " << Ring << " RingX = " << RingX << std::endl;
      break;
    }
  }

  int iSide = (posZ > 0) ? 1 : -1; 
  if((RingX - int(Ring) == 2 || RingX - int(Ring) == 3) && Ring-int(Ring) > 0.9 ) Ring += 1.;

  int iRing = iSide * int(Ring);


  float correctPhi = (posZ > 0) ? posPhi : (pigreco - posPhi);
  if(posPhi < 0 && posZ > 0) correctPhi += 2*pigreco;
  iphi_m = correctPhi/cellSide + 1.;

  //  std::cout << " computing ring_m = " << int(Ring) << " iphi_m = " << iphi_m  << " layer_m = " << layer_m << " iSide = " << iSide << std::endl;

  myDet_m = HGCScintillatorDetId(type_m, layer_m, iRing, iphi_m, false, 0);
  iphi_m = myDet_m.iphi();
  ring_m = myDet_m.ring();

  // std::cout << " post computing ring_m = " << myDet_m.ring() << " iradius = " << myDet_m.iradius() << " iphi_m = " << myDet_m.iphi()  
  // 	    << " layer_m = " << myDet_m.layer() 
  // 	    << " zSide = " << myDet_m.zside() << std::endl;
}



int SupportScintillator::get_Ring(){
  return myDet_m.ring();
}

int SupportScintillator::get_iPhi(){
  return myDet_m.iphi();
}


void SupportScintillator::get_ScintID(HGCScintillatorDetId& detId){
  detId = myDet_m;
  return;
}
