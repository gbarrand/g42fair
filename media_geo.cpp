// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file exlib.license for terms.

//exlib_build_use geant4 inlib

#include <G4Material.hh>
#include <G4UnitsTable.hh>

//////////////////////////////////////////////////////////////////////////////
/// code from inlib : ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//#include <inlib/sjust>
//#include <inlib/mnmx>

#include <vector>
#include <string>

namespace inlib {

enum what { leading, trailing, both };

inline bool strip(std::string& a_string,what a_type = both,char a_char = ' '){
  //return true = some stripping had been done.

  std::string::size_type l = a_string.length();
  if(l==0) return false;

  switch ( a_type ) {
  case leading:{
    char* pos = (char*)a_string.c_str();
    for(std::string::size_type i=0;i<l;i++,pos++) {
      if(*pos!=a_char) {
        a_string = a_string.substr(i,l-i);
        return (i?true:false); //i=0 : same string.
      }
    }
    // all chars are a_char :
    a_string.clear();
    return true;
    }break;
  case trailing:{
    char* pos = (char*)a_string.c_str();
    pos += (l-1);
    std::string::size_type i = l-1;
    std::string::const_reverse_iterator it;
    for(it=a_string.rbegin();it!=a_string.rend();++it,i--,pos--) {
      if(*pos!=a_char) {
        a_string = a_string.substr(0,i+1);
        return (i==(l-1)?false:true); //i==l-1 : same string.
      }
    }
    // all chars are a_char :
    a_string.clear();
    return true;
    }break;
  case both:{
    bool stat_lead = strip(a_string,leading,a_char);
    bool stat_trail = strip(a_string,trailing,a_char);
    if(stat_lead) return true;
    if(stat_trail) return true;
    }break;
  //default:break;
  }
  return false; //nothing done.
}


enum side { side_left, side_right, side_middle }; //have side_ for SWIG.

inline void justify(std::string& a_string,size_t a_size,side a_side = side_left,char a_c = ' '){
  // a_size is the final string length.
  strip(a_string);
  if(a_size<=a_string.size()) {
    a_string.resize(a_size);
  } else {
    if(a_side==side_left) {
      a_string = a_string + std::string(a_size-a_string.size(),a_c);
    } else if(a_side==side_right) {
      a_string = std::string(a_size-a_string.size(),a_c) + a_string;
    } else if(a_side==side_middle) {
      size_t l = a_size - a_string.size();
      size_t h = l/2;
      if(h*2==l) { //even number of spaces :
        a_string = std::string(h,a_c) + a_string + std::string(h,a_c);
      } else { // odd number of spaces :
        a_string = std::string(h,a_c) + a_string + std::string(h+1,a_c);
      }
    }
  }
}

//template <class T>
//inline T mn(const T& a,const T& b) {return (a<b?a:b);}    

template <class T>
inline T mx(const T& a,const T& b) {return (a>b?a:b);}    

}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>

inline void write_media_geo(const G4Material& a_material,
                            size_t a_npckov,
                            const G4double a_ENERGY[],       //ppckov
                            const G4double a_ABSORPTION[],   //absco
                            const G4double a_EFFICIENCY[],   //effic
                            const G4double a_RINDEX1[]       //rindex
			    ) {

  // first line :
  // <string19:name> <int:ncomponent> <read:ca[ncomponent]> <read:cz[ncomponent] <real:density> <real:cw[ncomponent]>
  // with :
  //   ca : mass number
  //   cz : atomic number
  //   cw : weight of the component in a mixture.
  //   density in g/cm3.
  std::string name = a_material.GetName();
  inlib::justify(name,19);
  std::cout << name;
  size_t nelem = a_material.GetNumberOfElements();
  std::cout << nelem << " ";
 {for(size_t ielem=0;ielem<nelem;ielem++) {
    const G4Element* elem = a_material.GetElement(ielem);
    std::cout << elem->GetN() << " ";
  }}
 {for(size_t ielem=0;ielem<nelem;ielem++) {
    const G4Element* elem = a_material.GetElement(ielem);
    std::cout << elem->GetZ() << " ";
  }}
  std::cout << a_material.GetDensity()/(CLHEP::g/CLHEP::cm3) << " ";
  const G4double* fracs = a_material.GetFractionVector();
 {for(size_t ielem=0;ielem<nelem;ielem++) {
    std::cout << fracs[ielem] << " ";
  }}
  std::cout << std::endl;

  // second line :
  // <string19:spaces> <int:sensFlag> <int:fldFlag> <real:fld> <real:epsil>
  // with :
  //   sensFlag : sensitivity flag.
  //   fldFlag : field flag.
  //   fld : Maximum field value in kilogauss
  //   epsil : Boundary crossing precision (unit ?).
  std::string empty19;
  inlib::justify(empty19,19);
  std::cout << empty19 << "0 0 20. .001" << std::endl;

  // third line :
  // <npckov>
  // following <npckov> lines :
  // for i in npckov
  //   <string19:spaces> <real:ppckov[i]:eV> <real:absco[i]:mm> <real:effic[i]:mm> <real:rindex[i]>
  // first pass to get max string length for numbers :
  size_t smx = 0;
  std::ostringstream oss;
 {for(size_t ipckov=0;ipckov<a_npckov;ipckov++) {
    oss.str("");     
    oss << a_ENERGY[ipckov]/CLHEP::eV;
    smx = inlib::mx<size_t>(smx,oss.str().size());
    
    oss.str("");     
    oss << a_ABSORPTION[ipckov]/CLHEP::mm;
    smx = inlib::mx<size_t>(smx,oss.str().size());
    
    oss.str("");     
    oss << a_EFFICIENCY[ipckov]/CLHEP::mm;
    smx = inlib::mx<size_t>(smx,oss.str().size());
    
    oss.str("");     
    oss << a_RINDEX1[ipckov];
    smx = inlib::mx<size_t>(smx,oss.str().size());
  }}
  smx += 2;
  
  std::cout << empty19 << a_npckov << std::endl;
 {std::string svalue;
  for(size_t ipckov=0;ipckov<a_npckov;ipckov++) {
    std::cout << empty19;

   {oss.str("");     
    oss << a_ENERGY[ipckov]/CLHEP::eV;
    svalue = oss.str();
    inlib::justify(svalue,smx);
    std::cout << svalue;}
    
   {oss.str("");     
    oss << a_ABSORPTION[ipckov]/CLHEP::mm;
    svalue = oss.str();
    inlib::justify(svalue,smx);
    std::cout << svalue;}
    
   {oss.str("");     
    oss << a_EFFICIENCY[ipckov]/CLHEP::mm;
    svalue = oss.str();
    inlib::justify(svalue,smx);
    std::cout << svalue;}
    
   {oss.str("");     
    oss << a_RINDEX1[ipckov];
    svalue = oss.str();
    inlib::justify(svalue,smx);
    std::cout << svalue;}
    
   std::cout << std::endl;
  }}
  
}

using namespace CLHEP;

inline void MEMPHYS_DetectorConstruction_ConstructMaterials_Water() {

  G4double density,a;
  
  //---Water
  
  a = 1.01*g/mole;
  G4Element* elH = new G4Element("Hydrgen","H", 1,a);
  
  a = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O", 8,a);
  
  density = 1.00*g/cm3;
  G4Material* Water = new G4Material("MEMPHYS_Water",density,2);
  Water->AddElement(elH, 2);
  Water->AddElement(elO, 1);

  //---Black sheet
  /*
  a = 12.01*g/mole;
  G4Element* elC = new G4Element("Carbon","C", 6,a);
  
  density = 0.95*g/cm3;
  G4Material* Blacksheet = new G4Material("MEMPHYS_Blacksheet",density,2);
  Blacksheet->AddElement(elC, 1);
  Blacksheet->AddElement(elH, 2);
  */
  
  // -------------------------------------------------------------
  // Generate & Add Material Properties Table
  // -------------------------------------------------------------
  
   //From SFDETSIM water absorption
  const G4int NUMENTRIES_water=60;
  
  G4double ENERGY_water[NUMENTRIES_water] =  //G.Barrand : 4*15 = 60. Ok.
    { 1.56962e-09*GeV, 1.58974e-09*GeV, 1.61039e-09*GeV, 1.63157e-09*GeV, 
      1.65333e-09*GeV, 1.67567e-09*GeV, 1.69863e-09*GeV, 1.72222e-09*GeV, 
      1.74647e-09*GeV, 1.77142e-09*GeV, 1.7971e-09*GeV,  1.82352e-09*GeV, 
      1.85074e-09*GeV, 1.87878e-09*GeV, 1.90769e-09*GeV, 1.93749e-09*GeV, 
      1.96825e-09*GeV, 1.99999e-09*GeV, 2.03278e-09*GeV, 2.06666e-09*GeV,
      2.10169e-09*GeV, 2.13793e-09*GeV, 2.17543e-09*GeV, 2.21428e-09*GeV, 
      2.25454e-09*GeV, 2.29629e-09*GeV, 2.33962e-09*GeV, 2.38461e-09*GeV, 
      2.43137e-09*GeV, 2.47999e-09*GeV, 2.53061e-09*GeV, 2.58333e-09*GeV, 
      2.63829e-09*GeV, 2.69565e-09*GeV, 2.75555e-09*GeV, 2.81817e-09*GeV, 
      2.88371e-09*GeV, 2.95237e-09*GeV, 3.02438e-09*GeV, 3.09999e-09*GeV,
      3.17948e-09*GeV, 3.26315e-09*GeV, 3.35134e-09*GeV, 3.44444e-09*GeV, 
      3.54285e-09*GeV, 3.64705e-09*GeV, 3.75757e-09*GeV, 3.87499e-09*GeV, 
      3.99999e-09*GeV, 4.13332e-09*GeV, 4.27585e-09*GeV, 4.42856e-09*GeV, 
      4.59258e-09*GeV, 4.76922e-09*GeV, 4.95999e-09*GeV, 5.16665e-09*GeV, 
      5.39129e-09*GeV, 5.63635e-09*GeV, 5.90475e-09*GeV, 6.19998e-09*GeV };

  // M Fechner : new ; define the water refraction index using refsg.F 
  //from skdetsim using the whole range.   
  G4double RINDEX1[NUMENTRIES_water] =  //G.Barrand : 7*7+6+5 = 49+11 = 60. Ok.
    { 1.32885, 1.32906, 1.32927, 1.32948, 1.3297,  1.32992, 1.33014, 
      1.33037, 1.3306,  1.33084, 1.33109, 1.33134, 1.3316,  1.33186, 1.33213,
      1.33241, 1.3327,  1.33299, 1.33329, 1.33361, 1.33393, 1.33427, 1.33462,
      1.33498, 1.33536, 1.33576, 1.33617, 1.3366,  1.33705, 1.33753, 1.33803,
      1.33855, 1.33911, 1.3397,  1.34033, 1.341,   1.34172, 1.34248, 1.34331,
      1.34419, 1.34515, 1.3462,  1.34733, 1.34858, 1.34994, 1.35145, 1.35312,
      1.35498, 1.35707, 1.35943, 1.36211, 1.36518, 1.36872, 1.37287, 1.37776,
      1.38362, 1.39074, 1.39956, 1.41075, 1.42535 };
  
  G4double ABSORPTION_water[NUMENTRIES_water] =  //G.Barrand : 5*11+4+1 = 55+5 = 60. Ok.
    {25.3504*cm, 31.7938*cm, 39.9915*cm, 50.454*cm,  63.85*cm, 
     81.0584*cm, 103.24*cm,  131.93*cm,  169.172*cm, 217.694*cm, 
     224.921*cm, 249.688*cm, 262.674*cm, 273*cm,     321.13*cm, 339.789*cm,
     351.617*cm, 363.108*cm, 385.802*cm, 461.042*cm, 707.714*cm, 
     1038.42*cm, 1383.7*cm,  1558.36*cm, 1722.65*cm, 1939.11*cm, 
     2092.49*cm, 2240.14*cm, 2962.96*cm, 4967.03*cm, 6368.58*cm, 
     8207.56*cm, 10634.2*cm, 13855.3*cm, 18157.3*cm, 23940.2*cm, 
     31766*cm,   42431.6*cm, 57074.9*cm, 77335.9*cm, 105598*cm, 
     145361*cm,  192434*cm,  183898*cm,  176087*cm,  168913*cm, 162301*cm, 
     156187*cm,  150516*cm,  145243*cm,  140327*cm,  135733*cm, 131430*cm, 
     127392*cm,  123594*cm,  120016*cm,  116640*cm,  113448*cm, 110426*cm, 
     107562*cm};
  
  // M Fechner: Rayleigh scattering -- as of version 4.6.2 of GEANT,
  // one may use one's own Rayleigh scattering lengths (the buffer is no
  // longer overwritten for "water", see 4.6.2 release notes)
  
  // RAYFF = 1/ARAS, for those of you who know SKdetsim...
  // actually that's not quite right because the scattering models
  // are different; in G4 there is no scattering depolarization
  // std value at SK = 0.6. But Mie scattering is implemented
  // in SKdetsim and not in G4
  
  // april 2005 : reduced reflections, let's increase scattering...
  //   G4double RAYFF = 1.0/1.65;
  G4double RAYFF = 1.0/1.5;
  
  G4double RAYLEIGH_water[NUMENTRIES_water] = {  //G.Barrand : 3*20 = 60. Ok.
    167024.4*cm*RAYFF, 158726.7*cm*RAYFF, 150742*cm*RAYFF,
    143062.5*cm*RAYFF, 135680.2*cm*RAYFF, 128587.4*cm*RAYFF,
    121776.3*cm*RAYFF, 115239.5*cm*RAYFF, 108969.5*cm*RAYFF,
    102958.8*cm*RAYFF, 97200.35*cm*RAYFF, 91686.86*cm*RAYFF,
    86411.33*cm*RAYFF, 81366.79*cm*RAYFF, 76546.42*cm*RAYFF,
    71943.46*cm*RAYFF, 67551.29*cm*RAYFF, 63363.36*cm*RAYFF,
    59373.25*cm*RAYFF, 55574.61*cm*RAYFF, 51961.24*cm*RAYFF,
    48527.00*cm*RAYFF, 45265.87*cm*RAYFF, 42171.94*cm*RAYFF,
    39239.39*cm*RAYFF, 36462.50*cm*RAYFF, 33835.68*cm*RAYFF,
    31353.41*cm*RAYFF, 29010.30*cm*RAYFF, 26801.03*cm*RAYFF,
    24720.42*cm*RAYFF, 22763.36*cm*RAYFF, 20924.88*cm*RAYFF,
    19200.07*cm*RAYFF, 17584.16*cm*RAYFF, 16072.45*cm*RAYFF,
    14660.38*cm*RAYFF, 13343.46*cm*RAYFF, 12117.33*cm*RAYFF,
    10977.70*cm*RAYFF, 9920.416*cm*RAYFF, 8941.407*cm*RAYFF,
    8036.711*cm*RAYFF, 7202.470*cm*RAYFF, 6434.927*cm*RAYFF,
    5730.429*cm*RAYFF, 5085.425*cm*RAYFF, 4496.467*cm*RAYFF,
    3960.210*cm*RAYFF, 3473.413*cm*RAYFF, 3032.937*cm*RAYFF,
    2635.746*cm*RAYFF, 2278.907*cm*RAYFF, 1959.588*cm*RAYFF,
    1675.064*cm*RAYFF, 1422.710*cm*RAYFF, 1200.004*cm*RAYFF,
    1004.528*cm*RAYFF, 833.9666*cm*RAYFF, 686.1063*cm*RAYFF
  };

  G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX", ENERGY_water, RINDEX1, NUMENTRIES_water);
  myMPT1->AddProperty("ABSLENGTH",ENERGY_water, ABSORPTION_water, NUMENTRIES_water);
  myMPT1->AddProperty("RAYLEIGH",ENERGY_water,RAYLEIGH_water,NUMENTRIES_water);
  Water->SetMaterialPropertiesTable(myMPT1);

  //G.Barrand : WARNING : below numbers, needed for media.geo, are of my own. (Not found in MEMPHYS).
  G4double EFFICIENCY[NUMENTRIES_water] = {  //G.Barrand : 10*6 = 60.Ok.
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m
  };
  write_media_geo(*Water,NUMENTRIES_water,ENERGY_water,ABSORPTION_water,EFFICIENCY,RINDEX1);
  
  /*
  G4double BLACKABS_blacksheet[NUMENTRIES_water] =
     { 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 
       1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 
       1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm};
   
  G4MaterialPropertiesTable *myMPT4 = new G4MaterialPropertiesTable();
  myMPT4->AddProperty("ABSLENGTH", ENERGY_water, BLACKABS_blacksheet, NUMENTRIES_water);
  Blacksheet->SetMaterialPropertiesTable(myMPT4);
  */
  
}

inline void WCSimDetectorConstruction_ConstructMaterials_Water() {

  G4double density;
  G4double a;

  //---Water
  
  a = 1.01*g/mole;
  G4Element* elH = new G4Element("Hydrogen","H", 1,a);
  
  a = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O", 8,a);
  
  density = 1.00*g/cm3;
  G4Material* Water = new G4Material("WCSim_Water",density,2);
  Water->AddElement(elH, 2);
  Water->AddElement(elO, 1);

  //---Gd doped Water

  a = 157.25*g/mole;
  G4Element* Gd = new G4Element("Gadolinium","Gd", 64,a);

  density = 1.00*g/cm3;
  G4Material* DopedWater = new G4Material("Doped Water",density,2);
  DopedWater->AddMaterial(Water,99.9*perCent);
  DopedWater->AddElement(Gd,0.1*perCent);

  //---Black sheet
  a = 12.01*g/mole;
  G4Element* elC = new G4Element("Carbon","C", 6,a);

  density = 0.95*g/cm3;
  G4Material* Blacksheet = new G4Material("WCSim_Blacksheet",density,2);
  Blacksheet->AddElement(elC, 1);
  Blacksheet->AddElement(elH, 2);

  //From SFDETSIM water absorption
  const G4int NUMENTRIES_water=60;

  G4double ENERGY_water[NUMENTRIES_water] =  //G.Barrand : 4*15 = 60. Ok.
     { 1.56962e-09*GeV, 1.58974e-09*GeV, 1.61039e-09*GeV, 1.63157e-09*GeV, 
       1.65333e-09*GeV, 1.67567e-09*GeV, 1.69863e-09*GeV, 1.72222e-09*GeV, 
       1.74647e-09*GeV, 1.77142e-09*GeV, 1.7971e-09*GeV,  1.82352e-09*GeV, 
       1.85074e-09*GeV, 1.87878e-09*GeV, 1.90769e-09*GeV, 1.93749e-09*GeV, 
       1.96825e-09*GeV, 1.99999e-09*GeV, 2.03278e-09*GeV, 2.06666e-09*GeV,
       2.10169e-09*GeV, 2.13793e-09*GeV, 2.17543e-09*GeV, 2.21428e-09*GeV, 
       2.25454e-09*GeV, 2.29629e-09*GeV, 2.33962e-09*GeV, 2.38461e-09*GeV, 
       2.43137e-09*GeV, 2.47999e-09*GeV, 2.53061e-09*GeV, 2.58333e-09*GeV, 
       2.63829e-09*GeV, 2.69565e-09*GeV, 2.75555e-09*GeV, 2.81817e-09*GeV, 
       2.88371e-09*GeV, 2.95237e-09*GeV, 3.02438e-09*GeV, 3.09999e-09*GeV,
       3.17948e-09*GeV, 3.26315e-09*GeV, 3.35134e-09*GeV, 3.44444e-09*GeV, 
       3.54285e-09*GeV, 3.64705e-09*GeV, 3.75757e-09*GeV, 3.87499e-09*GeV, 
       3.99999e-09*GeV, 4.13332e-09*GeV, 4.27585e-09*GeV, 4.42856e-09*GeV, 
       4.59258e-09*GeV, 4.76922e-09*GeV, 4.95999e-09*GeV, 5.16665e-09*GeV, 
       5.39129e-09*GeV, 5.63635e-09*GeV, 5.90475e-09*GeV, 6.19998e-09*GeV };


  // M Fechner : new ; define the water refraction index using refsg.F 
  //from skdetsim using the whole range.   
  G4double RINDEX1[NUMENTRIES_water] =  //G.Barrand : 7*7+6+5 = 49+11 = 60. Ok.
     {1.32885, 1.32906, 1.32927, 1.32948, 1.3297,  1.32992, 1.33014, 
      1.33037, 1.3306,  1.33084, 1.33109, 1.33134, 1.3316,  1.33186, 1.33213,
      1.33241, 1.3327,  1.33299, 1.33329, 1.33361, 1.33393, 1.33427, 1.33462,
      1.33498, 1.33536, 1.33576, 1.33617, 1.3366,  1.33705, 1.33753, 1.33803,
      1.33855, 1.33911, 1.3397,  1.34033, 1.341,   1.34172, 1.34248, 1.34331,
      1.34419, 1.34515, 1.3462,  1.34733, 1.34858, 1.34994, 1.35145, 1.35312,
      1.35498, 1.35707, 1.35943, 1.36211, 1.36518, 1.36872, 1.37287, 1.37776,
      1.38362, 1.39074, 1.39956, 1.41075, 1.42535};
   

  G4double ABWFF = 1.0;

  // Get from the tuning parameters
  //ABWFF = WCSimTuningParams->GetAbwff();
  ABWFF = 1.30;

  //T. Akiri: Values from Skdetsim 
  G4double ABSORPTION_water[NUMENTRIES_water] =  //G.Barrand : 5*12 = 60. Ok.
      {
        16.1419*cm*ABWFF,  18.278*cm*ABWFF, 21.0657*cm*ABWFF, 24.8568*cm*ABWFF, 30.3117*cm*ABWFF, 
	38.8341*cm*ABWFF, 54.0231*cm*ABWFF, 81.2306*cm*ABWFF, 120.909*cm*ABWFF, 160.238*cm*ABWFF, 
	193.771*cm*ABWFF, 215.017*cm*ABWFF, 227.747*cm*ABWFF,  243.85*cm*ABWFF, 294.036*cm*ABWFF, 
	321.647*cm*ABWFF,  342.81*cm*ABWFF, 362.827*cm*ABWFF, 378.041*cm*ABWFF, 449.378*cm*ABWFF,
        739.434*cm*ABWFF, 1114.23*cm*ABWFF, 1435.56*cm*ABWFF, 1611.06*cm*ABWFF, 1764.18*cm*ABWFF, 
	2100.95*cm*ABWFF,  2292.9*cm*ABWFF, 2431.33*cm*ABWFF,  3053.6*cm*ABWFF, 4838.23*cm*ABWFF, 
	6539.65*cm*ABWFF, 7682.63*cm*ABWFF, 9137.28*cm*ABWFF, 12220.9*cm*ABWFF, 15270.7*cm*ABWFF, 
	19051.5*cm*ABWFF, 23671.3*cm*ABWFF, 29191.1*cm*ABWFF, 35567.9*cm*ABWFF,   42583*cm*ABWFF,
        49779.6*cm*ABWFF, 56465.3*cm*ABWFF,   61830*cm*ABWFF, 65174.6*cm*ABWFF, 66143.7*cm*ABWFF,   
	  64820*cm*ABWFF,   61635*cm*ABWFF, 57176.2*cm*ABWFF, 52012.1*cm*ABWFF, 46595.7*cm*ABWFF, 
	41242.1*cm*ABWFF, 36146.3*cm*ABWFF, 31415.4*cm*ABWFF, 27097.8*cm*ABWFF, 23205.7*cm*ABWFF, 
	19730.3*cm*ABWFF, 16651.6*cm*ABWFF, 13943.6*cm*ABWFF, 11578.1*cm*ABWFF, 9526.13*cm*ABWFF
      };

   // M Fechner: Rayleigh scattering -- as of version 4.6.2 of GEANT,
   // one may use one's own Rayleigh scattering lengths (the buffer is no
   // longer overwritten for "water", see 4.6.2 release notes)

  // RAYFF = 1/ARAS, for those of you who know SKdetsim...
  // actually that's not quite right because the scattering models
  // are different; in G4 there is no scattering depolarization
  // std value at SK = 0.6. But Mie scattering is implemented
  // in SKdetsim and not in G4

   
  // april 2005 : reduced reflections, let's increase scattering...
  // sep 09: for the large detector like superK the old values are much better
  //G4double RAYFF = 1.0/1.65;  //old
  //G4double RAYFF = 1.0/1.5;  

  G4double RAYFF = 0.625;

  // Get from the tuning parameters
  //RAYFF = WCSimTuningParams->GetRayff();
  RAYFF = 0.75;

  //T. Akiri: Values from Skdetsim 
  G4double RAYLEIGH_water[NUMENTRIES_water] = {  //G.Barrand : 5*12 = 60. Ok.
      386929*cm*RAYFF,  366249*cm*RAYFF,  346398*cm*RAYFF,  327355*cm*RAYFF,  309097*cm*RAYFF,  
      291603*cm*RAYFF,  274853*cm*RAYFF,  258825*cm*RAYFF,  243500*cm*RAYFF,  228856*cm*RAYFF,  
      214873*cm*RAYFF,  201533*cm*RAYFF,  188816*cm*RAYFF,  176702*cm*RAYFF,  165173*cm*RAYFF,
      154210*cm*RAYFF,  143795*cm*RAYFF,  133910*cm*RAYFF,  124537*cm*RAYFF,  115659*cm*RAYFF,  
      107258*cm*RAYFF, 99318.2*cm*RAYFF, 91822.2*cm*RAYFF,   84754*cm*RAYFF, 78097.3*cm*RAYFF, 
     71836.5*cm*RAYFF,   65956*cm*RAYFF, 60440.6*cm*RAYFF, 55275.4*cm*RAYFF, 50445.6*cm*RAYFF,
       45937*cm*RAYFF, 41735.2*cm*RAYFF, 37826.6*cm*RAYFF, 34197.6*cm*RAYFF, 30834.9*cm*RAYFF, 
     27725.4*cm*RAYFF, 24856.6*cm*RAYFF, 22215.9*cm*RAYFF, 19791.3*cm*RAYFF, 17570.9*cm*RAYFF,   
       15543*cm*RAYFF, 13696.6*cm*RAYFF, 12020.5*cm*RAYFF, 10504.1*cm*RAYFF, 9137.15*cm*RAYFF,
     7909.45*cm*RAYFF,  6811.3*cm*RAYFF, 5833.25*cm*RAYFF,  4966.2*cm*RAYFF, 4201.36*cm*RAYFF, 
     3530.28*cm*RAYFF, 2944.84*cm*RAYFF, 2437.28*cm*RAYFF, 2000.18*cm*RAYFF,  1626.5*cm*RAYFF, 
     1309.55*cm*RAYFF, 1043.03*cm*RAYFF, 821.016*cm*RAYFF,  637.97*cm*RAYFF, 488.754*cm*RAYFF
  };

  G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX", ENERGY_water, RINDEX1, NUMENTRIES_water);
  myMPT1->AddProperty("ABSLENGTH",ENERGY_water, ABSORPTION_water, NUMENTRIES_water);
  myMPT1->AddProperty("RAYLEIGH",ENERGY_water,RAYLEIGH_water,NUMENTRIES_water);

  Water->SetMaterialPropertiesTable(myMPT1);

  //G.Barrand : WARNING : below numbers, needed for media.geo, are of my own. (Not found in MEMPHYS).
  G4double EFFICIENCY[NUMENTRIES_water] = {  //G.Barrand : 10*6 = 60.Ok.
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m,
    0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m, 0.001*m
  };
  write_media_geo(*Water,NUMENTRIES_water,ENERGY_water,ABSORPTION_water,EFFICIENCY,RINDEX1);
  
 /*
  //Gd doped water has the same optical properties as pure water
  DopedWater->SetMaterialPropertiesTable(myMPT1);
   
  G4double BLACKABS_blacksheet[NUMENTRIES_water] =
     { 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 
       1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 
       1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm};
   
  G4MaterialPropertiesTable *myMPT4 = new G4MaterialPropertiesTable();
  myMPT4->AddProperty("ABSLENGTH", ENERGY_water, BLACKABS_blacksheet, NUMENTRIES_water);
  Blacksheet->SetMaterialPropertiesTable(myMPT4);
 */
  
}

// See fairroot/geobase/FairGeoMedium.h and fairroot/geobase/FairGeoMedium.cxx/read() :
// <int:ncomponent> <read:ca[ncomponent]> <read:cz[ncomponent] <real:density> <real:cw[ncomponent]>
// ca : mass number
// cz : atomic number
// cw : weight of the component in a mixture.
// density in g/cm3.

// weightFac : sign of (ncomponent) /** Factor for weights (1: relative w., -1: w. by number of atoms)*/

// <int:sensFlag> <int:fldFlag> <real:fld> <real:epsil>
// sensFlag : sensitivity flag.
// fldFlag : field flag.
// fld : Maximum field value in kilogauss
// epsil : Boundary crossing precision (unit ?).

// <npckov>
// for i in npckov
//   <real:ppckov[i]> <real:absco[i]> <real:effic[i]> <real:rindex[i]>

// ppckov : //[npckov] /** Photon momentum*/
// absco  : //[npckov] /** Absoption length*/
// effic  : //[npckov] /** Detection efficiency*/
// rindex : //[npckov] /** Refraction index*/

int main() {
  MEMPHYS_DetectorConstruction_ConstructMaterials_Water();
  WCSimDetectorConstruction_ConstructMaterials_Water();
  return EXIT_SUCCESS;
}
