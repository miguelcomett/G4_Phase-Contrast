//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
/////////////////////////////////////////////////////
// Code by luca brombal, INFN -Trieste - 12.06.2009
/////////////////////////////////////////////////////

//
// $Id: PepiDetectorConstruction.cc 71079 2013-06-10 20:37:11Z ihrivnac $
//
/// \file PepiDetectorConstruction.cc
/// \brief Implementation of the PepiDetectorConstruction class

#include "PepiDetectorConstruction.hh"
#include "PepiDetectorMessenger.hh"
//#include "PepiImageQualityPhantomParam.hh"

#include "G4RunManager.hh"

#include "G4NistManager.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4PSCylinderSurfaceCurrent.hh"
#include "PepiPSPixiRad.hh"
//#include "PepiPSIoC.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDParticleWithEnergyFilter.hh"

#include "G4VisAttributes.hh"
#include "G4UIcommand.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <tuple>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4ThreadLocal G4bool PepiDetectorConstruction::fConstructedSDandField = false;

PepiDetectorConstruction::PepiDetectorConstruction()
: G4VUserDetectorConstruction(),
  fConstructed(false),
  fConstructedSDandField(false),
  fBidimensional(false),
  fWorldMaterial(0),
  fDetectorMaterial(0),
  fMaskMaterial(0),
  fObjectMaterial(0),
  fObject2Material(0),
  fObject3Material(0),
  fObject4Material(0),
  fSubMaterial(0),
  fSphereMaterial(0),
  fWorldSolid(0),
  fPixiRadSolid(0),
  fPixelSolid(0),
  fM2subSolid(0),
  fM1subSolid(0),
  fM2Solid(0),
  fM1Solid(0),
  fEnvelopeM2Solid(0),
  fEnvelopeM1Solid(0),
  fObjectSolid(0),
  fObject2Solid(0),
  fObject3Solid(0),
  fObject4Solid(0),
  fSphereSolid(0),
  fWorldLogical(0),
  fPixiRadLogical(0),
  fPixelLogical(0),
  fM2Logical(0),
  fM1Logical(0),
  fM1subLogical(0),
  fM2subLogical(0),  
  fEnvelopeM2Logical(0),
  fEnvelopeM1Logical(0),
  fObjectLogical(0),
  fObject2Logical(0),
  fObject3Logical(0),
  fObject4Logical(0),
  fSphereLogical(0),
  fWorldPhysical(0),
  fPixiRadPhysical(0),
  fPixelPhysical(0),
  fM2Physical(0),
  fM1Physical(0),
  fM1subPhysical(0),
  fM2subPhysical(0),  
  fEnvelopeM2Physical(0),
  fEnvelopeM1Physical(0),
  fObjectPhysical(0),
  fObject2Physical(0),
  fObject3Physical(0),
  fObject4Physical(0),
  fSpherePhysical(0),
  fRotAngle(0*deg),
  fTrans(0*um),
  fDith(0*um),
  fObjectDetDistance(10*cm),
  fSrcObjDistance(0.7*m),
  fWorldSizeX(0),
  fWorldSizeY(0),
  fWorldSizeZ(0),
  fOffset(50*cm),
  fnPixelsX(0),
  fnPixelsY(0),
  fPixiRadSizeX(0),
  fPixiRadSizeY(0),
  fPixiRadSizeZ(0),
  fMaskThickness(100*um),
  fM2Aperture(55*um),
  fM2Pitch(110*um),
  fSubThickness(500*um),
  fSourcePosZ(-85*cm),
  fThreshold1(3*keV),
  fThreshold2(3*keV),
  fThreshold3(3*keV),
  fThreshold4(3*keV),
  fThreshold5(3*keV),
  fThreshold6(3*keV),
  fThreshold7(3*keV),
  fThreshold8(3*keV),
  fAcquisitionType("doublemask"),
  fDetType("0COL"),
  fCheckOverlaps(false),
  fMessenger(0)
{
  // - All geometrical parameters depend on the object size
  // The object is defined as a full cylinder 
  // Inside the object there will be details of different 
  // materials 

  fWorldSizeX = 1*m; // Define tamaño de mundo en X
  fWorldSizeY = 1*m; // Define tamaño de mundo en Y
  fWorldSizeZ = 2.3*m; // Define tamaño de mundo en Z

  fObjSizeR = 5.1*cm;  
  fObjSizeY = 6.0*cm; 
  fTras = 100000*um;   // Translation of the sample-mask
  fDithe = 0*um;  // Steps of dithering-Translation of the sample

  fSkinThickness = 1.45*mm;

  fDetailSizeR = 2.5*cm;

  fPixelSizeX =  55*um; // Define tamaño de pixel en X, no definido al inicio
  fPixelSizeY =  55*um; // Define tamaño de pixel en Y, no definido al inicio
  fPixelSizeZ =  1000*um; // Define espesor de pixel en Z, no definido al inicio

  fMessenger = new PepiDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiDetectorConstruction::~PepiDetectorConstruction()
{ 
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* PepiDetectorConstruction::Construct()
{ 
  if (!fConstructed)
  {
    fConstructed = true;
    // - Define the Materials
    DefineMaterials();
    // - Define the Volumes
    DefineVolumes();
  } 

  return fWorldPhysical;
}

void PepiDetectorConstruction::ConstructSDandField()
{
  if(!fConstructedSDandField)
  {
    fConstructedSDandField = true;
    // - Define the Sensitive Detectors
    DefineDetectors();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nist = G4NistManager::Instance();

  // ========================================
  //              MATERIALS
  // ========================================

  G4Material* CdTe          = nist->FindOrBuildMaterial("G4_CADMIUM_TELLURIDE");
  G4Material* Air           = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* PlexiGlass    = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4Material* Water         = nist->FindOrBuildMaterial("G4_WATER");
  G4Material* Nylon  	    = nist->FindOrBuildMaterial("G4_NYLON-6-6");
  G4Material* PolyCarbonate = nist->FindOrBuildMaterial("G4_POLYCARBONATE");
  G4Material* Silicon       = nist->FindOrBuildMaterial("G4_Si");
  G4Material* Graphite      = nist->FindOrBuildMaterial("G4_GRAPHITE");
  G4Material* Gold          = nist->FindOrBuildMaterial("G4_Au");
  G4Material* Bone 	    = nist->FindOrBuildMaterial("G4_B-100_BONE");
  G4Material* Skin          = nist -> FindOrBuildMaterial("G4_SKIN_ICRP");
  
  // ========================================
  //              ADDITIONAL MATERIALS
  // ========================================
  //Additional elements
  G4double zH,aH,zC,aC,zO,aO,zN,aN,zAr,aAr,zI,aI,zNa,aNa,zP,aP,zS,aS,zCl,aCl,zK,aK,zFe,aFe,zAl,aAl,aCa,zCa,aSi,zSi,fractionmassA, densityA,fractionmassPI, densityPI,fractionmassB, densityB, densityAl, densityHPA,densityFGR,densityDBT1, fractionmassFGR,densityADP,fractionmassADP,densityDBT,fractionmassDBT,fractionmassDBT1,densityPQ595,fractionmassPQ595,fractionmassPQ3070,densityPQ3070,fractionmassPQ5050,densityPQ5050,densitySilicone,densityPMMA,fractionmassPMMA,fractionmassMuscle,densityMuscle,fractionmassAorta,densityAorta;
  
  G4String nameAll, symbolH,nameH, symbolC,nameC, symbolO,nameO, symbolN,nameN, symbolAr,nameAr, symbolI,nameI, symbolNa,nameNa, symbolP,nameP, symbolS,nameS, symbolCl,nameCl, symbolK,nameK, symbolFe,nameFe, symbolAl,nameSi,symbolSi,nameA, namePI, nameB, nameAl, nameHPA,nameCa,symbolCa,nameFGR,nameADP,nameDBT,nameDBT1,namePQ595,namePQ3070,namePQ5050,nameSilicone,namePMMA,nameMuscle,nameAorta;
  
  G4int ncomponentsA, ncomponentsPI, ncomponentsB, ncomponentsAl, natomsAl, ncomponentsHPA, natomsHPA,natomsSilicone, ncomponentsFGR,ncomponentsADP,ncomponentsDBT,ncomponentsDBT1,ncomponentsPQ595,ncomponentsPQ3070,ncomponentsPQ5050,ncomponentsSilicone,ncomponentsPMMA,ncomponentsMuscle,ncomponentsAorta;

  aN = 14.01*g/mole;
  G4Element* Nitrogen1 = new G4Element(nameN="Nitrogen",symbolN="N" , zN= 7., aN);
  aO = 16.00*g/mole;
  G4Element* Oxygen1 = new G4Element(nameO="Oxygen" ,symbolO="O" , zO= 8., aO);
  aH = 1.00*g/mole;
  G4Element* Hydrogen1 = new G4Element(nameH="Hydrogen",symbolH="H" , zH= 1., aH);
  aC = 12.00*g/mole;
  G4Element* Carbon1 = new G4Element(nameC="Carbon" ,symbolC="C" , zC= 6., aC);
  aAr = 40.00*g/mole;
  G4Element* Argon1 = new G4Element(nameAr="Argon",symbolAr="Ar" , zAr= 18., aAr);
  aI = 127.00*g/mole;
  G4Element* Iodine1 = new G4Element(nameI="Iodine" ,symbolI="I" , zI= 53., aI);
  aNa = 23.00*g/mole;
  G4Element* Sodium1 = new G4Element(nameNa="Sodium",symbolNa="Na" , zNa= 11., aNa);
  aP = 31.00*g/mole;
  G4Element* Phosphorus1 = new G4Element(nameP="Phosphorus" ,symbolP="P" , zP= 15., aP);
  aS = 32.00*g/mole;
  G4Element* Sulfur1 = new G4Element(nameS="Sulfur",symbolS="S" , zS= 16., aS);
  aCl = 35.40*g/mole;
  G4Element* Chlorine1 = new G4Element(nameCl="Chlorine" ,symbolCl="Cl" , zCl= 17., aCl);
  aK = 39.10*g/mole;
  G4Element* Potassium1 = new G4Element(nameK="Potassium",symbolK="K" , zK= 19., aK);
  aFe = 55.90*g/mole;
  G4Element* Iron1 = new G4Element(nameFe="Iron" ,symbolFe="Fe" , zFe= 26., aFe);
  aAl = 27.00*g/mole;
  G4Element* Aluminium1 = new G4Element(nameAll="Aluminium" ,symbolAl="Al" , zAl= 13., aAl);
  aCa = 40.10*g/mole;
  G4Element* Calcium1 = new G4Element(nameCa="Calcium" ,symbolCa="Ca" , zCa= 20., aCa);
  aSi = 28.10*g/mole;
  G4Element* Silicon1 = new G4Element(nameSi="Silicon" ,symbolSi="Si" , zSi= 14., aSi);
//  aMo = 95.95*g/mole;
  G4Material* Molybdenum1 = nist->FindOrBuildMaterial("G4_Mo");
//  aRh = 102.91*g/mole;
  G4Material* Rhodium1 = nist->FindOrBuildMaterial("G4_Rh");
//  aAg = 107.86*g/mole;
  G4Material* Silver1 = nist->FindOrBuildMaterial("G4_Ag");
  
  
  
  //Air
  densityA = 0.0013*g/cm3;
  G4Material* AAir = new G4Material(nameA="AAir",densityA,ncomponentsA=4);
  AAir ->AddElement(Carbon1, fractionmassA=0.0097*perCent);
  AAir ->AddElement(Oxygen1, fractionmassA=23.1801*perCent);
  AAir ->AddElement(Nitrogen1, fractionmassA=75.5301*perCent);
  AAir ->AddElement(Argon1, fractionmassA=1.2801*perCent);
  
  //P.Iodine
  densityPI = 0.4*g/cm3;
  G4Material* PIodine = new G4Material(namePI="PIodine",densityPI,ncomponentsPI=5);
  PIodine ->AddElement(Carbon1, fractionmassPI=1.981*perCent);
  PIodine ->AddElement(Oxygen1, fractionmassPI=80.371*perCent);
  PIodine ->AddElement(Nitrogen1, fractionmassPI=0.381*perCent);
  PIodine ->AddElement(Hydrogen1, fractionmassPI=10.321*perCent);
  PIodine ->AddElement(Iodine1, fractionmassPI=6.946*perCent);
  
  //Blood
  densityB = 1.16*g/cm3;
  G4Material* Blood = new G4Material(nameB="Blood",densityB,ncomponentsB=10);
  Blood ->AddElement(Carbon1, fractionmassB=11*perCent);
  Blood ->AddElement(Oxygen1, fractionmassB=74.5*perCent);
  Blood ->AddElement(Nitrogen1, fractionmassB=3.3*perCent);
  Blood ->AddElement(Hydrogen1, fractionmassB=10.2*perCent);
  Blood ->AddElement(Sodium1, fractionmassB=0.1*perCent);
  Blood ->AddElement(Phosphorus1, fractionmassB=0.1*perCent);
  Blood ->AddElement(Sulfur1, fractionmassB=0.2*perCent);
  Blood ->AddElement(Chlorine1, fractionmassB=0.3*perCent);
  Blood ->AddElement(Potassium1, fractionmassB=0.2*perCent);
  Blood ->AddElement(Iron1, fractionmassB=0.1*perCent);
  
  //Alumina
  densityAl = 3.95*g/cm3;
  G4Material* Alumina = new G4Material(nameAl="Alumina",densityAl,ncomponentsAl=2);
  Alumina->AddElement(Aluminium1, natomsAl=2);
  Alumina->AddElement(Oxygen1, natomsAl=3);
  
  //HP
  densityHPA = 3.16*g/cm3;
  G4Material* Hydroxyapatite = new G4Material(nameHPA="Hydroxyapatite",densityHPA,ncomponentsHPA=4);
  Hydroxyapatite->AddElement(Calcium1, natomsHPA=10);
  Hydroxyapatite->AddElement(Oxygen1, natomsHPA=26);
  Hydroxyapatite->AddElement(Phosphorus1, natomsHPA=6);
  Hydroxyapatite->AddElement(Hydrogen1, natomsHPA=2);
  
  //Fibroglandular
  densityFGR = 1.020*g/cm3;
  G4Material*  FibroGland= new G4Material(nameFGR="Fibroglandular",densityFGR,ncomponentsFGR=8);
  FibroGland->AddElement(Carbon1, fractionmassFGR=33.2*perCent);
  FibroGland->AddElement(Oxygen1, fractionmassFGR=52.7*perCent);
  FibroGland->AddElement(Phosphorus1, fractionmassFGR=0.1*perCent);
  FibroGland->AddElement(Hydrogen1, fractionmassFGR=10.6*perCent);
  FibroGland->AddElement(Nitrogen1, fractionmassFGR=3.0*perCent);
  FibroGland->AddElement(Sodium1, fractionmassFGR=0.1*perCent);
  FibroGland->AddElement(Sulfur1, fractionmassFGR=0.2*perCent);
  FibroGland->AddElement(Chlorine1, fractionmassFGR=0.1*perCent);
  
  //Aorta
  densityAorta = 1.05*g/cm3;
  G4Material*  Aorta= new G4Material(nameAorta="Aorta",densityAorta,ncomponentsAorta=9);
  Aorta->AddElement(Carbon1, fractionmassAorta=14.7*perCent);
  Aorta->AddElement(Oxygen1, fractionmassAorta=69.8*perCent);
  Aorta->AddElement(Phosphorus1, fractionmassAorta=0.4*perCent);
  Aorta->AddElement(Hydrogen1, fractionmassAorta=9.9*perCent);
  Aorta->AddElement(Nitrogen1, fractionmassAorta=4.2*perCent);
  Aorta->AddElement(Sodium1, fractionmassAorta=0.2*perCent);
  Aorta->AddElement(Sulfur1, fractionmassAorta=0.3*perCent);
  Aorta->AddElement(Potassium1, fractionmassAorta=0.1*perCent);
  Aorta->AddElement(Calcium1, fractionmassAorta=0.4*perCent);
  
  //Adipose
  densityADP = 0.95*g/cm3;
  G4Material*  Adipose= new G4Material(nameADP="Adipose",densityADP,ncomponentsADP=7);
  Adipose->AddElement(Carbon1, fractionmassADP=59.7798*perCent);
  Adipose->AddElement(Oxygen1, fractionmassADP=27.7916*perCent);
  Adipose->AddElement(Sodium1, fractionmassADP=0.1297*perCent);
  Adipose->AddElement(Hydrogen1, fractionmassADP=11.3975*perCent);
  Adipose->AddElement(Nitrogen1, fractionmassADP=0.711*perCent);
  Adipose->AddElement(Sulfur1, fractionmassADP=0.0904*perCent);
  Adipose->AddElement(Chlorine1, fractionmassADP=0.1*perCent);
  
   //PMMA
  densityPMMA = 1.20*g/cm3;
  G4Material*  Methacrylate= new G4Material(namePMMA="PMMA",densityPMMA,ncomponentsPMMA=3);
  Methacrylate->AddElement(Carbon1, fractionmassPMMA=59.9846*perCent);
  Methacrylate->AddElement(Oxygen1, fractionmassPMMA=31.9613*perCent);
  Methacrylate->AddElement(Hydrogen1, fractionmassPMMA=8.0541*perCent);
  
  //Muscle
  densityMuscle = 1.05*g/cm3;
  G4Material*  Muscle= new G4Material(nameMuscle="Muscle",densityMuscle,ncomponentsMuscle=9);
  Muscle->AddElement(Hydrogen1, fractionmassMuscle=10.2*perCent);
  Muscle->AddElement(Carbon1, fractionmassMuscle=14.3*perCent);
  Muscle->AddElement(Nitrogen1, fractionmassMuscle=3.4*perCent);
  Muscle->AddElement(Oxygen1, fractionmassMuscle=71.0*perCent);
  Muscle->AddElement(Sodium1, fractionmassMuscle=0.1*perCent);
  Muscle->AddElement(Phosphorus1, fractionmassMuscle=0.2*perCent);
  Muscle->AddElement(Sulfur1, fractionmassMuscle=0.3*perCent);
  Muscle->AddElement(Chlorine1, fractionmassMuscle=0.1*perCent);
  Muscle->AddElement(Potassium1, fractionmassMuscle=0.4*perCent);
  
  //Polydimethylsiloxane 
  densitySilicone = 0.965*g/cm3;
  G4Material* Silicone = new G4Material(nameSilicone="Silicone",densitySilicone,ncomponentsSilicone=4);
  Silicone->AddElement(Carbon1, natomsSilicone=2);
  Silicone->AddElement(Oxygen1, natomsSilicone=1);
  Silicone->AddElement(Silicon1, natomsSilicone=1);
  Silicone->AddElement(Hydrogen1, natomsSilicone=6);
  
  //DenseBreast
  densityDBT = 1.005*g/cm3;
  G4Material* DenseBreast = new G4Material(nameDBT="DenseBreast",densityDBT,ncomponentsDBT=2);
  DenseBreast ->AddMaterial(Adipose, fractionmassDBT=15*perCent);
  DenseBreast ->AddMaterial(FibroGland, fractionmassDBT=85*perCent);
  
  //StandardBreast
  densityDBT1 = 0.97*g/cm3;
  G4Material* StandardBreast = new G4Material(nameDBT1="StandardBreast",densityDBT1,ncomponentsDBT1=2);
  StandardBreast ->AddMaterial(Adipose, fractionmassDBT1=50*perCent);
  StandardBreast ->AddMaterial(FibroGland, fractionmassDBT1=50*perCent);
  
  //Plaque10-90
  densityPQ595 = 1.396*g/cm3;
  G4Material* PlaqueS = new G4Material(namePQ595="PlaqueS",densityPQ595,ncomponentsPQ595=2);
  PlaqueS ->AddMaterial(Hydroxyapatite, fractionmassPQ595=10*perCent);
  PlaqueS ->AddMaterial(Methacrylate, fractionmassPQ595=90*perCent);
  
  //Plaque30-70
  densityPQ3070 = 1.788*g/cm3;
  G4Material* PlaqueI = new G4Material(namePQ3070="PlaqueI",densityPQ3070,ncomponentsPQ3070=2);
  PlaqueI ->AddMaterial(Hydroxyapatite, fractionmassPQ3070=30*perCent);
  PlaqueI ->AddMaterial(Methacrylate, fractionmassPQ3070=70*perCent);
  
  //Plaque20-80
  densityPQ5050 = 1.592*g/cm3;
  G4Material* PlaqueC = new G4Material(namePQ5050="PlaqueI",densityPQ5050,ncomponentsPQ5050=2);
  PlaqueC ->AddMaterial(Hydroxyapatite, fractionmassPQ5050=20*perCent);
  PlaqueC ->AddMaterial(Methacrylate, fractionmassPQ5050=80*perCent);

  // ========================================
  //         REFRACTION COEFFICIENTS
  // ========================================

  
  // Load delta coefficients and the respective energy interval
  // Values taken from: http://ts-imaging.science.unimelb.edu.au/Services/Simple/ICUtilXdata.aspx
  std::vector<double> GraphiteDelta = LoadDelta("../data/Graphite_delta.txt");
  std::vector<double> SiliconDelta = LoadDelta("../data/Silicon_delta.txt");
  std::vector<double> PlexiGlassDelta = LoadDelta("../data/PMMA_delta.txt");
  std::vector<double> NylonDelta = LoadDelta("../data/Nylon_delta.txt");
  std::vector<double> PolyCarbonateDelta = LoadDelta("../data/Polycharbonate_delta.txt");
  std::vector<double> WaterDelta = LoadDelta("../data/Water_delta.txt");
  std::vector<double> energies = LoadDelta("../data/Energy.txt");
  std::vector<double> AAirDelta = LoadDelta("../data/Air_delta.txt");
  std::vector<double> PIodineDelta = LoadDelta("../data/P.Iodine_delta.txt");
  std::vector<double> BloodDelta = LoadDelta("../data/Blood_delta.txt");
  std::vector<double> AluminaDelta = LoadDelta("../data/Alumina_delta.txt");
  std::vector<double> HydroxyapatiteDelta = LoadDelta("../data/HA_delta.txt");
  std::vector<double> AortaDelta = LoadDelta("../data/Aorta_delta.txt");
  std::vector<double> Breast5050Delta = LoadDelta("../data/Breast5050_delta.txt");
  std::vector<double> Breast8515Delta = LoadDelta("../data/Breast8515_delta.txt");
  std::vector<double> Plaque9010Delta = LoadDelta("../data/Plaque9010_delta.txt");
  std::vector<double> Plaque8020Delta = LoadDelta("../data/Plaque8020_delta.txt");
  std::vector<double> Plaque7030Delta = LoadDelta("../data/Plaque7030_delta.txt");
  std::vector<double> GlandularDelta = LoadDelta("../data/Glandular_delta.txt");
  std::vector<double> AdiposeDelta = LoadDelta("../data/Adipose_delta.txt");
  std::vector<double> MuscleDelta = LoadDelta("../data/Muscle_delta.txt");
  std::vector<double> SiliconeDelta = LoadDelta("../data/Silicone_delta.txt");
  std::vector<double> BoneDelta = LoadDelta("../data/Bone_delta.txt");
  std::vector<double> SkinDelta = LoadDelta("../data/Skin_delta.txt");
  	
  
  
  G4int NumEntries = static_cast<int>(energies.size()); 
  std::vector<double> GraphiteRindex(NumEntries); 
  std::vector<double> SiliconRindex(NumEntries);
  std::vector<double> PlexiGlassRindex(NumEntries);
  std::vector<double> NylonRindex(NumEntries);
  std::vector<double> PolyCarbonateRindex(NumEntries);
  std::vector<double> WaterRindex(NumEntries);
  std::vector<double> CdTeRindex(NumEntries);
  std::vector<double> AirRindex(NumEntries);
  std::vector<double> GoldRindex(NumEntries);
  std::vector<double> AAirRindex(NumEntries);
  std::vector<double> PIodineRindex(NumEntries);
  std::vector<double> BloodRindex(NumEntries);
  std::vector<double> AluminaRindex(NumEntries);
  std::vector<double> HARindex(NumEntries);
  std::vector<double> SiliconeRindex(NumEntries);
  std::vector<double> MolybdenumRindex(NumEntries);
  std::vector<double> RhodiumRindex(NumEntries);
  std::vector<double> SilverRindex(NumEntries);
  std::vector<double> AortaRindex(NumEntries);
  std::vector<double> Breast5050Rindex(NumEntries);
  std::vector<double> Breast8515Rindex(NumEntries);
  std::vector<double> Plaque9010Rindex(NumEntries);
  std::vector<double> Plaque8020Rindex(NumEntries);
  std::vector<double> Plaque7030Rindex(NumEntries);
  std::vector<double> AdiposeRindex(NumEntries);
  std::vector<double> GlandularRindex(NumEntries);
  std::vector<double> MuscleRindex(NumEntries);
  std::vector<double> BoneRindex(NumEntries);
  std::vector<double> SkinRindex(NumEntries);
  
  for (G4int i = 0; i < NumEntries; ++i)
  {
    PlexiGlassRindex[i]     = 1 - PlexiGlassDelta[i];
    GraphiteRindex[i]       = 1 - GraphiteDelta[i];
    SiliconRindex[i]        = 1 - SiliconDelta[i];

    WaterRindex[i]          = 1 - WaterDelta[i];
    NylonRindex[i]	    = 1 - NylonDelta[i];
    PolyCarbonateRindex[i]  = 1 - PolyCarbonateDelta[i];
    energies[i] = energies[i]*keV;
    CdTeRindex[i]           = 1;
    AirRindex[i]            = 1;
    GoldRindex[i]           = 1;
    AAirRindex[i]          = 1 - AAirDelta[i];
    PIodineRindex[i]	    = 1 - PIodineDelta[i];
    BloodRindex[i]  = 1 - BloodDelta[i];
    AluminaRindex[i]          = 1 - AluminaDelta[i];
    HARindex[i]	    = 1 - HydroxyapatiteDelta[i];
    SiliconeRindex[i]	    = 1 - SiliconeDelta[i];
    MolybdenumRindex[i]	    = 1;
    RhodiumRindex[i]	    = 1;
    SilverRindex[i]	    = 1;
    AortaRindex[i]	    = 1 - AortaDelta[i];
    Breast5050Rindex[i]	    = 1 - Breast5050Delta[i];
    Breast8515Rindex[i]	    = 1 - Breast8515Delta[i];
    Plaque9010Rindex[i]	    = 1 - Plaque9010Delta[i];
    Plaque8020Rindex[i]	    = 1 - Plaque8020Delta[i];
    Plaque7030Rindex[i]	    = 1 - Plaque7030Delta[i];
    AdiposeRindex[i]	    = 1 - AdiposeDelta[i];
    GlandularRindex[i]	    = 1 - GlandularDelta[i];
    MuscleRindex[i]	    = 1 - MuscleDelta[i];
    BoneRindex[i]	    = 1 - BoneDelta[i];
    SkinRindex[i]	    = 1 - SkinDelta[i];
 }
 

  G4MaterialPropertiesTable* AAirMatPropTbl = new G4MaterialPropertiesTable();
  AAirMatPropTbl->AddProperty("RINDEX",energies.data(),AAirRindex.data(),NumEntries);
  AAir->SetMaterialPropertiesTable(AAirMatPropTbl);
  
  G4MaterialPropertiesTable* PIodineMatPropTbl = new G4MaterialPropertiesTable();
  PIodineMatPropTbl->AddProperty("RINDEX",energies.data(),PIodineRindex.data(),NumEntries);
  PIodine->SetMaterialPropertiesTable(PIodineMatPropTbl);
  
  G4MaterialPropertiesTable* BloodMatPropTbl = new G4MaterialPropertiesTable();
  BloodMatPropTbl->AddProperty("RINDEX",energies.data(),BloodRindex.data(),NumEntries);
  Blood->SetMaterialPropertiesTable(BloodMatPropTbl);
  
  G4MaterialPropertiesTable* AluminaMatPropTbl = new G4MaterialPropertiesTable();
  AluminaMatPropTbl->AddProperty("RINDEX",energies.data(),AluminaRindex.data(),NumEntries);
  Alumina->SetMaterialPropertiesTable(AluminaMatPropTbl);
  
  G4MaterialPropertiesTable* HAMatPropTbl = new G4MaterialPropertiesTable();
  HAMatPropTbl->AddProperty("RINDEX",energies.data(),HARindex.data(),NumEntries);
  Hydroxyapatite->SetMaterialPropertiesTable(HAMatPropTbl);
  
  G4MaterialPropertiesTable* DBTMatPropTbl = new G4MaterialPropertiesTable();
  DBTMatPropTbl->AddProperty("RINDEX",energies.data(),Breast8515Rindex.data(),NumEntries);
  DenseBreast->SetMaterialPropertiesTable(DBTMatPropTbl); 
  
  G4MaterialPropertiesTable* DBT1MatPropTbl = new G4MaterialPropertiesTable();
  DBT1MatPropTbl->AddProperty("RINDEX",energies.data(),Breast5050Rindex.data(),NumEntries);
  StandardBreast->SetMaterialPropertiesTable(DBT1MatPropTbl);
  
  G4MaterialPropertiesTable* AortaMatPropTbl = new G4MaterialPropertiesTable();
  AortaMatPropTbl->AddProperty("RINDEX",energies.data(),AortaRindex.data(),NumEntries);
  Aorta->SetMaterialPropertiesTable(AortaMatPropTbl);
  
  
  G4MaterialPropertiesTable* SiliconeMatPropTbl = new G4MaterialPropertiesTable();
  SiliconeMatPropTbl->AddProperty("RINDEX",energies.data(),SiliconeRindex.data(),NumEntries);
  Silicone->SetMaterialPropertiesTable(SiliconeMatPropTbl);
  
  G4MaterialPropertiesTable* PlSMatPropTbl = new G4MaterialPropertiesTable();
  PlSMatPropTbl->AddProperty("RINDEX",energies.data(),Plaque9010Rindex.data(),NumEntries);
  PlaqueS->SetMaterialPropertiesTable(PlSMatPropTbl);
  
  G4MaterialPropertiesTable* PlIMatPropTbl = new G4MaterialPropertiesTable();
  PlIMatPropTbl->AddProperty("RINDEX",energies.data(),Plaque7030Rindex.data(),NumEntries);
  PlaqueI->SetMaterialPropertiesTable(PlIMatPropTbl);
  
  G4MaterialPropertiesTable* PlCMatPropTbl = new G4MaterialPropertiesTable();
  PlCMatPropTbl->AddProperty("RINDEX",energies.data(),Plaque8020Rindex.data(),NumEntries);
  PlaqueC->SetMaterialPropertiesTable(PlCMatPropTbl); 

  G4MaterialPropertiesTable* GraphiteMatPropTbl = new G4MaterialPropertiesTable();
  GraphiteMatPropTbl->AddProperty("RINDEX",energies.data(),GraphiteRindex.data(),NumEntries);
  Graphite->SetMaterialPropertiesTable(GraphiteMatPropTbl);  

  G4MaterialPropertiesTable* SiliconMatPropTbl = new G4MaterialPropertiesTable();
  SiliconMatPropTbl->AddProperty("RINDEX",energies.data(),SiliconRindex.data(),NumEntries);
  Silicon->SetMaterialPropertiesTable(SiliconMatPropTbl);  

  G4MaterialPropertiesTable* WaterMatPropTbl = new G4MaterialPropertiesTable();
  WaterMatPropTbl->AddProperty("RINDEX",energies.data(),WaterRindex.data(),NumEntries);
  Water->SetMaterialPropertiesTable(WaterMatPropTbl);    

  G4MaterialPropertiesTable* PlexiGlassMatPropTbl = new G4MaterialPropertiesTable();
  PlexiGlassMatPropTbl->AddProperty("RINDEX",energies.data(),PlexiGlassRindex.data(),NumEntries);
  PlexiGlass->SetMaterialPropertiesTable(PlexiGlassMatPropTbl);

  G4MaterialPropertiesTable* NylonMatPropTbl = new G4MaterialPropertiesTable();
  NylonMatPropTbl->AddProperty("RINDEX",energies.data(),NylonRindex.data(),NumEntries);
  Nylon->SetMaterialPropertiesTable(NylonMatPropTbl);

  G4MaterialPropertiesTable* PolyCarbonateMatPropTbl = new G4MaterialPropertiesTable();
  PolyCarbonateMatPropTbl->AddProperty("RINDEX",energies.data(),PolyCarbonateRindex.data(),NumEntries);
  PolyCarbonate->SetMaterialPropertiesTable(PolyCarbonateMatPropTbl);

  G4MaterialPropertiesTable* CdTeMatPropTbl = new G4MaterialPropertiesTable();
  CdTeMatPropTbl->AddProperty("RINDEX",energies.data(),CdTeRindex.data(),NumEntries);
  CdTe->SetMaterialPropertiesTable(CdTeMatPropTbl);

  G4MaterialPropertiesTable* AirMatPropTbl = new G4MaterialPropertiesTable();
  AirMatPropTbl->AddProperty("RINDEX",energies.data(),AirRindex.data(),NumEntries);
  Air->SetMaterialPropertiesTable(AirMatPropTbl);

  G4MaterialPropertiesTable* GoldMatPropTbl = new G4MaterialPropertiesTable();
  GoldMatPropTbl->AddProperty("RINDEX",energies.data(),GoldRindex.data(),NumEntries);
  Gold->SetMaterialPropertiesTable(GoldMatPropTbl);
  
  G4MaterialPropertiesTable* MoMatPropTbl = new G4MaterialPropertiesTable();
  MoMatPropTbl->AddProperty("RINDEX",energies.data(),MolybdenumRindex.data(),NumEntries);
  Molybdenum1->SetMaterialPropertiesTable(MoMatPropTbl);
  
  G4MaterialPropertiesTable* RhMatPropTbl = new G4MaterialPropertiesTable();
  RhMatPropTbl->AddProperty("RINDEX",energies.data(),RhodiumRindex.data(),NumEntries);
  Rhodium1->SetMaterialPropertiesTable(RhMatPropTbl);
  
  G4MaterialPropertiesTable* AgMatPropTbl = new G4MaterialPropertiesTable();
  AgMatPropTbl->AddProperty("RINDEX",energies.data(),SilverRindex.data(),NumEntries);
  Silver1->SetMaterialPropertiesTable(AgMatPropTbl);
  
  G4MaterialPropertiesTable* PMMAMatPropTbl = new G4MaterialPropertiesTable();
  PMMAMatPropTbl->AddProperty("RINDEX",energies.data(),PlexiGlassRindex.data(),NumEntries);
  Methacrylate->SetMaterialPropertiesTable(PMMAMatPropTbl);

  G4MaterialPropertiesTable* AdiposeMatPropTbl = new G4MaterialPropertiesTable();
  AdiposeMatPropTbl->AddProperty("RINDEX",energies.data(),AdiposeRindex.data(),NumEntries);
  Adipose->SetMaterialPropertiesTable(AdiposeMatPropTbl);
  
  G4MaterialPropertiesTable* MuscleMatPropTbl = new G4MaterialPropertiesTable();
  MuscleMatPropTbl->AddProperty("RINDEX",energies.data(),MuscleRindex.data(),NumEntries);
  Muscle->SetMaterialPropertiesTable(MuscleMatPropTbl);

  G4MaterialPropertiesTable* BoneMatPropTbl = new G4MaterialPropertiesTable();
  BoneMatPropTbl->AddProperty("RINDEX",energies.data(),BoneRindex.data(),NumEntries);
  Bone->SetMaterialPropertiesTable(BoneMatPropTbl);

  G4MaterialPropertiesTable* SkinMatPropTbl = new G4MaterialPropertiesTable();
  SkinMatPropTbl->AddProperty("RINDEX",energies.data(),SkinRindex.data(),NumEntries);
  Skin->SetMaterialPropertiesTable(SkinMatPropTbl); 
  
  // ========================================
  //           DEFAULT MATERIALS
  // ========================================

  fWorldMaterial    = Air;// Vacuum1;
  fIonCMaterial     = Air;// Vacuum1;
  fDetectorMaterial = CdTe;
  fMaskMaterial     = Gold;
  fObjectMaterial   = Skin;
  fObject2Material  = Adipose;
  fObject3Material  = Muscle;
  fObject4Material  = Bone;
  fSubMaterial	    = Graphite;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::DefineVolumes()
{
  
  // - PIXIRAD Parameters

  fnPixelsX = 256;

  // - The y-dimension of PixiRad changes if we want
  // a bidimensional acquisition or not
  if (fBidimensional)
  {
    fnPixelsY = 256;
  }
  else
  {
    fnPixelsY = 1;
  }

  fPixiRadSizeX = fPixelSizeX * fnPixelsX; // - Tamaño detector X=Tamaño Pixel en X*n.pixeles en X
  fPixiRadSizeY = fPixelSizeY * fnPixelsY; // - Tamaño detector Y=Tamaño Pixel en Y*n.pixeles en Y
  fPixiRadSizeZ = fPixelSizeZ; // - Espesor de sensor

  // ========================================
  //                  WORLD
  // ========================================
  
  // - Build the WORLD as an unrotated Box in (0,0,0)
  fWorldSolid = new G4Box("World",                         //its name
                          fWorldSizeX/2,                   //its size
                          fWorldSizeY/2,
                          fWorldSizeZ/2);          
   
  fWorldLogical = new G4LogicalVolume(fWorldSolid,         //its solid
                                      fWorldMaterial,      //its material
                                      "World");            //its name
                       
  fWorldPhysical =  new G4PVPlacement(0,                   //no rotation
                                      G4ThreeVector(),     //at (0,0,0)
                                      fWorldLogical,       //its logical volume
                                      "World",             //its name
                                      0,                   //its mother  volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps 


  // ========================================
  //                DETECTOR
  // ========================================

  // - Build the DETECTOR PIXEL UNIT as a Box
  fPixelSolid = new G4Box("Pixel",                          //its name
                          fPixelSizeX/2,                    //its size
                          fPixelSizeY/2,
                          fPixelSizeZ/2);
                     
  fPixelLogical = new G4LogicalVolume(fPixelSolid,          //its solid
                                      fDetectorMaterial,    //its material
                                      "PixelLV");           //its name

  // Build the Detector Envelope 
  fPixiRadSolid =  new G4Box("PixiRad",                    //its name                 
                             fPixiRadSizeX/2,              //its size
                             fPixiRadSizeY/2,
                             fPixiRadSizeZ/2);        
      
  fPixiRadLogical =  new G4LogicalVolume(fPixiRadSolid,    //its solid
                                         fWorldMaterial,   //its material
                                         "PixiRad");       //its name
                    
  // - Place the physical copies of the pixel in a x-y matrix
  // The full detector is built from the top-left corner in
  // [*][*][*][*][*][*][*][*][*][*][*][*] #1 row (fnPixelsX long)
  // [*][*][*][*]........................ #2 row (fnPixelsX long)
  // ....................................
  // .................................... # fnPixelsX * fnPixels Y
  G4int copy_no=0;  

  for (G4int iY = 0; iY < fnPixelsY ; iY++)
  {
    for (G4int iX = 0; iX < fnPixelsX ; iX++)
    {
      G4double x = + iX*fPixelSizeX - fPixiRadSizeX/2;// + fPixelSizeX/2;
      G4double y = - iY*fPixelSizeY + fPixiRadSizeY/2;// - fPixelSizeY/2;
      
      G4ThreeVector position = G4ThreeVector(x, y, 0);
      G4String  name = "Pixel_" + G4UIcommand::ConvertToString(copy_no);

      fPixelPhysical =  new G4PVPlacement(0,                           //its rotation
                                          position,                    //its position
                                          fPixelLogical,               //its logical volume
                                          name,                        //its name
                                          fPixiRadLogical,             //its mother volume
                                          false,                       //no boolean operation
                                          copy_no,                     //copy number
                                          fCheckOverlaps);             //checking overlaps 

      copy_no++;                              
    }
  }

  // - Place the Detector Envelope in the World
  G4ThreeVector positionPixirad = G4ThreeVector(0, 0, fSourcePosZ+fSrcObjDistance+fObjectDetDistance);
  fPixiRadPhysical = new G4PVPlacement(0,                                                  //its rotation
                                       positionPixirad,					   //its position
                                       fPixiRadLogical,                                    //its logical volume
                                       "PixiRad",                                          //its name
                                       fWorldLogical,                                      //its mother volume
                                       false,                                              //no boolean operation
                                       0,                                                  //copy number
                                       fCheckOverlaps);                                    //checking overlaps


  G4Region* aRegion = new G4Region("PixiRad");
  fPixiRadLogical->SetRegion(aRegion);
  aRegion->AddRootLogicalVolume(fPixiRadLogical);
  // ========================================
  //                DETECTOR MASK M2 and SUBSTRATE
  // ========================================

  if (fAcquisitionType=="doublemask")
  {
  G4double mag_M2 = (fSrcObjDistance+fObjectDetDistance)/(fSrcObjDistance+fObjectDetDistance-(fMaskThickness/2+fPixiRadSizeZ/2)); // magnification of the mask M2
  G4ThreeVector M2Position = positionPixirad-G4ThreeVector(0*cm,0,(fMaskThickness+fPixiRadSizeZ)/2+1*nm);
  // - Build the MASK APERTURE UNIT as a Box
  CreateMask("M2", mag_M2,fM2Pitch, fM2Aperture, M2Position, fMaskThickness, fMaskMaterial, fEnvelopeM2Logical, fEnvelopeM2Physical);
  CreateSubstrate("M2sub", mag_M2, M2Position-G4ThreeVector(0,0,fSubThickness/2+fMaskThickness/2), fSubThickness, fSubMaterial, fM2subLogical, fM2subPhysical);
  }
  // ========================================
  //                 Objects
  // ========================================
  // - Build the Object 1 as a Tube
  G4double Rmin=4.45*cm;
  G4double Rmid=4.95*cm;
  G4double Rs=2.25*cm;
                                       
  fObjectSolid = new G4Tubs("Tube",                         //its name
                            Rmid,                     //its inner radius
                            fObjSizeR,                     //its outer radius
                            fObjSizeY,                 //its half lenght Z
                            0,               //its initial phi angle
                            2.*pi);                //its angle of the segment

  fObjectLogical = new G4LogicalVolume(fObjectSolid,       //its solid
                                       fObjectMaterial,    //its material
                                       "TubeLV");          //its name                                     
  
  fScoringVolume=fObjectLogical;
  G4ThreeVector objectPosition = G4ThreeVector(0.*mm+fDithe, 0, fSourcePosZ+fSrcObjDistance);
  
  
 fObject2Solid = new G4Tubs("Tube2",                         //its name
                            Rmin,                     //its inner radius
                            Rmid,                     //its outer radius
                            fObjSizeY,                 //its half lenght Z
                            0,               //its initial phi angle
                            2.*pi);                //its angle of the segment

  fObject2Logical = new G4LogicalVolume(fObject2Solid,       //its solid
                                       fObject2Material,    //its material
                                       "Tube2LV");          //its name
  G4ThreeVector objectPosition2 = G4ThreeVector(0.*mm+fDithe, 0, fSourcePosZ+fSrcObjDistance);
  

  
  fObject3Solid = new G4Tubs("Tube3",                         //its name
                            Rs,                     //its inner radius
                            Rmin,                     //its outer radius
                            fObjSizeY,                 //its half lenght Z
                            0,               //its initial phi angle
                            2.*pi);                //its angle of the segment

  fObject3Logical = new G4LogicalVolume(fObject3Solid,       //its solid
                                       fObject3Material,    //its material
                                       "Tube3LV");          //its name

  G4ThreeVector objectPosition3 = G4ThreeVector(0.*mm+fDithe, 0,fSourcePosZ+fSrcObjDistance);
  
  fObject4Solid = new G4Tubs("Tube4",                         //its name
                            0,                     //its inner radius
                            Rs,                     //its outer radius
                            fObjSizeY,                 //its half lenght Z
                            0,               //its initial phi angle
                            2.*pi);                //its angle of the segment

  PorousBone = fObject4Solid;
  G4int numPores = 20;
  G4double poreRadius = 0.30 * cm;

  Pore = new G4Sphere("Pore", 0, poreRadius, 0, 2.*pi, 0, 1.*pi);

  for (int i = 1; i <= numPores; i++)
  {
        G4double r      = G4RandGauss::shoot(0.6, 0.25) * (Rs - 0);
        G4double theta  = G4UniformRand() * 2.*pi;
        G4double z      = (G4UniformRand() * (fObjSizeY - 0 - poreRadius)) + 0;
      
        G4double x = r * std::cos(theta);
        G4double y = r * std::sin(theta);

        G4ThreeVector porePosition = G4ThreeVector(x, y, -z);
        PorousBone = new G4SubtractionSolid("PorousBone", PorousBone, Pore, 0, porePosition);
  }

  fObject4Logical = new G4LogicalVolume(PorousBone,       //its solid
                                       fObject4Material,    //its material
                                       "Tube4LV");          //its name

  G4ThreeVector objectPosition4 = G4ThreeVector(0.*mm+fDithe, 0,fSourcePosZ+fSrcObjDistance);
  
  G4RotationMatrix* rotMat =  new G4RotationMatrix();
  rotMat->rotateX(fRotAngle);
                                      
  fObject4Physical = new G4PVPlacement(rotMat,              //its rotation
                                      objectPosition4,       //its translation
                                      fObject4Logical,      //its logical volume
                                      "Tube4",              //its name
                                      fWorldLogical,       //its mother volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps                                     
                                      

  fObject3Physical = new G4PVPlacement(rotMat,              //its rotation
                                      objectPosition3,       //its translation
                                      fObject3Logical,      //its logical volume
                                      "Tube3",              //its name
                                      fWorldLogical,       //its mother volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps


  fObject2Physical = new G4PVPlacement(rotMat,              //its rotation
                                      objectPosition2,       //its translation
                                      fObject2Logical,      //its logical volume
                                      "Tube2",              //its name
                                      fWorldLogical,       //its mother volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps

  fObjectPhysical = new G4PVPlacement(rotMat,              //its rotation
                                      objectPosition,       //its translation
                                      fObjectLogical,      //its logical volume
                                      "Tube",              //its name
                                      fWorldLogical,       //its mother volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking over


  // ========================================
  //       SAMPLE MASK M1 and substrate
  // ========================================
  if (fAcquisitionType=="doublemask"|| fAcquisitionType=="singlemask")
  {
  G4double mag_M1 = (fSrcObjDistance+fObjectDetDistance)/(fSrcObjDistance-(fMaskThickness/2+fObjSizeR)); // magnification of the mask M1 
  G4double rel_magg = fSrcObjDistance/(fSrcObjDistance-(fMaskThickness+2*fObjSizeR)/2); 
  G4ThreeVector M1Position = G4ThreeVector(fTras/rel_magg,0,fSourcePosZ+fSrcObjDistance)-G4ThreeVector(0,0,(fMaskThickness+2*fObjSizeR)/2);
  
  std::tie(fEnvelopeM1Logical,fEnvelopeM1Physical) = CreateMask("M1", mag_M1,fM2Pitch, fM2Aperture, M1Position, fMaskThickness, fMaskMaterial, fEnvelopeM1Logical, fEnvelopeM1Physical);
  
  std::tie(fM1subLogical,fM1subPhysical) = CreateSubstrate("M1sub", mag_M1, M1Position-G4ThreeVector(0,0,fSubThickness/2+fMaskThickness/2), fSubThickness, fSubMaterial, fM1subLogical,fM1subPhysical);
  }

  // ========================================
  //               ION CHAMBER
  // ========================================

  // - Build the ION CHAMBER as an unrotated Box
  fIonCSolid = new G4Box("IonChamber",                     //its name
                         fPixiRadSizeX/2,                            //its size
                         fPixiRadSizeY/2,
                         2*mm);          
   
  fIonCLogical = new G4LogicalVolume(fIonCSolid,          //its solid
                                     fIonCMaterial,       //its material
                                     "IonChamberLV");     //its name
                       
  G4ThreeVector IOCposition = positionPixirad - G4ThreeVector(0, 0, 10*mm);
  fIonCPhysical =  new G4PVPlacement(0,                   //no rotation
                                     IOCposition,            //at position
                                     fIonCLogical,        //its logical volume
                                     "IonChamber",        //its name
                                     fWorldLogical,       //its mother  volume
                                     false,               //no boolean operation
                                     0,                   //copy number
                                     fCheckOverlaps);     //checking overlaps 


  
  // ========================================
  //              VISUALIZATION
  // ========================================

  G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  fWorldLogical->SetVisAttributes(worldVisAtt);  
  // logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());  

  G4VisAttributes* objectVisAtt = new G4VisAttributes(G4Colour(0.94, 0.5, 0.5));
//  objectVisAtt->SetForceWireframe(true);//Gray
  objectVisAtt->SetForceSolid(true);
  // objectVisAtt->SetForceAuxEdgeVisible(true);     
  fObjectLogical->SetVisAttributes(objectVisAtt);

  G4VisAttributes* objectVisAtt2 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
//  objectVisAtt2->SetForceWireframe(true); //Yellow
  objectVisAtt2->SetForceSolid(true);
 //  objectVisAtt->SetForceAuxEdgeVisible(true);     
  fObject2Logical->SetVisAttributes(objectVisAtt2);

  G4VisAttributes* objectVisAtt3 = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
//  objectVisAtt3->SetForceWireframe(true); //Red
  objectVisAtt3->SetForceSolid(true);
  // objectVisAtt->SetForceAuxEdgeVisible(true);     
  fObject3Logical->SetVisAttributes(objectVisAtt3);

  G4VisAttributes* objectVisAtt4 = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
//  objectVisAtt3->SetForceWireframe(true); //White
  objectVisAtt4->SetForceSolid(true);
  // objectVisAtt->SetForceAuxEdgeVisible(true);     
  fObject4Logical->SetVisAttributes(objectVisAtt4);

  G4VisAttributes* ionCVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  ionCVisAtt->SetVisibility(false);
  ionCVisAtt->SetForceWireframe(true);
  fIonCLogical->SetVisAttributes(ionCVisAtt);  
 
  G4VisAttributes* pixelVisAtt = new G4VisAttributes(G4Colour(0.0,0,1.0));
  pixelVisAtt->SetForceWireframe(true);
  //pixelVisAtt->SetForceSolid(true);
  pixelVisAtt->SetVisibility(false);
  // pixelVisAtt->SetForceAuxEdgeVisible(true);
  fPixelLogical->SetVisAttributes(pixelVisAtt);

  G4VisAttributes* pixiradVisAtt = new G4VisAttributes(G4Colour(0.0,0,1.0));
  pixiradVisAtt->SetForceWireframe(true);
  pixiradVisAtt->SetVisibility(true);
  pixiradVisAtt->SetForceSolid(true);
  // pixelVisAtt->SetForceAuxEdgeVisible(true);
  fPixiRadLogical->SetVisAttributes(pixiradVisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::DefineDetectors()
{
  // ========================================
  //                 SCORING
  // ========================================

  G4SDManager* pSDMan = G4SDManager::GetSDMpointer();
  pSDMan->SetVerboseLevel(1);
  
  G4double eKMin = 1*keV;
  G4double eKMax = 150*keV;

  G4SDParticleWithEnergyFilter* gammaEKinFilter = new G4SDParticleWithEnergyFilter("gammaEKinFilter",eKMin,eKMax);
  gammaEKinFilter->add("gamma");

  // ========================================
  // - Define a Multifunctional detector 
  // ----> Ideal photon counter 
  // ========================================

  G4MultiFunctionalDetector* pixiRadSD = new G4MultiFunctionalDetector("PixiRadSD");
  pSDMan->AddNewDetector(pixiRadSD);
  SetSensitiveDetector("PixelLV",pixiRadSD);

  // - Ideal photon counter  scores the number of gammas that hit its -Z surface from outside
  // Surface is defined at the -Z surface.
  // Direction                  -Z   +Z
  //   0  IN || OUT            ->|<-  |
  //   1  IN                   ->|    |
  //   2  OUT                    |<-  |
  G4PSFlatSurfaceCurrent* sTotSurfCurrent = new G4PSFlatSurfaceCurrent("TotalSurfCurrent", 1);
  sTotSurfCurrent->SetFilter(gammaEKinFilter);
  sTotSurfCurrent->DivideByArea(false);
  sTotSurfCurrent->Weighted(false);


  pixiRadSD->RegisterPrimitive(sTotSurfCurrent);


  // -----------------------------------------------
  // - Photon counter (Pixirad) with realistic energy response and up to 3 thresholds per pixel
  if(fDetType=="1COL"|| fDetType=="2COL"||fDetType=="3COL"||fDetType=="4COL"|| fDetType=="5COL"||fDetType=="6COL"|| fDetType=="7COL"||fDetType=="8COL" )
  {
      PepiPSPixiRad* sPixiRad = new PepiPSPixiRad("Threshold1", fThreshold1, 
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

      pixiRadSD->RegisterPrimitive(sPixiRad);
      
      if(fDetType=="2COL"||fDetType=="3COL"||fDetType=="4COL"|| fDetType=="5COL"||fDetType=="6COL"|| fDetType=="7COL"|| fDetType=="8COL")
      {
       PepiPSPixiRad* sPixiRad2 = new PepiPSPixiRad("Threshold2", fThreshold2, 
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

       pixiRadSD->RegisterPrimitive(sPixiRad2);

       if(fDetType=="3COL"||fDetType=="4COL"|| fDetType=="5COL"||fDetType=="6COL"|| fDetType=="7COL"|| fDetType=="8COL")
       {
        PepiPSPixiRad* sPixiRad3 = new PepiPSPixiRad("Threshold3", fThreshold3,
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

        pixiRadSD->RegisterPrimitive(sPixiRad3);

        if(fDetType=="4COL"|| fDetType=="5COL"||fDetType=="6COL"|| fDetType=="7COL"|| fDetType=="8COL")
        {
         PepiPSPixiRad* sPixiRad4 = new PepiPSPixiRad("Threshold4", fThreshold4,
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

         pixiRadSD->RegisterPrimitive(sPixiRad4);

         if(fDetType=="5COL"||fDetType=="6COL"|| fDetType=="7COL"|| fDetType=="8COL")
         {
          PepiPSPixiRad* sPixiRad5 = new PepiPSPixiRad("Threshold5", fThreshold5,
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

          pixiRadSD->RegisterPrimitive(sPixiRad5);

          if(fDetType=="6COL"|| fDetType=="7COL"|| fDetType=="8COL")
          {
           PepiPSPixiRad* sPixiRad6 = new PepiPSPixiRad("Threshold6", fThreshold6,
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

           pixiRadSD->RegisterPrimitive(sPixiRad6);

           if(fDetType=="7COL"|| fDetType=="8COL")
           {
            PepiPSPixiRad* sPixiRad7 = new PepiPSPixiRad("Threshold7", fThreshold7,
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

            pixiRadSD->RegisterPrimitive(sPixiRad7);

            if(fDetType=="8COL")
            {
             PepiPSPixiRad* sPixiRad8 = new PepiPSPixiRad("Threshold8", fThreshold8,
                                                fnPixelsX, fnPixelsY,
                                                fPixelSizeX, fPixelSizeY, fPixelSizeZ, "keV");

             pixiRadSD->RegisterPrimitive(sPixiRad8);
            } 

           }
          }
         }
        }
       }
      
      }
  }
  
  // ========================================  
  // - Define a Multifunctional detector 
  // ----> ION CHAMBER to check if there is flux in front of the detector 
  // ========================================

  G4MultiFunctionalDetector* ionChamberSD = new G4MultiFunctionalDetector("IonChamberSD");
  pSDMan->AddNewDetector(ionChamberSD);
  SetSensitiveDetector("IonChamberLV",ionChamberSD);

  // - Ion Chamber scores the number of photons that enter its surface
  // Surface is defined at the -Z surface.
  // Direction                  -Z   +Z
  //   0  IN || OUT            ->|<-  |
  //   1  IN                   ->|    |
  //   2  OUT                    |<-  |
  G4PSFlatSurfaceCurrent* sCurrentIoC = new G4PSFlatSurfaceCurrent("CurrentIoC", 1, "permm2");
  sCurrentIoC->SetFilter(gammaEKinFilter);
  
  ionChamberSD->RegisterPrimitive(sCurrentIoC);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::SetObjectDetDistance(G4double objectDetDistance)
{
  if(!fConstructed)
  {
    if(objectDetDistance < 10*cm || objectDetDistance > 2*m)
    {
      G4cerr << objectDetDistance << " is out of bounds (Must be > 0.1 m AND < 2 m. - Command is ignored." << G4endl;
    }
    else
    {
      fObjectDetDistance = objectDetDistance;
    }
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// NEW for PEPI
void PepiDetectorConstruction::SetSourcePosZ(G4double sourcePosZ)
{
  if(!fConstructed)
  {
    if(sourcePosZ < -2.3/2*m)
    {
      G4cerr << sourcePosZ << "Source is probably outside the World volume - Command is ignored." << G4endl;
    }
    else
    {
      fSourcePosZ = sourcePosZ;
    }
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// NEW for PEPI
void PepiDetectorConstruction::SetSrcObjDistance(G4double srcObjDistance)
{
  if(!fConstructed)
  {
    if(srcObjDistance < 0*cm || srcObjDistance > 2.3*m)
    {
      G4cerr << srcObjDistance << "Source object distance is out of bounds (Must be > 0 m AND < World size. - Command is ignored." << G4endl;
    }
    else
    {
      fSrcObjDistance = srcObjDistance;
    }
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetMaskThickness(G4double maskThickness)
{
  if(!fConstructed)
  {
    if(maskThickness < 0*um || maskThickness > 1000*um)
    {
      G4cerr << maskThickness << "Mask thickness is out of bounds (Must be > 0 um AND < 1000 um (default 300 um).- Command is ignored." << G4endl;
    }
    else
    {
      fMaskThickness = maskThickness;
    }
    G4cout <<"The mask thickness is "<< fMaskThickness/um <<"um"<< G4endl;
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetM2Pitch(G4double m2Pitch)
{
  if(!fConstructed)
  {
    if(m2Pitch < 0*um || m2Pitch > 1000*um)
    {
      G4cerr << m2Pitch << "Mask pitch is out of bounds (Must be > 0 um AND < 1000 um (default 62 um).- Command is ignored." << G4endl;
    }
    else
    {
      fM2Pitch = m2Pitch;
    }
    G4cout <<"The detector mask pitch is "<< fM2Pitch/um <<"um"<< G4endl;
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetM2Aperture(G4double m2Aperture)
{
  if(!fConstructed)
  {
    if(m2Aperture < 0*um || m2Aperture > fM2Pitch)
    {
      G4cerr << m2Aperture << "Mask aperture is out of bounds (Must be > 0 um AND < pitch (default 15 um).- Command is ignored." << G4endl;
    }
    else
    {
      fM2Aperture = m2Aperture;
    }
    G4cout <<"The detector mask aperture is "<< fM2Aperture/um <<"um"<< G4endl;
  }
  else
  {
    G4cerr << "Cannot change geometry after initialization. - Command is ignored." << G4endl;
    return;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetThreshold(G4double threshold1, G4double threshold2, G4double threshold3, G4double threshold4, G4double threshold5, G4double threshold6, G4double threshold7, G4double threshold8)
{  
  if(!fConstructed)
  {
    fThreshold1=threshold1;
    fThreshold2=threshold2;
    fThreshold3=threshold3;
    fThreshold4=threshold4;
    fThreshold5=threshold5;
    fThreshold6=threshold6;
    fThreshold7=threshold7;
    fThreshold8=threshold8;
  }
  else
  {
    G4cerr << "Cannot change threshold after inizialization. - Command is ignored." << G4endl;
    return;
  }  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::SetBidimensionalAcquisition(G4bool bidimensional)
{
  if(fBidimensional) return;

  fBidimensional = bidimensional;
  //G4cout<< bidimensional << G4endl;
  if(!fConstructed) return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::SetEIMovements(G4double trans, G4double dith, G4double rotAngle)
{
  fTrans    = trans;
  fDith     = dith;
  fRotAngle = rotAngle;
  if(!fConstructed) return;
// Stepping and Dithering of the sample

 // object 1
  G4ThreeVector position1 = G4ThreeVector(0.*mm+fTrans+fDith, 0, fSourcePosZ+fSrcObjDistance);
  Move("Tube", fObjectLogical, fObjectPhysical, position1,fWorldLogical,fRotAngle);

 // object 2
  G4ThreeVector position2 = G4ThreeVector(0.*mm+fTrans+fDith, 0, fSourcePosZ+fSrcObjDistance);
  Move("Tube2", fObject2Logical, fObject2Physical, position2,fWorldLogical,fRotAngle);

 // object 3
  G4ThreeVector position3 = G4ThreeVector(0.*mm+fTrans+fDith, 0, fSourcePosZ+fSrcObjDistance);
  Move("Tube3", fObject3Logical, fObject3Physical, position3,fWorldLogical,fRotAngle);

 // object 4
  G4ThreeVector position4 = G4ThreeVector(0.*mm+fTrans+fDith, 0, fSourcePosZ+fSrcObjDistance);
  Move("Tube4", fObject4Logical, fObject4Physical, position4,fWorldLogical,fRotAngle);

  G4cout<<"Sample translated to " << (fTrans+fDith)/um << " um" <<G4endl;
  
  if(fAcquisitionType=="singlemask" || fAcquisitionType=="doublemask")
  {
 // Translation Sample Mask M1 and substrate
  G4double rel_mag = fSrcObjDistance/(fSrcObjDistance-(fMaskThickness+2*fObjSizeR)/2);
  G4ThreeVector position = G4ThreeVector(fTrans/rel_mag, 0, fSourcePosZ+fSrcObjDistance-(fMaskThickness+2*fObjSizeR)/2);
  Move("EnvelopeM1", fEnvelopeM1Logical, fEnvelopeM1Physical, position, fWorldLogical,0);
  Move("M1sub", fM1subLogical, fM1subPhysical, position-G4ThreeVector(0,0,fSubThickness/2+fMaskThickness/2),fWorldLogical,0);
 // print  
  G4cout<<"Sample Mask translated to " << fTrans/rel_mag/um << " um" <<G4endl;
  }

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiDetectorConstruction::SetObjectMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);

  if(pttoMaterial){
    fObjectMaterial = pttoMaterial;
    if(fConstructed) fObjectLogical->SetMaterial(fObjectMaterial);

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
  else G4cerr << materialChoice << " is not defined. - Command is ignored." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetAcquisitionType(G4String acquisitionType)
{
  fAcquisitionType = acquisitionType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetDetType(G4String detType)
{
  fDetType = detType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<double> PepiDetectorConstruction::LoadDelta(G4String FileName) 
{
G4cout<<"reading delta values from "<< FileName << G4endl;

  std::ifstream myData(FileName, std::ios::binary);
  std::vector<G4double> deltas;
  if(!myData.is_open())//file not open
    {
        G4cout<< FileName << "file does not exist in the specified path" << G4endl;
	return deltas;
    }
  double num = 0.0;
  while (myData >> num){
      deltas.push_back(num);
  }
  return deltas;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::tuple<G4LogicalVolume*,G4VPhysicalVolume*> PepiDetectorConstruction::CreateMask(G4String name, G4double mag, G4double pitch, G4double aperture, G4ThreeVector position, G4double thickness, G4Material* material, G4LogicalVolume* EnvelopeLogical, G4VPhysicalVolume* EnvelopePhysical)
{

  // - Build the MASK APERTURE UNIT as a Box
  G4Box* MSolid =     new G4Box(name,                               //its name
                          ((pitch-aperture)/mag)/2,   //its size
                          1.1*fPixiRadSizeY/2,
                          thickness/2);

  G4String lvname = name+"LV";                   
  G4LogicalVolume* MLogical = new G4LogicalVolume(MSolid,       //its solid
                                   material,  //its material
                                   lvname);        //its name

  // Build the MASK ENVELOPE 
  G4String envname = "Envelope"+name;                     
  G4Box* EnvelopeSolid =  new G4Box(envname,                  //its name                 
                        (1.1*fPixiRadSizeX/mag)/2,              //its size
                        1.1*fPixiRadSizeY/2,
                        thickness/2);        
  G4String lvenvname = "Envelope"+name+"LV";                         
  EnvelopeLogical =  new G4LogicalVolume(EnvelopeSolid,    //its solid
                                            fWorldMaterial,      //its material
                                            lvenvname);     //its name

  // - Place the physical copies of the mask aperture unit
   G4int copy_no=0;  

    for (G4int iX = 0; iX < int(fnPixelsX*fPixelSizeX/pitch)+2; iX++)
    {
          
      G4double x = + iX*pitch/mag - ((fPixiRadSizeX/mag)/2 + (pitch)/mag + ((pitch)/mag)/2);
      G4double y = 0;
 // G4cout << "position \n" << x <<G4endl;
      
      G4ThreeVector px_position = G4ThreeVector(x, y, 0);
      G4String  name1 = name + "_" + G4UIcommand::ConvertToString(copy_no);

      /*G4VPhysicalVolume* MPhysical =*/  new G4PVPlacement(0,                           //its rotation
                                       px_position,                    //its position
                                       MLogical,                  //its logical volume
                                       name1,                        //its name
                                       EnvelopeLogical,          //its mother volume
                                       false,                       //no boolean operation
                                       copy_no,                     //copy number
                                       fCheckOverlaps);             //checking overlaps 
      copy_no++;                              
    }

  // - Place the Sample Mask Envelope in the World
  EnvelopePhysical = new G4PVPlacement(0,                                                  //its rotation
                                          position,   				      //its position
                                          EnvelopeLogical,                                    //its logical volume
                                          envname,                                       //its name
                                          fWorldLogical,                                      //its mother volume
                                          false,                                              //no boolean operation
                                          0,                                                  //copy number
                                          fCheckOverlaps);                                    //checking overlaps


  G4Region* aRegion = new G4Region(name);
  EnvelopeLogical->SetRegion(aRegion);
  aRegion->AddRootLogicalVolume(EnvelopeLogical);
  
  

  G4VisAttributes* MVisAtt = new G4VisAttributes(G4Colour(0.8,0.6,0.));
//  M1VisAtt->SetForceWireframe(false);
  MVisAtt->SetForceSolid(true);
  MVisAtt->SetVisibility(true);
  // pixelVisAtt->SetForceAuxEdgeVisible(true);
  MLogical->SetVisAttributes(MVisAtt);

  G4VisAttributes* envelopeMVisAtt = new G4VisAttributes(G4Colour(1.0,0.,0.));
//  envelopeM1VisAtt->SetForceWireframe(false);
  envelopeMVisAtt->SetForceSolid(true);
  envelopeMVisAtt->SetVisibility(false);
  envelopeMVisAtt->SetForceAuxEdgeVisible(false);
  EnvelopeLogical->SetVisAttributes(envelopeMVisAtt);
  
  return std::make_tuple(EnvelopeLogical, EnvelopePhysical);  
  
}//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::tuple<G4LogicalVolume*,G4VPhysicalVolume*> PepiDetectorConstruction::CreateSubstrate(G4String name, G4double mag, G4ThreeVector position, G4double thickness, G4Material* material, G4LogicalVolume* SubLogical, G4VPhysicalVolume* SubPhysical)
{

  // - Build the substrate as a box
  G4Box* SubSolid = new G4Box(name,                               //its name
                          (1.2*fPixiRadSizeX/mag)/2,   //its size
                          (1.2*fPixiRadSizeY)/2,
                          thickness/2);
  G4String lvname = name + "LV";
     
  SubLogical = new G4LogicalVolume(SubSolid,       //its solid
                                       material,    //its material
                                       lvname);          //its name
 

  SubPhysical = new G4PVPlacement(0,              //its rotation
                                      position,       //its translation
                                      SubLogical,      //its logical volume
                                      name,              //its name
                                      fWorldLogical,       //its mother volume
                                      false,               //no boolean operation
                                      0,                   //copy number
                                      fCheckOverlaps);     //checking overlaps

  G4VisAttributes* SubVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
//  M1subVisAtt->SetForceWireframe(false);
  SubVisAtt->SetForceSolid(true);
  SubVisAtt->SetVisibility(true);
  SubVisAtt->SetForceAuxEdgeVisible(false);
  SubLogical->SetVisAttributes(SubVisAtt);
  
 return std::make_tuple(SubLogical, SubPhysical);
}  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiDetectorConstruction::Move(G4String name, G4LogicalVolume* Logical, G4VPhysicalVolume* Physical, G4ThreeVector position,G4LogicalVolume* WLogical,G4double rot)
{
// for lateral translation of masks and sample
G4RotationMatrix* rotMatt =  new G4RotationMatrix();
  rotMatt->rotateX(rot);

  Logical->RemoveDaughter(Physical);
  delete Physical;
  Physical = new G4PVPlacement(rotMatt,                                                  //its rotation
                                          position,   				      //its position
                                          Logical,                                    //its logical volume
                                          name,                                       //its name
                                          WLogical,                                      //its mother volume
                                          false,                                              //no boolean operation
                                          0,                                                  //copy number
                                          fCheckOverlaps);                                    //checking overlaps                                
}
void PepiDetectorConstruction::Move1(G4String name, G4LogicalVolume* Logical, G4VPhysicalVolume* Physical, G4ThreeVector position,G4LogicalVolume* WLogical,G4double rot)
{
// for lateral translation of masks and sample
G4RotationMatrix* rotMatt =  new G4RotationMatrix();
  rotMatt->rotateZ(rot);

  Logical->RemoveDaughter(Physical);
  delete Physical;
  Physical = new G4PVPlacement(rotMatt,                                                  //its rotation
                                          position,   				      //its position
                                          Logical,                                    //its logical volume
                                          name,                                       //its name
                                          WLogical,                                      //its mother volume
                                          false,                                              //no boolean operation
                                          0,                                                  //copy number
                                          fCheckOverlaps);                                    //checking overlaps                                
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
