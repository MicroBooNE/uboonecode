/** Check masses of materials in MicroBooNE geometry.

    Defines mass analytic calculation functions that keep
    track of how much mass there is of each material.

    For best performance, load this in ROOT using 
      .L microboone_masses.C+
 */

#include <iostream>  // for cout
#include <cmath>
#include <vector>
#include <TROOT.h>
#include <TSystem.h> // for gSystem
#include <TList.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TString.h>
#include <TRandom.h>
#include <TError.h>
#include <TGeoBBox.h>
#include <TH2F.h>

/** This function does analytic mass calculations.
    It relies on TGeoVolume::Capacity().
    It is similar to ROOT's TGeoVolume::WeightA(), except
     - it finds mass of each material inside the volume;
     - it doesn't change TGeoManager's "top volume";
     - it includes mass from materials with any density, including gasses.
*/
void FindMassAnalytic 
(TGeoVolume* vol,              /**< volume to scan */
 std::vector<double> & masses, /**< array of masses indexed by material index */
 double & total_mass,          /**< destination for total mass */
 bool update = false           /**< false=calculate mass; true=add to mass */) 
{
  double capacity = vol->Capacity();
  if (!update) {
    int nmat = vol->GetGeoManager()->GetListOfMaterials()->GetSize();
    masses.resize(nmat);
    masses.clear();
    total_mass = 0.0;
  }
  int nd = vol->GetNdaughters();
  for (int i=0; i<nd; i++) {
    TGeoVolume* dvol = vol->GetNode(i)->GetVolume();
    capacity -= dvol->Capacity();
    FindMassAnalytic(dvol, masses, total_mass, true);
  }
  double density = 0.0;
  int imat = -1;
  if (! vol->IsAssembly() ) {
    TGeoMaterial* mat = vol->GetMaterial();
    if (mat) {
      density = mat->GetDensity();
      imat = mat->GetIndex();
      //-- Note:  ROOT's "TGeoVolume::WeightA()" has a hard-coded
      //   cut on density of 0.01, below which the density is zeroed out.
      //   The code has a short explanatory comment: "do not weight gasses".
      //   This code does not have that cut. =GAHS.
    }
  }
  double dm = density * capacity * 0.001; // 0.001 is to convert to kg
  total_mass += dm;
  if (imat >= 0)
    masses[imat] += dm;
}


/** This function does a MC integration to estimate the mass of the specified
    volume to the given fractional precision, and also estimates the mass
    of each material in the volume (not necessarily to that precision).

    Optionally, it makes histograms showing distribution of the mass.

    It relies on TGeoManager::GetFindNode(x,y,z)

    It is based on ROOT's TGeoChecker::Weight(), and is very similar, except
     - it finds mass of each material inside the volume;
     - it includes mass from materials with any density, including gasses;
     - it has the optional feature of histogramming the mass location.

    Note unlike FindMassAnalytic(), this function does temporarily
    change TGeoManager's "top volume", and therefore might not be thread-safe.

    The mass location histograms are orthographic projections in xy, zx, and zy,
    one of each for each material. They are saved in TObjArray whose total
    length is 3*nmat, where nmat is the number of materials.  The index is
    3*imat + ixyz, where imat is the material index and 
    ixyz = 0 for xy, 1 for zx, and 2 for zy.
*/
void FindMassMC
(TGeoVolume* vol,              /**< volume to scan */
 std::vector<double> & masses, /**< array of masses indexed by material index */
 std::vector<double> & sigma_m, /**< array of mass uncertainty indexed by material index */
 double & total_mass,          /**< destination for total mass */
 double & sigma_total_mass,    /**< destination for total mass uncertainty */
 double precision=0.01,        /**< fractional precision desired */
 TObjArray* histos=0,          /**< pointer to array to receive histograms;
				    set to NULL if no histograms wanted */
 TGeoBBox *userbbox=0          /**< optional pointer to a different bounding
                                    box to use. Useful for zooming in to
                                    a particular sub volume of the volume */
 )
{
  TGeoManager* geoManager = vol->GetGeoManager();
  TList *matlist = geoManager->GetListOfMaterials();
  Int_t nmat = matlist ? matlist->GetSize() : 0;
  if (!nmat) {
    Error("FindMassMC", "no materials defined");
    return;
  }
  Int_t *nin = new Int_t[nmat];
  memset(nin, 0, nmat*sizeof(Int_t));
  TGeoBBox *box;
  if (userbbox) {
    box = userbbox;
  }
  else {
    TGeoShape *shape = vol->GetShape();
    if (!shape->InheritsFrom(TGeoBBox::Class())) {
      Error("FindMassMC", "passed a volume whose shape isn't a TGeoBBox");
      return;
    }
    box = (TGeoBBox *)shape;
  }
  Double_t dx = box->GetDX();
  Double_t dy = box->GetDY();
  Double_t dz = box->GetDZ();
  Double_t ox = (box->GetOrigin())[0];
  Double_t oy = (box->GetOrigin())[1];
  Double_t oz = (box->GetOrigin())[2];
  Double_t x,y,z;
  TGeoNode *node;
  TGeoMaterial *mat;
  Double_t vbox = 0.000008*dx*dy*dz; // m3
  Bool_t end = kFALSE;
  Double_t weight=0, sigma, eps, dens;
  Double_t eps0=1.;
  Int_t indmat;
  Int_t igen=0;
  Int_t iin = 0;
  if (histos) {
    //-- extend histogram array if needed
    if (histos->Capacity() < 3*nmat)
      histos->Expand(3*nmat);
    histos->SetLast(3*nmat-1);
    for (indmat=0; indmat<nmat; indmat++) {
      mat = (TGeoMaterial*)matlist->At(indmat);
      (*histos)[3*indmat] = new TH2F
	(Form("hxy_%s",mat->GetName()),
	 Form("XY (end view) distribution of %s in %s;x [cm];y [cm]", 
	      mat->GetName(), vol->GetName()),
	 100, ox-dx, ox+dx, 100, oy-dy, oy+dy);
      (*histos)[3*indmat+1] = new TH2F
	(Form("hzx_%s",mat->GetName()),
	 Form("ZX (top view) distribution of %s in %s;z [cm];x [cm]", 
	      mat->GetName(), vol->GetName()),
	 100, oz-dz, oz+dz, 100, ox-dx, ox+dx);
      (*histos)[3*indmat+2] = new TH2F
	(Form("hzy_%s",mat->GetName()),
	 Form("ZY (side view) distribution of %s in %s;z [cm];y [cm]", 
	      mat->GetName(), vol->GetName()),
	 100, oz-dz, oz+dz, 100, oy-dy, oy+dy);
    }
  }
  //-- temporarily set this volume as the top in its geoManager
  TGeoVolume * oldtop = geoManager->GetTopVolume();
  if (oldtop != vol)
    geoManager->SetTopVolume(vol);
  //-- *** be sure not to return from function without resetting top volume ***
  while (!end) {
    x = ox-dx+2*dx*gRandom->Rndm();
    y = oy-dy+2*dy*gRandom->Rndm();
    z = oz-dz+2*dz*gRandom->Rndm();
    node = geoManager->FindNode(x,y,z);
    igen++;
    if (!node) continue;
    mat = node->GetVolume()->GetMedium()->GetMaterial();
    indmat = mat->GetIndex();
    if (indmat<0) continue;
    nin[indmat]++;
    iin++;
    if (histos) {
      ((TH2F*)(histos->At(indmat*3)))->Fill(x,y);
      ((TH2F*)(histos->At(indmat*3+1)))->Fill(z,x);
      ((TH2F*)(histos->At(indmat*3+2)))->Fill(z,y);
    }
    if ((iin%100000)==0 || igen>1E8) {
      weight = 0;
      sigma = 0;
      for (indmat=0; indmat<nmat; indmat++) {
	mat = (TGeoMaterial*)matlist->At(indmat);
	dens = mat->GetDensity(); //  [g/cm3]
	// if (dens<1E-2) dens=0;  <<<--- Commented out so gasses are included.
	dens *= 1000.;            // [kg/m3]
	weight += dens*Double_t(nin[indmat]);
	sigma  += dens*dens*nin[indmat];
      }
      sigma = sqrt(sigma);
      eps = sigma/weight;
      weight *= vbox/Double_t(igen);
      sigma *= vbox/Double_t(igen);
      if (eps<precision || igen>1E8) {
	end = kTRUE;                      
      }
    }
  }
  //-- reset top node
  if (oldtop) 
    geoManager->SetTopVolume(oldtop);
  //-- clean up and store return values
  delete [] nin;
  total_mass = weight;
  sigma_total_mass = sigma;
  masses.resize(nmat);
  sigma_m.resize(nmat);
  for (indmat=0; indmat<nmat; indmat++) {
    mat = (TGeoMaterial*)matlist->At(indmat);
    dens = mat->GetDensity(); //  [g/cm3]
    dens *= 1000.;            // [kg/m3]
    masses[indmat] = dens*vbox*Double_t(nin[indmat])/Double_t(igen);
    sigma_m[indmat] = dens*vbox*sqrt(nin[indmat])/Double_t(igen);
  }
}

/** This function loads the geometry, does the checks, and
    prints the results.  Optional volName argument picks what
    volume to check. */
void microboone_masses(TString volName="volCryostat")
{
  using namespace std;

  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  TGeoManager::Import("microboonevX.gdml");

  TList* matlist = gGeoManager->GetListOfMaterials();
  matlist->Print();

  gGeoManager->GetTopNode();

  TGeoVolume *vol = gGeoManager->FindVolumeFast(volName);

  vector<double> masses;
  double total_mass = 0.0;

  FindMassAnalytic(vol, masses, total_mass);

  cout << "\nMasses in kg of each material inside " << vol->GetName() << ":\n";
  TIter matit(matlist);
  for (TGeoMaterial* mat=0; (mat=(TGeoMaterial*)matit()) != NULL; ) {
    int imat = mat->GetIndex();
    cout << " " << mat->GetName() << "\t" << masses[imat] << endl;
  }

  cout << "Total mass inside " << vol->GetName() << " = " << total_mass 
       << " kg.\n";

  vector<double> mc_masses;
  vector<double> mc_sigma_m;
  double mc_total_mass = 0.0;
  double sigma_mc_total_mass = 0.0;
  TObjArray* histos = new TObjArray();
  TFile *tf = new TFile("microboone_masses.root", "RECREATE");

  FindMassMC(vol, mc_masses, mc_sigma_m, mc_total_mass, sigma_mc_total_mass,
	     0.5e-3, histos);
  
  cout << "\nMC estimate of masses in kg of each material inside " << vol->GetName() << ":\n";
  matit.Reset();
  for (TGeoMaterial* mat=0; (mat=(TGeoMaterial*)matit()) != NULL; ) {
    int imat = mat->GetIndex();
    cout << " " << mat->GetName() << "\t" << mc_masses[imat] 
	 << " +/- " << mc_sigma_m[imat] << endl;
  }

  cout << "MC estimate of total mass inside " << vol->GetName() << " = " 
       << mc_total_mass << " +/- " << sigma_mc_total_mass << " kg.\n";

  histos->Write();
  gGeoManager->Write();

  tf->Close();
}
