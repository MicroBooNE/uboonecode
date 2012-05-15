typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;

microboone_geo(TString volName="volTPC")
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  TGeoManager::Import("microboone.gdml");

  drawopt optuboone[] = {
    {"volGround",       kOrange-7},
    {"volOverburden",       kOrange-7},
    {"volConcreteEnclosure", kGray},
    {"volConcreteEnclosureBottom", kGray},
    {0, 0}
  };

  for (int i=0;; ++i) {
     if (optuboone[i].volume==0) break;
       gGeoManager->FindVolumeFast(optuboone[i].volume)->SetLineColor(optuboone[i].color);
  }
  TList* mat = gGeoManager->GetListOfMaterials();
  TIter next(mat);
  TObject *obj;
  while (obj = next()) {
    obj->Print();
  }

  gGeoManager->GetTopNode();
  gGeoManager->CheckOverlaps(10e-24);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  gGeoManager->FindVolumeFast(volName)->Draw("ogl");

  cout << "LAr weight in TPC = "<< endl;
  float m_tpc = gGeoManager->FindVolumeFast("volTPCActive")->Weight();

//   TFile *tf = new TFile("microboone.root", "RECREATE");
 
//   gGeoManager->Write();

//   tf->Close();
}
