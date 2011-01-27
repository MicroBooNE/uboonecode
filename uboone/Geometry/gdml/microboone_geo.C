typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;

microboone_geo(TString volName="")
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
  gGeoManager->CheckOverlaps(0.01);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  //gGeoManager->GetTopVolume()->Draw();
  //if ( ! volName.IsNull() ) gGeoManager->FindVolumeFast(volName)->Draw("ogl");
  gGeoManager->FindVolumeFast("volWorld")->Draw("ogl");

  TFile *tf = new TFile("microboone.root", "RECREATE");
 
  gGeoManager->Write();

  tf->Close();
}
