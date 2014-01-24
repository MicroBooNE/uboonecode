typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;

microboone_granite_geo(TString volName="")
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  TGeoManager::Import("microboone-granite.gdml");

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
  //gGeoManager->CheckOverlaps(0.0000001);
  gGeoManager->CheckOverlaps(10e-24);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  //gGeoManager->GetTopVolume()->Draw();
  //if ( ! volName.IsNull() ) gGeoManager->FindVolumeFast(volName)->Draw("ogl");
  //gGeoManager->FindVolumeFast("volWorld")->Draw("ogl");
  gGeoManager->FindVolumeFast("volTPC")->Draw("ogl");
  //gGeoManager->FindVolumeFast("volTPCPlane")->Draw("ogl");

  TGeoVolume *TPC = gGeoManager->FindVolumeFast("volTPC");
  float m_tpc = TPC->Weight();
  TGeoVolume *Cathode = gGeoManager->FindVolumeFast("volCathodePlate");
  float m_cathode = Cathode->Weight();
  TGeoVolume *Ground = gGeoManager->FindVolumeFast("volGroundPlate");
  float m_ground = Ground->Weight();
  TGeoVolume *UVPlane = gGeoManager->FindVolumeFast("volTPCPlane");
  float m_uvplane = UVPlane->Weight();
  TGeoVolume *YPlane = gGeoManager->FindVolumeFast("volTPCPlaneVert");
  float m_yplane = YPlane->Weight();
  TGeoVolume *FieldCageH = gGeoManager->FindVolumeFast("volFieldCageTubeTop");
  float m_fchoriz = FieldCageH->Weight();
  TGeoVolume *FieldCageV = gGeoManager->FindVolumeFast("volFieldCageTubeFront");
  float m_fcvert = FieldCageV->Weight();

  float m_tpc_argon = m_tpc - ( m_cathode + m_ground + 2*m_uvplane + m_yplane + 50*(m_fchoriz + m_fcvert));
  //float m_tpc_argon = m_tpc - m_yplane;
  cout << "LAr weight in TPC = " << m_tpc_argon << " kg\n" <<endl;

  TFile *tf = new TFile("microboone-granite.root", "RECREATE");
 
  gGeoManager->Write();

  tf->Close();
}
