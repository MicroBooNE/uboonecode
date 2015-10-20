// E. Cheu
// Test macro to run in gmbrowser. 
//
// Divides canvas in two and plots a data histogram in one pad
//   and a reference histogram in the other pad.
void test(TCanvas *c1, TFile *data, TFile *ref)
{
  // Get histograms.
  TH1* h = data->Get("h_electrons_N");
  TH1* r = ref->Get("h_electrons_N");
  if (!h && !r) return;

  // Divide screen in two.
  c1->Divide(1, 2);

  // Plot first histogram.
  c1->cd(1);
  h->Draw();

  // Plot second histogram.
  c1->cd(2);
  r->SetLineColor(2);
  r->Draw();
}
