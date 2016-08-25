{

  gStyle->SetOptStat(0);
  
  TFile input_file("~wketchum/larsoft_development/PatternMatching/ac_trakana_output_chain.root");
  TTree *anatree;
  input_file.GetObject("ac_trkana_tree",anatree);

  TCanvas c2D("c2D","2D histo");
  TCanvas ceff("ceff","Efficiency/Mistag Rate");
  TCanvas croc("croc","ROC curve");
  
  anatree->SetMarkerStyle(6);
  TH2F* h2D = new TH2F("h2D","Track #Deltax vs. Matched Fraction; matched fraction; #Deltax (cm)",
		       101,0.0,1.1,100.,0.0,300.);
  h2D->SetMarkerStyle(6);
  anatree->Project("h2D","abs(track_end_x-track_start_x):pm_fraction");
  c2D.cd();
  h2D->Draw();

  const int NPOINTS = 101;
  const double stepsize = 1./(double)(NPOINTS-1);
  double frac[NPOINTS],frac_err[NPOINTS];
  double eff[NPOINTS],eff_err[NPOINTS];
  double mtag[NPOINTS],mtag_err[NPOINTS];
  double pur[NPOINTS],pur_err[NPOINTS];
  
  
  const TString min_dx_value("210.");
  const TString dx_var("abs(track_end_x-track_start_x)");

  const TString pass_sel = dx_var + ">" + min_dx_value;
  const TString fail_sel = dx_var + "<" + min_dx_value;
  const double total_pass = anatree->GetEntries(pass_sel);
  const double total_fail = anatree->GetEntries(fail_sel);
  const double total = total_pass+total_fail;
  
  for(int i=0; i<NPOINTS; ++i){
    frac_err[i] = 0;
    frac[i] = (double)i*stepsize;

    TString pmCut; pmCut.Form(" && pm_fraction>%f",(float)frac[i]);
    double pass = anatree->GetEntries(pass_sel+pmCut);
    double fail = anatree->GetEntries(fail_sel+pmCut);

    printf("frac=%f, pass=%f / %f, fail=%f / %f\n",frac[i],pass,total_pass,fail,total_fail);
    
    eff[i] = pass/total_pass; mtag[i] = fail/total_fail; pur[i] = pass/(pass+fail);
    eff_err[i] = 0.0; mtag_err[i] = 0.0; pur_err[i] = 0.0;
  }

  ceff.cd();
  
  TGraphErrors* gr_eff = new TGraphErrors(NPOINTS,frac,eff,frac_err,eff_err);
  gr_eff->SetName("gr_eff");
  gr_eff->SetTitle(";matched fraction;");
  gr_eff->SetMarkerColor(kBlue);
  gr_eff->SetMarkerStyle(21);
  gr_eff->SetMarkerSize(0.8);
  gr_eff->Draw("ALP");

  TGraphErrors* gr_mtag = new TGraphErrors(NPOINTS,frac,mtag,frac_err,mtag_err);
  gr_mtag->SetName("gr_mtag");
  gr_mtag->SetTitle(";matched fraction;");
  gr_mtag->SetMarkerColor(kRed);
  gr_mtag->SetMarkerStyle(21);
  gr_mtag->SetMarkerSize(0.8);
  gr_mtag->Draw("LP");

  leg_eff = new TLegend(0.2,0.2,0.5,0.4);
  leg_eff->AddEntry("gr_eff","Efficiency","lep");
  leg_eff->AddEntry("gr_mtag","Mistag Rate","lep");
  leg_eff->Draw();

  croc.cd();

  TGraphErrors* gr_roc = new TGraphErrors(NPOINTS-1,pur,eff,pur_err,eff_err);
  gr_roc->SetName("gr_roc");
  gr_roc->SetTitle("ROC Curve;purity;efficiency");
  gr_roc->SetMarkerColor(kRed);
  gr_roc->SetMarkerStyle(21);
  gr_roc->SetMarkerSize(0.8);
  gr_roc->Draw("ALP");


  TString pmCut_pass; pmCut_pass.Form(" && pm_fraction>%f",0.85);
  TString pmCut_fail; pmCut_fail.Form(" && pm_fraction<%f",0.85);

  const int NBINS = 20;

  TCanvas c_startx("c_startx");
  c_startx.cd();

  TH1F* h_start_x_pass = new TH1F("h_start_x_pass","Track start x; x (cm); Entries (Normalized)",
				  NBINS,-70.0,370.);
  h_start_x_pass->SetLineColor(kBlue);
  TH1F* h_start_x_fail = new TH1F("h_start_x_fail","Track start x; x (cm); Entries (Normalized)",
				  NBINS,-70.0,370.);
  h_start_x_fail->SetLineColor(kRed);
  anatree->Project("h_start_x_pass","track_start_x",pass_sel+pmCut_pass);
  anatree->Project("h_start_x_fail","track_start_x",pass_sel+pmCut_fail);
  
  h_start_x_pass->DrawNormalized();
  h_start_x_fail->DrawNormalized("same");
  
  leg_startx = new TLegend(0.55,0.7,0.85,0.85);
  leg_startx->AddEntry("h_start_x_pass","match > 85%","l");
  leg_startx->AddEntry("h_start_x_fail","match < 85%","l");
  leg_startx->Draw();

  TCanvas c_endx("c_endx");
  c_endx.cd();

  TH1F* h_end_x_pass = new TH1F("h_end_x_pass","Track end x; x (cm); Entries (Normalized)",
				  NBINS,-70.0,370.);
  h_end_x_pass->SetLineColor(kBlue);
  TH1F* h_end_x_fail = new TH1F("h_end_x_fail","Track end x; x (cm); Entries (Normalized)",
				  NBINS,-70.0,370.);
  h_end_x_fail->SetLineColor(kRed);
  anatree->Project("h_end_x_pass","track_end_x",pass_sel+pmCut_pass);
  anatree->Project("h_end_x_fail","track_end_x",pass_sel+pmCut_fail);
  
  h_end_x_pass->DrawNormalized();
  h_end_x_fail->DrawNormalized("same");

  leg_endx = new TLegend(0.55,0.7,0.85,0.85);
  leg_endx->AddEntry("h_end_x_pass","match > 85%","l");
  leg_endx->AddEntry("h_end_x_fail","match < 85%","l");
  leg_endx->Draw();


  TCanvas c_starty("c_starty");
  c_starty.cd();

  TH1F* h_start_y_pass = new TH1F("h_start_y_pass","Track start y; y (cm); Entries (Normalized)",
				  NBINS,-200.0,200.);
  h_start_y_pass->SetLineColor(kBlue);
  TH1F* h_start_y_fail = new TH1F("h_start_y_fail","Track start y; y (cm); Entries (Normalized)",
				  NBINS,-200.0,200.);
  h_start_y_fail->SetLineColor(kRed);
  anatree->Project("h_start_y_pass","track_start_y",pass_sel+pmCut_pass);
  anatree->Project("h_start_y_fail","track_start_y",pass_sel+pmCut_fail);
  
  h_start_y_pass->DrawNormalized();
  h_start_y_fail->DrawNormalized("same");
  
  leg_starty = new TLegend(0.55,0.7,0.85,0.85);
  leg_starty->AddEntry("h_start_y_pass","match > 85%","l");
  leg_starty->AddEntry("h_start_y_fail","match < 85%","l");
  leg_starty->Draw();

  TCanvas c_endy("c_endy");
  c_endy.cd();

  TH1F* h_end_y_pass = new TH1F("h_end_y_pass","Track end y; y (cm); Entries (Normalized)",
				  NBINS,-200.0,200.);
  h_end_y_pass->SetLineColor(kBlue);
  TH1F* h_end_y_fail = new TH1F("h_end_y_fail","Track end y; y (cm); Entries (Normalized)",
				  NBINS,-200.0,200.);
  h_end_y_fail->SetLineColor(kRed);
  anatree->Project("h_end_y_pass","track_end_y",pass_sel+pmCut_pass);
  anatree->Project("h_end_y_fail","track_end_y",pass_sel+pmCut_fail);
  
  h_end_y_pass->DrawNormalized();
  h_end_y_fail->DrawNormalized("same");

  leg_endy = new TLegend(0.55,0.7,0.85,0.85);
  leg_endy->AddEntry("h_end_y_pass","match > 85%","l");
  leg_endy->AddEntry("h_end_y_fail","match < 85%","l");
  leg_endy->Draw();

  TCanvas c_startz("c_startz");
  c_startz.cd();

  TH1F* h_start_z_pass = new TH1F("h_start_z_pass","Track start z; z (cm); Entries (Normalized)",
				  NBINS,-100.0,1100.);
  h_start_z_pass->SetLineColor(kBlue);
  TH1F* h_start_z_fail = new TH1F("h_start_z_fail","Track start z; z (cm); Entries (Normalized)",
				  NBINS,-100.0,1100.);
  h_start_z_fail->SetLineColor(kRed);
  anatree->Project("h_start_z_pass","track_start_z",pass_sel+pmCut_pass);
  anatree->Project("h_start_z_fail","track_start_z",pass_sel+pmCut_fail);
  
  h_start_z_pass->DrawNormalized();
  h_start_z_fail->DrawNormalized("same");
  
  leg_startz = new TLegend(0.55,0.7,0.85,0.85);
  leg_startz->AddEntry("h_start_z_pass","match > 85%","l");
  leg_startz->AddEntry("h_start_z_fail","match < 85%","l");
  leg_startz->Draw();

  TCanvas c_endz("c_endz");
  c_endz.cd();

  TH1F* h_end_z_pass = new TH1F("h_end_z_pass","Track end z; z (cm); Entries (Normalized)",
				  NBINS,-100.0,1100.);
  h_end_z_pass->SetLineColor(kBlue);
  TH1F* h_end_z_fail = new TH1F("h_end_z_fail","Track end z; z (cm); Entries (Normalized)",
				  NBINS,-100.0,1100.);
  h_end_z_fail->SetLineColor(kRed);
  anatree->Project("h_end_z_pass","track_end_z",pass_sel+pmCut_pass);
  anatree->Project("h_end_z_fail","track_end_z",pass_sel+pmCut_fail);
  
  h_end_z_pass->DrawNormalized();
  h_end_z_fail->DrawNormalized("same");

  leg_endz = new TLegend(0.55,0.7,0.85,0.85);
  leg_endz->AddEntry("h_end_z_pass","match > 85%","l");
  leg_endz->AddEntry("h_end_z_fail","match < 85%","l");
  leg_endz->Draw();

  TCanvas c_startz_trg("c_startz_trg");
  c_startz_trg.cd();

  TH1F* h_start_z_trg_pass = new TH1F("h_start_z_trg_pass","Track start z; z (cm); Entries (Normalized)",
				  NBINS,-100.0,1100.);
  h_start_z_trg_pass->SetLineColor(kBlue);
  TH1F* h_start_z_trg_fail = new TH1F("h_start_z_trg_fail","Track start z; z (cm); Entries (Normalized)",
				  NBINS,-100.0,1100.);
  h_start_z_trg_fail->SetLineColor(kRed);
  anatree->Project("h_start_z_trg_pass","track_start_z",pass_sel+pmCut_pass);
  anatree->Project("h_start_z_trg_fail","track_start_z",fail_sel+pmCut_pass);
  
  h_start_z_trg_pass->DrawNormalized();
  h_start_z_trg_fail->DrawNormalized("same");
  
  leg_startz_trg = new TLegend(0.55,0.7,0.85,0.85);
  leg_startz_trg->AddEntry("h_start_z_trg_pass","match > 85%, #DeltaX > 210","l");
  leg_startz_trg->AddEntry("h_start_z_trg_fail","match > 85%, #DeltaX < 210","l");
  leg_startz_trg->Draw();
  /*
  TCanvas c_endz("c_endz");
  c_endz.cd();

  TH1F* h_end_z_pass = new TH1F("h_end_z_pass","Track end z; z (cm); Entries (Normalized)",
				  NBINS,-100.0,1100.);
  h_end_z_pass->SetLineColor(kBlue);
  TH1F* h_end_z_fail = new TH1F("h_end_z_fail","Track end z; z (cm); Entries (Normalized)",
				  NBINS,-100.0,1100.);
  h_end_z_fail->SetLineColor(kRed);
  anatree->Project("h_end_z_pass","track_end_z",pass_sel+pmCut_pass);
  anatree->Project("h_end_z_fail","track_end_z",pass_sel+pmCut_fail);
  
  h_end_z_pass->DrawNormalized();
  h_end_z_fail->DrawNormalized("same");

  leg_endz = new TLegend(0.55,0.7,0.85,0.85);
  leg_endz->AddEntry("h_end_z_pass","match > 85%","l");
  leg_endz->AddEntry("h_end_z_fail","match < 85%","l");
  leg_endz->Draw();
*/
}
