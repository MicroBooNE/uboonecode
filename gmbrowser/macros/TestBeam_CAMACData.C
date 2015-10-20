void Veto_System_Monitor(TCanvas *c, TFile *data, TFile *ref)
{    
    // There are 16 Graphs
    TGraph** gr_array = new TGraph*[16];
    
    // Get Graphs from File
    for(int i = 0; i < 16; i++ ){
        gr_array[i] = (TGraph*)data->Get(Form("CAMAC_Data/Veto_System_Monitor_Scalar%d",i+1));    
    }

    // Draw Graphs - All on Single Page
    c->Divide(4,4);

    for(int i = 0; i < 16; i++ ){
        c->cd(i+1);
        gr_array[i]->Draw("AP");
    }
}

void Fire_TOF(TCanvas *c, TFile *data, TFile *ref)
{    
    TGraph* gr1 = (TGraph*)data->Get("CAMAC_Data/Time_of_flight_tdc_fire");

    c->cd(1); 
    gr1->SetFillColor(38);
    gr1->Draw("AB");
}

void Fire_Veto(TCanvas *c, TFile *data, TFile *ref)
{
    TGraph* gr1 = (TGraph*)data->Get("CAMAC_Data/Veto_System_Fire_Scalar");
    
    c->cd(1); 
    gr1->SetFillColor(38);
    gr1->Draw("AB");
}

