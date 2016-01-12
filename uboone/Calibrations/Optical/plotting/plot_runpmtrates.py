import os,sys
import ROOT as rt
#import numpy as np

input = rt.TFile(sys.argv[1],"OPEN")
eventtree = input.Get( "specalib/eventtree" )

out = rt.TFile("output_rateplots.root", "RECREATE")

pois = rt.TF1("pois","[0]*TMath::Poisson(x,[1])",0,50)

NCHANS = 32
nsamples = 1500.0
hchrate_run = {}
hchrate_tot = {}
hrate = rt.TH1D( "hrate", ";channel;Rate", 32, 0, 32 )
c = rt.TCanvas("c","",800,500)
c.Draw()
for ch in range(0,NCHANS):
    hname = "hrate_run_ch%d"%(ch)
    hchrate_run[ch] = rt.TH2D( hname, ";hours after start of run;npulses in beam window", 24, 0, 24, 50, 0, 50 )
    eventtree.Draw( "nchfires[%d]:event_unixtime*1.0e-6/(3600.0)>>%s"%(ch,hname), "chmaxamp < 2048+500" )
    hchrate_tot[ch] = hchrate_run[ch].ProjectionY( hname.replace("run","tot") )
    pois.SetParameter( 0, hchrate_tot[ch].GetMaximum() )
    pois.SetParameter( 1, hchrate_tot[ch].GetMean() )
    hchrate_tot[ch].Draw()
    hchrate_tot[ch].Fit( pois, "R", "", 0, 20 )
    fitmean = pois.GetParameter(1)
    fiterr  = pois.GetParError(1)
    eventrate = (fitmean/float(15.625e-9*nsamples))*1.0e-3
    rate_err = (fiterr/float(15.625e-9*nsamples))*1.0e-3
    hrate.SetBinContent( ch+1, eventrate )
    hrate.SetBinError( ch+1, rate_err )
    c.Update()
    c.SaveAs("figs/npulses_ch%d.png"%(ch))

    

out.Write()
