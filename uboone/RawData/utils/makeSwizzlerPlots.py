import ROOT
import array
import sys

#ROOT.gStyle.SetOptStat(0)

if len(sys.argv) < 2:
  ### Andy's temporary run number inputs
 runNumbers = [2316,2317,2319,2364,2365, 2398]
 runNumbers = [2365]
 runNumbers = [2704]
 runNumbers = [2952]
else:
  ### How to actually do it - USAGE: python makeSwizzlerPlots.py 100 101 102... etc (100..102 are run numbers)
  runNumbers = []
  runNumbers = [int(i) for i in sys.argv if i != sys.argv[0]]

for runNumber in runNumbers:  
  inFile = ROOT.TFile("daq_hist_"+str(runNumber)+".root")
  
  tree = inFile.Get("Debug/tMyTree")

  print "run ", runNumber, ", file = ", inFile
  
  save = False
  
  canv = ROOT.TCanvas()
  
  
  lastFrame = -1
  lastSample = -1
  lastTime = -1
  lastTimeBNB = -1
  lastTimeNuMI = -1
  lastTimeEXT = -1
  lastTimePMTBeam = -1
  lastTimePMTCosmic = -1
  
  nBNB = 0
  nNuMI = 0
  nEXT = 0
  nPMTBeam = 0
  nPMTCosmic = 0
  
  hDeltaTBetweenTriggers = ROOT.TH1D("Time between triggers","Time between triggers",300,1e-9,2)
  hDeltaTBetweenTriggersBNB = ROOT.TH1D("Time between triggers (BNB)","Time between triggers (BNB)",100,1e-9,3)
  hDeltaTBetweenTriggersNuMI = ROOT.TH1D("Time between triggers (NuMI)","Time between triggers (NuMI)",100,1e-9,10)
  hDeltaTBetweenTriggersEXT = ROOT.TH1D("Time between triggers (EXT)","Time between triggers (EXT)",100,1e-9,5)
  hDeltaTBetweenTriggersPMTBeam = ROOT.TH1D("Time between triggers (PMT Beam)","Time between triggers (PMT Beam)",100,1e-9,100)
  hDeltaTBetweenTriggersPMTCosmic = ROOT.TH1D("Time between triggers (PMT Cosmic)","Time between triggers (PMT Cosmic)",100,1e-9,5)
  hDeltaTBetweenTriggers.SetTitle("time between triggers;time / s; events")
  hDeltaTBetweenTriggersBNB.SetTitle("(BNB) time between triggers;time / s; events")
  hDeltaTBetweenTriggersNuMI.SetTitle("(NuMI) time between triggers;time / s; events")
  hDeltaTBetweenTriggersEXT.SetTitle("(EXT) time between triggers;time / s; events")
  hDeltaTBetweenTriggersPMTBeam.SetTitle("(PMT beam) time between triggers;time / s; events")
  hDeltaTBetweenTriggersPMTCosmic.SetTitle("(PMT cosmic) time between triggers;time / s; events")
  
  hDeltaTbetweenBeamAndActive = ROOT.TH1D("Time from beam to first active", "Time from beam to first active",100,1e-9,30)
  hDeltaTbetweenBeamAndActive.SetTitle("Time from beam to first active;time / s;")

  
  for event in tree:
    
    #print lastFrame
    if lastFrame < 0:
      firstEventTime = event.triggerTime/1e6
    if event.triggerBitBNB and lastTimeBNB < 0:
      nBNB+=1
      firstEventTimeBNB = event.triggerTime/1e6
    if event.triggerBitNuMI and lastTimeNuMI < 0:
      nNuMI+=1
      firstEventTimeNuMI = event.triggerTime/1e6
    if event.triggerBitEXT and lastTimeEXT < 0:
      nEXT+=1
      firstEventTimeEXT = event.triggerTime/1e6
    if event.triggerBitPMTBeam and lastTimePMTBeam < 0:
      nPMTBeam+=1
      firstEventTimePMTBeam = event.triggerTime/1e6
    if event.triggerBitPMTCosmic and lastTimePMTCosmic < 0:
      nPMTCosmic+=1
      firstEventTimePMTCosmic = event.triggerTime/1e6
    if lastFrame >0:
      frameDiff = event.triggerFrame - lastFrame
      sampleDiff = event.triggerSample - lastSample
      timeDiff = event.triggerTime/1e6 - lastTime
      hDeltaTBetweenTriggers.Fill(timeDiff)
      if event.triggerBitBNB and lastTimeBNB > 0:    
        nBNB+=1
        timeDiffBNB = event.triggerTime/1e6 - lastTimeBNB
        hDeltaTBetweenTriggersBNB.Fill(timeDiffBNB)
      if event.triggerBitNuMI and lastTimeNuMI > 0:
        nNuMI+=1
        timeDiffNuMI = event.triggerTime/1e6 - lastTimeNuMI
        hDeltaTBetweenTriggersNuMI.Fill(timeDiffNuMI)
      if event.triggerBitEXT and lastTimeEXT > 0:    
        nEXT+=1
        timeDiffEXT = event.triggerTime/1e6 - lastTimeEXT
        hDeltaTBetweenTriggersEXT.Fill(timeDiffEXT)
      if event.triggerBitPMTBeam and lastTimePMTBeam > 0:    
        nPMTBeam+=1
        timeDiffPMTBeam = event.triggerTime/1e6 - lastTimePMTBeam
        hDeltaTBetweenTriggersPMTBeam.Fill(timeDiffPMTBeam)
      if event.triggerBitPMTCosmic and lastTimePMTCosmic > 0:    
        nPMTCosmic+=1
        timeDiffPMTCosmic = event.triggerTime/1e6 - lastTimePMTCosmic
        hDeltaTBetweenTriggersPMTCosmic.Fill(timeDiffPMTCosmic)
      if event.triggerActive:
        hDeltaTbetweenBeamAndActive.Fill(event.triggerTime - max(lastTimeBNB,lastTimeNuMI))
  #    print frameDiff, sampleDiff, timeDiff
    lastFrame = event.triggerFrame
    lastSample = event.triggerSample
    lastTime = event.triggerTime/1e6
    if event.triggerBitBNB:    
      lastTimeBNB = event.triggerTime/1e6
    if event.triggerBitNuMI:    
      lastTimeNuMI = event.triggerTime/1e6
    if event.triggerBitEXT:    
      lastTimeEXT = event.triggerTime/1e6
    if event.triggerBitPMTBeam:    
      lastTimePMTBeam = event.triggerTime/1e6
    if event.triggerBitPMTCosmic:    
      lastTimePMTCosmic = event.triggerTime/1e6
  
  averageTrigRate = tree.GetEntries() / (lastTime - firstEventTime) 
  if nBNB:
    averageTrigRateBNB = nBNB / (lastTimeBNB - firstEventTimeBNB) 
  else:
    averageTrigRateBNB = 0
  if nNuMI:
    averageTrigRateNuMI = nNuMI / (lastTimeNuMI - firstEventTimeNuMI) 
  else:
    averageTrigRateNuMI = 0
  if nEXT:
    averageTrigRateEXT = nEXT / (lastTimeEXT - firstEventTimeEXT)
  else:
    averageTrigRateEXT = 0
  if nPMTBeam:
    averageTrigRatePMTBeam = nPMTBeam / (lastTimePMTBeam - firstEventTimePMTBeam)
  else:
    averageTrigRatePMTBeam = 0
  if nPMTCosmic and (lastTimePMTCosmic - firstEventTimePMTCosmic):
    averageTrigRatePMTCosmic = nPMTCosmic / (lastTimePMTCosmic - firstEventTimePMTCosmic)
  else:
    averageTrigRatePMTCosmic = 0
  print "rates: total, BNB, NuMI, EXT, PMTBeam, PMTCosmic"
  print averageTrigRate, averageTrigRateBNB, averageTrigRateNuMI, averageTrigRateEXT, averageTrigRatePMTBeam, averageTrigRatePMTCosmic
  print "diff = ", averageTrigRate - averageTrigRateBNB - averageTrigRateNuMI - averageTrigRateEXT - averageTrigRatePMTBeam - averageTrigRatePMTCosmic
  hTriggerRates = ROOT.TH1D("trigger rates","trigger rates",6,0,6)
  hTriggerRates.GetXaxis().SetBinLabel(1, "BNB")
  hTriggerRates.GetXaxis().SetBinLabel(2, "NuMI")
  hTriggerRates.GetXaxis().SetBinLabel(3, "EXT")
  hTriggerRates.GetXaxis().SetBinLabel(4, "PMT Beam")
  hTriggerRates.GetXaxis().SetBinLabel(5, "PMT Cosmic")
  hTriggerRates.GetXaxis().SetBinLabel(6, "TOTAL")
  hTriggerRates.SetBinContent(1, averageTrigRateBNB)
  hTriggerRates.SetBinContent(2, averageTrigRateNuMI)
  hTriggerRates.SetBinContent(3, averageTrigRateEXT)
  hTriggerRates.SetBinContent(4, averageTrigRatePMTBeam)
  hTriggerRates.SetBinContent(5, averageTrigRatePMTCosmic)
  hTriggerRates.SetBinContent(6, averageTrigRate)
  hTriggerRates.SetTitle("trigger rates;trigger type;rate / Hz")
  hTriggerRates.Draw()
  if save:
    canv.SaveAs("plots/triggerRates"+str(runNumber)+".eps")
  else:
    raw_input()
  canv.SetLogy(True)
  canv.SetLogx(True)
  hDeltaTBetweenTriggers.Draw()
  if save:
    canv.SaveAs("plots/dTbetweenTriggersAll"+str(runNumber)+".eps")
  else:
    raw_input()
  canv.SetLogx(False)
  hDeltaTBetweenTriggersBNB.Draw()
  if save:
    canv.SaveAs("plots/dTbetweenTriggersBNB"+str(runNumber)+".eps")
  else:
    raw_input()
  hDeltaTBetweenTriggersNuMI.Draw()
  if save:
    canv.SaveAs("plots/dTbetweenTriggersNuMI"+str(runNumber)+".eps")
  else:
    raw_input()
  hDeltaTBetweenTriggersEXT.Draw()
  if save:
    canv.SaveAs("plots/dTbetweenTriggersEXT"+str(runNumber)+".eps")
  else:
    raw_input()
  hDeltaTBetweenTriggersPMTBeam.Draw()
  if save:
    canv.SaveAs("plots/dTbetweenTriggersPMTBeam"+str(runNumber)+".eps")
  else:
    raw_input()
  hDeltaTBetweenTriggersPMTCosmic.Draw()
  if save:
    canv.SaveAs("plots/dTbetweenTriggersPMTCosmic"+str(runNumber)+".eps")
  else:
    raw_input()
  canv.SetLogy(False)
  
  hDeltaTbetweenBeamAndActive.Draw()
  if save:
    canv.SaveAs("plots/dTbetweenActiveAndLastBeam"+str(runNumber)+".eps")
  else:
    raw_input()
  
  #raw_input()
  #lastFrame = -1
  #lastSample = -1
  #print "trigger board"
  #for event in tree:
  #  if lastFrame >0:
  #    frameDiff = event.triggerFrame - lastFrame
  #    sampleDiff = event.triggerSample - lastSample
  #    print frameDiff, sampleDiff
  #  lastFrame = event.triggerFrame
  #  lastSample = event.triggerSample
  ###################################
  ## trigger board - FEM differences
  ###################################
  hFrameDiff = ROOT.TH1D("frame difference", "frame difference", 5, -2, 3)
  hSampleDiff = ROOT.TH1D("sample difference", "sample difference", 10, -5, 5)
  hFrameDiffBNB = ROOT.TH1D("frame difference BNB", "frame difference BNB", 5, -2, 3)
  hSampleDiffBNB = ROOT.TH1D("sample difference BNB", "sample difference BNB", 10, -5, 5)
  hFrameDiffNuMI = ROOT.TH1D("frame difference NuMI", "frame difference NuMI", 5, -2, 3)
  hSampleDiffNuMI = ROOT.TH1D("sample difference NuMI", "sample difference NuMI", 10, -5, 5)
  hFrameDiffEXT = ROOT.TH1D("frame difference EXT", "frame difference EXT", 5, -2, 3)
  hSampleDiffEXT = ROOT.TH1D("sample difference EXT", "sample difference EXT", 10, -5, 5)


  hFrameDiffbetweenROAndBNBTrig = ROOT.TH1D("frame diff from RO to trigger BNB","frame diff from RO to trigger BNB",10,-5,5)
  hFrameDiffbetweenROAndNuMITrig = ROOT.TH1D("frame diff from RO to trigger NuMI","frame diff from RO to trigger NuMI",10,-5,5)
  hFrameDiffbetweenROAndEXTTrig = ROOT.TH1D("frame diff from RO to trigger EXT","frame diff from RO to trigger EXT",10,-5,5)
  hSampleDiffbetweenROAndBNBTrig = ROOT.TH1D("sample diff from RO to trigger BNB","sample diff from RO to trigger BNB",40,-20,20)
  hSampleDiffbetweenROAndNuMITrig = ROOT.TH1D("sample diff from RO to trigger NuMI","sample diff from RO to trigger NuMI",40,-20,20)
  hSampleDiffbetweenROAndEXTTrig = ROOT.TH1D("sample diff from RO to trigger EXT","sample diff from RO to trigger EXT",40,-20,20)

  hSampleDiffbetweenRWMAndBNBTrig = ROOT.TH1D("sample diff from trigger BNB to RWM","sample diff from trigger BNB to RWM",1600,-96000,96000)
  hFrameDiffbetweenRWMAndBNBTrig = ROOT.TH1D("frame diff from trigger BNB to RWM","frame diff from trigger BNB to RWM",40,-20,20)

  hTimeDiffTriggerToWaveForm = ROOT.TH1D("waveform time - trigger time","waveform time - trigger time",64,-3200,4800)
  

  for event in tree:
    triggerFrame = event.triggerFrame
    triggerSample = event.triggerSample / 32
    FEM5triggerFrame = event.FEM5triggerFrame
    FEM5triggerSample = event.FEM5triggerSample
    FEM6triggerFrame = event.FEM6triggerFrame
    FEM6triggerSample = event.FEM6triggerSample

    for i in xrange(event.N_PMT_waveforms - 1):
      if event.PMT_waveform_times[i] > 1e-200 and event.PMT_waveform_times[i] < 1e50:
      #print event.PMT_waveform_times[i], event.triggerTime
        hTimeDiffTriggerToWaveForm.Fill( event.PMT_waveform_times[i] - event.triggerTime )
        if event.PMT_waveform_times[i] - event.triggerTime > 4800:
          print "Overflow!", event.PMT_waveform_times[i] - event.triggerTime, event.PMT_waveform_times[i]
        if event.PMT_waveform_times[i] - event.triggerTime < -3200:
          print "Underflow!", event.PMT_waveform_times[i] - event.triggerTime, event.PMT_waveform_times[i]

#    if event.triggerBitBNB:
#      print "trigger time - RO time", event.triggerTime - event.RO_BNBtriggerTime
#      print event.triggerTime, event.RO_BNBtriggerTime
  
    hFrameDiff.Fill(FEM5triggerFrame - triggerFrame)
    hSampleDiff.Fill(FEM5triggerSample - triggerSample)
    if event.triggerBitBNB:
      hFrameDiffBNB.Fill(FEM5triggerFrame - triggerFrame)
      hSampleDiffBNB.Fill(FEM5triggerSample - triggerSample)
      hFrameDiffbetweenROAndBNBTrig.Fill(triggerFrame - event.RO_BNBtriggerFrame)
      hSampleDiffbetweenROAndBNBTrig.Fill(event.triggerSample - event.RO_BNBtriggerSample)
      if event.RO_RWMtriggerFrame > 0:
        hFrameDiffbetweenRWMAndBNBTrig.Fill(event.RO_RWMtriggerFrame - event.triggerFrame) ## RWM signal
        hSampleDiffbetweenRWMAndBNBTrig.Fill(event.RO_RWMtriggerSample - event.triggerSample) ## RWM signal
    if event.triggerBitNuMI:
      hFrameDiffNuMI.Fill(FEM5triggerFrame - triggerFrame)
      hSampleDiffNuMI.Fill(FEM5triggerSample - triggerSample )
      hFrameDiffbetweenROAndNuMITrig.Fill(triggerFrame - event.RO_NuMItriggerFrame)
      hSampleDiffbetweenROAndNuMITrig.Fill(event.triggerSample - event.RO_NuMItriggerSample)
    if event.triggerBitEXT:
      hFrameDiffEXT.Fill(FEM5triggerFrame - triggerFrame)
      hSampleDiffEXT.Fill(FEM5triggerSample - triggerSample)
      hFrameDiffbetweenROAndEXTTrig.Fill(triggerFrame - event.RO_EXTtriggerFrame)
      hSampleDiffbetweenROAndEXTTrig.Fill(event.triggerSample - event.RO_EXTtriggerSample)
#      if event.RO_RWMtriggerFrame > 0:  ### TEMP, FOR THE EXT RWM check
#        hFrameDiffbetweenRWMAndBNBTrig.Fill(event.RO_RWMtriggerFrame - event.triggerFrame) ## RWM signal
#        hSampleDiffbetweenRWMAndBNBTrig.Fill(event.RO_RWMtriggerSample - event.triggerSample) ## RWM signal
#        if (abs(event.RO_RWMtriggerSample - event.triggerSample) > 780):
  
  #  hFrameDiff.Fill(FEM6triggerFrame - triggerFrame)
  #  hSampleDiff.Fill(FEM6triggerSample - triggerSample)
    referenceFrame = FEM5triggerFrame
    referenceSample = FEM5triggerSample
    if FEM5triggerFrame != referenceFrame:
      print "Frames don't agree! (PMT FEM5)", FEM5triggerFrame, referenceFrame
    if FEM6triggerFrame != referenceFrame:
      print "Frames don't agree! (PMT FEM6)", FEM6triggerFrame, referenceFrame
    if event.TPC1triggerFrame != referenceFrame:
      print "Frames don't agree! (TPC1)", event.TPC1triggerFrame, referenceFrame
    if event.TPC2triggerFrame != referenceFrame:
      print "Frames don't agree! (TPC2)", event.TPC2triggerFrame, referenceFrame
    if event.TPC3triggerFrame != referenceFrame:
      print "Frames don't agree! (TPC3)", event.TPC3triggerFrame, referenceFrame
    if event.TPC4triggerFrame != referenceFrame:
      print "Frames don't agree! (TPC4)", event.TPC4triggerFrame, referenceFrame
    if event.TPC5triggerFrame != referenceFrame:
      print "Frames don't agree! (TPC5)", event.TPC5triggerFrame, referenceFrame
    if event.TPC6triggerFrame != referenceFrame:
      print "Frames don't agree! (TPC6)", event.TPC6triggerFrame, referenceFrame
    if event.TPC7triggerFrame != referenceFrame:
      print "Frames don't agree! (TPC7)", event.TPC7triggerFrame, referenceFrame
    if event.TPC8triggerFrame != referenceFrame:
      print "Frames don't agree! (TPC8)", event.TPC8triggerFrame, referenceFrame
    if event.TPC9triggerFrame != referenceFrame:
      print "Frames don't agree! (TPC9)", event.TPC9triggerFrame, referenceFrame
    if FEM5triggerSample != referenceSample:
      print "Samples don't agree! (PMT FEM5)", FEM5triggerSample, referenceSample
    if FEM6triggerSample != referenceSample:
      print "Samples don't agree! (PMT FEM6)", FEM6triggerSample, referenceSample
    if event.TPC1triggerSample != referenceSample:
      print "Samples don't agree! (TPC1)", event.TPC1triggerSample, referenceSample
    if event.TPC2triggerSample != referenceSample:
      print "Samples don't agree! (TPC2)", event.TPC2triggerSample, referenceSample
    if event.TPC3triggerSample != referenceSample:
      print "Samples don't agree! (TPC3)", event.TPC3triggerSample, referenceSample
    if event.TPC4triggerSample != referenceSample:
      print "Samples don't agree! (TPC4)", event.TPC4triggerSample, referenceSample
    if event.TPC5triggerSample != referenceSample:
      print "Samples don't agree! (TPC5)", event.TPC5triggerSample, referenceSample
    if event.TPC6triggerSample != referenceSample:
      print "Samples don't agree! (TPC6)", event.TPC6triggerSample, referenceSample
    if event.TPC7triggerSample != referenceSample:
      print "Samples don't agree! (TPC7)", event.TPC7triggerSample, referenceSample
    if event.TPC8triggerSample != referenceSample:
      print "Samples don't agree! (TPC8)", event.TPC8triggerSample, referenceSample
    if event.TPC9triggerSample != referenceSample:
      print "Samples don't agree! (TPC9)", event.TPC9triggerSample, referenceSample


#    if FEM5triggerFrame != FEM6triggerFrame:
#      print "FEM trigger frames don't agree!", FEM5triggerFrame, FEM6triggerFrame
#      print "triggerSample = ", triggerSample
#    if FEM5triggerSample != FEM6triggerSample:
#      print "FEM trigger samples don't agree!", FEM5triggerSample, FEM6triggerSample
#      print "triggerSample = ", triggerSample

  hTimeDiffTriggerToWaveForm.Draw()
  if save:
    canv.SaveAs("plots/PMT_flash_dtFromTrigger"+str(runNumber)+".eps")
  else:
    raw_input()
  
  hFrameDiffbetweenROAndBNBTrig.Draw()
  if save:
    canv.SaveAs("plots/ROtoTrigFrameDiffBNB"+str(runNumber)+".eps")
  else:
    raw_input()
  hSampleDiffbetweenROAndBNBTrig.Draw()
  if save:
    canv.SaveAs("plots/ROtoTrigSampleDiffBNB"+str(runNumber)+".eps")
  else:
    raw_input()
  hFrameDiffbetweenROAndNuMITrig.Draw()
  if save:
    canv.SaveAs("plots/ROtoTrigFrameDiffNuMI"+str(runNumber)+".eps")
  else:
    raw_input()
  hSampleDiffbetweenROAndNuMITrig.Draw()
  if save:
    canv.SaveAs("plots/ROtoTrigSampleDiffNuMI"+str(runNumber)+".eps")
  else:
    raw_input()

  hFrameDiffbetweenROAndEXTTrig.Draw()
  if save:
    canv.SaveAs("plots/ROtoTrigFrameDiffEXT"+str(runNumber)+".eps")
  else:
    raw_input()
  hSampleDiffbetweenROAndEXTTrig.Draw()
  if save:
    canv.SaveAs("plots/ROtoTrigSampleDiffEXT"+str(runNumber)+".eps")
  else:
    raw_input()
  
  hFrameDiffbetweenRWMAndBNBTrig.Draw()
  if save:
    canv.SaveAs("plots/BNBTrigToRWMFrameDiff"+str(runNumber)+".eps")
  else:
    raw_input()
  hSampleDiffbetweenRWMAndBNBTrig.Draw()
  if save:
    canv.SaveAs("plots/BNBTrigToRWMSampleDiff"+str(runNumber)+".eps")
  else:
    raw_input()


 
  hFrameDiff.Draw()
  if save:
    canv.SaveAs("plots/TriggerFrameDifference"+str(runNumber)+".eps")
  else:
    raw_input()
  hSampleDiff.Draw()
  if save:
    canv.SaveAs("plots/TriggerSampleDifference"+str(runNumber)+".eps")
  else:
    raw_input()
  
  hFrameDiffBNB.Draw()
  if save:
    canv.SaveAs("plots/TriggerFrameDifferenceBNB"+str(runNumber)+".eps")
  else:
    raw_input()
  hSampleDiffBNB.Draw()
  if save:
    canv.SaveAs("plots/TriggerSampleDifferenceBNB"+str(runNumber)+".eps")
  else:
    raw_input()
  
  hFrameDiffNuMI.Draw()
  if save:
    canv.SaveAs("plots/TriggerFrameDifferenceNuMI"+str(runNumber)+".eps")
  else:
    raw_input()
  hSampleDiffNuMI.Draw()
  if save:
    canv.SaveAs("plots/TriggerSampleDifferenceNuMI"+str(runNumber)+".eps")
  else:
    raw_input()
  
  hFrameDiffEXT.Draw()
  if save:
    canv.SaveAs("plots/TriggerFrameDifferenceEXT"+str(runNumber)+".eps")
  else:
    raw_input()
  hSampleDiffEXT.Draw()
  if save:
    canv.SaveAs("plots/TriggerSampleDifferenceEXT"+str(runNumber)+".eps")
  else:
    raw_input()
  
  ################################
  ## compression
  ################################
  hCrateByCrateCompression = ROOT.TH1D("crate by crate compression", "crate by crate compression",10,0,10)
  hMeanCompression = ROOT.TH1D("compression", "compression",100,0,10)
  hEventByEventCompression = ROOT.TH1D("event by event compression", "event by event compression",tree.GetEntries(),0,tree.GetEntries())
  hCrateByCrateCompression.SetTitle(";crate number;compression factor")
  hMeanCompression.SetTitle(";compression factor;number of events")
  hEventByEventCompression.SetTitle(";event number;compression factor")
  ADCwords0 = 0
  ADCwords1 = 0
  ADCwords2 = 0
  ADCwords3 = 0
  ADCwords4 = 0
  ADCwords5 = 0
  ADCwords6 = 0
  ADCwords7 = 0
  ADCwords8 = 0
  ADCwords9 = 0
  NumWords0 = 0
  NumWords1 = 0
  NumWords2 = 0
  NumWords3 = 0
  NumWords4 = 0
  NumWords5 = 0
  NumWords6 = 0
  NumWords7 = 0
  NumWords8 = 0
  NumWords9 = 0
  i = 0
  for event in tree:
    i += 1
    ADCwordsEvent = event.ADCwords_crate0
    ADCwordsEvent += event.ADCwords_crate1
    ADCwordsEvent += event.ADCwords_crate2
    ADCwordsEvent += event.ADCwords_crate3
    ADCwordsEvent += event.ADCwords_crate4
    ADCwordsEvent += event.ADCwords_crate5
    ADCwordsEvent += event.ADCwords_crate6
    ADCwordsEvent += event.ADCwords_crate7
    ADCwordsEvent += event.ADCwords_crate8
    ADCwordsEvent += event.ADCwords_crate9
    NumWordsEvent = event.NumWords_crate0
    NumWordsEvent += event.NumWords_crate1
    NumWordsEvent += event.NumWords_crate2
    NumWordsEvent += event.NumWords_crate3
    NumWordsEvent += event.NumWords_crate4
    NumWordsEvent += event.NumWords_crate5
    NumWordsEvent += event.NumWords_crate6
    NumWordsEvent += event.NumWords_crate7
    NumWordsEvent += event.NumWords_crate8
    NumWordsEvent += event.NumWords_crate9
    if NumWordsEvent:
      hMeanCompression.Fill(ADCwordsEvent / float(NumWordsEvent))
      hEventByEventCompression.SetBinContent(i, ADCwordsEvent / float(NumWordsEvent))
  
    ADCwords0 += event.ADCwords_crate0
    ADCwords1 += event.ADCwords_crate1
    ADCwords2 += event.ADCwords_crate2
    ADCwords3 += event.ADCwords_crate3
    ADCwords4 += event.ADCwords_crate4
    ADCwords5 += event.ADCwords_crate5
    ADCwords6 += event.ADCwords_crate6
    ADCwords7 += event.ADCwords_crate7
    ADCwords8 += event.ADCwords_crate8
    ADCwords9 += event.ADCwords_crate9
    NumWords0 += event.NumWords_crate0
    NumWords1 += event.NumWords_crate1
    NumWords2 += event.NumWords_crate2
    NumWords3 += event.NumWords_crate3
    NumWords4 += event.NumWords_crate4
    NumWords5 += event.NumWords_crate5
    NumWords6 += event.NumWords_crate6
    NumWords7 += event.NumWords_crate7
    NumWords8 += event.NumWords_crate8
    NumWords9 += event.NumWords_crate9
  
  if NumWords1:
    hCrateByCrateCompression.SetBinContent(1,ADCwords1/float(NumWords1))
  if NumWords2:
    hCrateByCrateCompression.SetBinContent(2,ADCwords2/float(NumWords2))
  if NumWords3:
    hCrateByCrateCompression.SetBinContent(3,ADCwords3/float(NumWords3))
  if NumWords4:
    hCrateByCrateCompression.SetBinContent(4,ADCwords4/float(NumWords4))
  if NumWords5:
    hCrateByCrateCompression.SetBinContent(5,ADCwords5/float(NumWords5))
  if NumWords6:
    hCrateByCrateCompression.SetBinContent(6,ADCwords6/float(NumWords6))
  if NumWords7:
    hCrateByCrateCompression.SetBinContent(7,ADCwords7/float(NumWords7))
  if NumWords8:
    hCrateByCrateCompression.SetBinContent(8,ADCwords8/float(NumWords8))
  if NumWords9:
    hCrateByCrateCompression.SetBinContent(9,ADCwords9/float(NumWords9))
  
  hMeanCompression.Draw()
  if save:
    if NumWordsEvent:
      canv.SaveAs("plots/compression"+str(runNumber)+".eps")
  else:
    raw_input()
  
  hEventByEventCompression.SetMinimum(0)
  hEventByEventCompression.Draw()
  if save:
    if NumWordsEvent:
      canv.SaveAs("plots/compressionByEvent"+str(runNumber)+".eps")
  else:
    raw_input()
  hCrateByCrateCompression.Draw()
  if save:
    if NumWordsEvent:
      canv.SaveAs("plots/compressionByCrate"+str(runNumber)+".eps")
  else:
    raw_input()
  
#  continue # instead of exit...
  
  ################################
  ## simple N_discriminators plot
  ################################
  
  N = dict()
  for event in tree:
    for channel in xrange(40):
      if channel not in N:
        N[channel]=0
      N[channel] += event.N_discriminators[channel]
  
  h = ROOT.TH1D("N_discriminators","N_discriminators",40,0,40)
  for bin in xrange(40):
    h.SetBinContent(bin+1, N[bin])
  
  h.SetTitle("Number of discriminators fired in run;channel number;number of discriminators")
  
  h.Draw("hist")
  if save:
    canv.SaveAs("plots/N_discriminators.eps")
  else:
    raw_input()
  
  
  ################################
  ## get trigger + discriminators plot
  ################################
  h = dict()
  for channel in xrange(40):
    h[channel] = ROOT.TH1D("discriminators"+str(channel),"discriminators"+str(channel),100,-3200*2,3200*4)
  for event in tree:
    tf = event.triggerFrame
    ts = event.triggerSample
    for i in xrange(100):
      for channelNumber in xrange(40):
        f = event.discriminatorFrame[i + (channelNumber)*100]
        s = event.discriminatorSample[i + (channelNumber)*100]
    
        if f > 0 and s > 0:
          entry = (f - tf)*3200 + s - ts
          if event.discriminatorType[i + (channelNumber)*100] == 1:
  #          print "discriminator 2", channelNumber, "frame = ", event.discriminatorFrame[i + (channelNumber)*100], ", sample =",event.discriminatorSample[i + (channelNumber)*100]
            h[channelNumber].Fill(entry)
  
  for channel in xrange(40):  
    #print "channel", channel
    h[channel].SetTitle("channel "+str(channel)+";2MHz sample since trigger;number of discriminators")
    h[channel].Draw("hist")
    trigLine = ROOT.TArrow(0,h[channel].GetMaximum(),0,0)
    trigLine.SetLineWidth(2)
    trigLine.Draw()
    if save:
      canv.SaveAs("plots/sampleSinceTrigger_channel"+str(channel)+".eps")
    else:
      raw_input()
  
  
  
