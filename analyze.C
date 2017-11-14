
{
    TFile * f = TFile::Open("g4out.root");
    TTree * tree = (TTree*)f->Get("ntuple/ntuple");
 
    gStyle->SetOptStat(1111111);
   
    bool drawDelta = false;

    TCanvas * el_canvas = new TCanvas();
    el_canvas->Divide(2);
    el_canvas->cd(1);
    if(drawDelta) tree->Draw("timingVector-timingVector[0] >> el1(1000,0,1000)","processVector == 6");
    else tree->Draw("timingVector >> el1(1000,0,1000)","processVector == 6");
    gPad->SetLogy();
    el_canvas->cd(2);
    if(drawDelta) tree->Draw("timingVector-timingVector[0] >> el2(100,0,1e22)","processVector == 6");
    else tree->Draw("timingVector >> el2(100)","processVector == 6 && timingVector > 1000");
    
    TCanvas * inel_canvas = new TCanvas();
    inel_canvas->Divide(2);
    inel_canvas->cd(1);
    if(drawDelta) tree->Draw("timingVector-timingVector[0] >> inel1(1000,0,1000)","processVector == 7");
    else tree->Draw("timingVector >> inel1(1000,0,1000)","processVector == 7");
    gPad->SetLogy();
    inel_canvas->cd(2);
    if(drawDelta) tree->Draw("timingVector-timingVector[0] >> inel2(100,0,1e22)","processVector == 7");
    else tree->Draw("timingVector >> inel2(100)","processVector == 7 && timingVector > 1000");

    TCanvas * cap_canvas = new TCanvas();
    cap_canvas->Divide(2);
    cap_canvas->cd(1);
    if(drawDelta) tree->Draw("timingVector-timingVector[0] >> cap1(1000,0,1000)","processVector == 8");
    else tree->Draw("timingVector >> cap1(1000,0,1000)","processVector == 8");
    gPad->SetLogy();
    cap_canvas->cd(2);
    if(drawDelta) tree->Draw("timingVector-timingVector[0] >> cap2(100,0,1e22)","processVector == 8");
    else tree->Draw("timingVector >> cap2(100)","processVector == 8 && timingVector > 1000");


}
