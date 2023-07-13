R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("gallery/Event.h")
#include "tools.h"

// -----------------
// Particle PDGcode:
// https://pdg.lbl.gov/2019/reviews/rpp2018-rev-monte-carlo-numbering.pdf
// -----------------
// #e– 	11
// νe 	12
// μ– 	13
// νμ 	14
// τ– 	15
// ντ 	16
// τ'– 	17
// ντ' 	18 
// p+   2212
// n    2112

bool debug = false;

void PlotSimEnergyDeposit(int nfiles=0, int nevent=42001, int pdgcode=11, string simtag = "largeant:LArG4DetectorServicevolTPCActive", string prod="HDVDPRod2_nu_numu2nue_nue2nutau_30deg_NoRawDigits_reco_v09_53_00d02", bool Debug = false){

  debug = Debug;
  SetDebug(debug);

  string filelist = "./list/"+prod+".list";
  vector<string> filenames = ReadFileList(nfiles,filelist);

  art::InputTag sim_tag(simtag);

  // X-axis
  vector<string> xtitle = {"Z (cm)","Y (cm)","Z (cm)"};
  vector<float>  xmin   = {0   ,-650,0   };
  vector<float>  xmax   = {900 , 650,900 };
  // Y-axis limits
  vector<string> ytitle = {"X (cm)","X (cm)","Y (cm)"};
  vector<float>  ymin   = {-325,-325,-650};
  vector<float>  ymax   = { 325, 325, 650};

  // Setup graphs for each particle species (e+/e- mu+/mu- and p)
  vector<string> pdgnames= {"e^{-}"  ,"e^{+}","#mu^{-}","#mu^{+}","p^{+}" };
  vector<int> pdgcodes   = {11       ,-11    ,-13      ,13       ,2212    };
  Color_t colors[5]      = {kAzure+10,kBlue  ,kOrange  ,kRed     ,kGreen+2};

  vector<vector<TGraph*> > gXYZ(pdgcodes.size());
  for(int i=0;i<gXYZ.size();i++){
    gXYZ[i].resize(3);
    for(int j=0;j<gXYZ[i].size();j++){
      gXYZ[i][j] = new TGraph(); 
      gXYZ[i][j]->SetName(Form("gXYZi_%d_%d",pdgcodes[i],j)); 
      gXYZ[i][j]->SetMarkerColor(colors[i]);                  gXYZ[i][j]->SetLineColor(colors[i]);
      gXYZ[i][j]->GetXaxis()->SetTitle(xtitle[j].c_str());    gXYZ[i][j]->GetYaxis()->SetTitle(ytitle[j].c_str());
    } // end of XYZ loop
  } // end of particle loop

  vector<TH2F*> hXYZ(3);
  for(int i=0;i<hXYZ.size();i++) 
    hXYZ[i] = new TH2F(Form("hXYZ_%d",i),Form(";%s;%s",xtitle[i].c_str(),ytitle[i].c_str()),int(xmax[i]-xmin[i]),xmin[i],xmax[i],int(ymax[i]-ymin[i]),ymin[i],ymax[i]);

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    if(ev.eventAuxiliary().event() < nevent) continue;
    if(ev.eventAuxiliary().event() > nevent) break;

    cout << "Event: " << ev.eventAuxiliary().event() << endl;

    // Retrieve list of energy deposits per event
    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag);
    size_t ndepos = depolist->size();
    cout << "  # deposition: " << ndepos << endl;
  
    // loop over list of depos in this event
    for (size_t depo_i = 0; depo_i < ndepos; depo_i++) {
      const sim::SimEnergyDeposit& depo = depolist->at(depo_i);

      // Retrieve position of energy deposit
      geo::Point_t xyz = depo.MidPoint();

      // Check if particle is within the pdgcodes list - otherwise skip it
      vector<int>::iterator it;
      it = std::find(pdgcodes.begin(), pdgcodes.end(), depo.PdgCode());

      if(it != pdgcodes.end()){
        gXYZ[it - pdgcodes.begin()][0]->SetPoint(gXYZ[it - pdgcodes.begin()][0]->GetN(),xyz.Z(),xyz.X());
        gXYZ[it - pdgcodes.begin()][1]->SetPoint(gXYZ[it - pdgcodes.begin()][1]->GetN(),xyz.Y(),xyz.X());
        gXYZ[it - pdgcodes.begin()][2]->SetPoint(gXYZ[it - pdgcodes.begin()][2]->GetN(),xyz.Z(),xyz.Y());
      }

      // skip particles that does not have the correct pdgcode
      if(depo.PdgCode() != pdgcode)          continue;

      hXYZ[0]->Fill(xyz.Z(),xyz.X(),depo.Energy());
      hXYZ[1]->Fill(xyz.Y(),xyz.X(),depo.Energy());
      hXYZ[2]->Fill(xyz.Z(),xyz.Y(),depo.Energy());

    } // end of depo loop
  } // end of event loop

  // Plot energy deposition color maps for the selected pdgcode particle
  TCanvas* cXYZE = new TCanvas("cXYZE","Geant4 spatial energy deposition distributions");
  cXYZE->Divide(2,2);
  for(int i=0;i<hXYZ.size();i++){ cXYZE->cd(i+1); gPad->SetLogz(); hXYZ[i]->SetStats(0); hXYZ[i]->Draw("colz");}

  // Build the legend to display the particle names that belong to each graph
  TLegend* leg = new TLegend(0.2,0.2,0.8,0.8);
  leg->SetLineWidth(0); leg->SetHeader("Particle list:");
  for(int i=0;i<gXYZ.size();i++){ 
    bool legSet = false;
    for(int j=0;j<gXYZ[i].size();j++){ 
      if(gXYZ[i][j]->GetN() > 0 && !legSet){ 
        legSet = true;
        leg->AddEntry(gXYZ[i][j],pdgnames[i].c_str(),"L");
      }
    }
  }

  cXYZE->SaveAs("simenergydeposit/cXYZE.pdf"); 
  cXYZE->SaveAs("simenergydeposit/cXYZE.root"); 

  // Plot position of energy deposition particle-wise with color coding
  TCanvas* cXYZ = new TCanvas("cXYZ","Geant4 spatial energy deposition distributions");
  cXYZ->Divide(2,2);
  for(int j=0;j<gXYZ[0].size();j++){
    cXYZ->cd(j+1);
    // Draw the legend only on the 1st pad
    bool gSet = false;
    for(int i=0;i<gXYZ.size();i++){
      if(!gSet && gXYZ[i][j]->GetN() > 0){    gXYZ[i][j]->Draw("AP"); gSet=true;}
      else if(gSet && gXYZ[i][j]->GetN() > 0) gXYZ[i][j]->Draw("P");
    }
  }
  cXYZ->cd(4);
  leg->Draw();
 
  cXYZ->SaveAs("simenergydeposit/cXYZ.pdf"); 
  cXYZ->SaveAs("simenergydeposit/cXYZ.root"); 
}
