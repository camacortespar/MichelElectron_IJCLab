/*
  G4 criteria selection for Michel electrons
  Particle PDGcode:
  https://pdg.lbl.gov/2019/reviews/rpp2018-rev-monte-carlo-numbering.pdf
*/

R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "tools.h"

//---PlotTopology---//
void PlotTopology(int nfiles, int event_i, int event_f, string simtag1, string simtag2, string listname, bool muon)
{
  //X-axis
  vector<string> xtitle = {"X (cm)", "Z (cm)", "Z (cm)"};
  vector<float>  xmin   = {-300, -50, -50};
  vector<float>  xmax   = {300, 325, 325};
  //Y-axis
  vector<string> ytitle = {"Y (cm)", "Y (cm)", "X (cm)"};
  vector<float>  ymin   = {-300, -300, -305};
  vector<float>  ymax   = {300, 300, 300};
  //Setup graphs for each particle species (e+/e- mu+/mu- and p)
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes    = {11, -11, 13, -13, 2212};
  Color_t colors[6]       = {kBlue+2, kAzure+8, kOrange+7, kRed+2, kGreen+2, kMagenta+2};

  //Space distributions
  vector<vector<TGraph*>> gXYZ(pdgcodes.size());
  for(int i = 0; i < gXYZ.size(); i++){
    gXYZ[i].resize(3);    //xyz
    for(int j = 0; j < gXYZ[i].size(); j++){
      gXYZ[i][j] = new TGraph(); 
      gXYZ[i][j]->SetName(Form("gXYZ_%d_%d", pdgcodes[i], j)); 
      gXYZ[i][j]->SetMarkerColor(colors[i]);
      gXYZ[i][j]->SetLineColor(colors[i]);
      gXYZ[i][j]->GetXaxis()->SetTitle(xtitle[j].c_str());
      gXYZ[i][j]->GetYaxis()->SetTitle(ytitle[j].c_str());
    } //end of XYZ loop
  }   //end of particle loop

  //Energy depositions
  vector<TH2F*> hXYZ(3);
  for(int i = 0; i < hXYZ.size(); i++){ 
    hXYZ[i] = new TH2F(Form("hXYZ_%d", i), Form(";%s;%s", xtitle[i].c_str(), ytitle[i].c_str()), int(xmax[i]-xmin[i]), xmin[i], xmax[i], int(ymax[i]-ymin[i]), ymin[i], ymax[i]);  
  }

  //To distinguish between muon/antimuon generated
  int pdg;
  string generated;
  if(muon) {pdg = pdgcodes[2]; generated = "muon";}
  if(!muon){pdg = pdgcodes[3]; generated = "antimuon";}

  //All tags that you need
  art::InputTag sim_tag1(simtag1);
  art::InputTag sim_tag2(simtag2);  
  //Search your files
  string filelist          = "./list/"+listname+".list";
  vector<string> filenames = ReadFileList(nfiles, filelist);

  //Go for events!
  for(gallery::Event ev(filenames); !ev.atEnd(); ev.next()){
    if(ev.eventAuxiliary().event() < event_i) continue;
    if(ev.eventAuxiliary().event() > event_f) break; 
    if(ev.eventAuxiliary().event() == 817) continue;    //unsolved problem with this event!  
      
    cout << "Event: " << ev.eventAuxiliary().event() << endl;
      
    //Retrieve list of energy deposits per event
    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag1);
    size_t ndepos       = depolist->size();
    cout << "# Depositions: " << ndepos << endl;

    //Depo loop for spatial and energy distribution
    for(size_t depo_i = 0; depo_i < ndepos; depo_i++){
      const sim::SimEnergyDeposit& depo = depolist->at(depo_i);
  	  
	    //Retrieve position of energy deposit
	    geo::Point_t xyz = depo.MidPoint();
      	  
	    //Check if particle is within the pdgcodes list - otherwise skip it
	    vector<int>::iterator it;
	    it = std::find(pdgcodes.begin(), pdgcodes.end(), depo.PdgCode());
      
	    if(it != pdgcodes.end()){
	      gXYZ[it - pdgcodes.begin()][0]->SetPoint(gXYZ[it - pdgcodes.begin()][0]->GetN(), xyz.X(), xyz.Y());
	      gXYZ[it - pdgcodes.begin()][1]->SetPoint(gXYZ[it - pdgcodes.begin()][1]->GetN(), xyz.Z(), xyz.Y());
	      gXYZ[it - pdgcodes.begin()][2]->SetPoint(gXYZ[it - pdgcodes.begin()][2]->GetN(), xyz.Z(), xyz.X());
	    }

      //Muon/antimuon map energy deposit
      if(depo.PdgCode() == pdg){
        hXYZ[0]->Fill(xyz.X(), xyz.Y(), depo.Energy());
	      hXYZ[1]->Fill(xyz.Z(), xyz.Y(), depo.Energy());
	      hXYZ[2]->Fill(xyz.Z(), xyz.X(), depo.Energy());
      }
	  }   //end of depo loop
  }     //end of event loop

  //Plot position of energy deposition particle-wise with color coding
  TCanvas* cXYZ = new TCanvas("cXYZ", "G4 Spatial Distributions");
  cXYZ->Divide(2,2); 
  for(int j = 0; j < gXYZ[0].size(); j++){
      cXYZ->cd(j+1);
      //Draw every pad  
      bool gSet = false;
      for(int i = 0; i < gXYZ.size(); i++){
	        if(!gSet && gXYZ[i][j]->GetN() > 0){gXYZ[i][j]->Draw("AP"); gSet = true;}
	        else if(gSet && gXYZ[i][j]->GetN() > 0){gXYZ[i][j]->Draw("P");}
	    }
  }

  //Plot energy deposition color maps for muon
  TCanvas* cXYZE = new TCanvas("cXYZE", "G4 Muon/AntiMuon Energy Deposition");
  cXYZE->Divide(2,2);
  for(int i = 0; i < hXYZ.size(); i++){
      cXYZE->cd(i+1); gPad->SetLogz(); hXYZ[i]->SetStats(0); hXYZ[i]->Draw("colz");
  }

  //Build the legend to display the particle names that belong to each graph
  TLegend* leg = new TLegend(0.2, 0.2, 0.8, 0.8);
  leg->SetLineWidth(0); leg->SetHeader("Particle list:");
  for(int i = 0; i < gXYZ.size(); i++){ 
      bool legSet = false;
      for(int j = 0; j < gXYZ[i].size(); j++){ 
	      if(gXYZ[i][j]->GetN() > 0 && !legSet){legSet = true; leg->AddEntry(gXYZ[i][j], pdgnames[i].c_str(), "L");}
	    }
  }
  
  //Choose where you put legend
  cXYZ->cd(4);
  leg->Draw();

  //Save it!
  cXYZ->SaveAs(("../images/simmichelanalysis/"+generated+"_cXYZ_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cXYZ->SaveAs(("../images/simmichelanalysis/"+generated+"_cXYZ_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cXYZE->SaveAs(("../images/simmichelanalysis/"+generated+"_cXYZE_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cXYZE->SaveAs(("../images/simmichelanalysis/"+generated+"_cXYZE_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());
} //end of PlotTopology function

//---MuonAnalysis---//
vector<pair<vector<double>, int>> MuonAnalysis(int nfiles, int event_i, int event_f, string simtag1, string simtag2, string listname, bool muon)
{
  //Particle options
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes    = {11, -11, 13, -13, 2212};
  Color_t colors[6]       = {kBlue+2, kAzure+8, kOrange+7, kRed+2, kGreen+2, kMagenta+2};
  
  //To distinguish between muon/antimuon generated
  int pdg;
  string partname;
  string generated;
  if(muon) {pdg = pdgcodes[2]; partname = pdgnames[2]; generated = "muon";}
  if(!muon){pdg = pdgcodes[3]; partname = pdgnames[3]; generated = "antimuon";}
/*
  vector<string> ytitle_muon = {("Deposited "+partname+" Energy (MeV)").c_str(), ("Remaining "+partname+" Energy (MeV)").c_str()};
  vector<TGraph*> gmuon(2);
  for(int i = 0; i < gmuon.size(); i++){
    gmuon[i] = new TGraph();
    gmuon[i]->SetName(Form("gmuon_%d", i));
    gmuon[i]->SetMarkerColor(colors[i+3]);
    gmuon[i]->SetLineColor(colors[i+3]);
    gmuon[i]->SetMarkerStyle(5);
    gmuon[i]->SetMarkerSize(0.5);
    gmuon[i]->GetXaxis()->SetTitle("Track Length (cm)");
    gmuon[i]->GetYaxis()->SetTitle(ytitle_muon[i].c_str());
  }
*/
  //All tags that you need
  art::InputTag sim_tag1(simtag1);
  art::InputTag sim_tag2(simtag2);  
  //Search your files
  string filelist          = "./list/"+listname+".list";
  vector<string> filenames = ReadFileList(nfiles, filelist);

  int muon_pc = 0;    //to account for all points from all events
  int ev_c    = 0;
  vector<pair<vector<double>, int>> muon_Info(event_f-(event_i-1));

  //Go for events!
  for(gallery::Event ev(filenames); !ev.atEnd(); ev.next()){
    //Physical quantities
    vector<double> muPos;       //hit position from start hit
    vector<double> muDepoE;     //energy deposition in each hit
    vector<double> muRemE;      //remaining energy in each hit

    if(ev.eventAuxiliary().event() < event_i) continue;
    if(ev.eventAuxiliary().event() > event_f) break; 

    cout << "Event: " << ev.eventAuxiliary().event() << endl;
      
    //Retrieve list of energy deposits per event
    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag1);
    size_t ndepos       = depolist->size();
    //Retrieve generated information   
    auto const montecarlo        = ev.getValidHandle<vector<simb::MCTruth>>(sim_tag2);
    const simb::MCTruth& gen     = montecarlo->at(0);
    const simb::MCParticle& muon = gen.GetParticle(0);

    double muon_E0 = (1000)*muon.E();   //in MeV
    cout << "Muon/AntiMuon Energy: " << muon_E0 << " MeV" << endl;
    double muon_AllDepoE   = 0;
    int    muon_LastDepoIt = -1;
    geo::Point_t muon_FirstHit;
    geo::Point_t muon_LastHit;
    TVector3 muon_Eje;
    vector<double> muAxis;

    //Depo loop for muon/antimuon information
    for(size_t depo_i = 0; depo_i < ndepos; depo_i++){
      const sim::SimEnergyDeposit& depo = depolist->at(depo_i);
  	  
	    //Retrieve position of energy deposit
	    geo::Point_t xyz = depo.MidPoint();

      if(depo.PdgCode() == pdg){        
        if(depo_i == 0){muon_FirstHit = depo.MidPoint();}

        double muon_FirstPosition  = sqrt(pow(muon_FirstHit.X(), 2)+pow(muon_FirstHit.Y(), 2)+pow(muon_FirstHit.Z(), 2));
        double muon_ActualPosition = sqrt(pow(xyz.X(), 2)+pow(xyz.Y(), 2)+pow(xyz.Z(), 2));
        double muon_TrackLength    = muon_ActualPosition - muon_FirstPosition;
        double muon_RemEnergy      = muon_E0 - depo.Energy();
        muon_AllDepoE             += depo.Energy();
        
        //Energy deposit and position for muon/antimuon
        muPos.push_back(abs(muon_TrackLength));
        muDepoE.push_back(depo.Energy());
        muRemE.push_back(muon_RemEnergy);

        muon_E0         = muon_RemEnergy;      //update actual energy
        muon_LastDepoIt = depo_i;              //save last hit
      } //end of muon/anitmuon condition
	  }   //end of depo loop 
    
    if(muPos.size() == 0){cout << "There is no muon/antimuon information here!" << endl; ev_c++; continue;}    //go back to begining if you don't have muon information
    
    const sim::SimEnergyDeposit& muon_LastDepo = depolist->at(muon_LastDepoIt);
    muon_LastHit = muon_LastDepo.MidPoint();
    cout << "Muon/AntiMuon Total Deposit Energy: " << muon_AllDepoE << " MeV" << endl;
    cout << "Muon/AntiMuon Last Hit: " << muon_LastHit << endl;
    //Direction vector  
    muon_Eje(0) = muon_FirstHit.X()-muon_LastHit.X();
    muon_Eje(1) = muon_FirstHit.Y()-muon_LastHit.Y();
    muon_Eje(2) = muon_FirstHit.Z()-muon_LastHit.Z();
    muon_Eje    = muon_Eje.Unit();

    for(int i = 0; i < 3; i++){muAxis.push_back(muon_Eje(i));}
  /*
    for (int i = 0; i < muPos.size(); i++){   
        gmuon[0]->SetPoint(muon_pc, muPos[i], muDepoE[i]);
        gmuon[1]->SetPoint(muon_pc, muPos[i], muRemE[i]);
        muon_pc++;
    }   //end of muon/antimuon information loop
  */
    //pair<vector<TVector3>, geo::Point_t> p(muon_Axis, muon_LastHit);
    pair<vector<double>, int> p(muAxis, muon_LastDepoIt);
    muon_Info[ev_c] = p;
    ev_c++;
  } //end of event loop
/*
  //Plot muon/antimuon information
  TCanvas* cmuon_1 = new TCanvas("cmuon_1", "Muon/AntiMuon DepoE");
  cmuon_1->cd();
  gmuon[0]->Draw("AP");
  cmuon_1->Update();

  TCanvas* cmuon_2 = new TCanvas("cmuon_2", "Muon/AntiMuon RemE");
  cmuon_2->cd();
  gmuon[1]->Draw("AP");
  cmuon_2->Update();

  //Save it!
  cmuon_1->SaveAs(("../images/simmichelanalysis/"+generated+"_depoE_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cmuon_1->SaveAs(("../images/simmichelanalysis/"+generated+"_depoE_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cmuon_2->SaveAs(("../images/simmichelanalysis/"+generated+"_remE_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cmuon_2->SaveAs(("../images/simmichelanalysis/"+generated+"_remE_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());
*/
  return muon_Info;
}

//---ElectronAnalysis---//
void ElectronAnalysis(int nfiles, int event_i, int event_f, string simtag1, string simtag2, string simtag3, string listname, vector<pair<vector<double>, int>> muon_Info, bool muon)
{
  //X-axis
  vector<string> xtitle = {"#theta (degrees)", "Energy (MeV)", "Energy (MeV)", "Energy (MeV)", "Energy (MeV)"};
  vector<int> bins      = {180, 90, 90, 90, 90};
  vector<float> xmin    = {-10, 0, 0, 0, 0};
  vector<float> xmax    = {200, 360, 100, 60, 60};
  //Y-axis
  vector<string> ytitle_elec = {"Total Deposited e^{-}/e^{+} Energy (MeV)", "Completeness"};
  //Particle options
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes    = {11, -11, 13, -13, 2212};
  Color_t colors[7]       = {kBlue+2, kAzure+8, kOrange+7, kRed+2, kGreen-3, kMagenta+2, kBlue-6};

  //To distinguish between muon/antimuon generated
  int pdg;
  string partname;
  string generated;
  if(muon) {pdg = pdgcodes[2]; partname = pdgnames[2]; generated = "muon";}
  if(!muon){pdg = pdgcodes[3]; partname = pdgnames[3]; generated = "antimuon";}

  vector<TH1D*> helec(5);
  for(int i = 0; i < helec.size(); i++){
    helec[i] = new TH1D(Form("helec_%d", i), Form(";%s;Counts", xtitle[i].c_str()), bins[i], xmin[i], xmax[i]);
    if(muon) {helec[i]->SetLineColor(colors[6]);  helec[i]->SetFillColor(colors[6]);}
    if(!muon){helec[i]->SetLineColor(colors[1]);  helec[i]->SetFillColor(colors[1]);} 
  }
  helec[0]->SetLineColor(colors[5]);  helec[0]->SetFillColor(colors[5]);

  vector<TGraph*> gene_sel(2);
  for(int i = 0; i < gene_sel.size(); i++){
    gene_sel[i] = new TGraph();
    gene_sel[i]->SetName(Form("gene_sel_%d", i));
    if(muon) {gene_sel[i]->SetMarkerColor(colors[6]);  gene_sel[i]->SetLineColor(colors[6]); gene_sel[i]->SetFillColor(colors[6]);}
    if(!muon){gene_sel[i]->SetMarkerColor(colors[1]);  gene_sel[i]->SetLineColor(colors[1]); gene_sel[i]->SetFillColor(colors[1]);} 
    gene_sel[i]->SetMarkerStyle(3);
    gene_sel[i]->SetMarkerSize(0.25);
    gene_sel[i]->GetXaxis()->SetTitle("#theta (degrees)");
    gene_sel[i]->GetYaxis()->SetTitle(ytitle_elec[i].c_str());
  }

  TH1D* hposi = new TH1D("", Form("Decay MC e^{+};%s;Counts", xtitle[4].c_str()), bins[4], xmin[4], xmax[4]);
  hposi->SetLineColor(colors[4]);  hposi->SetFillColor(colors[4]);

  TH1D* hdiff = new TH1D("", Form("E_{e^{+}, MC} - E_{e^{+}, Selected};%s;Counts", xtitle[4].c_str()), bins[4], -xmax[4], xmax[4]);
  hdiff->SetLineColor(colors[3]);  //hdiff->SetFillColor(colors[3]);

  //All tags that you need
  art::InputTag sim_tag1(simtag1);
  art::InputTag sim_tag2(simtag2);
  art::InputTag sim_tag3(simtag3); 
  //Search your files
  string filelist          = "./list/"+listname+".list";
  vector<string> filenames = ReadFileList(nfiles, filelist);

  vector<double> Posi_Energy(event_f-(event_i-1));           //to store each decay positron energy from mc
  const int Nbins = 180;
  vector<float> SelectedEnergy(Nbins);
  vector<float> SelectedCompleteness(Nbins);
  int elec_pc = 0;                                          //to account for all points from all events
  int ev_c    = 0;
  
  //Go for events!
  for(gallery::Event ev(filenames); !ev.atEnd(); ev.next()){
    //Physical quantities
    vector<double> elecDist;    //distance of electron hit from last muon hit
    vector<double> elecDepoE;   //deposited energy of selected electrons

    if(ev.eventAuxiliary().event() < event_i) continue;
    if(ev.eventAuxiliary().event() > event_f) break; 
      
    //Retrieve list of energy deposits per event
    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag1);
    size_t ndepos       = depolist->size();

    //Decay mc positron information ~ only for antimuon events
    if(!muon){
      auto const posicarlo = ev.getValidHandle<vector<simb::MCParticle>>(sim_tag3);
      size_t nposi         = posicarlo->size();

      for(size_t pos_i = 0; pos_i < nposi; pos_i++){
        const simb::MCParticle& positron = posicarlo->at(pos_i);
        if(positron.PdgCode() == -11 && positron.Process() == "Decay"){Posi_Energy[ev_c] = (1000)*positron.E(); hposi->Fill(Posi_Energy[ev_c]);}
      }    
    } 

    if(muon_Info[ev_c].first.size() == 0){cout << "There is no muon information here!" << endl; ev_c++; continue;}    //go back to begining if you don't have muon information

    //Information from muon/anitmuon analysis
    TVector3 muon_Axis;        
    vector<double> eje = muon_Info[ev_c].first;      
    int muon_LastDepoIt = muon_Info[ev_c].second;
    for(int i = 0; i < eje.size(); i++){muon_Axis(i) = eje[i];}
    const sim::SimEnergyDeposit& muon_LastDepo = depolist->at(muon_LastDepoIt);
    geo::Point_t muon_LastHit = muon_LastDepo.MidPoint();

    //Electron information
    double Rmax  = 25.0;        //max. distance from last muon hit (in cm) ~ electron analysis
    double theta;               //cone opening angle
    int theta_bin; 
    double elec_DepoTot   = 0;
    double elec_DepoMask  = 0; 
    double elec_DepoSph   = 0;
    double elec_DepoSelec = 0;       

    TVector3 w_Hit;             //weighted electron hits
    TVector3 elec_Bary;         //coordinates of barycenter shower
    TVector3 r_Axis;            //axis cone vector
    TVector3 r_ElecHit;         //all electron hits vector
    TVector3 r_ElecHit_Sel;     //selected electron hits vector
    
    //Let's create our cone
    for(size_t edepo_i = 0; edepo_i < ndepos; edepo_i++){
      const sim::SimEnergyDeposit& edepo = depolist->at(edepo_i);
  	  
	    //Retrieve position of energy deposit
	    geo::Point_t e_xyz = edepo.MidPoint();
      r_ElecHit(0) = e_xyz.X()-muon_LastHit.X();
      r_ElecHit(1) = e_xyz.Y()-muon_LastHit.Y();
      r_ElecHit(2) = e_xyz.Z()-muon_LastHit.Z();
      //Distance of electron hit from muon last hit
      double r_aux       = sqrt(pow((e_xyz.X()-muon_LastHit.X()), 2)+pow((e_xyz.Y()-muon_LastHit.Y()), 2)+pow((e_xyz.Z()-muon_LastHit.Z()), 2));
      double theta_muon  = muon_Axis.Angle(r_ElecHit)*TMath::RadToDeg(); 
      if(abs(edepo.PdgCode()) == 11){  
        elec_DepoTot += edepo.Energy();  
        if(abs(theta_muon) < 20.0){continue;}      //muon/antimuon track masking 
        elec_DepoMask += edepo.Energy();
        if(r_aux <= Rmax){
          w_Hit(0)     += edepo.Energy()*e_xyz.X();
          w_Hit(1)     += edepo.Energy()*e_xyz.Y();
          w_Hit(2)     += edepo.Energy()*e_xyz.Z();
          elec_DepoSph += edepo.Energy();
        } //end of sphere condition
      }   //enf of electron condition
    }     //end of depo loop

    helec[1]->Fill(elec_DepoTot);
    helec[2]->Fill(elec_DepoMask);
    helec[3]->Fill(elec_DepoSph);

    if(elec_DepoSph <= 0){cout << "There is no Michel electron information here!" << endl; ev_c++; continue;}    //go back to begining if you don't have Michel electron information

    //Fill barycenter shower vector
    for(int i = 0; i < 3; i++){elec_Bary(i) = w_Hit(i)/elec_DepoSph;}
    cout << "Barycenter: (" << elec_Bary.X() << "," << elec_Bary.Y() << "," << elec_Bary.Z() << ")" << endl;

    //Find normalized axis cone vector  
    r_Axis(0) = elec_Bary.X()-muon_LastHit.X();
    r_Axis(1) = elec_Bary.Y()-muon_LastHit.Y();
    r_Axis(2) = elec_Bary.Z()-muon_LastHit.Z();
    r_Axis    = r_Axis.Unit();

    vector<float> elec_DistriE(Nbins);  

    //Depo loop for michel electron selection
    for(size_t edepo_i = 0; edepo_i < ndepos; edepo_i++){
      const sim::SimEnergyDeposit& edepo = depolist->at(edepo_i);
  	  
	    //Retrieve position of energy deposit
	    geo::Point_t e_xyz = edepo.MidPoint();
      
      if(abs(edepo.PdgCode()) == 11){        
        r_ElecHit_Sel(0)   = e_xyz.X()-muon_LastHit.X();
        r_ElecHit_Sel(1)   = e_xyz.Y()-muon_LastHit.Y();
        r_ElecHit_Sel(2)   = e_xyz.Z()-muon_LastHit.Z();
        double r_aux       = r_ElecHit_Sel.Mag();
        r_ElecHit_Sel      = r_ElecHit_Sel.Unit();
        double theta_muon  = muon_Axis.Angle(r_ElecHit_Sel)*TMath::RadToDeg();    //angle between muon direction and hit

        if(abs(theta_muon) < 20.0){continue;}                                     //muon/antimuon track masking condition

        theta     = r_Axis.Angle(r_ElecHit_Sel)*TMath::RadToDeg();                //angle between axis and hit in degrees
        theta_bin = floor(((theta*Nbins)/180));

        if(r_aux <= Rmax){                                                        //containment sphere condition                   
          elec_DistriE[theta_bin] += edepo.Energy();                      
          if(theta <= 70){elec_DepoSelec += edepo.Energy();}                      //selection cone condition       
          //Fill angle distribution for sphere events 
          helec[0]->Fill(theta);
        } //end of containtment sphere condition
      }   //end of electron hit condition
    }     //end of depo loop

    helec[4]->Fill(elec_DepoSelec);
    
    //double sumcomp   = 0.0;
    double sumcomp = 0.0;
    for(int i = 0; i < Nbins; i++){
      sumcomp += elec_DistriE[i]/elec_DepoSph;              //portion of e energy inside cs + tm zone
      SelectedEnergy[i]       += elec_DistriE[i]/elec_DepoSph;      
      SelectedCompleteness[i] += sumcomp;  
    }

    hdiff->Fill(Posi_Energy[ev_c]-elec_DepoSelec);          //energy difference between the mc and the selected one
    ev_c++;   
  } //end of event loop

  cout << "# Events: " << ev_c << endl;
  
  for(int i = 0; i < Nbins; i++){
    SelectedEnergy[i]       /= ev_c;
    SelectedCompleteness[i] /= ev_c;
    gene_sel[0]->SetPoint(i, (i*180)/Nbins, SelectedEnergy[i]);
    gene_sel[1]->SetPoint(i, (i*180)/Nbins, SelectedCompleteness[i]);
  }

  //Get some numbers
  double y_me;

  for(int i = 0; i < gene_sel[1]->GetN(); ++i){
    double x, y;
    gene_sel[1]->GetPoint(i, x, y);
    if(x == 70){y_me = y;}
  }
  cout << "Completeness (theta = 70°): " << y_me << endl;
  
  //Plot electron information
  TCanvas* ctheta = new TCanvas("ctheta", "Selected Theta Distribution");
  ctheta->cd();
  ctheta->SetLogy();
  helec[0]->SetStats(0);
  helec[0]->Draw(); 
  ctheta->Update();

  TCanvas* cene = new TCanvas("cene", "Total Energy Distribution");
  cene->cd();
  cene->SetLogy();  
  helec[1]->SetStats(0);
  helec[1]->Draw(); 
  cene->Update();

  TCanvas* cene_2 = new TCanvas("cene_2", "MTM");
  cene_2->cd();
  cene_2->SetLogy();
  helec[2]->SetStats(0);
  helec[2]->Draw(); 
  cene_2->Update();

  TCanvas* cene_3 = new TCanvas("cene_3", "MTM + CS");
  cene_3->cd();
  cene_3->SetLogy();
  helec[3]->SetStats(0);
  helec[3]->Draw();; 
  cene_3->Update();

  TCanvas* cene_4 = new TCanvas("cene_4", "MTM + CS + SC");
  cene_4->cd();
  cene_4->SetLogy();
  helec[4]->SetStats(0);
  helec[4]->Draw(); 
  cene_4->Update();

  TCanvas* celec_sel = new TCanvas("celec_sel", "Selected Electron Information");
  celec_sel->Divide(1,2);
  for(int i = 0; i < gene_sel.size(); i++){
    celec_sel->cd(i+1); gene_sel[i]->Draw("AB");
  } 
  celec_sel->Update();

  //Decay mc positron information
  if(!muon){
    THStack* stack = new THStack("stack", "");
    stack->Add(hposi);
    stack->Add(helec[4]);

    TCanvas* cposi = new TCanvas("cposi", "Decay MC e^{+}");
    cposi->cd();
    cposi->SetLogy();
    stack->Draw("hist");

    //Legend
    TLegend* legend = new TLegend(0.125, 0.8, 0.225, 0.88);
    legend->AddEntry(hposi, " Decay MC Positron", "l");
    legend->AddEntry(helec[4], " Selected e^{+}/e^{-}", "l");
    legend->SetBorderSize(0);
    legend->Draw();

    TCanvas* cdiff = new TCanvas("cdiff", "Positron Energy Diff");
    cdiff->cd();
    cdiff->SetLogy();
    hdiff->SetStats(0);
    hdiff->Draw(); 
    cdiff->Update();

    cposi->SaveAs(("../images/simmichelanalysis/"+generated+"_posi_energy_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
    cposi->SaveAs(("../images/simmichelanalysis/"+generated+"_posi_energy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

    cdiff->SaveAs(("../images/simmichelanalysis/"+generated+"_diff_posi_energy_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
    cdiff->SaveAs(("../images/simmichelanalysis/"+generated+"_diff_posi_energy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());
  }

  //Save it!
  ctheta->SaveAs(("../images/simmichelanalysis/"+generated+"_theta_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  ctheta->SaveAs(("../images/simmichelanalysis/"+generated+"_theta_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cene->SaveAs(("../images/simmichelanalysis/"+generated+"_Eenergy_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cene->SaveAs(("../images/simmichelanalysis/"+generated+"_Eenergy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cene_2->SaveAs(("../images/simmichelanalysis/"+generated+"_MTMenergy_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cene_2->SaveAs(("../images/simmichelanalysis/"+generated+"_MTMenergy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cene_3->SaveAs(("../images/simmichelanalysis/"+generated+"_CSenergy_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cene_3->SaveAs(("../images/simmichelanalysis/"+generated+"_CSenergy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cene_4->SaveAs(("../images/simmichelanalysis/"+generated+"_Selenergy_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cene_4->SaveAs(("../images/simmichelanalysis/"+generated+"_Selenergy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  celec_sel->SaveAs(("../images/simmichelanalysis/"+generated+"_sel_elec_energy_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  celec_sel->SaveAs(("../images/simmichelanalysis/"+generated+"_sel_elec_energy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());
}

bool debug = false;

// ----- M A I N ----- //

void SimMichelAnalysis
(
  int nfiles      = 1,            //# files in your list
  int event_i     = 1,            //starting event
  int event_f     = 1,            //final event
  string simtag1  = "branch1",
  string simtag2  = "branch2",
  string simtag3  = "branch3", 
  string listname = "name",
  bool muon       = true,         //choose muon (true) or anti muon (false) ~ be careful with the listname!
  bool debug      = false
)
{
  debug = Debug;
  SetDebug(debug);

  //PlotTopology(nfiles, event_i, event_f, simtag1, simtag2, listname, muon);
  vector<pair<vector<double>, int>> muon_Info = MuonAnalysis(nfiles, event_i, event_f, simtag1, simtag2, listname, muon);
  ElectronAnalysis(nfiles, event_i, event_f, simtag1, simtag2, simtag3, listname, muon_Info, muon);
}