R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "tools.h"

// -----------------
// Particle PDGcode:
// https://pdg.lbl.gov/2019/reviews/rpp2018-rev-monte-carlo-numbering.pdf
// -----------------
// e– 	11
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

//-----Main-----//

void MuonGroove
(
  int nfiles     = 1,
  int event_i    = 1,
  int event_f    = 1,
  string simtag1 = "branch1",
  string simtag2 = "branch2", 
  string name    = "list owner",
  int mpdg       = 1,
  bool Debug     = false
)
{
  debug = Debug;
  SetDebug(debug);

  //All tags that you need
  art::InputTag sim_tag1(simtag1);
  art::InputTag sim_tag2(simtag2);  
  //Search your files
  string filelist          = "./list/camilo_"+name+".list";
  vector<string> filenames = ReadFileList(nfiles, filelist);
  
  //X-axis
  vector<string> xtitle = {"X (cm)", "Z (cm)", "Z (cm)"};
  vector<float>  xmin   = {-325, -350, -350};
  vector<float>  xmax   = {325, 350, 350};
  //Y-axis
  vector<string> ytitle = {"Y (cm)", "Y (cm)", "X (cm)"};
  vector<float>  ymin   = {-650, -650, -325};
  vector<float>  ymax   = {650, 650, 325};
  vector<string> ytitle_muon = {"Deposited #mu^{+} Energy (MeV)", "Remaining #mu^{+} Energy (MeV)"};
  vector<string> ytitle_elec = {"Total Deposited e^{-}/e^{+} Energy (MeV)", "Completeness"};

  //Setup graphs for each particle species (e+/e- mu+/mu- and p)
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes    = {11, -11, 13, -13, 2212};
  Color_t colors[7]       = {kBlue+2, kAzure+8, kOrange-3, kRed+2, kGreen+2, kViolet-5, kBlue-6};
  
  //---Plots---///
  
  //Space distributions
  vector<vector<TGraph*> > gXYZ(pdgcodes.size());
  for(int i = 0; i < gXYZ.size(); i++){
    gXYZ[i].resize(3);  //xyz
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

  //Muon information
  vector<TGraph*> gmuon(2);
  for(int i = 0; i < gmuon.size(); i++){
    gmuon[i] = new TGraph();
    gmuon[i]->SetName(Form("gmuon_%d", i));
    gmuon[i]->SetMarkerColor(colors[i+3]);
    gmuon[i]->SetLineColor(colors[i+3]);
    gmuon[i]->SetMarkerStyle(6);
    gmuon[i]->GetXaxis()->SetTitle("Track Length (cm)");
    gmuon[i]->GetYaxis()->SetTitle(ytitle_muon[i].c_str());
  }

  TH1D* hmuon = new TH1D("", "Simulated Cosmic Muon Energy Distribution; Energy (MeV); Counts", 90, 0, 1400);
  hmuon->SetLineColor(colors[2]);
  hmuon->SetFillColor(colors[2]);

  //Electron information
  TH1D* htheta = new TH1D("htheta", "Selected Angle Distribution; #theta (degrees); Counts", 180, -10, 200);
  htheta->SetLineColor(colors[5]);
  htheta->SetFillColor(colors[5]);

  TH1D* hene = new TH1D("hene", "Energy Distribution for PDGCode = 11, -11; Energy (MeV); Counts", 90, 0, 300);
  hene->SetLineColor(colors[6]);  //6
  hene->SetFillColor(colors[6]);

  TH1D* hene_2 = new TH1D("hene_2", "Muon Track Masking; Energy (MeV); Counts", 90, 0, 100);
  hene_2->SetLineColor(colors[6]);
  hene_2->SetFillColor(colors[6]);

  TH1D* hene_3 = new TH1D("hene_3", "Muon Track Masking + Containment Sphere; Energy (MeV); Counts", 90, 0, 60);
  hene_3->SetLineColor(colors[6]);
  hene_3->SetFillColor(colors[6]);

  TH1D* hene_4 = new TH1D("hene_4", "Muon Track Masking + Containment Sphere + Selection Cone; Energy (MeV); Counts", 90, 0, 60);
  hene_4->SetLineColor(colors[6]);
  hene_4->SetFillColor(colors[6]);

  vector<TGraph*> gene(2);
  for(int i = 0; i < gene.size(); i++){
    gene[i] = new TGraph();
    gene[i]->SetName(Form("gene_%d", i));
    gene[i]->SetMarkerColor(colors[i]);
    gene[i]->SetLineColor(colors[i]);
    gene[i]->SetFillColor(colors[i]); 
    gene[i]->SetMarkerStyle(3);
    gene[i]->SetMarkerSize(0.25);
    gene[i]->GetXaxis()->SetTitle("#theta (degrees)");
    gene[i]->GetYaxis()->SetTitle(ytitle_elec[i].c_str());
  }

  vector<TGraph*> gene_sel(2);
  for(int i = 0; i < gene_sel.size(); i++){
    gene_sel[i] = new TGraph();
    gene_sel[i]->SetName(Form("gene_sel_%d", i));
    gene_sel[i]->SetMarkerColor(colors[6]);
    gene_sel[i]->SetLineColor(colors[6]);
    gene_sel[i]->SetFillColor(colors[6]); 
    gene_sel[i]->SetMarkerStyle(3);
    gene_sel[i]->SetMarkerSize(0.25);
    gene_sel[i]->GetXaxis()->SetTitle("#theta (degrees)");
    gene_sel[i]->GetYaxis()->SetTitle(ytitle_elec[i].c_str());
  }

  int muon_pc = 0, elec_pc = 0;   //to account for all points from all events

  //---Yeah Science!---//

  const int Nbins = 180;
  vector<float> EnergyPerAngleEvent(Nbins);
  vector<float> SelectedEnergyPerAngleEvent(Nbins);
  vector<float> CompletenessOnlyAngleEvent(Nbins);
  vector<float> CompletenessRadiusAndAngleEvent(Nbins);
  int ev_counter = 0;
  
  //Go for events!
  for(gallery::Event ev(filenames); !ev.atEnd(); ev.next()){
    //Physical quantities
    vector<double> muPos;       //muon hit position from start hit
    vector<double> muDepoE;     //muon energy deposition in each hit
    vector<double> muRemE;      //muon remaining energy in each hit
    vector<double> elecDist;    //distance of electron hit from last muon hit
    vector<double> elecDepoE;   //deposited energy of selected electrons

    if(ev.eventAuxiliary().event() < event_i) continue;
    if(ev.eventAuxiliary().event() > event_f) break; 

    if(ev.eventAuxiliary().event() == 817){continue;}  
      
    cout << "Event: " << ev.eventAuxiliary().event() << endl;
      
    //Retrieve list of energy deposits per event
    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag1);
    size_t ndepos       = depolist->size();
    cout << "# Depositions " << ndepos << endl;
    //Retrieve generated information   
    auto const montcarlo         = ev.getValidHandle<vector<simb::MCTruth>>(sim_tag2);
    const simb::MCTruth& gen     = montcarlo->at(0);
    const simb::MCParticle& muon = gen.GetParticle(0);

    double muon_E0 = (1000)*muon.E();   //in MeV
    cout << "Muon Energy: " << muon_E0 << " MeV" << endl;
    double muon_AllDepoE   = 0;
    int    muon_LastDepoIt = -1;
    geo::Point_t muon_FirstHit;
    TVector3 Muon_Axis;
    hmuon->Fill(muon_E0);

    //Depo loop for spatial and energy distribution + muon information
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

      //Muon information
      if(depo.PdgCode() == mpdg){
        //Map energy deposit
        hXYZ[0]->Fill(xyz.X(), xyz.Y(), depo.Energy());
	      hXYZ[1]->Fill(xyz.Z(), xyz.Y(), depo.Energy());
	      hXYZ[2]->Fill(xyz.Z(), xyz.X(), depo.Energy());
        
        if(depo_i == 0){muon_FirstHit = depo.MidPoint();}

        double muon_FirstPosition  = sqrt(pow(muon_FirstHit.X(), 2)+pow(muon_FirstHit.Y(), 2)+pow(muon_FirstHit.Z(), 2));
        double muon_ActualPosition = sqrt(pow(xyz.X(), 2)+pow(xyz.Y(), 2)+pow(xyz.Z(), 2));
        double muon_TrackLength    = muon_ActualPosition - muon_FirstPosition;
        double muon_RemEnergy      = muon_E0 - depo.Energy();
        muon_AllDepoE             += depo.Energy();
        
        //Energy deposit and position for muon
        muPos.push_back(abs(muon_TrackLength));
        muDepoE.push_back(depo.Energy());
        muRemE.push_back(muon_RemEnergy);

        muon_E0         = muon_RemEnergy;    //update the actual muon energy
        muon_LastDepoIt = depo_i;              //save muon last hit
      } //end of muon condition
	  }   //end of depo loop 
    
    if(muPos.size() == 0){cout << "There is no muon information here!" << endl; continue;}    //go back to begining if you don't have muon information

    const sim::SimEnergyDeposit& muon_LastDepo = depolist->at(muon_LastDepoIt);
    geo::Point_t muon_LastHit                  = muon_LastDepo.MidPoint();
    cout << "Muon Total Deposit Energy: " << muon_AllDepoE << " MeV" << endl;
    cout << "Muon Last Hit: " << muon_LastHit << " || " << muon_LastDepoIt << endl;
    //Muon direction vector  
    Muon_Axis(0) = muon_FirstHit.X()-muon_LastHit.X();
    Muon_Axis(1) = muon_FirstHit.Y()-muon_LastHit.Y();
    Muon_Axis(2) = muon_FirstHit.Z()-muon_LastHit.Z();
    Muon_Axis    = Muon_Axis.Unit();
    /*
    //Fill muon information
    for (int i = 0; i < muPos.size(); i++){
        //Get muon position, deposited and remaining energy for each hit 
        double muon_Posicion = muPos[i];
        double muon_Deposito = muDepoE[i];
        double muon_Restante = muRemE[i];       
  
        gmuon[0]->SetPoint(muon_pc, muon_Posicion, muon_Deposito);
        gmuon[1]->SetPoint(muon_pc, muon_Posicion, muon_Restante);
        muon_pc++;
    }   //end of muon information loop
    */
    //Electron information
    double Rmax  = 25.0;        //max. distance from last muon hit (in cm) ~ electron analysis
    double theta;               //cone opening angle
    int theta_bin;
    double elec_DepoTot  = 0;
    double elec_DepoMask = 0; 
    double elec_DepoSph = 0;
    double elec_DepoSelec = 0;

    double sum_elec     = 0;
    double sum_selec    = 0;        

    TVector3 w_Hit;             //weighted electron hit position vector
    TVector3 elec_Bary;         //coordinates of barycenter shower
    TVector3 r_Axis;            //axis cone vector
    TVector3 r_ElecHit;
    TVector3 r_ElecHit_sel;
    
    //Let's create our cone
    for(size_t edepo_i = 0; edepo_i < ndepos; edepo_i++){
      const sim::SimEnergyDeposit& edepo = depolist->at(edepo_i);
  	  
	    //Retrieve position of energy deposit
	    geo::Point_t e_xyz = edepo.MidPoint();
      r_ElecHit(0) = e_xyz.X()-muon_LastHit.X();
      r_ElecHit(1) = e_xyz.Y()-muon_LastHit.Y();
      r_ElecHit(2) = e_xyz.Z()-muon_LastHit.Z();
      //Distance of electron hit from muon last hit
      double r_aux = sqrt(pow((e_xyz.X()-muon_LastHit.X()), 2)+pow((e_xyz.Y()-muon_LastHit.Y()), 2)+pow((e_xyz.Z()-muon_LastHit.Z()), 2));
      double theta_muon  = Muon_Axis.Angle(r_ElecHit)*TMath::RadToDeg(); 
      if(abs(edepo.PdgCode()) == 11){  
        elec_DepoTot+= edepo.Energy();  
        if(abs(theta_muon) < 20.0){continue;}   //No haga nada si el hit está cerca al muon track 
        elec_DepoMask += edepo.Energy();        //Solo e hits intentando evitar delta rays
        if(r_aux <= Rmax){
          w_Hit(0) += edepo.Energy()*e_xyz.X();
          w_Hit(1) += edepo.Energy()*e_xyz.Y();
          w_Hit(2) += edepo.Energy()*e_xyz.Z();
          elec_DepoSph += edepo.Energy();
        }
       }
    }   //end of depo loop

    hene->Fill(elec_DepoTot);
    hene_2->Fill(elec_DepoMask);
    hene_3->Fill(elec_DepoSph);

    if(elec_DepoSph <= 0) continue;
    //Fill barycenter shower vector
    for(int i = 0; i < 3; i++){
      elec_Bary(i) = w_Hit(i)/elec_DepoSph;
    }
    cout << "Barycenter: (" << elec_Bary.X() << "," << elec_Bary.Y() << "," << elec_Bary.Z() << ")" << endl;

    //Find normalized axis cone vector  
    r_Axis(0) = elec_Bary.X()-muon_LastHit.X();
    r_Axis(1) = elec_Bary.Y()-muon_LastHit.Y();
    r_Axis(2) = elec_Bary.Z()-muon_LastHit.Z();
    r_Axis    = r_Axis.Unit();

    vector<float> elec_DistriE(Nbins);
    vector<float> elec_DistriE_raw(Nbins);
  
    //Depo loop for electron information inside and outside containtment sphere
    for(size_t edepo_i = 0; edepo_i < ndepos; edepo_i++){
      const sim::SimEnergyDeposit& edepo = depolist->at(edepo_i);
  	  
	    //Retrieve position of energy deposit
	    geo::Point_t e_xyz = edepo.MidPoint();
      
      if(abs(edepo.PdgCode()) == 11){        
        r_ElecHit_sel(0) = e_xyz.X()-muon_LastHit.X();
        r_ElecHit_sel(1) = e_xyz.Y()-muon_LastHit.Y();
        r_ElecHit_sel(2) = e_xyz.Z()-muon_LastHit.Z();
        double r_aux = r_ElecHit_sel.Mag();
        r_ElecHit_sel = r_ElecHit_sel.Unit();
        double theta_muon  = Muon_Axis.Angle(r_ElecHit_sel)*TMath::RadToDeg();

        if(abs(theta_muon) < 20.0){continue;}
        sum_elec += 1;
        theta     = r_Axis.Angle(r_ElecHit_sel)*TMath::RadToDeg();    //angle between axis and hit in degrees
        theta_bin = floor(((theta*Nbins)/180));
        if(r_aux <= Rmax){  
          elec_DistriE[theta_bin] += edepo.Energy();
          if(theta <= 40){
            elec_DepoSelec += edepo.Energy();
          }
          //Fill angle distribution for selected events 
          htheta->Fill(theta);
          sum_selec += 1;
        } //end of containtment sphere condition
      }   //end of electron hit condition
    }     //end of depo loop

    hene_4->Fill(elec_DepoSelec);    
    
    double memory = 0.0;
    double memory_2 = 0.0;
    for(int i = 0; i < Nbins; i++){
      memory   += elec_DistriE[i]/elec_DepoTot;              //cuanta energía hay en el cono de selección con respecto a toda la energía e del evento
      memory_2 += elec_DistriE[i]/elec_DepoSph;              //cuanta energía hay en el cono de selección con respecto a toda la energía e de la esfera + track masking

      EnergyPerAngleEvent[i]             += elec_DistriE[i]/elec_DepoTot;
      SelectedEnergyPerAngleEvent[i]     += elec_DistriE[i]/elec_DepoSph;
      CompletenessOnlyAngleEvent[i]      += memory;      
      CompletenessRadiusAndAngleEvent[i] += memory_2;
    }
     
    cout << "# e Depositions: " << sum_elec << ", # Selected e Depositions: " << sum_selec << endl;
    cout << "Electron Total Deposit Energy: "<< elec_DepoTot <<" MeV || in Sphere: " << elec_DepoSph <<" MeV || in Cone: " << elec_DepoSelec << " MeV" << endl;
    ev_counter++;
    
  } //end of event loop

  cout << "# Events: " << ev_counter << endl;
  
  for(int i = 0; i < Nbins; i++){
    EnergyPerAngleEvent[i]        /= ev_counter;
    CompletenessOnlyAngleEvent[i] /= ev_counter;

    gene[0]->SetPoint(i, (i*180)/Nbins, EnergyPerAngleEvent[i]);
    gene[1]->SetPoint(i, (i*180)/Nbins, CompletenessOnlyAngleEvent[i]);

    SelectedEnergyPerAngleEvent[i]     /= ev_counter;
    CompletenessRadiusAndAngleEvent[i] /= ev_counter;

    gene_sel[0]->SetPoint(i, (i*180)/Nbins, SelectedEnergyPerAngleEvent[i]);
    gene_sel[1]->SetPoint(i, (i*180)/Nbins, CompletenessRadiusAndAngleEvent[i]);
  }

  double x_me;
  double y_dune;

  for(int i = 0; i < gene_sel[1]->GetN(); ++i){
    double x, y;
    gene_sel[1]->GetPoint(i, x, y);
    if(abs(y - 0.8) < 0.001){x_me = x;}
    if(x == 40){y_dune = y;}
  }

  cout << "Theta (for Completeness = 0.8): " << x_me << endl;
  cout << "Completeness (for theta = 40°): " << y_dune << endl;
  
  //---Plots---//
  
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
  TCanvas* cXYZE = new TCanvas("cXYZE", "G4 Muon Energy Deposition");
  cXYZE->Divide(2,2);
  for(int i = 0; i < hXYZ.size(); i++){
      cXYZE->cd(i+1); gPad->SetLogz(); hXYZ[i]->SetStats(0); hXYZ[i]->Draw("colz");
  }
  /*
  //Plot muon information
  TCanvas* cmuon = new TCanvas("cmuon", "Muon Information");
  cmuon->Divide(1,2);
  for(int i = 0; i < gmuon.size(); i++){
    cmuon->cd(i+1); gmuon[i]->Draw("AP");
  }  
  cmuon->Update();
  */
  TCanvas* chmuon = new TCanvas("chmuon", "Muon Distribution");
  chmuon->cd();
  hmuon->Draw(); 
  chmuon->Update();

  //Plot electron information
  TCanvas* ctheta = new TCanvas("ctheta", "Selected Theta Distribution");
  ctheta->cd();
  ctheta->SetLogy();
  htheta->SetStats(0);
  htheta->Draw(); 
  ctheta->Update();

  TCanvas* cene = new TCanvas("cene", "Total Energy Distribution");
  cene->cd();
  cene->SetLogy();  
  hene->SetStats(0);
  hene->Draw(); 
  cene->Update();

  TCanvas* cene_2 = new TCanvas("cene_2", "MTM");
  cene_2->cd();
  cene_2->SetLogy();
  hene_2->SetStats(0);
  hene_2->Draw(); 
  cene_2->Update();

  TCanvas* cene_3 = new TCanvas("cene_3", "MTM + CS");
  cene_3->cd();
  cene_3->SetLogy();
  hene_3->SetStats(0);
  hene_3->Draw(); 
  cene_3->Update();

  TCanvas* cene_4 = new TCanvas("cene_4", "MTM + CS + SC");
  cene_4->cd();
  cene_4->SetLogy();
  //hene_4->SetStats(0);
  hene_4->Draw(); 
  cene_4->Update();
  
  TCanvas* celec = new TCanvas("celec", "Electron Information");
  celec->Divide(1,2);
  for(int i = 0; i < gene.size(); i++){
    celec->cd(i+1); gene[i]->Draw("AP");
  }  
  celec->Update();

  TCanvas* celec_sel = new TCanvas("celec_sel", "Selected Electron Information");
  celec_sel->Divide(1,2);
  for(int i = 0; i < gene_sel.size(); i++){
    celec_sel->cd(i+1); gene_sel[i]->Draw("AP");
  }  
  celec_sel->Update();
  
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
  cXYZ->SaveAs(("../images/muongroove/"+name+"_cXYZ_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  cXYZ->SaveAs(("../images/muongroove/"+name+"cXYZ_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cXYZE->SaveAs(("../images/muongroove/"+name+"_cXYZE_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  cXYZE->SaveAs(("../images/muongroove/"+name+"_cXYZE_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str()); 
  /*
  cmuon->SaveAs(("../images/muongroove/muon_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cmuon->SaveAs(("../images/muongroove/muon_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());
  */
  chmuon->SaveAs(("../images/muongroove/"+name+"_energy_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  chmuon->SaveAs(("../images/muongroove/"+name+"_energy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());
 
  ctheta->SaveAs(("../images/muongroove/"+name+"_theta_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  ctheta->SaveAs(("../images/muongroove/"+name+"_theta_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cene->SaveAs(("../images/muongroove/"+name+"_Eenergy_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  cene->SaveAs(("../images/muongroove/"+name+"_Eenergy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cene_2->SaveAs(("../images/muongroove/"+name+"_MTMenergy_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  cene_2->SaveAs(("../images/muongroove/"+name+"_MTMenergy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cene_3->SaveAs(("../images/muongroove/"+name+"_CSenergy_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  cene_3->SaveAs(("../images/muongroove/"+name+"_CSenergy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  cene_4->SaveAs(("../images/muongroove/"+name+"_Selenergy_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  cene_4->SaveAs(("../images/muongroove/"+name+"_Selenergy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());
  
  celec->SaveAs(("../images/muongroove/"+name+"elec_energy_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  celec->SaveAs(("../images/muongroove/"+name+"elec_energy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  celec_sel->SaveAs(("../images/muongroove/"+name+"sel_elec_energy_"+to_string(event_i)+"_"+to_string(event_f)+".png").c_str()); 
  celec_sel->SaveAs(("../images/muongroove/"+name+"sel_elec_energy_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());
  
}