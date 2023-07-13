R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
R__ADD_INCLUDE_PATH("gallery/Event.h")
#include "tools.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include <iterator>

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

//----Main-----//

void MuonGroove
(
  int nfiles = 1,
  int event_i = 1,
  int event_f = 1,
  string simtag1 = "branch1",
  string simtag2 = "branch2", 
  string prod = "list owner",
  bool Debug = false
)
 
{
  debug = Debug;
  SetDebug(debug);

  //All tags that u need
  art::InputTag sim_tag1(simtag1);
  art::InputTag sim_tag2(simtag2);  
  //Search ur files
  string filelist = "./list/"+prod+"_muon.list";
  vector<string> filenames = ReadFileList(nfiles, filelist);
  
  //X-axis
  vector<string> xtitle = {"Z (cm)", "X (cm)", "Z (cm)"};
  vector<float>  xmin   = {0, -325, 0};
  vector<float>  xmax   = {900, 325, 900};
  //Y-axis
  vector<string> ytitle = {"X (cm)", "Y (cm)", "Y (cm)"};
  vector<float>  ymin   = {-325, -650, -650};
  vector<float>  ymax   = { 325, 650, 650};
  vector<string> ytitle_muon = {"Deposited #mu^{-} Energy (MeV)", "Remaining #mu^{-} Energy (MeV)"};

  //Setup graphs for each particle species (e+/e- mu+/mu- and p)
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes    = {11, -11, 13, -13, 2212};
  Color_t colors[5]       = {kBlue+2, kAzure+2, kOrange-3, kRed, kGreen+2};

  //Physical quantities
  vector<double> muPos;
  vector<double> muDepoE;
  vector<double> muRemE;
  vector<double> elecDist;    //distance of electron hit from last muon hit
  vector<double> elecDepoE;   //deposited energy of selected electrons
  
  //---Plots---///
  
  //Tracks
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

  //Energy deposition
  vector<TH2F*> hXYZ(3);
  for(int i = 0; i < hXYZ.size(); i++){ 
    hXYZ[i] = new TH2F(Form("hXYZ_%d", i), Form(";%s;%s", xtitle[i].c_str(), ytitle[i].c_str()), int(xmax[i]-xmin[i]), xmin[i], xmax[i], int(ymax[i]-ymin[i]), ymin[i], ymax[i]);
  }

  //Muon information
  vector<TGraph*> hmuon(2);
  for(int i = 0; i < hmuon.size(); i++){
    hmuon[i] = new TGraph();
    hmuon[i]->SetName(Form("hmuon_%d", i));
    hmuon[i]->SetMarkerColor(colors[i]);
    hmuon[i]->SetLineColor(colors[i]);
    //hmuon[i]->SetMarkerSize(2.0);
    hmuon[i]->GetXaxis()->SetTitle("Track Length (cm)");
    hmuon[i]->GetYaxis()->SetTitle(ytitle_muon[i].c_str());
  }

  //Electron information
  TGraph* helec = new TGraph();
  helec->SetMarkerColor(colors[3]);
  helec->SetLineColor(colors[3]);
  helec->SetMarkerSize(2.0);
  helec->GetXaxis()->SetTitle("R (cm)");
  helec->GetYaxis()->SetTitle("Total Deposited e Energy (MeV)");
  
  
  //Go for events!
  for(gallery::Event ev(filenames); !ev.atEnd(); ev.next()){

    if(ev.eventAuxiliary().event() < event_i) continue;
    if(ev.eventAuxiliary().event() > event_f) break;   
      
    cout << "Event: " << ev.eventAuxiliary().event() << endl;
      
    //Retrieve list of energy deposits per event
    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag1);
    size_t ndepos = depolist->size();
    cout << "# Deposition: " << ndepos << endl;
    //Retrieve generated muon information   
    auto const montcarlo = ev.getValidHandle<vector<simb::MCTruth>>(sim_tag2);
    const simb::MCTruth& gen = montcarlo->at(0);
    const simb::MCParticle& muon = gen.GetParticle(0);

    geo::Point_t muon_FirstHit;
    double muon_E0 = (1000)*muon.E();   //in MeV
    double muon_AllDepoE;
    int    muon_LastDepoIt = -1;
    cout << "Muon Energy: " << muon_E0 << " MeV" << endl;

    //Loop over list of depos in this event
    for(size_t depo_i = 0; depo_i < ndepos; depo_i++){

      const sim::SimEnergyDeposit& depo = depolist->at(depo_i);
  	  
	    //Retrieve position of energy deposit
	    geo::Point_t xyz = depo.MidPoint();
      	  
	    //Check if particle is within the pdgcodes list - otherwise skip it
	    vector<int>::iterator it;
	    it = std::find(pdgcodes.begin(), pdgcodes.end(), depo.PdgCode());
      
	    if(it != pdgcodes.end()){
	      gXYZ[it - pdgcodes.begin()][0]->SetPoint(gXYZ[it - pdgcodes.begin()][0]->GetN(), xyz.Z(), xyz.X());
	      gXYZ[it - pdgcodes.begin()][1]->SetPoint(gXYZ[it - pdgcodes.begin()][1]->GetN(), xyz.X(), xyz.Y());
	      gXYZ[it - pdgcodes.begin()][2]->SetPoint(gXYZ[it - pdgcodes.begin()][2]->GetN(), xyz.Z(), xyz.Y());
	    }

      //Muon information
      if(depo.PdgCode() == pdgcodes[2]){
        //Map energy deposit
        hXYZ[0]->Fill(xyz.Z(), xyz.X(), depo.Energy());
	      hXYZ[1]->Fill(xyz.X(), xyz.Y(), depo.Energy());
	      hXYZ[2]->Fill(xyz.Z(), xyz.Y(), depo.Energy());
        
        if(depo_i == 0){muon_FirstHit = depo.MidPoint();}

        double muon_ActualPosition = sqrt(pow(xyz.X(), 2)+pow(xyz.Y(), 2)+pow(xyz.Z(), 2));
        double muon_FirstPosition  = sqrt(pow(muon_FirstHit.X(), 2)+pow(muon_FirstHit.Y(), 2)+pow(muon_FirstHit.Z(), 2));
        double muon_TrackLength    = muon_FirstPosition - muon_ActualPosition;
        double muon_RemEnergy      = muon_E0 - depo.Energy();
        muon_AllDepoE             += depo.Energy();
        
        //Energy deposit and position for muon
        muPos.push_back(abs(muon_TrackLength));
        muDepoE.push_back(depo.Energy());
        muRemE.push_back(muon_RemEnergy);
        
        //cout << depo_i << " " << xyz.X() << " "<< depo.Energy() << " " << depo.PdgCode() << " " << muon_AllDepoE << endl;

        muon_E0 = muon_RemEnergy;    //update the actual muon energy
        muon_LastDepoIt = depo_i;    //save muon last hit
      }
	    //Skip particles that does not have the correct pdgcode 
	    else {continue;}
	  } //end of depo loop 
    
    const sim::SimEnergyDeposit& muon_LastDepo = depolist->at(muon_LastDepoIt);
    geo::Point_t muon_LastHit = muon_LastDepo.MidPoint();
    cout << "Muon Total Deposit Energy: " << muon_AllDepoE << " MeV" << endl;
    cout << "Muon Last Hit: " << muon_LastDepoIt << " || " << muon_LastHit << endl;
    //Fill muon information
    for (int i = 0; i < muPos.size(); i++){
        //Get muon position, deposited and remaining energy for each hit 
        double muon_Posicion = muPos[i];
        double muon_Deposito = muDepoE[i];
        double muon_Restante = muRemE[i];       
  
        hmuon[0]->SetPoint(i, muon_Posicion, muon_Deposito);
        hmuon[1]->SetPoint(i, muon_Posicion, muon_Restante);
    }

    //Electron information

    double delta = 10;                       //variation from last muon hit (in cm)
    double Rmax  = sqrt(3*pow(delta, 2));    //max. distance from last muon hit (in cm)
    
    //double elec_Total = 0;

    for(double r = 0.0; r <= Rmax; r += 0.001){

      double elec_DepoAux = 0;
      for(size_t edepo_i = 0; edepo_i < ndepos; edepo_i++){

        const sim::SimEnergyDeposit& edepo = depolist->at(edepo_i);
  	  
	      //Retrieve position of energy deposit
	      geo::Point_t e_xyz = edepo.MidPoint();
        //Distance from muon last hit
        double r_aux = sqrt(pow((muon_LastHit.X()-e_xyz.X()), 2)+pow((muon_LastHit.Y()-e_xyz.Y()), 2)+pow((muon_LastHit.Z()-e_xyz.Z()), 2));
      
        if(edepo.PdgCode() == pdgcodes[0] && r_aux <= r){
          elec_DepoAux += edepo.Energy();
          //cout << "Selected electrons:" << edepo_i << " " << edepo.PdgCode() << endl;
          //elec_Total += 1;    //selected electrons counting 
        }
      }   //end depo loop

      elecDist.push_back(r);
      elecDepoE.push_back(elec_DepoAux);
      
      //cout << r << " " << elec_DepoAux << endl;

    }
    
    //Fill electron information
    for (int i = 0; i < elecDist.size(); i++){
        // 
        double elec_Distancia = elecDist[i];
        double elec_Deposito = elecDepoE[i];  

        //cout << "--------" << endl;
        //cout << elec_Distancia << " " << elec_Deposito << endl;    
  
        helec->SetPoint(i, elec_Distancia, elec_Deposito);
    }

  } //end of event loop

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

  //Plot muon information
  TCanvas* cmuon = new TCanvas("cmuon", "Muon Information");
  cmuon->Divide(1,2);
  for(int i = 0; i < hmuon.size(); i++){
    cmuon->cd(i+1); hmuon[i]->Draw("AP");
  }  
  cmuon->Update();

  //Plot electron information
  TCanvas* celec = new TCanvas("celec", "Elec Information");
  celec->cd();
  helec->Draw("AP"); 
  celec->Update();

  //Build the legend to display the particle names that belong to each graph
  TLegend* leg = new TLegend(0.2, 0.2, 0.8, 0.8);
  leg->SetLineWidth(0); leg->SetHeader("Particle list:");
  for(int i = 0; i < gXYZ.size(); i++){ 
      bool legSet = false;
      for(int j = 0; j < gXYZ[i].size(); j++){ 
	      if(gXYZ[i][j]->GetN() > 0 && !legSet){legSet = true; leg->AddEntry(gXYZ[i][j],pdgnames[i].c_str(), "L");}
	    }
  }
  
  //Choose where you put legend
  cXYZ->cd(4);
  leg->Draw();

  //Save it!

  cXYZE->SaveAs("../images/muongroove/cXYZE.pdf"); 
  cXYZE->SaveAs("../images/muongroove/cXYZE.root"); 
  
  cXYZ->SaveAs("../images/muongroove/cXYZ.pdf"); 
  cXYZ->SaveAs("../images/muongroove/cXYZ.root"); 

  cmuon->SaveAs("../images/muongroove/muon_info.pdf"); 
  cmuon->SaveAs("../images/muongroove/muon_info.root");

  celec->SaveAs("../images/muongroove/elec_info.pdf"); 
  celec->SaveAs("../images/muongroove/elec_info.root");
}
