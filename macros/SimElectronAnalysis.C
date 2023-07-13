/*
  G4 analysis for electrons in LAr
*/
R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
R__ADD_INCLUDE_PATH("gallery/Event.h")
#include "tools.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include <iterator>

bool debug = false;

// ----- M A I N ----- //

void SimElectronAnalysis
(
  int nfiles      = 1,               //# files in your list                           
  int event_i     = 1,               //starting event
  int event_f     = 1,               //final event
  string simtag1  = "branch1",
  string simtag2  = "branch2", 
  string listname = "name",    //# files in your list
  bool Debug      = false
) 
{
  debug = Debug;
  SetDebug(debug);

  //All tags that u need
  art::InputTag sim_tag1(simtag1);
  art::InputTag sim_tag2(simtag2);  
  //Search ur files
  string filelist = "./list/"+listname+".list";
  vector<string> filenames = ReadFileList(nfiles, filelist);
  
  //X-axis
  vector<string> xtitle = {"X (cm)", "Z (cm)", "Z (cm)"};
  vector<float>  xmin   = {-325, -350, -350};
  vector<float>  xmax   = {325, 350, 350};
  //Y-axis
  vector<string> ytitle = {"Y (cm)", "Y (cm)", "X (cm)"};
  vector<float>  ymin   = {-650, -650, -325};
  vector<float>  ymax   = {650, 650, 325};
  vector<string> ytitle_elec = {"Deposited e^{-} Energy (MeV)", "Completeness"};

  //Setup graphs for each particle species (e+/e- mu+/mu- and p)
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes    = {11, -11, 13, -13, 2212};
  Color_t colors[5]       = {kRed+2, kOrange-3, kGreen+3, kBlue-7, kViolet-6};
 
  vector<vector<double>> meanComp(event_f-(event_i-1));
  int ev_i = 0;
  double Rmax = 100.0;
  double Rstep = 0.5;
  
  //---Plots---///

  //MC electron information
  TH1D* helec = new TH1D("helec", "Simulated Electron Energy Distribution; Energy (MeV); Counts", 180, 0, 100);
  helec->SetLineColor(colors[3]);
  //helec->SetFillColor(colors[3]);

  //Selected electron information
  vector<TGraph*> gelec(2);
  for(int i = 0; i < gelec.size(); i++){
    gelec[i] = new TGraph();
    gelec[i]->SetName(Form("gelec_%d", i));
    gelec[i]->SetMarkerColor(colors[i+2]);
    gelec[i]->SetLineColor(colors[i+2]);
    gelec[i]->SetMarkerStyle(3);
    gelec[i]->SetMarkerSize(0.05);
    gelec[i]->GetXaxis()->SetTitle("R (cm)");
    gelec[i]->GetYaxis()->SetTitle(ytitle_elec[i].c_str());
  }
  
  TGraph* gcomp = new TGraph();
  gcomp->SetMarkerColor(colors[4]);
  gcomp->SetLineColor(colors[4]);
  gcomp->SetMarkerStyle(3);
  gcomp->SetMarkerSize(0.5);
  gcomp->GetXaxis()->SetTitle("R (cm)");
  gcomp->GetYaxis()->SetTitle(ytitle_elec[1].c_str());
            

  int elec_pc = 0;   //to account for all points from all events

  //Go for events!
  for(gallery::Event ev(filenames); !ev.atEnd(); ev.next()){
    //Physical quantities
    vector<double> elecDist;                       //distance of electron hit from "last muon hit"
    vector<double> elecDepoE;                      //deposited energy of selected electrons
    vector<double> elecComp;                       //completeness vector

    if(ev.eventAuxiliary().event() < event_i) continue;
    if(ev.eventAuxiliary().event() > event_f) break;   
      
    cout << "Event: " << ev.eventAuxiliary().event() << endl;
      
    //Retrieve list of energy deposits per event
    auto const depolist = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(sim_tag1);
    size_t ndepos       = depolist->size();
    cout << "# Deposition: " << ndepos << endl;
    //Retrieve generated electron information   
    auto const montecarlo        = ev.getValidHandle<vector<simb::MCTruth>>(sim_tag2);
    const simb::MCTruth& gen     = montecarlo->at(0);
    const simb::MCParticle& elec = gen.GetParticle(0);
    
    //Electron information
    geo::Point_t elec_FirstHit;
    double elec_E0 = (1000)*elec.E();                           //in MeV
    cout << "Electron Energy: " << elec_E0 << " MeV" << endl;

    helec->Fill(elec_E0);
    
    for(double r = 0.0; r <= Rmax; r += Rstep){
      double elec_AllDepoE = 0;

      for(size_t depo_i = 0; depo_i < ndepos; depo_i++){
        const sim::SimEnergyDeposit& depo = depolist->at(depo_i);

        if(depo_i == 0){elec_FirstHit = depo.MidPoint();}
  	  
	      //Retrieve position of energy deposit
	      geo::Point_t xyz = depo.MidPoint();

        //Distance from electron first hit
        double r_aux = sqrt(pow((elec_FirstHit.X()-xyz.X()), 2)+pow((elec_FirstHit.Y()-xyz.Y()), 2)+pow((elec_FirstHit.Z()-xyz.Z()), 2));
      
        if(depo.PdgCode() == pdgcodes[0] && r_aux <= r){elec_AllDepoE += depo.Energy();}
      } //end of depo loop

      elecDist.push_back(r);
      elecDepoE.push_back(elec_AllDepoE);
      elecComp.push_back(elec_AllDepoE/elec_E0);
    } //end of radii loop
  
    //Fill selected electron information
    for(int i = 0; i < elecDist.size(); i++){
        double elec_Distancia = elecDist[i];
        double elec_Deposito  = elecDepoE[i]; 
        double elec_Completo  = elecComp[i];      
  
        gelec[0]->SetPoint(elec_pc, elec_Distancia, elec_Deposito);
        gelec[1]->SetPoint(elec_pc, elec_Distancia, elec_Completo);
        elec_pc++;        
    }

    //Event information
    meanComp[ev_i] = elecComp;
    ev_i++;
  } //end of event loop

  //Averaging completeness for each event
  vector<double> sums(meanComp[0].size(), 0);
  double sum = 0;
  for(size_t i = 0; i < meanComp[0].size(); i++){
    for(const auto & event : meanComp){sums[i] += event[i];}
    }

  for(size_t i = 0; i < sums.size(); i++) {
    sums[i] /= meanComp.size();
  }

  //Let's get back radii information
  vector<double> radiosto;
  for(double r = 0.0; r <= Rmax; r += Rstep){radiosto.push_back(r);}

  //Fill event information
  for(int i = 0; i < sums.size(); i++){
    double Radio         = radiosto[i];
    double comp_Promedio = sums[i];      
  
    gcomp->SetPoint(i, Radio, comp_Promedio);       
  }

  double y_me;
  double y_dune;

  for(int i = 0; i < gcomp->GetN(); ++i) {
    double x, y;
    gcomp->GetPoint(i, x, y);
    if(x == 25){y_me = y;}
    if(x == 20){y_dune = y;}
  }

  cout << "Completeness (R = 25 cm): " << y_me << endl;
  cout << "Completeness (R = 20 cm): " << y_dune << endl;

  //Plot mc electron information
  TCanvas* cene = new TCanvas("cene", "Energy Distribution");
  cene->cd();
  helec->Draw(); 
  cene->Update();

  //Plot selected electron information
  TCanvas* celec = new TCanvas("celec", "Electron Information");
  celec->Divide(1,2);
  for(int i = 0; i < gelec.size(); i++){
    celec->cd(i+1); gelec[i]->Draw("AP");
  }  
  celec->Update();

  TCanvas* comp = new TCanvas("comp", "Comp Information");
  comp->cd();
  gcomp->Draw("AP");
  comp->Update();

  //Save it! 
  cene->SaveAs(("../images/simelectronanalysis/mcEdist_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  cene->SaveAs(("../images/simelectronanalysis/mcEdist_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  celec->SaveAs(("../images/simelectronanalysis/gen_elec_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  celec->SaveAs(("../images/simelectronanalysis/gen_elec_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());

  comp->SaveAs(("../images/simelectronanalysis/comp_"+to_string(event_i)+"_"+to_string(event_f)+".pdf").c_str()); 
  comp->SaveAs(("../images/simelectronanalysis/comp_"+to_string(event_i)+"_"+to_string(event_f)+".root").c_str());
}