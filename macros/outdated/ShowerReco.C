R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/AnalysisBase/Calorimetry.h")
R__ADD_INCLUDE_PATH("lardataobj/RawData/RawDigit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Cluster.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Hit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/PFParticle.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Track.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Wire.h")

#include "tools.h"
#include "TGraph.h"
#include "TVector3.h"
#include <cmath>
#include <iostream>
#include <string>
#include <time.h>       //for clock_t, clock(), CLOCKS_PER_SEC
#include <utility>
#include <vector>

using namespace art;

//---List[Event][Track/Shower]([Point][X,Y,Z], "pdg code")---//
vector<vector<pair<vector<vector<double>>, int const>>> ListXYZ(string rootfilelabel, string recotag1, string recotag2)
{
  art::InputTag reco_tag1(recotag1);
  art::InputTag reco_tag2(recotag2);
  vector<string> filename(1, rootfilelabel);
    
  vector<vector<pair<vector<vector<double>>, int const>>> EvShXYZ;

  //Go for events!    
  for(gallery::Event ev(filename); !ev.atEnd(); ev.next()){

    cout << "Event: " << ev.eventAuxiliary().event() << endl;

    //Retrieve list of reconstructed particles        
    auto particles = ev.getValidHandle<vector<recob::PFParticle>>(reco_tag1);
    size_t npart   = particles->size();        
    const art::FindManyP<recob::SpacePoint> findSpacePoints(particles, ev, reco_tag1);     //here are the space point of each particle 
    vector<pair<vector<vector<double>>, int const>> ShXYZ;                                //position vector related with pdgcode
    
    //Retrieve space points for each particle
    for(size_t i = 0, sz = npart; i != sz; ++i){

      vector<const recob::SpacePoint*> SPs;
      const vector<art::Ptr<recob::SpacePoint>> pfpSpacePoints = findSpacePoints.at(particles->at(i).Self());       //pointer vector
      int const pdg = particles->at(i).PdgCode();
      //auto reco3d_SPs = ev.getValidHandle<vector<recob::SpacePoint>>(recotag2);           //Get the list of reconstructed SpacePoints (SPs)

      cout << "Particle: " << i << " ||  PDGCode: " << pdg << endl;

      for(auto pointer : pfpSpacePoints){

          SPs.push_back(pointer.get());

      }

      cout << "Pandora SP Size: " << SPs.size() << endl;

      art::FindManyP<recob::Hit> fmHits(SPs, ev, reco_tag1);   //Get the list of hits from each plane associated to each SP

      cout << "Pandora Hits Size: " << fmHits.size() << endl;

      vector<vector<double>> Pts(SPs.size());                 //position vector
      double Q = 0;

      for(size_t j = 0, sz = SPs.size(); j != sz; ++j){

        const vector<art::Ptr<recob::Hit>> Hits = fmHits.at(j);

        for (size_t k = 0, sz = Hits.size(); k != sz; ++k){
        
            cout << sz << " " << Hits[k]->WireID().Plane << endl;
        
          //if(Hits->v == 2){

            //Condition to avoid twice counted SPs
            double newQ = Hits[0]->Integral();

            if(newQ != Q){

			        //SpacePoints 
              double X = SPs.at(j)->XYZ()[0];   Pts.at(j).push_back(X);
              double Y = SPs.at(j)->XYZ()[1];   Pts.at(j).push_back(Y);
              double Z = SPs.at(j)->XYZ()[2];   Pts.at(j).push_back(Z);

              //Charge
              Pts.at(j).push_back(newQ);
				
              Q = newQ;

            }

          //}
        
        }

      cout << "4Vector: (" << Pts.at(j).at(0) << ", " << Pts.at(j).at(1) << ", " << Pts.at(j).at(2) << ", " << Pts.at(j).at(3) << ")" << endl; 

      }
      
      pair<vector<vector<double>>, int const> p(Pts , pdg);      //space points linked with pdgcode
      ShXYZ.push_back(p);

    }  //particle loop

    EvShXYZ.push_back(ShXYZ);

  }       

  return EvShXYZ;

}

//---Analysis---//
void Analysis(vector<pair<vector<vector<double>>, int const>> ShXYZ, int event_i)
{
  cout << "Event: " << event_i << endl;
  int npart = ShXYZ.size();
  cout << "# Particles: " << npart << endl;

  //Muon information
  int muIndex;

  //Electron information
  vector<double> elecIndex;

  for(size_t i = 0; i != npart; ++i){

    pair<vector<vector<double>>, int const> p = ShXYZ.at(i);        //paired space points - pdgcode
    vector<vector<double>> PtsXYZ = p.first;
    int const pdg = p.second;
    int npoints = PtsXYZ.size();

    if(pdg == 13){

      double muon_FirstPosition  = sqrt(pow(PtsXYZ.at(0).at(0), 2)+pow(PtsXYZ.at(0).at(1), 2)+pow(PtsXYZ.at(0).at(2), 2));
      double muon_LastPosition   = sqrt(pow(PtsXYZ.at(npoints-1).at(0), 2)+pow(PtsXYZ.at(npoints-1).at(1), 2)+pow(PtsXYZ.at(npoints-1).at(2), 2)); 
      double muon_TrackLength    = muon_LastPosition - muon_FirstPosition;

      if(abs(muon_TrackLength) >= 70){muIndex = i;}

    }

    if(pdg == 11){

      elecIndex.push_back(i);

    }

    cout << "Particle: " << i << " || PDGCode: " << pdg << " || # Space Points: " << npoints << endl;
  }  //particles loop

  pair<vector<vector<double>>, int const> muon_Selected = ShXYZ.at(muIndex);        //selected muon information
  vector<vector<double>> muon_XYZ = muon_Selected.first;
  int muon_Points = muon_XYZ.size();
  vector<TVector3> muPoints(muon_Points);

  for(size_t i = 0; i != muon_Points; ++i){
    
    //XYZs for selected muon
    muPoints[i](0) = muon_XYZ.at(i).at(0);  muPoints[i](1) = muon_XYZ.at(i).at(1);  muPoints[i](2) = muon_XYZ.at(i).at(2);
    cout << "Selected Muon Positions: (" << muPoints[i].X() << "," << muPoints[i].Y() << "," << muPoints[i].Z() << ")" << endl;    

  }
  
  cout << "Selected Muon Index: " << muIndex << " || # Points of Selected Muon: " << muon_Points << endl;
  cout << "Selected Muon Last Position: (" << muPoints[muon_Points-1].X() << "," << muPoints[muon_Points-1].Y() << "," << muPoints[muon_Points-1].Z() << ")" << endl;  

  for(size_t i = 0, sz = elecIndex.size(); i != sz; i++){

    pair<vector<vector<double>>, int const> elec_Selected = ShXYZ.at(elecIndex[i]);        //selected electron information
    vector<vector<double>> elec_XYZ = elec_Selected.first;
    int elec_Points = elec_XYZ.size();
    vector<TVector3> elecPoints(elec_Points);

    double muon_LastPosition  = muPoints[muon_Points-1].Mag();
    double elec_FirstPosition = sqrt(pow(elec_XYZ.at(0).at(0), 2)+pow(elec_XYZ.at(0).at(1), 2)+pow(elec_XYZ.at(0).at(2), 2));
    double elec_LastPosition  = sqrt(pow(elec_XYZ.at(elec_Points-1).at(0), 2)+pow(elec_XYZ.at(elec_Points-1).at(1), 2)+pow(elec_XYZ.at(elec_Points-1).at(2), 2));
    double elec_ShowerLength  = elec_LastPosition - elec_FirstPosition;
    double closeness          = elec_FirstPosition - muon_LastPosition;

    cout << "Electron ShowerLength: " << abs(elec_ShowerLength) << endl;

    if(abs(closeness) <= 10){

      for(size_t i = 0; i != elec_Points; ++i){
    
        //XYZs for selected muon
        elecPoints[i](0) = elec_XYZ.at(i).at(0);  elecPoints[i](1) = elec_XYZ.at(i).at(1);  elecPoints[i](2) = elec_XYZ.at(i).at(2);
        cout << "Selected Electron Positions: (" << elecPoints[i].X() << "," << elecPoints[i].Y() << "," << elecPoints[i].Z() << ")" << endl;
        

      }   

    cout << "Selected Electron Index: " << elecIndex[i] << " || # Points of Selected Electron: " << elec_Points << endl;
    cout << "Selected Electron First Position: (" << elecPoints[0].X() << "," << elecPoints[0].Y() << "," << elecPoints[0].Z() << ")" << endl;   
    
    }
    else{continue;}

  }

}

//---Plots---//
void Plotting(vector<pair<vector<vector<double>>, int const>> ShXYZ, string config = "ShowerReco")
{
  //Setup graphs for each particle species (e+/e- mu+/mu- and p)
  vector<string> pdgnames = {"e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "p^{+}"};
  vector<int> pdgcodes = {11, -11, 13, -13, 2212};
  Color_t colors[6]    = {kBlue+2, kAzure+2, kOrange-3, kRed+2, kGreen+2, kViolet-5};

  //Tracks - showers
  vector<TMultiGraph*> mgXYZ(3);
  for(int i = 0; i < mgXYZ.size(); i++){
    mgXYZ[i] = new TMultiGraph();
  }

  TLegend* leg = new TLegend(0.2, 0.2, 0.8, 0.8);
  leg->SetLineWidth(0); leg->SetHeader("Identified particles:");
    
  for(size_t i = 0, sz = ShXYZ.size(); i != sz; ++i){

    pair<vector<vector<double>>, int const> p = ShXYZ.at(i);        //paired space points - pdgcode
    vector<vector<double>> PtsXYZ = p.first;
    int npoints = PtsXYZ.size();
    double X[npoints], Y[npoints], Z[npoints];
    int const pdg = p.second;

    for(size_t j = 0; j != npoints; ++j){
      
      X[j] = PtsXYZ.at(j).at(0);  Y[j] = PtsXYZ.at(j).at(1);  Z[j] = PtsXYZ.at(j).at(2);        //XYZs for each type of particle

    }

    TGraph* gXY = new TGraph(npoints, X, Y);
    if(pdg == pdgcodes[0]){gXY->SetMarkerColor(colors[0]); gXY->SetMarkerStyle(7); mgXYZ[0]->Add(gXY, "P"); mgXYZ[0]->SetTitle("Z Projection; X (cm); Y (cm)");}
    else if(pdg == pdgcodes[2]){gXY->SetMarkerColor(colors[2]); mgXYZ[0]->Add(gXY, "P");}
    else{mgXYZ[0]->Add(gXY, "C");}      //draw a line for everything else     
        
    TGraph* gYZ = new TGraph(npoints, Y, Z);
    if(pdg == pdgcodes[0]){gYZ->SetMarkerColor(colors[0]); gYZ->SetMarkerStyle(7); mgXYZ[1]->Add(gYZ, "P"); mgXYZ[1]->SetTitle("X Projection; Y (cm); Z (cm)");}
    else if(pdg == pdgcodes[2]){gYZ->SetMarkerColor(colors[2]); mgXYZ[1]->Add(gYZ, "P");}
    else{mgXYZ[1]->Add(gYZ, "C");}

    TGraph* gXZ = new TGraph(npoints, X, Z);
    if(pdg == pdgcodes[0]){gXZ->SetMarkerColor(colors[0]); gXZ->SetMarkerStyle(7); mgXYZ[2]->Add(gXZ, "P"); mgXYZ[2]->SetTitle("Y Projection; X (cm); Z (cm)");}
    else if(pdg == pdgcodes[2]){gXZ->SetMarkerColor(colors[2]); mgXYZ[2]->Add(gXZ, "P");}
    else{mgXYZ[2]->Add(gXZ, "C");}

    //Appropiate legend    
    if(pdg == pdgcodes[0]){leg->AddEntry(gXY, ("Particle: " + to_string(i) + " || " + pdgnames[0]).c_str(), "L");}
    else if(pdg == pdgcodes[2]){leg->AddEntry(gXY, ("Particle: " + to_string(i) + " || " + pdgnames[2]).c_str(), "L");}
    else{leg->AddEntry(gXY, ("Particle: " + to_string(i) + " || PDG: " + to_string(pdg)).c_str(), "L");}

  }

  TCanvas *cXYZ = new TCanvas("cXYZ", "Reco Spatial Distributions");
  cXYZ->Divide(2,2);
  for(int i = 0; i < mgXYZ.size(); i++){
    cXYZ->cd(i+1); mgXYZ[i]->Draw("APL");
  }
  cXYZ->cd(4);
  leg->Draw();

  cXYZ->SaveAs(("../images/reco/cXYZ_" + config + ".pdf").c_str());
  cXYZ->SaveAs(("../images/reco/cXYZ_" + config + ".root").c_str());
    
}

//-----Main-----//

void ShowerReco
(
  int nfiles = 1,
  int event_i = 1,
  string recotag1 = "branch1",
  string recotag2 = "branch2",
  string prod = "list owner"
)
{
  //Search ur files
  string filelist = "./list/"+prod+"_reco.list";
  vector<string> filenames = ReadFileList(nfiles, filelist); 

  ListXYZ(filenames[nfiles-1], recotag1, recotag2); 

  //Analysis(ListXYZ(filenames[0], recotag1, recotag2).at(event_i-1), event_i);
  //Plotting(ListXYZ(filenames[0]).at(event_i-1), "Cosmics");
}

