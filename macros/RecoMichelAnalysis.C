/*
  Linking space points from gaushit with space points from pandora with pdfcode
*/
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
vector<vector<pair<vector<double>, int const>>> ListXYZ(string rootfilelabel, string recotag1, string recotag2)
{
  art::InputTag reco_tag1(recotag1);
  art::InputTag reco_tag2(recotag2);
  vector<string> filename(1, rootfilelabel);
        
  vector<vector<pair<vector<double>, int const>>> EvXYZ_Hit;

  for(gallery::Event ev(filename); !ev.atEnd(); ev.next()){
    cout << "Event: " << ev.eventAuxiliary().event() << endl;

    //---Gaushit---//
    auto const hitHandle = ev.getValidHandle<std::vector<recob::Hit>>(reco_tag2);   //access to gaushit information
    vector<art::Ptr<recob::Hit>> hits;                                              //pointers vector
    art::fill_ptr_vector(hits, hitHandle);                                          //fills vector with pointers

    const art::FindManyP<recob::SpacePoint> fmsph(hitHandle, ev, reco_tag1);

    vector<vector<double>> Pts(hits.size());
    vector<int> sp_ID(hits.size(), -1);
    double Q = 0;

    for(size_t i = 0, sz_i = hits.size(); i != sz_i; i++){
      for(int j = 0; j < 5; j++) Pts.at(i).push_back(-999);

      vector<art::Ptr<recob::SpacePoint>> sps = fmsph.at(i);                        //SP associated with this hit

      if(sps.size() == 1){                                                          //if there is a SP related with the hit
        art::Ptr<recob::SpacePoint> sp = sps.front();                               //retrieve SP

        //We need just the collective plane information
        if(hits[i]->WireID().Plane == 2){
          double newQ = hits[i]->SummedADC();                                        //to avoid twice counted SPs                                   

          if(newQ != Q){
            //SpacePoints 
            Pts.at(i).at(0) = sp->XYZ()[0];
            Pts.at(i).at(1) = sp->XYZ()[1];
            Pts.at(i).at(2) = sp->XYZ()[2];
            //Charge
            Pts.at(i).at(3) = newQ;
            //Plane ID
            Pts.at(i).at(4) = hits[i]->WireID().Plane;
            //SP ID
            sp_ID[i] = sp->ID();

            Q = newQ;
          }   //end Q cond
        }     //end plane cond
      }       //end sp cond
    }         //end hits loop
    
    vector<pair<vector<double>, int const>> XYZ_Hit;                                  //SP + charge + plane ID vector related to ID
    for(size_t i = 0, sz_i = Pts.size(); i != sz_i; i++){
      pair<vector<double>, int const> p_hits(Pts.at(i), sp_ID[i]);
      XYZ_Hit.push_back(p_hits);
    } 

    //---PFParticle---//
    //Retrieve list of reconstructed particles        
    auto particles = ev.getValidHandle<vector<recob::PFParticle>>(reco_tag1);
    size_t npart   = particles->size();        
    const art::FindManyP<recob::SpacePoint> findSpacePoints(particles, ev, reco_tag1);     //here are the space point of each particle 
    vector<pair<vector<vector<double>>, int const>> ShXYZ;

    vector<vector<vector<double>>> partPts(npart);
    vector<vector<double>> pfpTotalPts;
    vector<vector<int>> SPID(npart);
    vector<double> pfpSPID;

    //Retrieve space points for each particle
    for(size_t i = 0, sz_i = npart; i != sz_i; ++i){

      vector<const recob::SpacePoint*> pfpSPs;
      const vector<art::Ptr<recob::SpacePoint>> pfpSpacePoints = findSpacePoints.at(particles->at(i).Self());       //pointer vector
      int const pdg = particles->at(i).PdgCode();


      for(auto pointer : pfpSpacePoints){pfpSPs.push_back(pointer.get());}

      vector<vector<double>> pfpPts(pfpSPs.size());                 

      for(size_t j = 0, sz_j = pfpSPs.size(); j != sz_j; ++j){
        //SpacePoints
        pfpPts.at(j).push_back(pfpSPs.at(j)->XYZ()[0]);
        pfpPts.at(j).push_back(pfpSPs.at(j)->XYZ()[1]);
        pfpPts.at(j).push_back(pfpSPs.at(j)->XYZ()[2]);
        //PDGCode
        pfpPts.at(j).push_back(pdg);
        //SP ID
        SPID[i].push_back(pfpSPs.at(j)->ID());
      } //end sp loop
      partPts[i] = pfpPts;
    }   //end particle loop

    vector<pair<vector<double>, int const>> pfpXYZ;         //SP + pdgcode vector related to ID
    for(size_t i = 0, sz_i = npart; i != sz_i; i++){
        int npts = partPts[i].size();
        for(size_t j = 0, sz_j = npts; j != sz_j; j++){          
          pair<vector<double>, int const> p(partPts[i].at(j), SPID[i].at(j));
          pfpXYZ.push_back(p);
        }
    }   //end fill pfpXYZ loop

    cout << "# SP's from Hits: " << XYZ_Hit.size() << " || # SP's from PFP:" << pfpXYZ.size()  << endl;

    //Looking for matches
    size_t size1 = XYZ_Hit.size();
    size_t size2 = pfpXYZ.size();
    int sum = 0;

    for(size_t i = 0; i < size1; i++){
      for(size_t j = 0; j < size2; j++){
        if(XYZ_Hit.at(i).second == pfpXYZ.at(j).second){sum += 1;}
      }
    }

    cout << "# of matches: " << sum << endl;
    EvXYZ_Hit.push_back(XYZ_Hit);
  } //end event loop    

  return EvXYZ_Hit;
}

// ----- M A I N ----- //

void RecoMichelAnalysis
(
  int nfiles      = 1,               //# files in your list                           
  //int event_i   = 1,               //starting event
  //int event_f   = 1,               //final event
  string recotag1 = "branch1",       //"pandora" for PFParticle
  string recotag2 = "branch2",       //"gaushit" for hit information
  string listname = "name"
)
{
  //Search ur files
  string filelist = "./list/"+listname+".list";
  vector<string> filenames = ReadFileList(nfiles, filelist); 

  ListXYZ(filenames[nfiles-1], recotag1,  recotag2); 
}