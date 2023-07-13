R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/AnalysisBase/Calorimetry.h")
R__ADD_INCLUDE_PATH("lardataobj/RawData/RawDigit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Cluster.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Hit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/PFParticle.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Track.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Wire.h")
//R__ADD_INCLUDE_PATH("protoduneana/Utilities/ProtoDUNEPFParticleUtils.h")

#include "tools.h"
#include "TGraph.h"
#include "TVector3.h"
#include <cmath>
#include <iostream>
#include <string>
#include <time.h>       //for clock_t, clock(), CLOCKS_PER_SEC
#include <utility>
#include <vector>
//#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

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

        art::InputTag cluster_tag("pandora");
        art::InputTag hit_tag("gaushit");
        auto const hitHandle = ev.getValidHandle<std::vector<recob::Hit>>(hit_tag);
        std::vector<art::Ptr<recob::Hit>> hits;
        art::fill_ptr_vector(hits, hitHandle);
        const art::FindManyP<recob::SpacePoint> fmsph(hitHandle, ev, cluster_tag);
        

       for (auto hit : hits) {
            std::vector<art::Ptr<recob::SpacePoint>> sps = fmsph.at(hit.key());
            if (sps.size() == 1) {
                art::Ptr<recob::SpacePoint> sp = sps.front();
                cout<<hit.key()<<" "<<hit->WireID().Plane<<" "<< sp->XYZ()[0]<<"  " <<sp->ID()<<endl;
            }

        }
    }
    return EvShXYZ;
}
void TestReco()
{
  ListXYZ("/silver/DUNE/andres-cortes/output/protodunevd_reco2.root", "pandora", "reco3d"); 

}


