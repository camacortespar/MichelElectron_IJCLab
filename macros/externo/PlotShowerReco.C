/*

*/
R__ADD_INCLUDE_PATH("lardataobj/RawData/RawDigit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Wire.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Hit.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Cluster.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/Track.h")
R__ADD_INCLUDE_PATH("lardataobj/RecoBase/PFParticle.h")
R__ADD_INCLUDE_PATH("lardataobj/AnalysisBase/Calorimetry.h")
R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("protoduneana/Utilities/ProtoDUNEShowerUtils.h")
#include <vector>
#include <cmath>
#include <iostream>
#include "TGraph.h"
#include <string>
#include <utility>
#include "TVector3.h"
#include <time.h>       // for clock_t, clock(), CLOCKS_PER_SEC

using namespace art;


std::vector<std::vector< pair< std::vector<std::vector<double>> , int const> >> ListXYZ(
    std::string rootfilelabel) // List[Event][Shower]([Point][X,Y,Z] , "code pdg")
{
    art::InputTag Module("pandora");
    vector<string> filename(1,rootfilelabel);
    
    std::vector<std::vector< pair< std::vector<std::vector<double>> , int const> >> PtsXYZShEv;
    
    for (gallery::Event ev(filename); !ev.atEnd(); ev.next())
    {
        auto const& particles = *ev.getValidHandle<vector<recob::PFParticle>>(Module);
        auto particless = ev.getValidHandle<vector<recob::PFParticle>>(Module);        
        const art::FindManyP<recob::SpacePoint> findSpacePoints(particless, ev, Module);
        std::vector< pair< std::vector<std::vector<double>> , int const> > PtsXYZSh;
        for (size_t i = 0, sz = particles.size(); i != sz; ++i)
        {
            std::vector<const recob::SpacePoint*> SPs;
            const std::vector<art::Ptr<recob::SpacePoint>> pfpSpacePoints = findSpacePoints.at(particles.at(i).Self());
            int const pdg = particles.at(i).PdgCode();
            
            for(auto pointer : pfpSpacePoints)
            {
                SPs.push_back(pointer.get());
            }
            
            std::vector<std::vector<double>> PtsXYZ(SPs.size());
            
            for (unsigned j=0, sz = SPs.size(); j!=sz;++j)
            {
                double X = SPs.at(j)->XYZ()[0];
                PtsXYZ.at(j).push_back(X);
                double Y = SPs.at(j)->XYZ()[1];
                PtsXYZ.at(j).push_back(Y);
                double Z = SPs.at(j)->XYZ()[2];
                PtsXYZ.at(j).push_back(Z);
            }
            pair< std::vector<std::vector<double>> , int const> p(PtsXYZ , pdg);
            PtsXYZSh.push_back(p);
        }
        PtsXYZShEv.push_back(PtsXYZSh);
    }
    return PtsXYZShEv;
}

void cout_ListXYZ(std::vector<std::vector< pair< std::vector<std::vector<double>> , int const> >> PtsXYZShEv)
{
    for (unsigned e=0, sz=PtsXYZShEv.size(); e!=sz; e++)
    {
        cout << endl << "*********************** EVENT #" << e <<"***********************" << endl;
        std::vector< pair< std::vector<std::vector<double>> , int const> > PtsXYZSh = PtsXYZShEv.at(e);
        for (size_t i = 0, sz = PtsXYZSh.size(); i != sz; ++i)
        {
            pair< std::vector<std::vector<double>> , int const> p = PtsXYZSh.at(i);
            std::vector<std::vector<double>> PtsXYZ = p.first;
            int pdg = p.second;        
            
            cout << endl << "Particle #" << i << "------  Code Pdg : " << pdg << "--------------" << endl;
            for (unsigned j=0, sz = PtsXYZ.size(); j!=sz;++j)
            {
                double X = PtsXYZ.at(j).at(0);
                double Y = PtsXYZ.at(j).at(1);
                double Z = PtsXYZ.at(j).at(2);
                cout << "SpacePoint n:"<< j <<"    (X: "<< X <<"   Y: "<< Y <<"  Z: "<< Z <<") "<<endl;
            }
        }
    }
}

void plotproj_1e_ShowerReco(std::vector< pair< std::vector<std::vector<double>> , int const> > PtsXYZSh, std::string title = "ShowerReco")
{
    TCanvas *c1 = new TCanvas("c1",title.c_str());
    c1->cd();
    c1->Divide(2,2);
    TMultiGraph* mgXY = new TMultiGraph();
    TMultiGraph* mgYZ = new TMultiGraph();
    TMultiGraph* mgXZ = new TMultiGraph();
    auto legend = new TLegend(0,0,1,1);
    legend->SetHeader("-Identified particles-","C");
    
    for (size_t i = 0, szSh = PtsXYZSh.size(); i != szSh; ++i)
    {
        pair< std::vector<std::vector<double>> , int const> p = PtsXYZSh.at(i);
        std::vector<std::vector<double>> PtsXYZ = p.first;
        int const pdg = p.second;
        int sz = PtsXYZ.size();
        double X[sz], Y[sz], Z[sz];
        for (unsigned j=0; j!=sz;++j)
        {
            X[j] = PtsXYZ.at(j).at(0);
            Y[j] = PtsXYZ.at(j).at(1);
            Z[j] = PtsXYZ.at(j).at(2);
        }
        int color = 99-(4*i);
        int marker = 6;
        TGraph* gXY = new TGraph(sz,X,Y);
        gXY->SetMarkerColor(color);
        gXY->SetMarkerStyle(marker);
        if (pdg==11)
        {
            mgXY->Add(gXY,"P");
        }
        else
        {
            mgXY->Add(gXY,"C");
        }
        
        TGraph* gYZ = new TGraph(sz,Y,Z);
        gYZ->SetMarkerColor(color);
        gYZ->SetMarkerStyle(marker);
        if (pdg==11)
        {
            mgYZ->Add(gYZ,"P");
        }
        else
        {
            mgYZ->Add(gYZ,"C");
        }
        
        TGraph* gXZ = new TGraph(sz,X,Z);
        gXZ->SetMarkerColor(color);
        gXZ->SetMarkerStyle(marker);
        if (pdg==11)
        {
            mgXZ->Add(gXZ,"P");
        }
        else
        {
            mgXZ->Add(gXZ,"C");
        }
        std::string label = "Particle " + to_string(i) + "  -> type : " + to_string(pdg);
        if (pdg==11)
        {
            gXY->SetFillStyle(3002);
            TLegendEntry *le = legend->AddEntry(gXY,label.c_str(),"f");
            le->SetFillStyle(3002);
            le->SetTextSize(0.05);
        }
        else
        {
            TLegendEntry *le = legend->AddEntry(gXY,label.c_str(),"L");
            le->SetTextSize(0.05);
        }
        
        //le->SetMarkerStyle(21);
    }
    c1->cd(4);
    legend->SetNColumns(2);
    legend->Draw();
    
    c1->cd(1);
    mgXY->SetTitle("Z Projection;X(cm);Y(cm)");
    mgXY->Draw("A p l");
    
    c1->cd(2);
    mgYZ->SetTitle("X Projection;Y(cm);Z(cm)");
    mgYZ->Draw("A p l");
    
    c1->cd(3);
    mgXZ->SetTitle("Y Projection;X(cm);Z(cm)");
    mgXZ->Draw("A p l");
    
}



// Test possible :
/*
    std::vector<std::vector<std::vector<std::vector<double>>>> PtsXYZShEv = ListXYZ("../output/10elec_1GeV_shower_reco.root");
    std::vector<std::vector<std::vector<double>>> PtsXYZSh = PtsXYZShEv.at(0);
    plotproj_1e_ShowerReco(PtsXYZSh);
*/



/*
void plotproj_ShowerReco(std::vector<std::vector<std::vector<std::vector<double>>>> PtsXYZShEv, std::string title = "ShowerReco")
{
    int NbEv = PtsXYZShEv.size();
    for (int e=0; e!=NbEv; e++)
    {
        TCanvas *c1 = new TCanvas("c1",title.c_str());
        c1->cd();
        c1->Divide(2,2);
        TMultiGraph* mgXY = new TMultiGraph();
        TMultiGraph* mgYZ = new TMultiGraph();
        TMultiGraph* mgXZ = new TMultiGraph();
        auto legend = new TLegend(0,0,1,1);
        legend->SetHeader("-Identified particles-","C");
        std::vector<std::vector<std::vector<double>>> PtsXYZSh = PtsXYZShEv.at(e);
        
        for (size_t i = 0, szSh = PtsXYZSh.size(); i != szSh; ++i)
        {
            pair< std::vector<std::vector<double>> , int const> p = PtsXYZSh.at(i);
            std::vector<std::vector<double>> PtsXYZ = p.first;
            int const pdg = p.second;
            int sz = PtsXYZ.size();
            double X[sz], Y[sz], Z[sz];
            for (unsigned j=0; j!=sz;++j)
            {
                X[j] = PtsXYZ.at(j).at(0);
                Y[j] = PtsXYZ.at(j).at(1);
                Z[j] = PtsXYZ.at(j).at(2);
            }
            TGraph* gXY = new TGraph(sz,X,Y);
            gXY->SetMarkerColor(99-(4*i));
            gXY->SetMarkerStyle(6);
            if (pdg==11)
            {
                mgXY->Add(gXY,"P");
            }
            else
            {
                mgXY->Add(gXY,"C");
            }
            
            TGraph* gYZ = new TGraph(sz,Y,Z);
            gYZ->SetMarkerColor(99-(4*i));
            gYZ->SetMarkerStyle(6);
            if (pdg==11)
            {
                mgYZ->Add(gYZ,"P");
            }
            else
            {
                mgYZ->Add(gYZ,"C");
            }
            
            TGraph* gXZ = new TGraph(sz,X,Z);
            gXZ->SetMarkerColor(99-(4*i));
            gXZ->SetMarkerStyle(6);
            if (pdg==11)
            {
                mgXZ->Add(gXZ,"P");
            }
            else
            {
                mgXZ->Add(gXZ,"C");
            }
            std::string label = "Particle " + to_string(i) + "  -> type : " + to_string(pdg);
            legend->AddEntry(gXY,label.c_str(),"P");
            legend->GetEntry()->SetMarkerSize(100);
        }
        
        c1->cd(1);
        mgXY->SetTitle("Z Projection;X(cm);Y(cm)");
        mgXY->Draw("A pmc plc");
        
        c1->cd(2);
        mgYZ->SetTitle("X Projection;Y(cm);Z(cm)");
        mgYZ->Draw("A pmc plc");
        
        c1->cd(3);
        mgXZ->SetTitle("Y Projection;X(cm);Z(cm)");
        mgXZ->Draw("A pmc plc");
        
        c1->cd(4);
        legend->SetNColumns(2);
        legend->Draw();
        
        c1->Draw();
        c1->Upgrade();
    }
}*/

//*************************************************************************************************************

// void plotproj_ShowerReco(std::vector<std::vector<std::vector<std::vector<double>>>> PtsXYZShEv)
// {
//     TCanvas *c1 = new TCanvas("c1","ShowerReco all Events");
//     int NbEv = PtsXYZShEv.size();
//     c1->Divide(3,NbEv);
//     
//     for (int e=0; e!=NbEv; e++)
//     {
//         auto legend = new TLegend();
//         legend->SetHeader("Showers Legend","C");
//         std::vector<std::vector<std::vector<double>>> PtsXYZSh = PtsXYZShEv.at(e);
//         
//         TMultiGraph* mgXY = new TMultiGraph();
//         TMultiGraph* mgYZ = new TMultiGraph();
//         TMultiGraph* mgXZ = new TMultiGraph();
//         
//         for (size_t i = 0, szSh = PtsXYZSh.size(); i != szSh; ++i)
//             {
//                 std::vector<std::vector<double>> PtsXYZ = PtsXYZSh.at(i);
//                 int sz = PtsXYZ.size();
//                 double X[sz], Y[sz], Z[sz];
//                 for (unsigned j=0; j!=sz;++j)
//                 {
//                     X[j] = PtsXYZ.at(j).at(0);
//                     Y[j] = PtsXYZ.at(j).at(1);
//                     Z[j] = PtsXYZ.at(j).at(2);
//                 }
//                 TGraph* gXY = new TGraph(sz,X,Y);
//                 gXY->SetMarkerColor(99-(4*i));
//                 gXY->SetMarkerStyle(6);
//                 mgXY->Add(gXY);
//                 
//                 TGraph* gYZ = new TGraph(sz,Y,Z);
//                 gYZ->SetMarkerColor(99-(4*i));
//                 gYZ->SetMarkerStyle(6);
//                 mgYZ->Add(gYZ);
//                 
//                 TGraph* gXZ = new TGraph(sz,X,Z);
//                 gXZ->SetMarkerColor(99-(4*i));
//                 gXZ->SetMarkerStyle(6);
//                 mgXZ->Add(gXZ);
//                 std::string label = "Shower #" + to_string(i);
//                 legend->AddEntry(gXY,label.c_str(),"P");
//             }
//         
//         c1->cd(1+e);
//         mgXY->SetTitle("Z Projection;X;Y");
//         mgXY->Draw("AP");
//         legend->Draw();
//         
//         c1->cd(2+e);
//         mgYZ->SetTitle("X Projection;Y;Z");
//         mgYZ->Draw("AP");
//         
//         c1->cd(3+e);
//         mgXZ->SetTitle("Y Projection;X;Z");
//         mgXZ->Draw("AP");
//     }
// }


//Prototype
/*
void plotproj_ShowerReco(std::vector<std::vector<std::vector<std::vector<double>>>> PtsXYZShEv, std::string Title = "ShowerReco")
{
    int NbEv = PtsXYZShEv.size();
    for (int e=0; e!=NbEv; e++)
    {
        cout << "Event #" << e << "---------------" << endl;
        std::vector<std::vector<std::vector<double>>> PtsXYZSh = PtsXYZShEv.at(e);
        std::string title = Title + " Event#" + to_string(e);
        plotproj_1e_ShowerReco(PtsXYZSh,title);
        cout << "OK" << endl << endl;
    }
}
*/


void PlotShowerReco()
{

    //filename = "/silver/DUNE/queau/output/root_10elec_1GeV/.root";
//     std::vector<std::vector< pair< std::vector<std::vector<double>> , int const> >> PtsXYZShEv = ListXYZ(filename); //10elec_1GeV_shower_reco.root
    //std::vector< pair< std::vector<std::vector<double>> , int const> > PtsXYZSh = PtsXYZShEv.at(0);
    
    plotproj_1e_ShowerReco(ListXYZ("../output/root_10elec_1GeV/1elec1Gev_C.root").at(0), "Reco cosmics");
    //plotproj_1e_ShowerReco(ListXYZ("../output/root_10elec_1GeV/1elec1Gev_TBHD.root").at(0), "Reco testbeamHD");
    //plotproj_1e_ShowerReco(ListXYZ("../output/root_10elec_1GeV/1elec1Gev_Nu.root").at(0), "Reco neutrino");
}
