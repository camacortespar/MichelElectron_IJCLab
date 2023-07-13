bool Debug = false;

double degTOrad = TMath::Pi()/180.;

void SetDebug(int DEBUG){Debug=DEBUG;}

double GetMax(TH1D* hist){
  int binmax = hist->GetMaximumBin();
  return hist->GetBinContent(binmax);
/*
  TF1* fitmax = new TF1("fitmax","[0]*TMath::Gaus([1],[2])",hist->GetBinCenter(binmax-2),hist->GetBinCenter(binmax+2));
  fitmax->SetParameter(0,hist->GetBinContent(binmax));
  fitmax->SetParameter(1,hist->GetBinCenter(binmax));
  fitmax->SetParameter(2,1);

  hist->Fit(fitmax,"QR");

  return fitmax->GetParameter(0);
*/
}

double GetWidth(TH1D* hist, double thrs = 0.5){
  int binmax = hist->GetMaximumBin();
  double max = GetMax(hist); 

  double lowamp1 = 1, highamp1 = 1;
  double lowamp2 = 1, highamp2 = 1;
  int    lowbin1 = 0, highbin1 = 0;
  int    lowbin2 = 0, highbin2 = 0;

  for(int bin = binmax; bin > 0; bin--){
    if(hist->GetBinContent(bin) / max >= thrs) lowbin2 = bin;
    if(hist->GetBinContent(bin) / max < thrs) {lowbin1 = bin; break;}
  }
  for(int bin = binmax; bin < hist->GetNbinsX(); bin++){
    if(hist->GetBinContent(bin) / max >= thrs)  highbin1 = bin-1;
    if(hist->GetBinContent(bin) / max <  thrs) {highbin2 = bin-1; break;}
  }

  double lowedge = hist->GetBinCenter(lowbin1) + (hist->GetBinCenter(lowbin2)-hist->GetBinCenter(lowbin1))/(hist->GetBinContent(lowbin2)-hist->GetBinContent(lowbin1))*(thrs*max-hist->GetBinContent(lowbin1));

  double highedge = hist->GetBinCenter(highbin1) + (hist->GetBinCenter(highbin2)-hist->GetBinCenter(highbin1))/(hist->GetBinContent(highbin2)-hist->GetBinContent(highbin1))*(thrs*max-hist->GetBinContent(highbin1));

  return highedge-lowedge;
}

void GetWidthEdges(TH1D* hist, double &lowedge, double &highedge, double thrs = 0.5){
  int binmax = hist->GetMaximumBin();
  double max = GetMax(hist); 

  double lowamp1 = 1, highamp1 = 1;
  double lowamp2 = 1, highamp2 = 1;
  int    lowbin1 = 0, highbin1 = 0;
  int    lowbin2 = 0, highbin2 = 0;

  for(int bin = binmax; bin > 0; bin--){
    if(hist->GetBinContent(bin) / max >= thrs) lowbin2 = bin;
    if(hist->GetBinContent(bin) / max < thrs) {lowbin1 = bin; break;}
  }
  for(int bin = binmax; bin < hist->GetNbinsX(); bin++){
    if(hist->GetBinContent(bin) / max >= thrs)  highbin1 = bin;
    if(hist->GetBinContent(bin) / max <  thrs) {highbin2 = bin; break;}
  }

  lowedge = hist->GetBinCenter(lowbin1) + (hist->GetBinCenter(lowbin2)-hist->GetBinCenter(lowbin1))/(hist->GetBinContent(lowbin2)-hist->GetBinContent(lowbin1))*(thrs*max-hist->GetBinContent(lowbin1));

  highedge = hist->GetBinCenter(highbin1) + (hist->GetBinCenter(highbin2)-hist->GetBinCenter(highbin1))/(hist->GetBinContent(highbin2)-hist->GetBinContent(highbin1))*(thrs*max-hist->GetBinContent(highbin1));

  return 0;
}

void GetAngles(const recob::Track& track, double &theta, double &phi){

  double x_fst=-1, y_fst=-1, z_fst=-1;
  double x_lst=-1, y_lst=-1, z_lst=-1;

  for(size_t hit=0; hit<track.NPoints(); ++hit){
    if(track.LocationAtPoint(hit).X() == -999) continue;
    if(x_fst == -1){
      x_fst = track.LocationAtPoint(hit).X();
      y_fst = track.LocationAtPoint(hit).Y();
      z_fst = track.LocationAtPoint(hit).Z();
    }
    x_lst = track.LocationAtPoint(hit).X();
    y_lst = track.LocationAtPoint(hit).Y();
    z_lst = track.LocationAtPoint(hit).Z();
  }

  phi   = atan((y_lst-y_fst)/(z_lst-z_fst)) / degTOrad;
  theta = atan((x_fst-x_lst)/sqrt(pow(y_lst-y_fst,2)+pow(z_lst-z_fst,2))) / degTOrad;
}

bool HitIsValid(const art::Ptr<recob::Hit> hit,
                const recob::TrackHitMeta* thm,
                const recob::Track& track)
{ 
  if(!thm)                                     return false;
  if (thm->Index() == 2147483647)              return false;
  if (!track.HasValidPoint(thm->Index()))      return false;
  return true;
}


std::vector<std::vector<unsigned>> OrganizeHitsSnippets(
      const std::vector<art::Ptr<recob::Hit>>& hits,
      const std::vector<const recob::TrackHitMeta*>& thms,
      const recob::Track& track,
      unsigned nplanes);

std::vector<std::vector<unsigned>> OrganizeHitsSnippets(
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const std::vector<const recob::TrackHitMeta*>& thms,
  const recob::Track& track,
  unsigned nplanes)
{
  // In this case, we need to only accept one hit in each snippet
  // Snippets are counted by the Start, End, and Wire. If all these are the same for a hit, then they are on the same snippet.
  //
  // If there are multiple valid hits on the same snippet, we need a way to pick the best one.
  // (TODO: find a good way). The current method is to take the one with the highest charge integral.
  struct HitIdentifier {
    int startTick;
    int endTick;
    int wire;
    float integral;

    // construct
    explicit HitIdentifier(const recob::Hit& hit)
      : startTick(hit.StartTick())
      , endTick(hit.EndTick())
      , wire(hit.WireID().Wire)
      , integral(hit.Integral())
    {}

    // Defines whether two hits are on the same snippet
    inline bool operator==(const HitIdentifier& rhs) const
    {
      return startTick == rhs.startTick && endTick == rhs.endTick && wire == rhs.wire;
    }

    // Defines which hit to pick between two both on the same snippet
    inline bool operator>(const HitIdentifier& rhs) const { return integral > rhs.integral; }
  };

  std::vector<std::vector<unsigned>> ret(nplanes);
  std::vector<std::vector<HitIdentifier>> hit_idents(nplanes);
  for (unsigned i = 0; i < hits.size(); i++) {
    if (HitIsValid(hits[i], thms[i], track)) {
      HitIdentifier this_ident(*hits[i]);

      // check if we have found a hit on this snippet before
      bool found_snippet = false;
      for (unsigned j = 0; j < ret[hits[i]->WireID().Plane].size(); j++) {
        if (this_ident == hit_idents[hits[i]->WireID().Plane][j]) {
          found_snippet = true;
          if (this_ident > hit_idents[hits[i]->WireID().Plane][j]) {
            ret[hits[i]->WireID().Plane][j] = i;
            hit_idents[hits[i]->WireID().Plane][j] = this_ident;
          }
          break;
        }
      }
      if (!found_snippet) {
        ret[hits[i]->WireID().Plane].push_back(i);
        hit_idents[hits[i]->WireID().Plane].push_back(this_ident);
      }
    }
  }
  return ret;
}

void WriteEventList(string filename, vector<string> file_list, vector<int> ev_list, vector<int> trk_list){
  if(Debug) std::cout << "    WriteEventList - start " << endl;

  ofstream File(filename.c_str());
  if(File){ 
    for(int i=0;i<file_list.size();i++) File << file_list[i] << " " << ev_list[i] << " " << trk_list[i] << "\n";
    File.close();
  }
  else cerr << "    ReadFileList - " << filename << " not found. Skipped." << endl;

  if(Debug) std::cout << "    WriteEventList - finished " << endl;
}

vector<string> ReadFileList(int nfiles, string FileName)
{
  if(Debug) std::cout << "    ReadFileList - start " << endl;
  vector<string> vec;
  string file;
  string filename = FileName;
  ifstream File(filename.c_str());

  if(File){
    cout << "    ReadFileList - " << filename << " opened"<< endl;
    while(File.good()){
      File >> file;
      vec.push_back(file);
      if(nfiles>0 && vec.size() == nfiles) break;
    }
  }
  else{
    cerr << "    ReadFileList - " << filename << " not found. Exit." << endl;
    exit(1);
  }

  File.close();

  if(nfiles == 0) vec.pop_back();
  cout << "    ReadFileList - " << vec.size() << " files found in " << FileName << endl;
  cout << "    finished" << endl;

  return vec;
}

vector<string> ReadFileEventList(int nfiles, string FileName)
{
  if(Debug) std::cout << "    ReadFileList - start " << endl;
  vector<string> vec;
  string file, evt, trk;
  string filename = FileName;
  ifstream File(filename.c_str());

  if(File){
    cout << "    ReadFileList - " << filename << " opened"<< endl;
    while(File.good()){
      File >> file >> evt >> trk;
      vec.push_back(file);
      if(nfiles>0 && vec.size() == nfiles) break;
    }
  }
  else{
    cerr << "    ReadFileList - " << filename << " not found. Exit." << endl;
    exit(1);
  }

  File.close();

  if(nfiles == 0) vec.pop_back();
  cout << "    ReadFileList - " << vec.size() << " files found in " << FileName << endl;
  cout << "    finished" << endl;

  return vec;
}

void ReadEventList(int nevents, string FileName, vector<string> &file_list, vector<int> &ev_list, vector<int> &trk_list)
{
  if(Debug) std::cout << "    ReadEventList - start " << endl;
  vector<int> vec;
  string file, evt, trk;
  string filename = FileName;
  ifstream File(filename.c_str());

  if(File){
    cout << "    ReadEventList - " << filename << " opened"<< endl;
    while(File.good()){
      File >> file >> evt >> trk;
      file_list.push_back(file);
      ev_list.push_back(atoi(evt.c_str()));
      trk_list.push_back(atoi(trk.c_str()));
      if(nevents>0 && ev_list.size() == nevents) break;
    }
  }
  else{
    cerr << "    ReadEventList - " << filename << " not found. Exit." << endl;
    exit(1);
  }

  File.close();

  if(nevents == 0){ file_list.pop_back(); ev_list.pop_back(); trk_list.pop_back();}
  cout << "    ReadEventList - " << ev_list.size() << " events found in " << FileName << endl;
  cout << "    finished" << endl;

  return;
}


void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;

   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;

   for (Int_t i=0;i<Nx;i++) {

      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }

      for (Int_t j=0;j<Ny;j++) {

         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }

         C->cd(0);

         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);

         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);

         pad->Draw();
      }
   }
}
