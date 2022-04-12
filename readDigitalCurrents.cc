
#include "readDigitalCurrents.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>           // for SubsysReco

//#include "SvtxEvaluator.h"
#include "tpc/TpcDefs.h"

#include <phool/PHCompositeNode.h>

#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>

#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>

#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>

using namespace std;

bool IsOverFrame(double r, double phi);

bool IsOverFrame(double r, double phi){
  //these parameters are taken from Feb 12 drawings of frames.
  double tpc_frame_side_gap=0.8;//mm //space between radial line and start of frame
  double tpc_frame_side_width=2.6;//mm //thickness of frame
  double tpc_margin=0.0;//mm // extra gap between edge of frame and start of GEM holes
  
  double tpc_frame_r3_outer=758.4;//mm inner edge of larger-r frame of r3
  double tpc_frame_r3_inner=583.5;//mm outer edge of smaller-r frame of r3
 
  double tpc_frame_r2_outer=574.9;//mm inner edge of larger-r frame of r2
  double tpc_frame_r2_inner=411.4;//mm outer edge of smaller-r frame of r2
 
  double tpc_frame_r1_outer=402.6;//mm inner edge of larger-r frame of r1
  double tpc_frame_r1_inner=221.0;//mm outer edge of smaller-r frame of r1
 
  //double tpc_sec0_phi=0.0;//get_double_param("tpc_sec0_phi");

  //if the coordinate is in the radial spaces of the frames, return true:
  if (r<tpc_frame_r1_inner+tpc_margin)
    return true;
  if (r>tpc_frame_r1_outer-tpc_margin  && r<tpc_frame_r2_inner+tpc_margin)
    return true;
  if (r>tpc_frame_r2_outer-tpc_margin  && r<tpc_frame_r3_inner+tpc_margin)
    return true;
  if (r>tpc_frame_r3_outer-tpc_margin)
    return true;

  //if the coordinate is within gap+width of a sector boundary, return true:
  //note that this is not a line of constant radius, but a linear distance from a radius.

  //find the two spokes we're between:
  double pi = 2 * acos(0.0);

  float sectorangle=(pi/6);
  float nsectors=phi/sectorangle;
  int nsec=floor(nsectors);
  float reduced_phi=phi-nsec*sectorangle; //between zero and sixty degrees.
  float dist_to_previous=r*sin(reduced_phi);
  float dist_to_next=r*sin(sectorangle-reduced_phi);
  if (dist_to_previous<tpc_frame_side_gap+tpc_frame_side_width+tpc_margin)
    return true;
  if (dist_to_next<tpc_frame_side_gap+tpc_frame_side_width+tpc_margin)
    return true;
  
  return false;
}

//____________________________________________________________________________..
readDigitalCurrents::readDigitalCurrents(const std::string &name, const std::string &filename):
 SubsysReco(name)
 , hm(nullptr)
 , _filename(filename)
 ,_ampIBFfrac(0.02)
 ,_collSyst(0)
// , outfile(nullptr)
 {
  std::cout << "readDigitalCurrents::readDigitalCurrents(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
readDigitalCurrents::~readDigitalCurrents()
{
  std::cout << "readDigitalCurrents::~readDigitalCurrents() Calling dtor" << std::endl;
  delete hm;

}

//____________________________________________________________________________..
int readDigitalCurrents::Init(PHCompositeNode *topNode)
{
  std::cout << "readDigitalCurrents::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  int nz=72;
  double z_rdo=108*cm;

  int nr=159;
  //const int nphi=128*3;
  //const int nz=62*2;
  //double z_rdo=105.5*cm;
  //double rmin=20*cm;
  double rmax=78*cm;

  hm = new Fun4AllHistoManager("HITHIST");
  const int r_bins_N = 51;
  double r_bins[r_bins_N+1] = {217.83, 
                              311.05,317.92,323.31,329.27,334.63,340.59,345.95,351.91,357.27,363.23,368.59,374.55,379.91,385.87,391.23,397.19,402.49,
                              411.53,421.70,431.90,442.11,452.32,462.52,472.73,482.94,493.14,503.35,513.56,523.76,533.97,544.18,554.39,564.59,574.76,
                              583.67,594.59,605.57,616.54,627.51,638.48,649.45,660.42,671.39,682.36,693.33,704.30,715.27,726.24,737.21,748.18,759.11};

  const int nphi=205;
  double phi_bins[nphi+1] = {0.,0.0068, 0.038675, 0.07055, 0.102425, 0.1343, 0.166175, 0.19805, 
                            0.229925, 0.2618, 0.293675, 0.32555, 0.357425, 0.3893, 0.421175, 0.45305, 0.484925, 
                            0.5168, 0.5304, 0.562275, 0.59415, 0.626025, 0.6579, 0.689775, 0.72165, 0.753525, 0.7854, 
                            0.817275, 0.84915, 0.881025, 0.9129, 0.944775, 0.97665, 1.008525, 1.0404, 1.054, 1.085875, 
                            1.11775, 1.149625, 1.1815, 1.213375, 1.24525, 1.277125, 1.309, 1.340875, 1.37275, 1.404625, 1.4365, 
                            1.468375, 1.50025, 1.532125, 1.564, 1.5776, 1.609475, 1.64135, 1.673225, 1.7051, 1.736975, 1.76885, 
                            1.800725, 1.8326, 1.864475, 1.89635, 1.928225, 1.9601, 1.991975, 2.02385, 2.055725, 2.0876, 2.1012, 
                            2.133075, 2.16495, 2.196825, 2.2287, 2.260575, 2.29245, 2.324325, 2.3562, 2.388075, 2.41995, 2.451825, 
                            2.4837, 2.515575, 2.54745, 2.579325, 2.6112, 2.6248, 2.656675, 2.68855, 2.720425, 2.7523, 2.784175, 2.81605, 
                            2.847925, 2.8798, 2.911675, 2.94355, 2.975425, 3.0073, 3.039175, 3.07105, 3.102925, 3.1348, 3.1484, 3.180275, 
                            3.21215, 3.244025, 3.2759, 3.307775, 3.33965, 3.371525, 3.4034, 3.435275, 3.46715, 3.499025, 3.5309, 3.562775, 
                            3.59465, 3.626525, 3.6584, 3.672, 3.703875, 3.73575, 3.767625, 3.7995, 3.831375, 3.86325, 3.895125, 3.927, 3.958875, 
                            3.99075, 4.022625, 4.0545, 4.086375, 4.11825, 4.150125, 4.182, 4.1956, 4.227475, 4.25935, 4.291225, 4.3231, 4.354975, 
                            4.38685, 4.418725, 4.4506, 4.482475, 4.51435, 4.546225, 4.5781, 4.609975, 4.64185, 4.673725, 4.7056, 4.7192, 4.751075, 
                            4.78295, 4.814825, 4.8467, 4.878575, 4.91045, 4.942325, 4.9742, 5.006075, 5.03795, 5.069825, 5.1017, 5.133575, 5.16545, 
                            5.197325, 5.2292, 5.2428, 5.274675, 5.30655, 5.338425, 5.3703, 5.402175, 5.43405, 5.465925, 5.4978, 5.529675, 5.56155, 
                            5.593425, 5.6253, 5.657175, 5.68905, 5.720925, 5.7528, 5.7664, 5.798275, 5.83015, 5.862025, 5.8939, 5.925775, 5.95765, 
                            5.989525, 6.0214, 6.053275, 6.08515, 6.117025, 6.1489, 6.180775, 6.21265, 6.244525, 6.2764,2*pi};

  
  //double phi_bins[nphi+1];
  //for (int p=0;p<=nphi;p++){
  //  phi_bins[p]=6.28319/nphi*p;
  //} 
  double z_bins[2*nz+1];
  for (int z=0;z<=2*nz;z++){
    z_bins[z]=-z_rdo+z_rdo/nz*z;
  } 

  _h_R  = new TH1F("_h_R" ,"_h_R;R, [m]"   ,r_bins_N ,r_bins);
  _h_hits  = new TH1F("_h_hits" ,"_h_hits;N, [hit]"   ,1e5,0-0.5,1e5-0.5);
  //_h_SC_ibf  = new TH3F("_h_SC_ibf" ,"_h_SC_ibf;#phi, [rad];R, [m];Z, [m]"   ,nphi,0,6.28319,nr,rmin,rmax,2*nz,-z_rdo,z_rdo);
  //_h_DC_SC = new TH3F("_h_DC_SC" ,"_h_DC_SC;#phi, [rad];R, [m];Z, [m]"   ,nphi,0,6.28319,nr,rmin,rmax,2*nz,-z_rdo,z_rdo);
  _h_hit_XY = new TH2F("_h_hit_XY" ,"_h_hit_XY;X, [m];Y, [m]"   ,4*nr,-1*rmax,rmax,4*nr,-1*rmax,rmax);
  //_h_DC_E = new TH2F("_h_DC_E" ,"_h_DC_E;ADC;E"   ,200,-100,2e3-100,500,-100,5e3-100);
  
  //double phi_bins[nphi+1];
  //for (int p=0;p<=nphi;p++){
  //  phi_bins[p]=6.28319/nphi*p;
  //} 
  //double z_bins[2*nz+1];
  //for (int z=0;z<=2*nz;z++){
  //  z_bins[z]=-z_rdo+z_rdo/nz*z;
  //} 
  _h_SC_ibf  = new TH3F("_h_SC_ibf" ,"_h_SC_ibf;#phi, [rad];R, [m];Z, [m]"   ,nphi,phi_bins,r_bins_N ,r_bins,2*nz,z_bins);
  _h_DC_SC = new TH3F("_h_DC_SC" ,"_h_DC_SC;#phi, [rad];R, [m];Z, [m]"   ,nphi,phi_bins,r_bins_N ,r_bins,2*nz,z_bins);
  _h_DC_SC_XY = new TH3F("_h_DC_SC_XY" ,"_h_DC_SC_XY;X, [m];Y, [m];Z, [m]"   ,4*nr,-1*rmax,rmax,4*nr,-1*rmax,rmax,2*nz,-z_rdo,z_rdo);
  _h_DC_E = new TH2F("_h_DC_E" ,"_h_DC_E;ADC;E"   ,200,-100,2e3-100,500,-100,5e3-100);
  hm->registerHisto(_h_R );
  hm->registerHisto(_h_hits );
  hm->registerHisto(_h_DC_SC );
  hm->registerHisto(_h_DC_SC_XY );
  hm->registerHisto(_h_hit_XY );
  hm->registerHisto(_h_DC_E );
  hm->registerHisto(_h_SC_ibf );
  //outfile = new TFile(_filename.c_str(), "RECREATE");
  _event_timestamp = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::InitRun(PHCompositeNode *topNode)
{
  std::cout << "readDigitalCurrents::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  std::string line;
  //AA collisions timestamps
  std::string txt_file = "/sphenix/user/shulga/Work/IBF/DistortionMap/timestamps_50kHz.txt";
  int start_line = 3;
  if(_collSyst==1){
    //pp collisions timestamps
    txt_file = "/phenix/u/hpereira/sphenix/work/g4simulations/timestamps_3MHz.txt";
    //txt_file = "/sphenix/user/shulga/Work/IBF/DistortionMap/timestamps_50kHz.txt";
    start_line = 2;
  }
  ifstream InputFile (txt_file);
  if (InputFile.is_open()){
    int n_line=0;
    while ( getline (InputFile,line) )
    {
      n_line++;
      if(n_line>start_line){
        std::istringstream is( line );
        double n[2] = {0,0};
        int i = 0;
        while( is >> n[i] ) {    
            i++;    
        }
        _timestamps[n[0]]=n[1];
        if(n_line<10){
          cout<<n[1]<<endl;
        }
        _keys.push_back(int(n[0]));
      }
    }
    InputFile.close();
  }

  else cout << "Unable to open file:"<<txt_file<<endl; 

  TFile *MapsFile; 
  //if(_fUseIBFMap){
    MapsFile = new TFile("/sphenix/user/shulga/Work/IBF/DistortionMap/IBF_Map.root","READ");
    if ( MapsFile->IsOpen() ) printf("Gain/IBF Maps File opened successfully\n");
    //_h_modules_anode       = (TH2F*)MapsFile ->Get("h_modules_anode")      ->Clone("_h_modules_anode");
    _h_modules_measuredibf = (TH2F*)MapsFile ->Get("h_modules_measuredibf")->Clone("_h_modules_measuredibf");
  //}
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::process_event(PHCompositeNode *topNode)
{
  double bX = _beamxing;
  //float bX = 1508071;

  //double z_bias_avg = 0;
  //if (_fAvg==1){ 
  //  z_bias_avg=1.05*(float) rand()/RAND_MAX;
  //}
  int bemxingsInFile = _keys.size();
  if (_evtstart>= bemxingsInFile) _evtstart=_evtstart-bemxingsInFile;
  int key = _keys.at(_evtstart);
  _event_timestamp = (float)_timestamps[key]*ns;//units in seconds
  _event_bunchXing = key;
  if(_evtstart%100==0) cout<<"_evtstart = "<<_evtstart<<endl;
  _evtstart++;

  //std::cout << "readDigitalCurrents::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  ostringstream nodename;
  set<std::string>::const_iterator iter;
  // //nodename << "G4HIT_TPC";
  nodename << "TRKR_HITSET";

  // //  SvtxEvaluator
  // SvtxEvaluator *hits = findNode::getClass<SvtxEvaluator>(topNode, nodename.str().c_str());
  // //int n_hits = 0;
  // if (hits){
  //   PHG4HitContainer::ConstRange hit_range = hits->getHits();
  // }
  //===================================
  // get node containing the digitized hits
  TrkrHitSetContainer* _hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, nodename.str().c_str());
  if (!_hitmap)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  ostringstream geo_nodename;
  geo_nodename << "CYLINDERCELLGEOM_SVTX";

  PHG4CylinderCellGeomContainer* _geom_container_ccgc = nullptr; 
  PHG4TpcCylinderGeomContainer* _geom_container_cgc = nullptr;
  if(_f_ccgc==1){
    _geom_container_ccgc = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, geo_nodename.str().c_str());
    if (!_geom_container_ccgc)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }else{
    _geom_container_cgc = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, geo_nodename.str().c_str());

    if (!_geom_container_cgc)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }


  // loop over all the hits
  // hits are stored in hitsets, so have to get the hitset first
  int n_hits = 0;
  double pi = 2 * acos(0.0);
  //float _event_bunchXing = 1508071;
  TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;iter != all_hitsets.second; ++iter){    
    //checking that the object is inside TPC
    if(TrkrDefs::getTrkrId(iter->first) == TrkrDefs::tpcId){
      TrkrHitSet::ConstRange range = iter->second->getHits();
      unsigned int layer = TrkrDefs::getLayer(iter->first);
      PHG4CylinderCellGeom *layergeom_ccgc = nullptr;
      PHG4TpcCylinderGeom *layergeom_cgc = nullptr;
      double radius = 0;
      if(_f_ccgc==1){
        layergeom_ccgc = _geom_container_ccgc->GetLayerCellGeom(layer);
        radius = layergeom_ccgc->get_radius()*cm;
      }else{
        layergeom_cgc = _geom_container_cgc->GetLayerCellGeom(layer);
        radius = layergeom_cgc->get_radius()*cm;
      }
      //PHG4TpcCylinderGeom *layergeom = _geom_container->GetLayerCellGeom(layer);
      //double radius = layergeom->get_radius()*cm;  // returns center of the layer    
      for(TrkrHitSet::ConstIterator hit_iter = range.first; hit_iter != range.second; ++hit_iter){
        int f_fill_ibf=0;

        //TrkrDefs::hitkey hit_key = hit_iter->first;
        unsigned short phibin = TpcDefs::getPad(hit_iter->first);
        unsigned short zbin = TpcDefs::getTBin(hit_iter->first); 
        double phi_center = 0;
        if(_f_ccgc==1){
          phi_center = layergeom_ccgc->get_phicenter(phibin);
        }else{
          phi_center = layergeom_cgc->get_phicenter(phibin);
        }
        if (phi_center<0) phi_center+=2*pi;

        float x = radius*cos(phi_center);
        float y = radius*sin(phi_center);
        
        _h_hit_XY->Fill( x, y);
        double z = 0;// layergeom->get_zcenter(zbin)*cm;  
        if(_f_ccgc==1){
          z = layergeom_ccgc->get_zcenter(zbin)*cm;
        }else{
          z = layergeom_cgc->get_zcenter(zbin)*cm;
        }
        TrkrHit *hit = hit_iter->second;
        unsigned short adc = hit->getAdc()-adc_pedestal;
        float E = hit->getEnergy();
        //double z = 0;
        //double z_prim = -1*1e10;
        double z_ibf =  -1*1e10;

        //if(!IsOverFrame(radius/mm,phi_center)){
          

          if(z>=0 && z<1.055*m){
            if(adc>=0)n_hits++;
            if(adc>=0)_h_DC_E->Fill(adc,E);

            //z_prim = z-(bX-_event_bunchXing)*106*vIon*ns;
            z_ibf = 1.055*m-(bX-_event_bunchXing)*106*vIon*ns;
            //if(n_hits%100==0)cout<<"z_ibf = "<<z_ibf<<"="<<"1.055*m-("<<bX<<"-"<<_event_bunchXing<<")*"<<106*vIon*ns<<endl;
            if( z_ibf>0 && z_ibf<1.055*m){
              f_fill_ibf=1;
            }
          }
          if(z<0 && z>-1.055*m){
            if(adc>=0)n_hits++;
            if(adc>=0)_h_DC_E->Fill(adc,E);

            //z_prim = z+(bX-_event_bunchXing)*106*vIon*ns;
            z_ibf = -1.055*m+(bX-_event_bunchXing)*106*vIon*ns;
            if( z_ibf<0 && z_ibf>-1.055*m){
              f_fill_ibf=1;
            }
          }        

        //Reading IBF and Gain weights according to X-Y position
        float w_ibf = 1.;
        //float w_gain = 1.;
        //if(_fUseIBFMap){
          int bin_x = _h_modules_measuredibf ->GetXaxis()->FindBin(x/mm);
          int bin_y = _h_modules_measuredibf ->GetYaxis()->FindBin(y/mm);
          w_ibf = _h_modules_measuredibf->GetBinContent(bin_x,bin_y);
          //w_gain = _h_modules_anode->GetBinContent(bin_x,bin_y);
          w_ibf = 1.;
        //}
          float w_adc = adc*w_ibf;
          _h_DC_SC->Fill(phi_center,radius,z,w_adc);        
          _h_DC_SC_XY->Fill(x,y,z,w_adc);
          if(f_fill_ibf==1){
            _h_SC_ibf  ->Fill(phi_center,radius,z_ibf,w_adc);
            _h_R->Fill(radius);
          }
        //}
        //if(n_hits%100==0) std::cout<<radius<<"|"<<phi_center<<"|"<<z<<std::endl;
      }
    }
  } 
    
  //TrkrHitSetContainer, TrkrHitSet, TrkrHit, and TrkrDefs objects used above in offline/packages/trackbase.
  _h_hits->Fill(n_hits);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::ResetEvent(PHCompositeNode *topNode)
{
  //std::cout << "readDigitalCurrents::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::EndRun(const int runnumber)
{
  std::cout << "readDigitalCurrents::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::End(PHCompositeNode *topNode)
{
  std::cout << "readDigitalCurrents::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  _h_R     ->Sumw2( false );
  _h_hits     ->Sumw2( false );
  _h_DC_E     ->Sumw2( false );
  _h_DC_SC    ->Sumw2( false );
  _h_hit_XY ->Sumw2( false );
  _h_DC_SC_XY ->Sumw2( false );
  _h_SC_ibf   ->Sumw2( false );
  hm->dumpHistos(_filename, "RECREATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::Reset(PHCompositeNode *topNode)
{
 std::cout << "readDigitalCurrents::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void readDigitalCurrents::Print(const std::string &what) const
{
  std::cout << "readDigitalCurrents::Print(const std::string &what) const Printing info for " << what << std::endl;
}

void readDigitalCurrents::SetEvtStart(int newEvtStart){
  _evtstart = newEvtStart;
  cout<<"Start event is set to: "<<newEvtStart<<endl;

}
void readDigitalCurrents::SetBeamXing(int newBeamXing){
  _beamxing = newBeamXing;
  cout<<"Initial BeamXing is set to: "<<newBeamXing<<endl;

}
void readDigitalCurrents::SetCollSyst(int coll_syst){
  _collSyst = coll_syst;
  std::string s_syst[2] = {"AA","pp"};
  cout<<"Collision system is set to: "<<s_syst[_collSyst]<<endl;

}

void readDigitalCurrents::SetIBF(float ampIBFfrac){
  _ampIBFfrac = ampIBFfrac;
  cout<<"IBF is set to: "<<_ampIBFfrac<<endl;
}

void readDigitalCurrents::SetCCGC(float f_ccgc){
  _f_ccgc = f_ccgc;
  cout<<"IBF is set to: "<<_f_ccgc<<endl;
}


