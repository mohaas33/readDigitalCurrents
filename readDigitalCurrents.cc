//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in readDigitalCurrents.h.
//
// readDigitalCurrents(const std::string &name = "readDigitalCurrents")
// everything is keyed to readDigitalCurrents, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// readDigitalCurrents::~readDigitalCurrents()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int readDigitalCurrents::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int readDigitalCurrents::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int readDigitalCurrents::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int readDigitalCurrents::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int readDigitalCurrents::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int readDigitalCurrents::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int readDigitalCurrents::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void readDigitalCurrents::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "readDigitalCurrents.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>           // for SubsysReco

//#include "SvtxEvaluator.h"
#include "tpc/TpcDefs.h"

#include <phool/PHCompositeNode.h>

#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

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
 
  double tpc_frame_r2_outer=574.9;//mm inner edge of larger-r frame of r3
  double tpc_frame_r2_inner=411.4;//mm outer edge of smaller-r frame of r3
 
  double tpc_frame_r1_outer=402.6;//mm inner edge of larger-r frame of r3
  double tpc_frame_r1_inner=221.0;//mm outer edge of smaller-r frame of r3
 
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

  int nr=159;
  int nphi=360;
  int nz=62*2;
  double z_rdo=105.5*cm;
  double rmin=20*cm;
  double rmax=78*cm;
  cout << "CalculateDistortions::Init(PHCompositeNode *topNode) Initializing" << endl;
  hm = new Fun4AllHistoManager("HITHIST");

  _h_hits  = new TH1F("_h_hits" ,"_h_hits;N, [hit]"   ,4000,0,1e6);
  _h_SC_ibf  = new TH3F("_h_SC_ibf" ,"_h_SC_ibf;#phi, [rad];R, [m];Z, [m]"   ,nphi,0,6.28319,nr,rmin,rmax,2*nz,-z_rdo,z_rdo);
  _h_DC_SC = new TH3F("_h_DC_SC" ,"_h_DC_SC;#phi, [rad];R, [m];Z, [m]"   ,nphi,0,6.28319,nr,rmin,rmax,2*nz,-z_rdo,z_rdo);
  _h_DC_SC_XY = new TH3F("_h_DC_SC_XY" ,"_h_DC_SC_XY;X, [m];Y, [m];Z, [m]"   ,4*nr,-1*rmax,rmax,4*nr,-1*rmax,rmax,2*nz,-z_rdo,z_rdo);
  _h_DC_E = new TH2F("_h_DC_E" ,"_h_DC_E;ADC;E"   ,200,-100,2e3-100,500,-100,5e3-100);
  hm->registerHisto(_h_hits );
  hm->registerHisto(_h_DC_SC );
  hm->registerHisto(_h_DC_SC_XY );
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

  PHG4CylinderCellGeomContainer* _geom_container =
        findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, geo_nodename.str().c_str());

  if (!_geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // loop over all the hits
  // hits are stored in hitsets, so have to get the hitset first
  int n_hits = 0;
  double pi = 2 * acos(0.0);
  //float _event_bunchXing = 1508071;
  TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;iter != all_hitsets.second; ++iter){    
    unsigned int layer = TrkrDefs::getLayer(iter->first);
    
    if(TrkrDefs::getTrkrId(iter->first) == TrkrDefs::tpcId){
      PHG4CylinderCellGeom *layergeom = _geom_container->GetLayerCellGeom(layer);
      double radius = layergeom->get_radius()*cm;  // returns center of the layer
      TrkrHitSet::ConstRange range = iter->second->getHits();
      for(TrkrHitSet::ConstIterator hit_iter = range.first; hit_iter != range.second; ++hit_iter){
        n_hits++;
        int f_fill_ibf=1;

        //TrkrDefs::hitkey hit_key = hit_iter->first;
        unsigned short phibin = TpcDefs::getPad(hit_iter->first);
        unsigned short zbin = TpcDefs::getTBin(hit_iter->first); 
        double phi_center = layergeom->get_phicenter(phibin);
        if (phi_center<0) phi_center+=2*pi;

        float x = radius*cos(phi_center);
        float y = radius*sin(phi_center);

        double z = layergeom->get_zcenter(zbin)*cm;  
        TrkrHit *hit = hit_iter->second;
        unsigned short adc = hit->getAdc();
        float E = hit->getEnergy();
        //double z = 0;
        //double z_prim = -1*1e10;
        double z_ibf =  -1*1e10;

        if(!IsOverFrame(radius/mm,phi_center)){
          if(z>=0){
            //z_prim = z-(bX-_event_bunchXing)*106*vIon*ns;
            z_ibf = 1.05-(bX-_event_bunchXing)*106*vIon*ns;
            //if(n_hits%100==0)cout<<"z_ibf = "<<z_ibf<<"="<<"1.05-("<<bX<<"-"<<_event_bunchXing<<")*"<<106*vIon*ns<<endl;
            if( z_ibf<=0){
              f_fill_ibf=0;
            }
          }
          if(z<0){
            //z_prim = z+(bX-_event_bunchXing)*106*vIon*ns;
            z_ibf = -1.05+(bX-_event_bunchXing)*106*vIon*ns;
            if( z_ibf>=0){
              f_fill_ibf=0;
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
        //}
          float w_adc = adc*w_ibf;
          _h_DC_E->Fill(adc,E);
          _h_DC_SC->Fill(phi_center,radius,z,w_adc);        
          _h_DC_SC_XY->Fill(x,y,z,w_adc);
          if(f_fill_ibf==1)_h_SC_ibf  ->Fill(phi_center,radius,z_ibf,w_adc);
        }
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
  _h_hits     ->Sumw2( false );
  _h_DC_E     ->Sumw2( false );
  _h_DC_SC    ->Sumw2( false );
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