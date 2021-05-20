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
#include <string>
#include <sstream>
#include <set>

using namespace std;

//____________________________________________________________________________..
readDigitalCurrents::readDigitalCurrents(const std::string &name, const std::string &filename):
 SubsysReco(name)
 , hm(nullptr)
 , _filename(filename)
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
  double cm=1e-2; //changed to make 'm' 1.0, for convenience.

  int nr=159;
  int nphi=360;
  int nz=62*2;
  double z_rdo=105.5*cm;
  double rmin=20*cm;
  double rmax=78*cm;
  cout << "CalculateDistortions::Init(PHCompositeNode *topNode) Initializing" << endl;
  hm = new Fun4AllHistoManager("HITHIST");

  _h_hits  = new TH1F("_h_hits" ,"_h_hits;N, [hit]"   ,4000,0,1e6);
  _h_DC_SC = new TH3F("_h_DC_SC" ,"_h_DC_SC;#phi, [rad];R, [m];Z, [m]"   ,nphi,0,6.28319,nr,rmin,rmax,2*nz,-z_rdo,z_rdo);
  _h_DC_E = new TH2F("_h_DC_E" ,"_h_DC_E;ADC;E"   ,200,-100,2e3-100,500,-100,5e3-100);
  hm->registerHisto(_h_hits );
  hm->registerHisto(_h_DC_SC );
  hm->registerHisto(_h_DC_E );
  //outfile = new TFile(_filename.c_str(), "RECREATE");
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::InitRun(PHCompositeNode *topNode)
{
  std::cout << "readDigitalCurrents::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::process_event(PHCompositeNode *topNode)
{
  std::cout << "readDigitalCurrents::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
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

  TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
       iter != all_hitsets.second;
       ++iter)
  {    
    unsigned int layer = TrkrDefs::getLayer(iter->first);
    PHG4CylinderCellGeom *layergeom = _geom_container->GetLayerCellGeom(layer);
    double radius = layergeom->get_radius();  // returns center of the layer
    
    TrkrHitSet::ConstRange range = iter->second->getHits();
    for(TrkrHitSet::ConstIterator hit_iter = range.first; hit_iter != range.second; ++hit_iter)
      {
        n_hits++;

       //TrkrDefs::hitkey hit_key = hit_iter->first;
        unsigned short phibin = TpcDefs::getPad(hit_iter->first);
        unsigned short zbin = TpcDefs::getTBin(hit_iter->first); 


        double phi_center = layergeom->get_phicenter(phibin);
        double z = layergeom->get_zcenter(zbin);  

        TrkrHit *hit = hit_iter->second;


       unsigned short adc = hit->getAdc();
       float E = hit->getEnergy();
       //double z = 0;

       _h_DC_E->Fill(adc,E);
       _h_DC_SC->Fill(radius,phi_center,z,adc);
       if(n_hits%100==0) std::cout<<radius<<"|"<<phi_center<<"|"<<z<<std::endl;
      }
    } 
  //TrkrHitSetContainer, TrkrHitSet, TrkrHit, and TrkrDefs objects used above in offline/packages/trackbase.
  _h_hits->Fill(n_hits);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::ResetEvent(PHCompositeNode *topNode)
{
  std::cout << "readDigitalCurrents::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
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
  _h_hits    ->Sumw2( false );
  _h_DC_E   ->Sumw2( false );
  hm->dumpHistos(_filename, "UPDATE");

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
