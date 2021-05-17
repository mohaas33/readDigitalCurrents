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

#include <phool/PHCompositeNode.h>

#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>

#include <phool/getClass.h>

#include <iostream>
#include <string>
#include <sstream>
#include <set>

using namespace std;

//____________________________________________________________________________..
readDigitalCurrents::readDigitalCurrents(const std::string &name):
 SubsysReco(name)
{
  std::cout << "readDigitalCurrents::readDigitalCurrents(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
readDigitalCurrents::~readDigitalCurrents()
{
  std::cout << "readDigitalCurrents::~readDigitalCurrents() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int readDigitalCurrents::Init(PHCompositeNode *topNode)
{
  std::cout << "readDigitalCurrents::Init(PHCompositeNode *topNode) Initializing" << std::endl;
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


  // loop over all the hits
  // hits are stored in hitsets, so have to get the hitset first

  TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
       iter != all_hitsets.second;
       ++iter)
  {    
    TrkrHitSet::ConstRange range = iter->second->getHits();
    for(TrkrHitSet::ConstIterator hitr = range.first; hitr != range.second; ++hitr)
      {
       //TrkrDefs::hitkey hit_key = hitr->first;
       TrkrHit *hit = hitr->second;
       unsigned short adc = hit->getAdc();
       std::cout<<adc<<std::endl;
      }
    } 
  //TrkrHitSetContainer, TrkrHitSet, TrkrHit, and TrkrDefs objects used above in offline/packages/trackbase.

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
