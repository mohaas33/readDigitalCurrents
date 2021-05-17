#pragma once
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <readDigitalCurrents.h>

#include <stdio.h>
//#include <sstream>

#include <string>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libreadDigitalCurrents.so)
R__LOAD_LIBRARY(libg4dst.so)

void Fun4All_FillDCMap(  const int nEvents = 10, const int eventsInFileStart = 0, const string &fname = "/sphenix/sim/sim01/sphnxpro/MDC1/pythia8_pp/PileUp/data/DST_TRKR_G4HIT_pythia8_mb-0000000001-01617.root" )
{
  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////

  Fun4AllServer *se = Fun4AllServer::instance();
  string cd_name = "readDigitalCurrents"+std::to_string(eventsInFileStart);
  //cout<<fname_tmp<<endl;
  //readDigitalCurrents *dist_calc = new readDigitalCurrents(cd_name, foutputname);
  readDigitalCurrents *dist_calc = new readDigitalCurrents();
  se->registerSubsystem(dist_calc);
  
  // this (DST) input manager just drives the event loop
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");
  in->fileopen(fname);
  se->registerInputManager(in);
  // events = 0 => run till end of input file
  if (nEvents <= 0)
  {
    return;
  }
  cout << endl << "Running over " << nEvents << " Events" << endl;
  se->run(nEvents);
  //}
  cout << endl << "Calling End in Fun4All_readDigitalCurrents.C" << endl;
  se->End();

  cout << endl << "All done, calling delete Fun4AllServer" << endl;
  delete se;

  cout << endl << "gSystem->Exit(0)" << endl;
  gSystem->Exit(0);
}