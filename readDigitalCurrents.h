// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef READDIGITALCURRENTS_H
#define READDIGITALCURRENTS_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>
#include <vector>

#pragma once 
#include <math.h>

class Fun4AllHistoManager;
class PHCompositeNode;

//class PHG4CylinderCellGeom;

//class TFile;
class TH1;
class TH2;
class TH3;


class readDigitalCurrents : public SubsysReco
{
 public:

  readDigitalCurrents(const std::string &name = "readDigitalCurrents", const std::string &filename = "DC_Hist_OnPlane_WIBF.root");

  virtual ~readDigitalCurrents();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void SetBeamXing(int newBeamXing);
  void SetEvtStart(int newEvtStart);
  void SetCollSyst(int coll_syst=0);
  void SetIBF(float ampIBFfrac=0.004);
  void SetCCGC(float f_ccgc=0);

  double pi = 2 * acos(0.0);

 protected:
   Fun4AllHistoManager *hm;
   std::string _filename;
   //TFile *outfile;
   std::map<int,int> _timestamps;
   std::vector<int> _keys;
   float _ampIBFfrac;
   int _collSyst;

 private:
   int _beamxing = 0;
   int _evtstart = 0;

    int _f_ccgc = 0;

    TH2*   _h_modules_measuredibf;

    TH1*   _h_R;
    TH1*   _h_hits;
    TH3*   _h_DC_SC;
    TH3*   _h_DC_SC_XY;
    TH2*   _h_hit_XY;
    TH2*   _h_DC_E;
    TH3*   _h_SC_ibf;
    float _event_timestamp;
    float _event_bunchXing;

    double adc_pedestal=74.4;
    double cm=1e1,m=1e3, mm=1; //changed to make 'm' 1.0, for convenience.
    float ns=1e-9,us=1e-6,ms=1e-3,s=1;
    float V=1;
    //float ionMobility=3.37*cm*cm/V/s;
    float ionMobility=1.65*cm*cm/V/s;
    float vIon=ionMobility*400*V/cm;

    float f=0.5;//for now, just pick the middle of the hit.  Do better later.
};

#endif // READDIGITALCURRENTS_H
