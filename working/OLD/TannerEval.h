#ifndef TANNEREVAL_H__
#define TANNEREVAL_H__

#include <fun4all/SubsysReco.h>
#include <set>
#include <string>

class PHCompositeNode;
class TFile;
class TTree;
class TH1;
class PHCompositeNode;
class PHHepMCGenEventMap;
class PHHepMCGenEvent;
class PHG4TruthInfoContainer;
class RawClusterContainer;
class RawCluster;
class SvtxTrackMap;
class GlobalVertex;
class SvtxTrackEval;
class SvtxEvalStack;
class CaloEvalStack;
class CaloTriggerInfo;

class TannerEval : public SubsysReco
{
  public:
    //Constructor
    TannerEval(const std::string &name = "TannerEval", const std::string &fname = "TannerEval.root");
    //Destructor
    virtual ~TannerEval();

    //Initialize SubsysReco
    int Init(PHCompositeNode *);
    int process_event(PHCompositeNode *);
    int End(PHCompositeNode *);
     //Event specific
    void doEventInfo(bool doEventInfo) { m_doEventInfo = doEventInfo; }
    void doGeo(bool doGeo) { m_doGeo = doGeo; }

    //void set_strict(bool set_strict) { m_set_strict = set_strict; }
    //Evaluation functions
    void doTracks(bool doTracks) { m_doTracks = doTracks; }
    void doHCALIN(bool doHCALIN) { m_doHCALIN = doHCALIN; }
    void doHCALOUT(bool doHCALOUT) { m_doHCALOUT = doHCALOUT; }
    void doCEMC(bool doCEMC) { m_doCEMC = doCEMC; }

   //Event specific

    //Particles
    void doHEPMCpart(bool doHEPMCpart) { m_doHEPMCpart = doHEPMCpart; }
    void doPHG4part(bool doPHG4part) { m_doPHG4part = doPHG4part; }

private:
  //output rootfile name
  std::string m_outname;

  //booleans from  public functions

  bool m_doEventInfo;
  bool m_doGeo;

  bool m_doTracks;
  bool m_doHCALIN;
  bool m_doHCALOUT;
  bool m_doCEMC;

  bool m_doHEPMCpart;
  bool m_doPHG4part;

  //VARIBLES TO BE OBTAINED FROM THE DST NODES

  //EVENT SPECIFIC information
  TFile *m_outfile;
  TFile *m_outfile_geometry;
  TTree *m_EventTree;
  TTree *m_GeoTree;
  TTree *m_TrackTree;
  TTree *m_HCALINTree;
  TTree *m_HCALOUTTree;
  TTree *m_CEMCTree;
  TTree *m_HEPMCTree;
  TTree *m_PHG4Tree;



  SvtxEvalStack *m_SVTX_EvalStack= nullptr;
  CaloEvalStack *m_HCALIN_EvalStack= nullptr;
  CaloEvalStack *m_HCALOUT_EvalStack= nullptr;
  CaloEvalStack *m_CEMC_EvalStack= nullptr;

  void getEventInfo(PHCompositeNode *topNode);
  void getTracks(PHCompositeNode *topNode);
  void getCEMC(PHCompositeNode *topNode);
  void getHCALIN(PHCompositeNode *topNode);
  void getHCALOUT(PHCompositeNode *topNode);
  void getHEPMC(PHCompositeNode *topNode);
  void getPHG4Truth(PHCompositeNode *topNode);

  void GrowTrees();
  //void TrimLeafs();
  void WaterSeeds();


  unsigned int m_iEvent;

  double m_cross_section;
  double m_event_weight;
  int m_n_generator_accepted;


  //HEPMC TTree Varibles
  int m_NHEPMCpartsinevent;
  int m_partid1;
  int m_partid2;
  double m_x1;
  double m_x2;
  int m_mpi;
  int m_process_id;
  double m_truthenergy;
  double m_trutheta;
  double m_truthphi;
  double m_truthpx;
  double m_truthpy;
  double m_truthpz;
  double m_truthpt;
  double m_truthp;
  int m_truthpid;


  //PH4G TTree variables
  int m_NPHG4partsinevent;
  int m_truthpid_PHG4;
  double m_truthpx_PHG4;
  double m_truthpy_PHG4;
  double m_truthpz_PHG4;
  double m_truthpt_PHG4;
  double m_truthp_PHG4;
  double m_truthenergy_PHG4;
  double m_trutheta_PHG4;
  double m_truthphi_PHG4;

  // Track TTree variables
  int m_Ntracksinevent;
  double m_tr_px;
  double m_tr_py;
  double m_tr_pz;
  double m_tr_p;
  double m_tr_pt;
  double m_tr_phi;
  double m_tr_eta;
  int m_charge;
  double m_chisq;
  int m_ndf;
  double m_dca;
  double m_tr_x;
  double m_tr_y;
  double m_tr_z;
  int m_truth_is_primary;
  double m_truthtrackpx;
  double m_truthtrackpy;
  double m_truthtrackpz;
  double m_truthtrackp;
  double m_truthtracke;
  double m_truthtrackpt;
  double m_truthtrackphi;
  double m_truthtracketa;
  int m_truthtrackpid;

  /// Cluster and tower variables
  //CEMC
  int m_Nclusters_CEMC;
  double m_clusenergy_CEMC;
  double m_cluseta_CEMC;
  double m_clustheta_CEMC;
  double m_clusphi_CEMC;
  double m_cluspt_CEMC;
  double m_cluspx_CEMC;
  double m_cluspy_CEMC;
  double m_cluspz_CEMC;
  double m_E_4x4_CEMC;
  int m_clus_trueID_CEMC;

  int m_Ntowers_CEMC;
  double m_tower_E_CEMC;
  int m_tower_iEta_CEMC;
  int m_tower_iPhi_CEMC;
  int m_tower_trueID_CEMC;

  //HCALIN
  int m_Nclusters_HCALIN;
  double m_clusenergy_HCALIN;
  double m_cluseta_HCALIN;
  double m_clustheta_HCALIN;
  double m_cluspt_HCALIN;
  double m_clusphi_HCALIN;
  double m_cluspx_HCALIN;
  double m_cluspy_HCALIN;
  double m_cluspz_HCALIN;
  double m_E_4x4_HCALIN;
  int m_clus_trueID_HCALIN;

  int m_Ntowers_HCALIN;
  double m_tower_E_HCALIN;
  int m_tower_iEta_HCALIN;
  int m_tower_iPhi_HCALIN;
  int m_tower_trueID_HCALIN;

  //HCALOUT
  int m_Nclusters_HCALOUT;
  double m_clusenergy_HCALOUT;
  double m_cluseta_HCALOUT;
  double m_clustheta_HCALOUT;
  double m_cluspt_HCALOUT;
  double m_clusphi_HCALOUT;
  double m_cluspx_HCALOUT;
  double m_cluspy_HCALOUT;
  double m_cluspz_HCALOUT;
  double m_E_4x4_HCALOUT;
  int m_clus_trueID_HCALOUT;

  int m_Ntowers_HCALOUT;
  double m_tower_E_HCALOUT;
  int m_tower_iPhi_HCALOUT;
  int m_tower_iEta_HCALOUT;
  int m_tower_trueID_HCALOUT;

  //specific geometric variables
  int m_ALLcalo_Ntowers;
  int m_ALLcalo_ID;
  int m_ALLcalo_towers_iEta;
  int m_ALLcalo_towers_iPhi;
  double m_ALLcalo_towers_Eta;
  double m_ALLcalo_towers_Phi;
  double m_ALLcalo_towers_x;
  double m_ALLcalo_towers_y;
  double m_ALLcalo_towers_z;


  enum calotype {
      kCEMC         = 0,
      kHCALIN       = 1,
      kHCALOUT       = 2
  };


  const int m_maxNtowers_ALLcalo = 5000000;


};



#endif
