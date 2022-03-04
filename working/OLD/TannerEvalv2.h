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
    ~TannerEval() override{};

    //Initialize SubsysReco
    int Init(PHCompositeNode *topNode) override;
    int process_event(PHCompositeNode *topNode) override;
    int End(PHCompositeNode *topNode) override;
     //Event specific
    void set_strict(bool set_strict) { m_set_strict = set_strict; }

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

  unsigned int m_iEvent;
  float m_cross_section;
  float m_event_weight;
  int m_n_generator_accepted;


  //HEPMC TTree Varibles
  int m_NHEPMCpartsinevent;
  int* m_partid1;
  int* m_partid2;
  float m_x1;
  float m_x2;
  int* m_mpi;
  int m_process_id;
  float* m_truthenergy;
  float* m_trutheta;
  float* m_truthphi;
  float* m_truthpx;
  float* m_truthpy;
  float* m_truthpz;
  float* m_truthpt;
  float* m_truthp;
  int* m_truthpid;


  //PH4G TTree variables
  int m_NPHG4partsinevent;
  int* m_truthpid_PHG4;
  float* m_truthpx_PHG4;
  float* m_truthpy_PHG4;
  float* m_truthpz_PHG4;
  float* m_truthpt_PHG4;
  float* m_truthp_PHG4;
  float* m_truthenergy_PHG4;
  float* m_trutheta_PHG4;
  float* m_truthphi_PHG4;

  // Track TTree variables
  int m_Ntracksinevent;
  float* m_tr_px;
  float* m_tr_py;
  float* m_tr_pz;
  float* m_tr_p;
  float* m_tr_pt;
  float* m_tr_phi;
  float* m_tr_eta;
  int* m_charge;
  float* m_chisq;
  int* m_ndf;
  float* m_dca;
  float* m_tr_x;
  float* m_tr_y;
  float* m_tr_z;
  int* m_truth_is_primary;
  float* m_truthtrackpx;
  float* m_truthtrackpy;
  float* m_truthtrackpz;
  float* m_truthtrackp;
  float* m_truthtracke;
  float* m_truthtrackpt;
  float* m_truthtrackphi;
  float* m_truthtracketa;
  int* m_truthtrackpid;

  /// Cluster and tower variables
  //CEMC
  int m_Nclusters_CEMC;
  float* m_clusenergy_CEMC;
  float* m_cluseta_CEMC;
  float* m_clustheta_CEMC;
  float* m_clusphi_CEMC;
  float* m_cluspt_CEMC;
  float* m_cluspx_CEMC;
  float* m_cluspy_CEMC;
  float* m_cluspz_CEMC;
  float* m_E_4x4_CEMC;
  int* m_clus_trueID_CEMC;

  int m_Ntowers_CEMC;
  float* m_tower_E_CEMC;
  int* m_tower_iEta_CEMC;
  int* m_tower_iPhi_CEMC;
  int* m_tower_trueID_CEMC;

  //HCALIN
  int m_Nclusters_HCALIN;
  float* m_clusenergy_HCALIN;
  float* m_cluseta_HCALIN;
  float* m_clustheta_HCALIN;
  float* m_cluspt_HCALIN;
  float* m_clusphi_HCALIN;
  float* m_cluspx_HCALIN;
  float* m_cluspy_HCALIN;
  float* m_cluspz_HCALIN;
  float* m_E_4x4_HCALIN;
  int* m_clus_trueID_HCALIN;

  int m_Ntowers_HCALIN;
  float* m_tower_E_HCALIN;
  int* m_tower_iEta_HCALIN;
  int* m_tower_iPhi_HCALIN;
  int* m_tower_trueID_HCALIN;

  //HCALOUT
  int m_Nclusters_HCALOUT;
  float* m_clusenergy_HCALOUT;
  float* m_cluseta_HCALOUT;
  float* m_clustheta_HCALOUT;
  float* m_cluspt_HCALOUT;
  float* m_clusphi_HCALOUT;
  float* m_cluspx_HCALOUT;
  float* m_cluspy_HCALOUT;
  float* m_cluspz_HCALOUT;
  float* m_E_4x4_HCALOUT;
  int* m_clus_trueID_HCALOUT;

  int m_Ntowers_HCALOUT;
  float* m_tower_E_HCALOUT;
  int* m_tower_iPhi_HCALOUT;
  int* m_tower_iEta_HCALOUT;
  int* m_tower_trueID_HCALOUT;

  //specific geometric variables
  int m_GEO_Ntowers;
  int m_GEO_ID;
  int* m_GEO_towers_iEta;
  int* m_GEO_towers_iPhi;
  float* m_GEO_towers_Eta;
  float* m_GEO_towers_Phi;
  float* m_GEO_towers_x;
  float* m_GEO_towers_y;
  float* m_GEO_towers_z;
  int* m_geometry_done;

  SvtxEvalStack* m_SVTX_EvalStack;
  CaloEvalStack* m_HCALIN_EvalStack;
  CaloEvalStack* m_HCALOUT_EvalStack;
  CaloEvalStack* m_CEMC_EvalStack;



  //output rootfile name
  bool m_set_strict;
  std::string m_outname;

  TFile* m_outfile;
  TFile* m_outfile_geometry;
  TTree* m_EventTree;
  TTree* m_GeoTree;



  void fillTTree(PHCompositeNode* topNode);
  void resetGeometryArrays();                             ///< reset the tree variables before filling for a new event
  void resetBuffer();                                     ///< reset the tree variables before filling for a new event


  const int m_maxNtowers_GEO = 5000000;
  const int m_maxNtowersCentral = 2000;
  const int m_maxNclusters = 2000;
  const int m_maxNtracks = 200;
  const int m_maxPHG4Part = 100000;
  const int m_maxNNEPMCPart = 1000;

  enum calotype {
      kCEMC         = 0,
      kHCALIN       = 1,
      kHCALOUT       = 2
  };



};



#endif
