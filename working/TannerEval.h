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
    TannerEval(const std::string& name = "TannerEval", const std::string& filename = "TannerEval.root");
    //Destructor
    ~TannerEval() override{};

    //Initialize SubsysReco
    int Init(PHCompositeNode *topNode) override;
    int process_event(PHCompositeNode *topNode) override;
    int End(PHCompositeNode *topNode) override;

    void set_strict(bool set_strict) { m_set_strict = set_strict; }

    void set_do_store_event_level_info(bool b) { m_doEventInfo = b; }
    void set_do_HCALIN(bool b) { m_doHCALIN = b; }
    void set_do_HCALOUT(bool b) { m_doHCALOUT = b; }
    void set_do_CEMC(bool b) { m_doCEMC = b; }
    void set_do_HITS(bool b) { m_doHits = b; }
    void set_do_TRACKS(bool b) { m_doTracks = b; }
    void set_do_CLUSTERS(bool b) { m_doClusters = b; }
    void set_do_VERTEX(bool b) { m_doVertex= b; }
    void set_do_PROJECTIONS(bool b) { m_doProjections = b; }
    void set_do_MCPARTICLES(bool b) { m_doMCparts = b; }
    void set_do_HEPMC(bool b) { m_doHEPMCparts = b; }
    void set_do_GEOMETRY(bool b) { m_doGeo= b; }

    // limit the tracing of towers and clusters back to the truth particles
    // to only those reconstructed objects above a particular energy
    // threshold (evaluation for objects above threshold unaffected)
    void set_reco_tracing_energy_threshold(float thresh)
    {
      m_reco_e_threshold = thresh;
    }
    void set_reco_tracing_energy_threshold_BECAL(float thresh)
    {
      m_reco_e_threshold_BECAL = thresh;
    }

    //! max depth/generation of the MC_particle/PHG4Particle that would be saved.
    void set_depth_MCstack(int d)
    {
      m_depth_MCstack = d;
    }


private:
//booleans from  public functions
  bool m_doEventInfo;
  bool m_doHCALIN;
  bool m_doHCALOUT;
  bool m_doCEMC;
  bool m_doHits;
  bool m_doTracks;
  bool m_doClusters;
  bool m_doVertex;
  bool m_doProjections;
  bool m_doMCparts;
  bool m_doHEPMCparts;
  bool m_doGeo;

  //VARIBLES TO BE OBTAINED FROM THE DST NODES

  unsigned int m_iEvent;
  float m_cross_section;
  float m_event_weight;
  int m_n_generator_accepted;

  //Track hits
  int m_nHitsLayers;
  int* m_hits_layerID;
  int* m_hits_trueID;
  float* m_hits_x;
  float* m_hits_y;
  float* m_hits_z;
  float* m_hits_t;

  //HEPMC TTree Varibles
  int m_nHEPMCparts;
  int m_process_id;
  float m_x1;
  float m_x2;
  int* m_partid1;
  int* m_partid2;
  int* m_mpi;
  float* m_truthenergy;
  float* m_trutheta;
  float* m_truthphi;
  float* m_truthpx;
  float* m_truthpy;
  float* m_truthpz;
  float* m_truthpt;
  float* m_truthp;
  int* m_truthpid;


  //MC TTree variables
  int m_nMCparts;
  int* m_MCpart_ID;
  int* m_MCpart_ID_parent;
  int* m_MCpart_PDG;
  float* m_MCpart_E;
  float* m_MCpart_px;
  float* m_MCpart_py;
  float* m_MCpart_pz;
  int* m_MCpart_BCID;

  // Track TTree variables
  int m_nTracks;
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
  float* m_tr_ID;
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
  float* m_truthtrackID;

  //Projections
  int m_nProjections;
  float* m_track_ProjTrackID;
  int* m_track_ProjLayer;
  float* m_track_TLP_x;
  float* m_track_TLP_y;
  float* m_track_TLP_z;
  float* m_track_TLP_t;
  float* m_track_TLP_true_x;
  float* m_track_TLP_true_y;
  float* m_track_TLP_true_z;
  float* m_track_TLP_true_t;



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
  //float* m_E_4x4_HCALIN;
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
  //float* m_E_4x4_HCALOUT;
  int* m_clus_trueID_HCALOUT;

  int m_Ntowers_HCALOUT;
  float* m_tower_E_HCALOUT;
  int* m_tower_iPhi_HCALOUT;
  int* m_tower_iEta_HCALOUT;
  int* m_tower_trueID_HCALOUT;

  // vertex
  float m_vertex_x;
  float m_vertex_y;
  float m_vertex_z;
  int m_vertex_NCont;
  float m_vertex_true_x;
  float m_vertex_true_y;
  float m_vertex_true_z;

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

  float m_reco_e_threshold;
  float m_reco_e_threshold_BECAL;
  int m_depth_MCstack;

  SvtxEvalStack* m_SVTX_EvalStack;
  CaloEvalStack* m_HCALIN_EvalStack;
  CaloEvalStack* m_HCALOUT_EvalStack;
  CaloEvalStack* m_CEMC_EvalStack;



  //output rootfile name
  bool m_set_strict;
  std::string m_filename;

  TFile* m_outfile;
  TFile* m_outfile_geometry;
  TTree* m_EventTree;
  TTree* m_GeoTree;

  // subroutines
  int GetProjectionIndex(std::string projname);           ///< return track projection index for given track projection layer
  std::string GetProjectionNameFromIndex(int projindex);  ///< return track projection layer name from projection index (see GetProjectionIndex)
  void fillOutputNtuples(PHCompositeNode* topNode);       ///< dump the evaluator information into ntuple for external analysis
  void resetGeometryArrays();                             ///< reset the tree variables before filling for a new event
  void resetBuffer();                                     ///< reset the tree variables before filling for a new event


  const int m_maxNHits = 10000;
  const int m_maxNProjections = 2000;
  const int m_maxNtowers_GEO = 5000000;
  const int m_maxNtowersCentral = 2000;
  const int m_maxNclusters = 2000;
  const int m_maxNtracks = 200;
  const int m_maxMCparts = 100000;
  const int m_maxHEPMCparts = 1000;

  enum calotype {
      kCEMC         = 0,
      kHCALIN       = 1,
      kHCALOUT       = 2
  };



};



#endif
