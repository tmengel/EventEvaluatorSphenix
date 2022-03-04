#include "TannerEval.h"

#include <g4eval/CaloEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>
#include <g4eval/CaloRawTowerEval.h>
#include <g4eval/CaloTruthEval.h>
#include <g4eval/SvtxEvalStack.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv2.h>

#include <calotrigger/CaloTriggerInfo.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>
#include <trackbase_historic/SvtxVertex.h>  // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>


#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phhepmc/PHGenIntegral.h>


#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator


#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <utility>
#include <cassert>
#include <sstream>
#include <string>

using namespace std;

TannerEval::TannerEval(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , m_doEventInfo(false)
  , m_doHCALIN(false)
  , m_doHCALOUT(false)
  , m_doCEMC(false)
  , m_doHits(false)
  , m_doTracks(false)
  , m_doClusters(false)
  , m_doVertex(false)
  , m_doProjections(false)
  , m_doMCparts(false)
  , m_doHEPMCparts(false)
  , m_doGeo(false)

  , m_iEvent(0)
  , m_cross_section(0)
  , m_event_weight(0)
  , m_n_generator_accepted(0)

  , m_nHitsLayers(0)
  , m_hits_layerID(0)
  , m_hits_trueID(0)
  , m_hits_x(0)
  , m_hits_y(0)
  , m_hits_z(0)
  , m_hits_t(0)

  , m_nHEPMCparts(0)
  , m_process_id(0)
  , m_x1(0)
  , m_x2(0)
  , m_partid1(0)
  , m_partid2(0)
  , m_mpi(0)
  , m_truthenergy(0)
  , m_trutheta(0)
  , m_truthphi(0)
  , m_truthpx(0)
  , m_truthpy(0)
  , m_truthpz(0)
  , m_truthpt(0)
  , m_truthp(0)
  , m_truthpid(0)

  , m_nMCparts(0)
  , m_MCpart_ID(0)
  , m_MCpart_ID_parent(0)
  , m_MCpart_PDG(0)
  , m_MCpart_E(0)
  , m_MCpart_px(0)
  , m_MCpart_py(0)
  , m_MCpart_pz(0)
  , m_MCpart_BCID(0)

  , m_nTracks(0)
  , m_tr_px(0)
  , m_tr_py(0)
  , m_tr_pz(0)
  , m_tr_p(0)
  , m_tr_pt(0)
  , m_tr_phi(0)
  , m_tr_eta(0)
  , m_charge(0)
  , m_chisq(0)
  , m_ndf(0)
  , m_dca(0)
  , m_tr_x(0)
  , m_tr_y(0)
  , m_tr_z(0)
  , m_tr_ID(0)
  , m_truth_is_primary(0)
  , m_truthtrackpx(0)
  , m_truthtrackpy(0)
  , m_truthtrackpz(0)
  , m_truthtrackp(0)
  , m_truthtracke(0)
  , m_truthtrackpt(0)
  , m_truthtrackphi(0)
  , m_truthtracketa(0)
  , m_truthtrackpid(0)
  , m_truthtrackID(0)

  , m_nProjections(0)
  , m_track_ProjTrackID(0)
  , m_track_ProjLayer(0)
  , m_track_TLP_x(0)
  , m_track_TLP_y(0)
  , m_track_TLP_z(0)
  , m_track_TLP_t(0)
  , m_track_TLP_true_x(0)
  , m_track_TLP_true_y(0)
  , m_track_TLP_true_z(0)
  , m_track_TLP_true_t(0)

  , m_Nclusters_CEMC(0)
  , m_clusenergy_CEMC(0)
  , m_cluseta_CEMC(0)
  , m_clustheta_CEMC(0)
  , m_clusphi_CEMC(0)
  , m_cluspt_CEMC(0)
  , m_cluspx_CEMC(0)
  , m_cluspy_CEMC(0)
  , m_cluspz_CEMC(0)
  , m_E_4x4_CEMC(0)
  , m_clus_trueID_CEMC(0)

  , m_Ntowers_CEMC(0)
  , m_tower_E_CEMC(0)
  , m_tower_iEta_CEMC(0)
  , m_tower_iPhi_CEMC(0)
  , m_tower_trueID_CEMC(0)

  , m_Nclusters_HCALIN(0)
  , m_clusenergy_HCALIN(0)
  , m_cluseta_HCALIN(0)
  , m_clustheta_HCALIN(0)
  , m_cluspt_HCALIN(0)
  , m_clusphi_HCALIN(0)
  , m_cluspx_HCALIN(0)
  , m_cluspy_HCALIN(0)
  , m_cluspz_HCALIN(0)

  , m_clus_trueID_HCALIN(0)

  , m_Ntowers_HCALIN(0)
  , m_tower_E_HCALIN(0)
  , m_tower_iEta_HCALIN(0)
  , m_tower_iPhi_HCALIN(0)
  , m_tower_trueID_HCALIN(0)

  , m_Nclusters_HCALOUT(0)
  , m_clusenergy_HCALOUT(0)
  , m_cluseta_HCALOUT(0)
  , m_clustheta_HCALOUT(0)
  , m_cluspt_HCALOUT(0)
  , m_clusphi_HCALOUT(0)
  , m_cluspx_HCALOUT(0)
  , m_cluspy_HCALOUT(0)
  , m_cluspz_HCALOUT(0)

  , m_clus_trueID_HCALOUT(0)

  , m_Ntowers_HCALOUT(0)
  , m_tower_E_HCALOUT(0)
  , m_tower_iPhi_HCALOUT(0)
  , m_tower_iEta_HCALOUT(0)
  , m_tower_trueID_HCALOUT(0)

  , m_vertex_x(0)
  , m_vertex_y(0)
  , m_vertex_z(0)
  , m_vertex_NCont(0)
  , m_vertex_true_x(0)
  , m_vertex_true_y(0)
  , m_vertex_true_z(0)

  , m_GEO_Ntowers(0)
  , m_GEO_ID(0)
  , m_GEO_towers_iEta(0)
  , m_GEO_towers_iPhi(0)
  , m_GEO_towers_Eta(0)
  , m_GEO_towers_Phi(0)
  , m_GEO_towers_x(0)
  , m_GEO_towers_y(0)
  , m_GEO_towers_z(0)
  , m_geometry_done(0)

  , m_reco_e_threshold(0.0)
  , m_reco_e_threshold_BECAL(0.0)
  , m_depth_MCstack(0)


  , m_SVTX_EvalStack(nullptr)
  , m_HCALIN_EvalStack(nullptr)
  , m_HCALOUT_EvalStack(nullptr)
  , m_CEMC_EvalStack(nullptr)

  , m_set_strict(false)
  , m_filename(filename)

  , m_outfile(nullptr)
  , m_outfile_geometry(nullptr)
  , m_EventTree(nullptr)
  , m_GeoTree(nullptr)


{
  m_hits_layerID = new int[m_maxNHits];
  m_hits_trueID = new int[m_maxNHits];
  m_hits_x = new float[m_maxNHits];
  m_hits_y = new float[m_maxNHits];
  m_hits_z = new float[m_maxNHits];
  m_hits_t = new float[m_maxNHits];

  m_partid1 = new int[m_maxHEPMCparts];
  m_partid2 = new int[m_maxHEPMCparts];
  m_mpi = new int[m_maxHEPMCparts];
  m_truthenergy = new float[m_maxHEPMCparts];
  m_trutheta = new float[m_maxHEPMCparts];
  m_truthphi = new float[m_maxHEPMCparts];
  m_truthpx  = new float[m_maxHEPMCparts];
  m_truthpy  = new float[m_maxHEPMCparts];
  m_truthpz  = new float[m_maxHEPMCparts];
  m_truthpt  = new float[m_maxHEPMCparts];
  m_truthp  = new float[m_maxHEPMCparts];
  m_truthpid  = new int[m_maxHEPMCparts];

  m_MCpart_ID = new int[m_maxMCparts];
  m_MCpart_ID_parent = new int[m_maxMCparts];
  m_MCpart_PDG = new int[m_maxMCparts];
  m_MCpart_E = new float[m_maxMCparts];
  m_MCpart_px = new float[m_maxMCparts];
  m_MCpart_py = new float[m_maxMCparts];
  m_MCpart_pz = new float[m_maxMCparts];
  m_MCpart_BCID = new int[m_maxMCparts];


  m_tr_px = new float[m_maxNtracks];
  m_tr_py = new float[m_maxNtracks];
  m_tr_pz = new float[m_maxNtracks];
  m_tr_p = new float[m_maxNtracks];
  m_tr_pt = new float[m_maxNtracks];
  m_tr_phi = new float[m_maxNtracks];
  m_tr_eta = new float[m_maxNtracks];
  m_charge = new int[m_maxNtracks];
  m_chisq = new float[m_maxNtracks];
  m_ndf = new int[m_maxNtracks];
  m_dca = new float[m_maxNtracks];
  m_tr_x = new float[m_maxNtracks];
  m_tr_y = new float[m_maxNtracks];
  m_tr_z = new float[m_maxNtracks];
  m_tr_ID = new float[m_maxNtracks];
  m_truth_is_primary = new int[m_maxNtracks];
  m_truthtrackpx = new float[m_maxNtracks];
  m_truthtrackpy = new float[m_maxNtracks];
  m_truthtrackpz = new float[m_maxNtracks];
  m_truthtrackp = new float[m_maxNtracks];
  m_truthtracke = new float[m_maxNtracks];
  m_truthtrackpt = new float[m_maxNtracks];
  m_truthtrackphi = new float[m_maxNtracks];
  m_truthtracketa = new float[m_maxNtracks];
  m_truthtrackpid = new int[m_maxNtracks];
  m_truthtrackID = new float[m_maxNtracks];

  m_track_ProjTrackID = new float[m_maxNProjections];
  m_track_ProjLayer = new int[m_maxNProjections];
  m_track_TLP_x = new float[m_maxNProjections];
  m_track_TLP_y = new float[m_maxNProjections];
  m_track_TLP_z = new float[m_maxNProjections];
  m_track_TLP_t = new float[m_maxNProjections];
  m_track_TLP_true_x = new float[m_maxNProjections];
  m_track_TLP_true_y = new float[m_maxNProjections];
  m_track_TLP_true_z = new float[m_maxNProjections];
  m_track_TLP_true_t = new float[m_maxNProjections];

  m_clusenergy_CEMC= new float[m_maxNclusters];
  m_cluseta_CEMC= new float[m_maxNclusters];
  m_clustheta_CEMC= new float[m_maxNclusters];
  m_clusphi_CEMC= new float[m_maxNclusters];
  m_cluspt_CEMC= new float[m_maxNclusters];
  m_cluspy_CEMC= new float[m_maxNclusters];
  m_cluspx_CEMC= new float[m_maxNclusters];
  m_cluspz_CEMC= new float[m_maxNclusters];
  m_E_4x4_CEMC= new float[m_maxNclusters];
  m_clus_trueID_CEMC= new int[m_maxNclusters];

  m_tower_E_CEMC= new float[m_maxNtowersCentral];
  m_tower_iEta_CEMC= new int[m_maxNtowersCentral];
  m_tower_iPhi_CEMC= new int[m_maxNtowersCentral];
  m_tower_trueID_CEMC= new int[m_maxNtowersCentral];

  m_clusenergy_HCALIN= new float[m_maxNclusters];
  m_cluseta_HCALIN= new float[m_maxNclusters];
  m_clustheta_HCALIN= new float[m_maxNclusters];
  m_clusphi_HCALIN= new float[m_maxNclusters];
  m_cluspt_HCALIN= new float[m_maxNclusters];
  m_cluspy_HCALIN= new float[m_maxNclusters];
  m_cluspx_HCALIN= new float[m_maxNclusters];
  m_cluspz_HCALIN= new float[m_maxNclusters];

  m_clus_trueID_HCALIN= new int[m_maxNclusters];

  m_tower_E_HCALIN= new float[m_maxNtowersCentral];
  m_tower_iEta_HCALIN= new int[m_maxNtowersCentral];
  m_tower_iPhi_HCALIN= new int[m_maxNtowersCentral];
  m_tower_trueID_HCALIN= new int[m_maxNtowersCentral];

  m_clusenergy_HCALOUT= new float[m_maxNclusters];
  m_cluseta_HCALOUT= new float[m_maxNclusters];
  m_clustheta_HCALOUT= new float[m_maxNclusters];
  m_clusphi_HCALOUT= new float[m_maxNclusters];
  m_cluspt_HCALOUT= new float[m_maxNclusters];
  m_cluspy_HCALOUT= new float[m_maxNclusters];
  m_cluspx_HCALOUT= new float[m_maxNclusters];
  m_cluspz_HCALOUT= new float[m_maxNclusters];

  m_clus_trueID_HCALOUT= new int[m_maxNclusters];

  m_tower_E_HCALOUT= new float[m_maxNtowersCentral];
  m_tower_iEta_HCALOUT= new int[m_maxNtowersCentral];
  m_tower_iPhi_HCALOUT= new int[m_maxNtowersCentral];
  m_tower_trueID_HCALOUT= new int[m_maxNtowersCentral];

  m_GEO_towers_iEta = new int[m_maxNtowers_GEO];
  m_GEO_towers_iPhi = new int[m_maxNtowers_GEO];
  m_GEO_towers_Eta = new float[m_maxNtowers_GEO];
  m_GEO_towers_Phi = new float[m_maxNtowers_GEO];
  m_GEO_towers_x = new float[m_maxNtowers_GEO];
  m_GEO_towers_y = new float[m_maxNtowers_GEO];
  m_GEO_towers_z = new float[m_maxNtowers_GEO];
  m_geometry_done = new int[20];

  for(int igem=0;igem<20;igem++){
      m_geometry_done[igem] = 0;
  }

}



int TannerEval::Init(PHCompositeNode* /*topNode*/){

    if (Verbosity() > 0 ) { cout << "Initializing TannerEval" <<endl;}
    //Create output TFLIES
    m_outfile = new TFile(m_filename.c_str(), "RECREATE");

    m_iEvent = 0;
    m_EventTree = new TTree("EventTree", "EventTree");

    if (m_doEventInfo){
       // Event level info. This isn't the most efficient way to store this info, but it's straightforward
       // within the structure of the class, so the size is small compared to the rest of the output.
       m_EventTree->Branch("m_cross_section", &m_cross_section, "m_cross_section/F");
       m_EventTree->Branch("m_event_weight", &m_event_weight, "m_event_weight/F");
       m_EventTree->Branch("m_n_generator_accepted", &m_n_generator_accepted, "m_n_generator_accepted/F");
     }
    if (m_doTracks){
       m_EventTree->Branch("m_nTracks", &m_nTracks ,"m_nTracks/I");
       m_EventTree->Branch("m_tr_px", m_tr_px, "m_tr_px[m_nTracks]/F");
       m_EventTree->Branch("m_tr_py", m_tr_py, "m_tr_py[m_nTracks]/F");
       m_EventTree->Branch("m_tr_pz", m_tr_pz, "m_tr_pz[m_nTracks]/F");
       m_EventTree->Branch("m_tr_p", m_tr_p, "m_tr_p[m_nTracks]/F");
       m_EventTree->Branch("m_tr_pt", m_tr_pt, "m_tr_pt[m_nTracks]/F");
       m_EventTree->Branch("m_tr_phi", m_tr_phi, "m_tr_phi[m_nTracks]/F");
       m_EventTree->Branch("m_tr_eta", m_tr_eta, "m_tr_eta[m_nTracks]/F");
       m_EventTree->Branch("m_charge", m_charge, "m_charge[m_nTracks]/I");
       m_EventTree->Branch("m_chisq", m_chisq, "m_chisq[m_nTracks]/F");
       m_EventTree->Branch("m_ndf", m_ndf, "m_ndf[m_nTracks]/I");
       m_EventTree->Branch("m_dca", m_dca, "m_dca[m_nTracks]/F");
       m_EventTree->Branch("m_tr_x", m_tr_x, "m_tr_x[m_nTracks]/F");
       m_EventTree->Branch("m_tr_y", m_tr_y, "m_tr_y[m_nTracks]/F");
       m_EventTree->Branch("m_tr_z", m_tr_z, "m_tr_z[m_nTracks]/F");
       m_EventTree->Branch("m_tr_ID", m_tr_ID, "m_tr_ID[m_nTracks]/F");
       m_EventTree->Branch("m_truth_is_primary", m_truth_is_primary, "m_truth_is_primary[m_nTracks]/I");
       m_EventTree->Branch("m_truthtrackpx", m_truthtrackpx, "m_truthtrackpx[m_nTracks]/F");
       m_EventTree->Branch("m_truthtrackpy", m_truthtrackpy, "m_truthtrackpy[m_nTracks]/F");
       m_EventTree->Branch("m_truthtrackpz", m_truthtrackpz, "m_truthtrackpz[m_nTracks]/F");
       m_EventTree->Branch("m_truthtrackp", m_truthtrackp, "m_truthtrackp[m_nTracks]/F");
       m_EventTree->Branch("m_truthtracke", m_truthtracke, "m_truthtracke[m_nTracks]/F");
       m_EventTree->Branch("m_truthtrackpt", m_truthtrackpt, "m_truthtrackpt[m_nTracks]/F");
       m_EventTree->Branch("m_truthtrackphi", m_truthtrackphi, "m_truthtrackphi[m_nTracks]/F");
       m_EventTree->Branch("m_truthtracketa", m_truthtracketa, "m_truthtracketa[m_nTracks]/F");
       m_EventTree->Branch("m_truthtrackpid", m_truthtrackpid, "m_truthtrackpid[m_nTracks]/I");
       m_EventTree->Branch("m_truthtrackID", m_truthtrackID, "m_truthtrackID[m_nTracks]/F");

     }
    if (m_doProjections){
      m_EventTree->Branch("m_nProjections", &m_nProjections, "m_nProjections/I");
      m_EventTree->Branch("m_track_ProjTrackID", m_track_ProjTrackID, "m_track_ProjTrackID[m_nProjections]/F");
      m_EventTree->Branch("m_track_ProjLayer", m_track_ProjLayer, "m_track_ProjLayer[m_nProjections]/I");
      m_EventTree->Branch("m_track_TLP_x", m_track_TLP_x, "m_track_TLP_x[m_nProjections]/F");
      m_EventTree->Branch("m_track_TLP_y", m_track_TLP_y, "m_track_TLP_y[m_nProjections]/F");
      m_EventTree->Branch("m_track_TLP_z", m_track_TLP_z, "m_track_TLP_z[m_nProjections]/F");
      m_EventTree->Branch("m_track_TLP_t", m_track_TLP_t, "m_track_TLP_t[m_nProjections]/F");
      m_EventTree->Branch("m_track_TLP_true_x", m_track_TLP_true_x, "m_track_TLP_true_x[m_nProjections]/F");
      m_EventTree->Branch("m_track_TLP_true_y", m_track_TLP_true_y, "m_track_TLP_true_y[m_nProjections]/F");
      m_EventTree->Branch("m_track_TLP_true_z", m_track_TLP_true_z, "m_track_TLP_true_z[m_nProjections]/F");
      m_EventTree->Branch("m_track_TLP_true_t", m_track_TLP_true_t, "m_track_TLP_true_t[m_nProjections]/F");
    }
    if (m_doHits){
       m_EventTree->Branch("m_nHits", &m_nHitsLayers, "m_nHits/I");
       m_EventTree->Branch("m_hits_layerID", m_hits_layerID, "m_hits_layerID[m_nHits]/I");
       m_EventTree->Branch("m_hits_trueID", m_hits_trueID, "m_hits_trueID[m_nHits]/I");
       m_EventTree->Branch("m_hits_x", m_hits_x, "m_hits_x[m_nHits]/F");
       m_EventTree->Branch("m_hits_y", m_hits_y, "m_hits_y[m_nHits]/F");
       m_EventTree->Branch("m_hits_z", m_hits_z, "m_hits_z[m_nHits]/F");
       m_EventTree->Branch("m_hits_t", m_hits_t, "m_hits_t[m_nHits]/F");

     }
    if (m_doHCALIN){

      if(m_doClusters){
     m_EventTree->Branch("m_Nclusters_HCALIN", &m_Nclusters_HCALIN,"m_Nclusters_HCALIN/I");
     m_EventTree->Branch("m_clusenergy_HCALIN", m_clusenergy_HCALIN, "m_clusenergy_HCALIN[m_Nclusters_HCALIN]/F");
     m_EventTree->Branch("m_cluseta_HCALIN", m_cluseta_HCALIN, "m_cluseta_HCALIN[m_Nclusters_HCALIN]/F");
     m_EventTree->Branch("m_clustheta_HCALIN", m_clustheta_HCALIN, "m_clustheta_HCALIN[m_Nclusters_HCALIN]/F");
     m_EventTree->Branch("m_cluspt_HCALIN", m_cluspt_HCALIN, "m_cluspt_HCALIN[m_Nclusters_HCALIN]/F");
     m_EventTree->Branch("m_clusphi_HCALIN", m_clusphi_HCALIN, "m_clusphi_HCALIN[m_Nclusters_HCALIN]/F");
     m_EventTree->Branch("m_cluspx_HCALIN", m_cluspx_HCALIN, "m_cluspx_HCALIN[m_Nclusters_HCALIN]/F");
     m_EventTree->Branch("m_cluspy_HCALIN", m_cluspy_HCALIN, "m_cluspy_HCALIN[m_Nclusters_HCALIN]/F");
     m_EventTree->Branch("m_cluspz_HCALIN", m_cluspz_HCALIN, "m_cluspz_HCALIN[m_Nclusters_HCALIN]/F");

     m_EventTree->Branch("m_clus_trueID_HCALIN", m_clus_trueID_HCALIN, "m_clus_trueID_HCALIN[m_Nclusters_HCALIN]/F");
   }

     m_EventTree->Branch("m_Ntowers_HCALIN",&m_Ntowers_HCALIN,"m_Ntowers_HCALIN/I");
     m_EventTree->Branch("m_tower_E_HCALIN",m_tower_E_HCALIN,"m_tower_E_HCALIN[m_Ntowers_HCALIN]/F");
     m_EventTree->Branch("m_tower_iEta_HCALIN",m_tower_iEta_HCALIN,"m_tower_iEta_HCALIN[m_Ntowers_HCALIN]/I");
     m_EventTree->Branch("m_tower_iPhi_HCALIN",m_tower_iPhi_HCALIN,"m_tower_iPhi_HCALIN[m_Ntowers_HCALIN]/I");
     m_EventTree->Branch("m_tower_trueID_HCALIN",m_tower_trueID_HCALIN,"m_tower_trueID_HCALIN[m_Ntowers_HCALIN]/I");


     }
    if (m_doCEMC){
      if(m_doClusters){
     m_EventTree->Branch("m_Nclusters_CEMC", &m_Nclusters_CEMC,"m_Nclusters_CEMC/I");
     m_EventTree->Branch("m_clusenergy_CEMC", m_clusenergy_CEMC, "m_clusenergy_CEMC[m_Nclusters_CEMC]/F");
     m_EventTree->Branch("m_cluseta_CEMC", m_cluseta_CEMC, "m_cluseta_CEMC[m_Nclusters_CEMC]/F");
     m_EventTree->Branch("m_clustheta_CEMC", m_clustheta_CEMC, "m_clustheta_CEMC[m_Nclusters_CEMC]/F");
     m_EventTree->Branch("m_cluspt_CEMC", m_cluspt_CEMC, "m_cluspt_CEMC[m_Nclusters_CEMC]/F");
     m_EventTree->Branch("m_clusphi_CEMC", m_clusphi_CEMC, "m_clusphi_CEMC[m_Nclusters_CEMC]/F");
     m_EventTree->Branch("m_cluspx_CEMC", m_cluspx_CEMC, "m_cluspx_CEMC[m_Nclusters_CEMC]/F");
     m_EventTree->Branch("m_cluspy_CEMC", m_cluspy_CEMC, "m_cluspy_CEMC[m_Nclusters_CEMC]/F");
     m_EventTree->Branch("m_cluspz_CEMC", m_cluspz_CEMC, "m_cluspz_CEMC[m_Nclusters_CEMC]/F");
     m_EventTree->Branch("m_E_4x4_CEMC", m_E_4x4_CEMC, "m_E_4x4_CEMC[m_Nclusters_CEMC]/F");
     m_EventTree->Branch("m_clus_trueID_CEMC", m_clus_trueID_CEMC, "m_clus_trueID_CEMC[m_Nclusters_CEMC]/F");
   }

     m_EventTree->Branch("m_Ntowers_CEMC",&m_Ntowers_CEMC,"m_Ntowers_CEMC/I");
     m_EventTree->Branch("m_tower_E_CEMC",m_tower_E_CEMC,"m_tower_E_CEMC[m_Ntowers_CEMC]/F");
     m_EventTree->Branch("m_tower_iEta_CEMC",m_tower_iEta_CEMC,"m_tower_iEta_CEMC[m_Ntowers_CEMC]/I");
     m_EventTree->Branch("m_tower_iPhi_CEMC",m_tower_iPhi_CEMC,"m_tower_iPhi_CEMC[m_Ntowers_CEMC]/I");
     m_EventTree->Branch("m_tower_trueID_CEMC",m_tower_trueID_CEMC,"m_tower_trueID_CEMC[m_Ntowers_CEMC]/I");
     }
    if (m_doHCALOUT){
      if(m_doClusters){
     m_EventTree->Branch("m_Nclusters_HCALOUT", &m_Nclusters_HCALOUT,"m_Nclusters_HCALOUT/I");
     m_EventTree->Branch("m_clusenergy_HCALOUT", m_clusenergy_HCALOUT, "m_clusenergy_HCALOUT[m_Nclusters_HCALOUT]/F");
     m_EventTree->Branch("m_cluseta_HCALOUT", m_cluseta_HCALOUT, "m_cluseta_HCALOUT[m_Nclusters_HCALOUT]/F");
     m_EventTree->Branch("m_clustheta_HCALOUT", m_clustheta_HCALOUT, "m_clustheta_HCALOUT[m_Nclusters_HCALOUT]/F");
     m_EventTree->Branch("m_cluspt_HCALOUT", m_cluspt_HCALOUT, "m_cluspt_HCALOUT[m_Nclusters_HCALOUT]/F");
     m_EventTree->Branch("m_clusphi_HCALOUT", m_clusphi_HCALOUT, "m_clusphi_HCALOUT[m_Nclusters_HCALOUT]/F");
     m_EventTree->Branch("m_cluspx_HCALOUT", m_cluspx_HCALOUT, "m_cluspx_HCALOUT[m_Nclusters_HCALOUT]/F");
     m_EventTree->Branch("m_cluspy_HCALOUT", m_cluspy_HCALOUT, "m_cluspy_HCALOUT[m_Nclusters_HCALOUT]/F");
     m_EventTree->Branch("m_cluspz_HCALOUT", m_cluspz_HCALOUT, "m_cluspz_HCALOUT[m_Nclusters_HCALOUT]/F");
     //m_EventTree->Branch("m_E_4x4_HCALOUT", m_E_4x4_HCALOUT, "m_E_4x4_HCALOUT[m_Nclusters_HCALOUT]/F");
     m_EventTree->Branch("m_clus_trueID_HCALOUT", m_clus_trueID_HCALOUT, "m_clus_trueID_HCALOUT[m_Nclusters_HCALOUT]/F");
}
     m_EventTree->Branch("m_Ntowers_HCALOUT",&m_Ntowers_HCALOUT,"m_Ntowers_HCALOUT/I");
     m_EventTree->Branch("m_tower_E_HCALOUT",m_tower_E_HCALOUT,"m_tower_E_HCALOUT[m_Ntowers_HCALOUT]/F");
     m_EventTree->Branch("m_tower_iEta_HCALOUT",m_tower_iEta_HCALOUT,"m_tower_iEta_HCALOUT[m_Ntowers_HCALOUT]/I");
     m_EventTree->Branch("m_tower_iPhi_HCALOUT",m_tower_iPhi_HCALOUT,"m_tower_iPhi_HCALOUT[m_Ntowers_HCALOUT]/I");
     m_EventTree->Branch("m_tower_trueID_HCALOUT",m_tower_trueID_HCALOUT,"m_tower_trueID_HCALOUT[m_Ntowers_HCALOUT]/I");
     }


    if (m_doHEPMCparts){
        m_EventTree->Branch("m_nHEPMCparts", &m_nHEPMCparts,"m_nHEPMCparts/I");
        m_EventTree->Branch("m_process_id", &m_process_id, "m_process_id/I");
        m_EventTree->Branch("m_x1", &m_x1, "m_x1/F");
        m_EventTree->Branch("m_x2", &m_x2, "m_x2/F");
        m_EventTree->Branch("m_partid1", m_partid1, "m_partid1[m_nHEPMCparts]/I");
        m_EventTree->Branch("m_partid2", m_partid2, "m_partid2[m_nHEPMCparts]/I");
        m_EventTree->Branch("m_mpi", m_mpi, "m_mpi[m_nHEPMCparts]/I");
        m_EventTree->Branch("m_truthenergy", m_truthenergy, "m_truthenergy[m_nHEPMCparts]/F");
        m_EventTree->Branch("m_trutheta", m_trutheta, "m_trutheta[m_nHEPMCparts]/F");
        m_EventTree->Branch("m_truthphi", m_truthphi, "m_truthphi[m_nHEPMCparts]/F");
        m_EventTree->Branch("m_truthpx", m_truthpx, "m_truthpx[m_nHEPMCparts]/F");
        m_EventTree->Branch("m_truthpy", m_truthpy, "m_truthpy[m_nHEPMCparts]/F");
        m_EventTree->Branch("m_truthpz", m_truthpz, "m_truthpz[m_nHEPMCparts]/F");
        m_EventTree->Branch("m_truthpt", m_truthpt, "m_truthpt[m_nHEPMCparts]/F");
        m_EventTree->Branch("m_truthpid", m_truthpid, "m_truthpid[m_nHEPMCparts]/I");

     }
    if (m_doMCparts){
       m_EventTree->Branch("m_nMCparts", &m_nMCparts,"m_nMCparts/I");
       m_EventTree->Branch("m_MCpart_ID", m_MCpart_ID, "m_MCpart_ID[m_nMCparts]/I");
       m_EventTree->Branch("m_MCpart_ID_parent", m_MCpart_ID_parent, "m_MCpart_ID_parent[m_nMCparts]/I");
       m_EventTree->Branch("m_MCpart_PDG", m_MCpart_PDG, "m_MCpart_PDG[m_nMCparts]/I");
       m_EventTree->Branch("m_MCpart_E", m_MCpart_E, "m_MCpart_E[m_nMCparts]/F");
       m_EventTree->Branch("m_MCpart_px", m_MCpart_px, "m_MCpart_px[m_nMCparts]/F");
       m_EventTree->Branch("m_MCpart_py", m_MCpart_py, "m_MCpart_py[m_nMCparts]/F");
       m_EventTree->Branch("m_MCpart_pz", m_MCpart_pz, "m_MCpart_pz[m_nMCparts]/F");
       m_EventTree->Branch("m_MCpart_BCID", m_MCpart_BCID, "m_MCpart_BCID[m_nMCparts]/I");

     }
    if (m_doVertex){
       // vertex
       m_EventTree->Branch("m_vertex_x", &m_vertex_x, "m_vertex_x/F");
       m_EventTree->Branch("m_vertex_y", &m_vertex_y, "m_vertex_y/F");
       m_EventTree->Branch("m_vertex_z", &m_vertex_z, "m_vertex_z/F");
       m_EventTree->Branch("m_vertex_NCont", &m_vertex_NCont, "m_vertex_NCont/I");
       m_EventTree->Branch("m_vertex_true_x", &m_vertex_true_x, "m_vertex_true_x/F");
       m_EventTree->Branch("m_vertex_true_y", &m_vertex_true_y, "m_vertex_true_y/F");
       m_EventTree->Branch("m_vertex_true_z", &m_vertex_true_z, "m_vertex_true_z/F");
     }
    if (m_doGeo){

       m_outfile_geometry = new TFile("geometry.root","RECREATE");

         m_GeoTree = new TTree("GeoTree", "GeoTree");
         m_GeoTree->Branch("m_GEO_ID", &m_GEO_ID, "nCaloHits/I");
         m_GeoTree->Branch("m_GEO_Ntowers", &m_GEO_Ntowers, "m_GEO_Ntowers/I");
         m_GeoTree->Branch("m_GEO_towers_iEta", m_GEO_towers_iEta,"calo_towers_iEta[m_GEO_Ntowers]/I");
         m_GeoTree->Branch("m_GEO_towers_iPhi", m_GEO_towers_iPhi,"calo_towers_iPhi[m_GEO_Ntowers]/I");
         m_GeoTree->Branch("m_GEO_towers_Eta", m_GEO_towers_Eta,"calo_towers_Eta[m_GEO_Ntowers]/F");
         m_GeoTree->Branch("m_GEO_towers_Phi", m_GEO_towers_Phi, "calo_towers_Phi[m_GEO_Ntowers]/F");
         m_GeoTree->Branch("m_GEO_towers_x", m_GEO_towers_x,"calo_towers_x[m_GEO_Ntowers]/F");
         m_GeoTree->Branch("m_GEO_towers_y", m_GEO_towers_y,"calo_towers_y[m_GEO_Ntowers]/F");
         m_GeoTree->Branch("m_GEO_towers_z", m_GEO_towers_z,"calo_towers_z[m_GEO_Ntowers]/F");


       }


    return Fun4AllReturnCodes::EVENT_OK;

}

int TannerEval::process_event(PHCompositeNode* topNode){


      if (Verbosity() > 0 ) cout << "Processing event via TannerEval" <<endl;



      if (m_doCEMC){
        //check if CEMC eval stack set_exist
        if(!m_CEMC_EvalStack){
          m_CEMC_EvalStack = new CaloEvalStack(topNode,"CEMC");
          m_CEMC_EvalStack->set_strict(m_set_strict);
          m_CEMC_EvalStack->set_verbosity(Verbosity()+1);
          }
        else { m_CEMC_EvalStack->next_event(topNode); }

        cout << "loaded EM Cal evalstack" << endl;



      }

      if (m_doHCALIN){
        if(!m_HCALIN_EvalStack){
          m_HCALIN_EvalStack = new CaloEvalStack(topNode,"HCALIN");
          m_HCALIN_EvalStack->set_strict(m_set_strict);
          m_HCALIN_EvalStack->set_verbosity(Verbosity()+1);
        }
        else {m_HCALIN_EvalStack->next_event(topNode); }

        cout << "loaded HCALIN evalstack" << endl;

        }

      if (m_doHCALOUT){

        if(!m_HCALOUT_EvalStack){
          m_HCALOUT_EvalStack = new CaloEvalStack(topNode,"HCALOUT");
          m_HCALOUT_EvalStack->set_strict(m_set_strict);
          m_HCALOUT_EvalStack->set_verbosity(Verbosity()+1);
        }

        else { m_HCALOUT_EvalStack->next_event(topNode);}

        cout << "loaded OUT HCal evalstack" << endl;


      }
      if(m_doTracks){
        if(!m_SVTX_EvalStack)
          {
            m_SVTX_EvalStack = new SvtxEvalStack(topNode);
            m_SVTX_EvalStack->set_verbosity(Verbosity());
          }

        else{ m_SVTX_EvalStack->next_event(topNode);}

      }

      fillOutputNtuples(topNode);
      //just in case
      m_iEvent++;
      cout << "Event " << m_iEvent << " complete" <<endl;

      return Fun4AllReturnCodes::EVENT_OK;


  }


void TannerEval::fillOutputNtuples(PHCompositeNode* topNode){

  if(m_doEventInfo){

      PHHepMCGenEventMap* hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
      if (hepmceventmap){
         cout << "saving event level info" << endl;

        for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
             eventIter != hepmceventmap->end();
             ++eventIter){

          PHHepMCGenEvent* hepmcevent = eventIter->second;

          if (hepmcevent){
            HepMC::GenEvent* truthevent = hepmcevent->getEvent();
            if (!truthevent){
              cout << PHWHERE << "no evt pointer under phhepmvgeneventmap found " << endl;
              return;
            }

            auto xsec = truthevent->cross_section();
            if (xsec) { m_cross_section = xsec->cross_section();}
            // Only fill the event weight if available.
            // The overall event weight will be stored in the last entry in the vector.
            auto weights = truthevent->weights();
            if (weights.size() > 0) {   m_event_weight = weights[weights.size() - 1]; }
          }
        }
      }

      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " PHHepMCGenEventMap node (for event level info) not found on node tree" << endl;
        }

      }

      // Retrieve the number of generator accepted events
      // Following how this was implemented in PHPythia8
      PHNodeIterator iter(topNode);
      PHCompositeNode *sumNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
      if (!sumNode)
      {
        cout << PHWHERE << "RUN Node missing doing nothing" << endl;
        return;
      }
      auto * integralNode = findNode::getClass<PHGenIntegral>(sumNode, "PHGenIntegral");
      if (integralNode)
      {
        m_n_generator_accepted = integralNode->get_N_Generator_Accepted_Event();
      }
      else
      {
        if (Verbosity() > 0){cout << PHWHERE << " PHGenIntegral node (for n generator accepted) not found on node tree. Continuing" << endl; }
      }

  }

  SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(m_doVertex){
    if (vertexmap){
      if (!vertexmap->empty()){
        if (Verbosity() > 0){
          cout << "saving vertex" << endl;
        }
        SvtxVertex* vertex = (vertexmap->begin())->second;

        m_vertex_x = vertex->get_x();
        m_vertex_y = vertex->get_y();
        m_vertex_z = vertex->get_z();
        m_vertex_NCont = vertex->size_tracks();
      }
      else {
        m_vertex_x = 0.;
        m_vertex_y = 0.;
        m_vertex_z = 0.;
        m_vertex_NCont = -1;
      }
    }

  }
  if(m_doHits){
    if (Verbosity() > 0)
    {
      cout << "saving hits" << endl;
    }
    m_nHitsLayers = 0;
    PHG4TruthInfoContainer* truthinfocontainerHits = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    for (int iIndex = 0; iIndex < 60; ++iIndex)
    {
      // you need to add your layer name here to be saved! This has to be done
      // as we do not want to save thousands of calorimeter hits!
      if ((GetProjectionNameFromIndex(iIndex).find("MVTX") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("INTT") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("TPC") != std::string::npos)  ||
          (GetProjectionNameFromIndex(iIndex).find("MICROMEGAS") != std::string::npos)
      ){
        string nodename = "G4HIT_" + GetProjectionNameFromIndex(iIndex);
        PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
        if (hits)
        {
          if (Verbosity() > 1)
          {
            cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;
          }
          PHG4HitContainer::ConstRange hit_range = hits->getHits();
          for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
          {
            if (Verbosity() > 1)
            {
              cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;
            }
            if(m_nHitsLayers > m_maxNHits){
              cout << __PRETTY_FUNCTION__ << " exceededed maximum hit array size! Please check where these hits come from!" << endl;
              break;
            }
            m_hits_x[m_nHitsLayers] = hit_iter->second->get_x(0);
            m_hits_y[m_nHitsLayers] = hit_iter->second->get_y(0);
            m_hits_z[m_nHitsLayers] = hit_iter->second->get_z(0);
            m_hits_t[m_nHitsLayers] = hit_iter->second->get_t(0);
            m_hits_layerID[m_nHitsLayers] = iIndex;
            if (truthinfocontainerHits)
            {
              PHG4Particle* particle = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());

              if (particle->get_parent_id() != 0)
              {
                PHG4Particle* g4particleMother = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());
                int mcSteps = 0;
                while (g4particleMother->get_parent_id() != 0)
                {
                  g4particleMother = truthinfocontainerHits->GetParticle(g4particleMother->get_parent_id());
                  if (g4particleMother == NULL) break;
                  mcSteps += 1;
                }
                if (mcSteps <= m_depth_MCstack)
                {
                  m_hits_trueID[m_nHitsLayers] = hit_iter->second->get_trkid();
                }
                else
                {
                  PHG4Particle* g4particleMother2 = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());
                  int mcSteps2 = 0;
                  while (g4particleMother2->get_parent_id() != 0 && (mcSteps2 < (mcSteps - m_depth_MCstack + 1)))
                  {
                    g4particleMother2 = truthinfocontainerHits->GetParticle(g4particleMother2->get_parent_id());
                    if (g4particleMother2 == NULL){
                      break;
                    } else {
                      m_hits_trueID[m_nHitsLayers] = g4particleMother2->get_parent_id();
                      mcSteps2 += 1;
                    }
                  }
                }
              }
              else
              {
                m_hits_trueID[m_nHitsLayers] = hit_iter->second->get_trkid();
              }
            }
            m_nHitsLayers++;

          }
          if (Verbosity() > 0)
          {
            cout << "saved\t" << m_nHitsLayers << "\thits for " << GetProjectionNameFromIndex(iIndex) << endl;
          }
        }
        else
        {
          if (Verbosity() > 0)
          {
            cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;
          }
          continue;
        }
      }
    }
  }
  if(m_doCEMC){

            //Do towers for CEMC from method in EventEvaluator.h
        if (Verbosity() > 0) { cout << "Grabbing CEMC event towers" << endl; }
        m_Ntowers_CEMC = 0;



        CaloRawTowerEval* towerevalCEMC = m_CEMC_EvalStack->get_rawtower_eval();
        string towernodeCEMC = "TOWER_CALIB_CEMC";
        RawTowerContainer *towersCEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeCEMC.c_str());

        if (towersCEMC){
          m_GEO_Ntowers = 0;
          if (Verbosity() > 0) { cout << "saving EMC towers" << endl; }
          string towergeomnodeCEMC = "TOWERGEOM_CEMC";
          RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeCEMC.c_str());

          if (towergeom){
          if(m_doGeo){
            cout<<"MADE IT TO GEO CEMC LOOP" <<endl;
            RawTowerGeomContainer::ConstRange all_towers = towergeom->get_tower_geometries();
            for (RawTowerGeomContainer::ConstIterator it = all_towers.first; it != all_towers.second; ++it){

              m_GEO_ID= kCEMC;
              m_GEO_towers_iEta[m_GEO_Ntowers] = it->second->get_bineta();
              m_GEO_towers_iPhi[m_GEO_Ntowers] = it->second->get_binphi();
              m_GEO_towers_Eta[m_GEO_Ntowers] = it->second->get_eta();
              m_GEO_towers_Phi[m_GEO_Ntowers] = it->second->get_phi();
              m_GEO_towers_x[m_GEO_Ntowers] = it->second->get_center_x();
              m_GEO_towers_y[m_GEO_Ntowers] = it->second->get_center_y();
              m_GEO_towers_z[m_GEO_Ntowers] = it->second->get_center_z();
              m_GEO_Ntowers++;

            }
            m_geometry_done[kCEMC]=1;
            m_GeoTree->Fill();
            resetGeometryArrays();
            cout<<"FILLING "<<endl;

            cout<<"CLEARING"<<endl;

          }
          cout <<"DID CEMC GEO" <<endl;


          RawTowerContainer::ConstRange begin_end = towersCEMC->getTowers();
          RawTowerContainer::ConstIterator rtiter;
          for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter){
            RawTower* tower = rtiter->second;
            if (tower){
              // min energy cut
              //if (tower->get_energy() < _reco_e_threshold) continue;

              m_tower_iEta_CEMC[m_Ntowers_CEMC] = tower->get_bineta();
              m_tower_iPhi_CEMC[m_Ntowers_CEMC] = tower->get_binphi();
              m_tower_E_CEMC[m_Ntowers_CEMC] = tower->get_energy();

              PHG4Particle* primary = towerevalCEMC->max_truth_primary_particle_by_energy(tower);
              if (primary){ m_tower_trueID_CEMC[m_Ntowers_CEMC]= primary->get_track_id();}
                // gflavor = primary->get_pid();
                // efromtruth = towerevalCEMC->get_energy_contribution(tower, primary);

              else{
                m_tower_trueID_CEMC[m_Ntowers_CEMC] = -10;
              }
              m_Ntowers_CEMC++;
            }
          }
          }
          else{
             cout << PHWHERE << " ERROR: Can't find " << towergeomnodeCEMC << endl;
          // return;
           }
          if (Verbosity() > 0) { cout << "saved\t" << m_Ntowers_CEMC << "\tCEMC towers" << endl; }
           }
        else {  if (Verbosity() > 0) { cout << PHWHERE << " ERROR: Can't find " << towernodeCEMC << endl; }
        // return;
        }
}
  if(m_doCEMC && m_doClusters){
        m_Nclusters_CEMC =0;
        /// Get the raw cluster container
        /// Note: other cluster containers exist as well. Check out the node tree when
        /// you run a simulation
        RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
        CaloRawClusterEval *clusterevalCEMC = m_CEMC_EvalStack->get_rawcluster_eval();

        if (!clusters)
        {
          cout << PHWHERE
               << "CEMC cluster node is missing, can't collect CEMC clusters"
               << endl;
          return;
        }

        /// Get the global vertex to determine the appropriate pseudorapidity of the clusters
        GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
        if (!vertexmap)
        {
          cout << "TannerEval::getEMClusters - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
          assert(vertexmap);  // force quit

          return;
        }

        if (vertexmap->empty())
        {
          cout << "TannerEval::getEMClusters - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
          return;
        }

        GlobalVertex *vtx = vertexmap->begin()->second;
        if (vtx == nullptr)   return;

        /// Trigger emulator
        //CaloTriggerInfo *trigger = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");

        /// Can obtain some trigger information if desired
        // if(trigger)
        //   {
        //     m_E_4x4_CEMC[m_Nclusters_CEMC] = trigger->get_best_EMCal_4x4_E();
        //   }

        RawClusterContainer::ConstRange begin_end = clusters->getClusters();
        RawClusterContainer::ConstIterator clusIter;
        cout <<"GRABBIN CLUSTERS ON CEMC" <<endl;
        /// Loop over the CEMC clusters
        for (clusIter = begin_end.first; clusIter != begin_end.second; ++clusIter)
        {
          /// Get this cluster
          RawCluster *cluster = clusIter->second;

          /// Get cluster characteristics
          /// This helper class determines the photon characteristics
          /// depending on the vertex position
          /// This is important for e.g. eta determination and E_T determination
          CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
          CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);
          m_clusenergy_CEMC[m_Nclusters_CEMC] = E_vec_cluster.mag();
          m_cluseta_CEMC[m_Nclusters_CEMC] = E_vec_cluster.pseudoRapidity();
          m_clustheta_CEMC[m_Nclusters_CEMC] = E_vec_cluster.getTheta();
          m_cluspt_CEMC[m_Nclusters_CEMC] = E_vec_cluster.perp();
          m_clusphi_CEMC[m_Nclusters_CEMC] = E_vec_cluster.getPhi();

          //If pt cut was wanted
          // if (m_cluspt_CEMC[m_Nclusters_CEMC] < m_mincluspt)
          //   continue;
          m_cluspx_CEMC[m_Nclusters_CEMC] =  (E_vec_cluster.perp())* cos(E_vec_cluster.getPhi());
          m_cluspy_CEMC[m_Nclusters_CEMC] = (E_vec_cluster.perp()) * sin(E_vec_cluster.getPhi());
          m_cluspz_CEMC[m_Nclusters_CEMC] = sqrt(((E_vec_cluster.mag()) *(E_vec_cluster.mag())) - (((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))) - (((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))));

          PHG4Particle* primary = clusterevalCEMC->max_truth_primary_particle_by_energy(cluster);
          if (primary){  m_clus_trueID_CEMC[m_Nclusters_CEMC] = primary->get_track_id(); }
          else { m_clus_trueID_CEMC[m_Nclusters_CEMC] = -10; }
          //fill the cluster tree with all CEMC clusters
          m_Nclusters_CEMC++;

    }

  }
  if(m_doHCALIN){


            //Do towers for HCALIN from method in EventEvaluator.h
        if (Verbosity() > 0) { cout << "Grabbing HCALIN event towers" << endl; }
        m_Ntowers_HCALIN = 0;



        CaloRawTowerEval* towerevalHCALIN = m_HCALIN_EvalStack->get_rawtower_eval();
        string towernodeHCALIN = "TOWER_CALIB_HCALIN";
        RawTowerContainer *towersHCALIN = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALIN.c_str());

        if (towersHCALIN){
          m_GEO_Ntowers = 0;
          if (Verbosity() > 0) { cout << "saving EMC towers" << endl; }
          string towergeomnodeHCALIN = "TOWERGEOM_HCALIN";
          RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALIN.c_str());

          if (towergeom){
          if(m_doGeo){
            cout<<"MADE IT TO GEO HCALIN LOOP" <<endl;
            RawTowerGeomContainer::ConstRange all_towers = towergeom->get_tower_geometries();
            for (RawTowerGeomContainer::ConstIterator it = all_towers.first; it != all_towers.second; ++it){

              m_GEO_ID= kHCALIN;
              m_GEO_towers_iEta[m_GEO_Ntowers] = it->second->get_bineta();
              m_GEO_towers_iPhi[m_GEO_Ntowers] = it->second->get_binphi();
              m_GEO_towers_Eta[m_GEO_Ntowers] = it->second->get_eta();
              m_GEO_towers_Phi[m_GEO_Ntowers] = it->second->get_phi();
              m_GEO_towers_x[m_GEO_Ntowers] = it->second->get_center_x();
              m_GEO_towers_y[m_GEO_Ntowers] = it->second->get_center_y();
              m_GEO_towers_z[m_GEO_Ntowers] = it->second->get_center_z();
              m_GEO_Ntowers++;

            }
            m_geometry_done[kHCALIN]=1;
              m_GeoTree->Fill();
              resetGeometryArrays();
            cout<<"FILLING "<<endl;

            cout<<"CLEARING"<<endl;

          }
          cout <<"DID HCALIN GEO" <<endl;


          RawTowerContainer::ConstRange begin_end = towersHCALIN->getTowers();
          RawTowerContainer::ConstIterator rtiter;
          for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter){
            RawTower* tower = rtiter->second;
            if (tower){
              // min energy cut
              //if (tower->get_energy() < _reco_e_threshold) continue;

              m_tower_iEta_HCALIN[m_Ntowers_HCALIN] = tower->get_bineta();
              m_tower_iPhi_HCALIN[m_Ntowers_HCALIN] = tower->get_binphi();
              m_tower_E_HCALIN[m_Ntowers_HCALIN] = tower->get_energy();

              PHG4Particle* primary = towerevalHCALIN->max_truth_primary_particle_by_energy(tower);
              if (primary){ m_tower_trueID_HCALIN[m_Ntowers_HCALIN]= primary->get_track_id();}
                // gflavor = primary->get_pid();
                // efromtruth = towerevalHCALIN->get_energy_contribution(tower, primary);

              else{ m_tower_trueID_HCALIN[m_Ntowers_HCALIN] = -10; }
              m_Ntowers_HCALIN++;
            }
          }
          }
          else{
             cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALIN << endl;
          // return;
           }
          if (Verbosity() > 0) { cout << "saved\t" << m_Ntowers_HCALIN << "\tHCALIN towers" << endl; }
           }
        else {  if (Verbosity() > 0) { cout << PHWHERE << " ERROR: Can't find " << towernodeHCALIN << endl; }
        // return;
        }
}
  if(m_doHCALIN && m_doClusters){
  m_Nclusters_HCALIN =0;
        /// Get the raw cluster container
        /// Note: other cluster containers exist as well. Check out the node tree when
        /// you run a simulation
        RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");
        CaloRawClusterEval *clusterevalHCALIN = m_HCALIN_EvalStack->get_rawcluster_eval();

        if (!clusters)
        {
          cout << PHWHERE
               << "HCALIN cluster node is missing, can't collect HCALIN clusters"
               << endl;
          return;
        }

        /// Get the global vertex to determine the appropriate pseudorapidity of the clusters
        GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
        if (!vertexmap)
        {
          cout << "TannerEval::getEMClusters - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
          assert(vertexmap);  // force quit

          return;
        }

        if (vertexmap->empty())
        {
          cout << "TannerEval::getEMClusters - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
          return;
        }

        GlobalVertex *vtx = vertexmap->begin()->second;
        if (vtx == nullptr)
          return;

        /// Trigger emulator
        // CaloTriggerInfo *trigger = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");
        //
        // /// Can obtain some trigger information if desired
        // if(trigger)
        //   {
        //   //  m_E_4x4_HCALIN[m_Nclusters_HCALIN] = trigger->get_best_EMCal_4x4_E();
        //   }

        RawClusterContainer::ConstRange begin_end = clusters->getClusters();
        RawClusterContainer::ConstIterator clusIter;
        cout <<"GRABBIN CLUSTERS ON HCALIN" <<endl;
        /// Loop over the HCALIN clusters
        for (clusIter = begin_end.first; clusIter != begin_end.second; ++clusIter)
        {
          /// Get this cluster
          RawCluster *cluster = clusIter->second;

          /// Get cluster characteristics
          /// This helper class determines the photon characteristics
          /// depending on the vertex position
          /// This is important for e.g. eta determination and E_T determination
          CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
          CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);
          m_clusenergy_HCALIN[m_Nclusters_HCALIN] = E_vec_cluster.mag();
          m_cluseta_HCALIN[m_Nclusters_HCALIN] = E_vec_cluster.pseudoRapidity();
          m_clustheta_HCALIN[m_Nclusters_HCALIN] = E_vec_cluster.getTheta();
          m_cluspt_HCALIN[m_Nclusters_HCALIN] = E_vec_cluster.perp();
          m_clusphi_HCALIN[m_Nclusters_HCALIN] = E_vec_cluster.getPhi();

          //If pt cut was wanted
          // if (m_cluspt_HCALIN[m_Nclusters_HCALIN] < m_mincluspt)
          //   continue;
          m_cluspx_HCALIN[m_Nclusters_HCALIN] =  (E_vec_cluster.perp())* cos(E_vec_cluster.getPhi());
          m_cluspy_HCALIN[m_Nclusters_HCALIN] = (E_vec_cluster.perp()) * sin(E_vec_cluster.getPhi());
          m_cluspz_HCALIN[m_Nclusters_HCALIN] = sqrt(((E_vec_cluster.mag()) *(E_vec_cluster.mag())) - (((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))) - (((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))));

          PHG4Particle* primary = clusterevalHCALIN->max_truth_primary_particle_by_energy(cluster);
          if (primary){  m_clus_trueID_HCALIN[m_Nclusters_HCALIN] = primary->get_track_id(); }
          else { m_clus_trueID_HCALIN[m_Nclusters_HCALIN] = -10; }
          //fill the cluster tree with all HCALIN clusters
          m_Nclusters_HCALIN++;

    }

  }
  if(m_doHCALOUT){
              //Do towers for HCALOUT from method in EventEvaluator.h
          if (Verbosity() > 0) { cout << "Grabbing HCALOUT event towers" << endl; }
          m_Ntowers_HCALOUT = 0;



          CaloRawTowerEval* towerevalHCALOUT = m_HCALOUT_EvalStack->get_rawtower_eval();
          string towernodeHCALOUT = "TOWER_CALIB_HCALOUT";
          RawTowerContainer *towersHCALOUT = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALOUT.c_str());

          if (towersHCALOUT){
            m_GEO_Ntowers = 0;
            if (Verbosity() > 0) { cout << "saving EMC towers" << endl; }
            string towergeomnodeHCALOUT = "TOWERGEOM_HCALOUT";
            RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALOUT.c_str());

            if (towergeom){
            if(m_doGeo){
              cout<<"MADE IT TO GEO HCALOUT LOOP" <<endl;
              RawTowerGeomContainer::ConstRange all_towers = towergeom->get_tower_geometries();
              for (RawTowerGeomContainer::ConstIterator it = all_towers.first; it != all_towers.second; ++it){

                m_GEO_ID= kHCALOUT;
                m_GEO_towers_iEta[m_GEO_Ntowers] = it->second->get_bineta();
                m_GEO_towers_iPhi[m_GEO_Ntowers] = it->second->get_binphi();
                m_GEO_towers_Eta[m_GEO_Ntowers] = it->second->get_eta();
                m_GEO_towers_Phi[m_GEO_Ntowers] = it->second->get_phi();
                m_GEO_towers_x[m_GEO_Ntowers] = it->second->get_center_x();
                m_GEO_towers_y[m_GEO_Ntowers] = it->second->get_center_y();
                m_GEO_towers_z[m_GEO_Ntowers] = it->second->get_center_z();
                m_GEO_Ntowers++;

              }
              m_geometry_done[kHCALOUT]=1;
                m_GeoTree->Fill();
                resetGeometryArrays();
              cout<<"FILLING "<<endl;

              cout<<"CLEARING"<<endl;

            }
            cout <<"DID HCALOUT GEO" <<endl;


            RawTowerContainer::ConstRange begin_end = towersHCALOUT->getTowers();
            RawTowerContainer::ConstIterator rtiter;
            for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter){
              RawTower* tower = rtiter->second;
              if (tower){
                // min energy cut
                //if (tower->get_energy() < _reco_e_threshold) continue;

                m_tower_iEta_HCALOUT[m_Ntowers_HCALOUT] = tower->get_bineta();
                m_tower_iPhi_HCALOUT[m_Ntowers_HCALOUT] = tower->get_binphi();
                m_tower_E_HCALOUT[m_Ntowers_HCALOUT] = tower->get_energy();

                PHG4Particle* primary = towerevalHCALOUT->max_truth_primary_particle_by_energy(tower);
                if (primary){ m_tower_trueID_HCALOUT[m_Ntowers_HCALOUT]= primary->get_track_id();}
                  // gflavor = primary->get_pid();
                  // efromtruth = towerevalHCALOUT->get_energy_contribution(tower, primary);

                else{ m_tower_trueID_HCALOUT[m_Ntowers_HCALOUT] = -10; }
                m_Ntowers_HCALOUT++;
              }
            }
            }
            else{
               cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALOUT << endl;
            // return;
             }
            if (Verbosity() > 0) { cout << "saved\t" << m_Ntowers_HCALOUT << "\tHCALOUT towers" << endl; }
             }
          else {  if (Verbosity() > 0) { cout << PHWHERE << " ERROR: Can't find " << towernodeHCALOUT << endl; }
          // return;
          }
}
  if(m_doHCALOUT && m_doClusters){
      m_Nclusters_HCALOUT =0;
          /// Get the raw cluster container
          /// Note: other cluster containers exist as well. Check out the node tree when
          /// you run a simulation
          RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");
          CaloRawClusterEval *clusterevalHCALOUT = m_HCALOUT_EvalStack->get_rawcluster_eval();

          if (!clusters)
          {
            cout << PHWHERE
                 << "HCALOUT cluster node is missing, can't collect HCALOUT clusters"
                 << endl;
            return;
          }

          /// Get the global vertex to determine the appropriate pseudorapidity of the clusters
          GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
          if (!vertexmap)
          {
            cout << "TannerEval::getEMClusters - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
            assert(vertexmap);  // force quit

            return;
          }

          if (vertexmap->empty())
          {
            cout << "TannerEval::getEMClusters - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
            return;
          }

          GlobalVertex *vtx = vertexmap->begin()->second;
          if (vtx == nullptr)
            return;

          /// Trigger emulator
          // CaloTriggerInfo *trigger = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");
          //
          // /// Can obtain some trigger information if desired
          // if(trigger)
          //   {
          //   //  m_E_4x4_HCALOUT[m_Nclusters_HCALOUT] = trigger->get_best_EMCal_4x4_E();
          //   }

          RawClusterContainer::ConstRange begin_end = clusters->getClusters();
          RawClusterContainer::ConstIterator clusIter;
          cout <<"GRABBIN CLUSTERS ON HCALOUT" <<endl;
          /// Loop over the HCALOUT clusters
          for (clusIter = begin_end.first; clusIter != begin_end.second; ++clusIter)
          {
            /// Get this cluster
            RawCluster *cluster = clusIter->second;

            /// Get cluster characteristics
            /// This helper class determines the photon characteristics
            /// depending on the vertex position
            /// This is important for e.g. eta determination and E_T determination
            CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
            CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);
            m_clusenergy_HCALOUT[m_Nclusters_HCALOUT] = E_vec_cluster.mag();
            m_cluseta_HCALOUT[m_Nclusters_HCALOUT] = E_vec_cluster.pseudoRapidity();
            m_clustheta_HCALOUT[m_Nclusters_HCALOUT] = E_vec_cluster.getTheta();
            m_cluspt_HCALOUT[m_Nclusters_HCALOUT] = E_vec_cluster.perp();
            m_clusphi_HCALOUT[m_Nclusters_HCALOUT] = E_vec_cluster.getPhi();

            //If pt cut was wanted
            // if (m_cluspt_HCALOUT[m_Nclusters_HCALOUT] < m_mincluspt)
            //   continue;
            m_cluspx_HCALOUT[m_Nclusters_HCALOUT] =  (E_vec_cluster.perp())* cos(E_vec_cluster.getPhi());
            m_cluspy_HCALOUT[m_Nclusters_HCALOUT] = (E_vec_cluster.perp()) * sin(E_vec_cluster.getPhi());
            m_cluspz_HCALOUT[m_Nclusters_HCALOUT] = sqrt(((E_vec_cluster.mag()) *(E_vec_cluster.mag())) - (((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))) - (((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))));

            PHG4Particle* primary = clusterevalHCALOUT->max_truth_primary_particle_by_energy(cluster);
            if (primary){  m_clus_trueID_HCALOUT[m_Nclusters_HCALOUT] = primary->get_track_id(); }
            else { m_clus_trueID_HCALOUT[m_Nclusters_HCALOUT] = -10; }
            //fill the cluster tree with all HCALOUT clusters
            m_Nclusters_HCALOUT++;

      }
  }


  if(m_doTracks){
      m_nTracks=0;
      m_nProjections = 0;
      /// SVTX tracks node
      SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

      if (!trackmap)
      {
        cout << PHWHERE
             << "SvtxTrackMap node is missing, can't collect tracks"
             << endl;
        return;
      }

      /// EvalStack for truth track matching
      /// Get the track evaluator
      SvtxTrackEval *trackeval = m_SVTX_EvalStack->get_track_eval();

      /// Get the range for primary tracks
      PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

      if (Verbosity() > 1)
      {
        cout << "Get the SVTX tracks" << endl;
      }

    //  m_nTracks = 0;
      for (SvtxTrackMap::Iter iter = trackmap->begin();
           iter != trackmap->end();
           ++iter)
      {
        SvtxTrack *track = iter->second;

        /// Get the reconstructed track info
        m_tr_px[m_nTracks]  = track->get_px();
        m_tr_py[m_nTracks]  = track->get_py();
        m_tr_pz[m_nTracks]  = track->get_pz();
        m_tr_p[m_nTracks]  = sqrt(m_tr_px[m_nTracks]  * m_tr_px[m_nTracks]  + m_tr_py[m_nTracks]  * m_tr_py[m_nTracks]  + m_tr_pz[m_nTracks]  * m_tr_pz[m_nTracks] );

        m_tr_pt[m_nTracks]  = sqrt(m_tr_px[m_nTracks]  * m_tr_px[m_nTracks]  + m_tr_py[m_nTracks]  * m_tr_py[m_nTracks] );

        // Make some cuts on the track to clean up sample
        if (m_tr_pt[m_nTracks]  < 0.5)   continue;

        m_tr_phi[m_nTracks]  = track->get_phi();
        m_tr_eta[m_nTracks]  = track->get_eta();

        m_charge[m_nTracks]  = track->get_charge();
        m_chisq[m_nTracks]  = track->get_chisq();
        m_ndf[m_nTracks]  = track->get_ndf();
        m_dca[m_nTracks]  = track->get_dca();
        m_tr_x[m_nTracks]  = track->get_x();
        m_tr_y[m_nTracks]  = track->get_y();
        m_tr_z[m_nTracks]  = track->get_z();
        m_tr_ID[m_nTracks] = track->get_id();

        /// Get truth track info that matches this reconstructed track
        PHG4Particle *truthtrack = trackeval->max_truth_particle_by_nclusters(track);
        m_truth_is_primary[m_nTracks] = truthinfo->is_primary(truthtrack);

        m_truthtrackpx[m_nTracks]  = truthtrack->get_px();
        m_truthtrackpy[m_nTracks]  = truthtrack->get_py();
        m_truthtrackpz[m_nTracks]  = truthtrack->get_pz();
        m_truthtrackp[m_nTracks]  = sqrt(m_truthtrackpx[m_nTracks] * m_truthtrackpx[m_nTracks] + m_truthtrackpy[m_nTracks] * m_truthtrackpy[m_nTracks] + m_truthtrackpz[m_nTracks] * m_truthtrackpz[m_nTracks]);

        m_truthtracke[m_nTracks]  = truthtrack->get_e();

        m_truthtrackpt[m_nTracks]  = sqrt(m_truthtrackpx[m_nTracks]  * m_truthtrackpx[m_nTracks]  + m_truthtrackpy[m_nTracks]  * m_truthtrackpy[m_nTracks] );
        m_truthtrackphi[m_nTracks]  = atan(m_truthtrackpy[m_nTracks]  / m_truthtrackpx[m_nTracks] );
        m_truthtracketa[m_nTracks]  = atanh(m_truthtrackpz[m_nTracks]  / m_truthtrackp[m_nTracks] );
        m_truthtrackpid[m_nTracks]  = truthtrack->get_pid();
        m_truthtrackID[m_nTracks]   = truthtrack->get_track_id();
        if(m_doProjections){

          for(SvtxTrack::ConstStateIter = trkstates = track->begin_states(); trkstates != track->end_states(); ++trkstates){

          }


        }






        m_nTracks++;

      }
  }

  m_nMCparts = 0;
  if(m_doMCparts){

    PHG4TruthInfoContainer* truthinfocontainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfocontainer)
    {
      if (Verbosity() > 0)
      {
        cout << "saving MC particles" << endl;
      }
      //GetParticleRange for all particles
      //GetPrimaryParticleRange for primary particles
      PHG4TruthInfoContainer::ConstRange range = truthinfocontainer->GetParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first; truth_itr != range.second; ++truth_itr)
      {
        PHG4Particle* g4particle = truth_itr->second;
        if (!g4particle) continue;

        int mcSteps = 0;
        PHG4Particle* g4particleMother = truth_itr->second;
        if (g4particle->get_parent_id() != 0)
        {
          while (g4particleMother->get_parent_id() != 0)
          {
            g4particleMother = truthinfocontainer->GetParticle(g4particleMother->get_parent_id());
            if (g4particleMother == NULL) break;
            mcSteps += 1;
          }
        }
        if (mcSteps > m_depth_MCstack) continue;

        // evaluating true primary vertex
        if (m_doVertex&& m_nMCparts == 0)
        {
          PHG4VtxPoint* vtx = truthinfocontainer->GetVtx(g4particle->get_vtx_id());
          if (vtx)
          {
            m_vertex_true_x = vtx->get_x();
            m_vertex_true_y = vtx->get_y();
            m_vertex_true_z = vtx->get_z();
          }
        }

        // in case of all MC particles, make restrictions on the secondary selection
        // if(g4particle->get_track_id()<0 && g4particle->get_e()<0.5) continue;
        // primary (g4particle->get_parent_id() == 0) selection via:
        // if(gtrackID < 0) continue;

        //using the e threshold also for the truth particles gets rid of all the low energy secondary particles
        if (g4particle->get_e() < m_reco_e_threshold) continue;

        m_MCpart_ID[m_nMCparts] = g4particle->get_track_id();
        m_MCpart_ID_parent[m_nMCparts] = g4particle->get_parent_id();
        m_MCpart_PDG[m_nMCparts] = g4particle->get_pid();
        m_MCpart_E[m_nMCparts] = g4particle->get_e();
        m_MCpart_px[m_nMCparts] = g4particle->get_px();
        m_MCpart_py[m_nMCparts] = g4particle->get_py();
        m_MCpart_pz[m_nMCparts] = g4particle->get_pz();
        //BCID added for G4Particle --  HEPMC particle matching
        m_MCpart_BCID[m_nMCparts] = g4particle->get_barcode();
        // TVector3 projvec(m_MCpart_px[0],m_MCpart_py[0],m_MCpart_pz[0]);
        // float projeta = projvec.Eta();
        m_nMCparts++;
      }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << m_nMCparts << "\tMC particles" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << endl;
      }
      return;
    }
}
  m_nHEPMCparts = 0;
  if(m_doHEPMCparts){

        /// Get the node from the node tree
        PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

        /// If the node was not properly put on the tree, return
        if (!hepmceventmap)
        {
          cout << PHWHERE
               << "HEPMC event map node is missing, can't collected HEPMC truth particles"
               << endl;
          return;
        }

        /// Could have some print statements for debugging with verbosity
        if (Verbosity() > 1)
        {
          cout << "Getting HEPMC truth particles " << endl;
        }

        /// You can iterate over the number of events in a hepmc event
        /// for pile up events where you have multiple hard scatterings per bunch crossing

        for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
             eventIter != hepmceventmap->end();
             ++eventIter)
        {
          /// Get the event
          PHHepMCGenEvent *hepmcevent = eventIter->second;

          if (hepmcevent)
          {
            /// Get the event characteristics, inherited from HepMC classes
            HepMC::GenEvent *truthevent = hepmcevent->getEvent();
            if (!truthevent)
            {
              cout << PHWHERE
                   << "no evt pointer under phhepmvgeneventmap found "
                   << endl;
              return;
            }

            /// Get the parton info
            HepMC::PdfInfo *pdfinfo = truthevent->pdf_info();

            /// Get the parton info as determined from HEPMC
            m_partid1[m_nHEPMCparts] = pdfinfo->id1();
            m_partid2[m_nHEPMCparts] = pdfinfo->id2();
            m_x1 = pdfinfo->x1();
            m_x2 = pdfinfo->x2();

            /// Are there multiple partonic intercations in a p+p event
            m_mpi[m_nHEPMCparts] = truthevent->mpi();

            /// Get the PYTHIA signal process id identifying the 2-to-2 hard process
            m_process_id = truthevent->signal_process_id();

            if (Verbosity() > 2)
            {
              cout << " Iterating over an event" << endl;
            }
            /// Loop over all the truth particles and get their information
            for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin();
                 iter != truthevent->particles_end();
                 ++iter)
            {
              /// Get each pythia particle characteristics
              m_truthenergy[m_nHEPMCparts] = (*iter)->momentum().e();
              m_truthpid[m_nHEPMCparts] = (*iter)->pdg_id();

              m_trutheta[m_nHEPMCparts] = (*iter)->momentum().pseudoRapidity();
              m_truthphi[m_nHEPMCparts] = (*iter)->momentum().phi();
              m_truthpx[m_nHEPMCparts] = (*iter)->momentum().px();
              m_truthpy[m_nHEPMCparts] = (*iter)->momentum().py();
              m_truthpz[m_nHEPMCparts] = (*iter)->momentum().pz();
              m_truthpt[m_nHEPMCparts] = sqrt(m_truthpx[m_nHEPMCparts] * m_truthpx[m_nHEPMCparts] + m_truthpy[m_nHEPMCparts] * m_truthpy[m_nHEPMCparts]);

              /// Fill the truth tree
              m_nHEPMCparts++;

            }
          }
        }

  }

  m_EventTree->Fill();
//  m_EventTree->Write();
  if (Verbosity() > 0){ cout << "Resetting buffer ..." << endl;}
  resetBuffer();
  if (Verbosity() > 0)
  {
    cout << "EventEvaluator buffer reset" << endl;
  }
  return;
}

int TannerEval::End(PHCompositeNode *topNode){
    if (Verbosity() > 0) { cout << "Ending TTANEVAL package" << endl; }

    m_outfile->cd();
    m_EventTree->Write();

    m_outfile->Close();
  //  delete m_outfile;

    if(m_doGeo){
      m_outfile_geometry->cd();
      m_GeoTree->Write();

      m_outfile_geometry->Close();
  //    delete m_outfile_geometry;
    }

    if (Verbosity() > 0)  {
      cout << "========================= " << Name() << "::End() ============================" << endl;
      cout << " " << m_iEvent << " events of output written to: " << m_filename<< endl;
      cout << "===========================================================================" << endl;
    }

   if(m_HCALIN_EvalStack) delete m_HCALIN_EvalStack;
   if(m_HCALOUT_EvalStack) delete m_HCALOUT_EvalStack;
   if(m_CEMC_EvalStack) delete m_CEMC_EvalStack;
   if(m_SVTX_EvalStack) delete m_SVTX_EvalStack;

  cout<< "=============================GOODBYE===================================="<<endl;
    return Fun4AllReturnCodes::EVENT_OK;

  }


int TannerEval::GetProjectionIndex(std::string projname)
  {
    if (projname.find("MVTX") != std::string::npos)
      return 1;
    else if (projname.find("INTT") != std::string::npos)
      return 2;
    else if (projname.find("TPC") != std::string::npos)
      return 3;
    else if (projname.find("MICROMEGAS") != std::string::npos)
      return 4;
    else if (projname.find("HCALIN") != std::string::npos)
      return 5;
    else if (projname.find("CEMC") != std::string::npos)
      return 6;
    else if (projname.find("HCALOUT") != std::string::npos)
      return 7;
    else
      return -1;

  }

std::string TannerEval::GetProjectionNameFromIndex(int projindex)
  {
    switch (projindex)
    {
      case 1:
        return "MVTX";
      case 2:
        return "INTT";
      case 3:
        return "TPC";
      case 4:
        return "MICROMEGAS";
      case 5:
        return "HCALIN";
      case 6:
        return "CEMC";
      case 7:
        return "HCALOUT";
      default:
        return "NOTHING";
    }
  }
void TannerEval::resetBuffer(){
if(m_doEventInfo){
  m_cross_section =0;
  m_event_weight = 0;
  m_n_generator_accepted = 0;
}
if(m_doHEPMCparts){
  m_nHEPMCparts=0;
  for(int spot=0; spot<m_maxHEPMCparts;spot++){
    m_partid1[spot]=0;
    m_partid2[spot]=0;
    m_mpi[spot]=0;
    m_truthenergy[spot]=0;
    m_trutheta[spot]=0;
    m_truthphi[spot]=0;
    m_truthpx[spot]=0;
    m_truthpy[spot]=0;
    m_truthpz[spot]=0;
    m_truthpt[spot]=0;
    m_truthp[spot]=0;
    m_truthpid[spot]=0;

  }
}
if(m_doMCparts){
  m_nMCparts=0;
  for(int spot=0; spot<m_maxMCparts; spot++){
    m_MCpart_ID[spot]=0;
    m_MCpart_ID_parent[spot]=0;
    m_MCpart_PDG[spot]=0;
    m_MCpart_E[spot]=0;
    m_MCpart_px[spot]=0;
    m_MCpart_py[spot]=0;
    m_MCpart_pz[spot]=0;
    m_MCpart_BCID[spot]=0;

  }
}
if(m_doTracks){

  m_nTracks=0;
  for(int spot=0; spot<m_maxNtracks;spot++){
    m_tr_px[spot]=0;
    m_tr_py[spot]=0;
    m_tr_pz[spot]=0;
    m_tr_p[spot]=0;
    m_tr_pt[spot]=0;
    m_tr_phi[spot]=0;
    m_tr_eta[spot]=0;
    m_charge[spot]=0;
    m_chisq[spot]=0;
    m_ndf[spot]=0;
    m_dca[spot]=0;
    m_tr_x[spot]=0;
    m_tr_y[spot]=0;
    m_tr_z[spot]=0;
    m_truth_is_primary[spot]=0;
    m_truthtrackpx[spot]=0;
    m_truthtrackpy[spot]=0;
    m_truthtrackpz[spot]=0;
    m_truthtrackp[spot]=0;
    m_truthtracke[spot]=0;
    m_truthtrackpt[spot]=0;
    m_truthtrackphi[spot]=0;
    m_truthtracketa[spot]=0;
    m_truthtrackpid[spot]=0;

  }

  if(m_doProjections){
    m_nProjections=0;
    for(int spot = 0; spot<m_maxNProjections; spot++){
      m_track_ProjLayer[spot] = -1;
      m_track_ProjTrackID[spot] = 0;
      m_track_TLP_x[spot] = 0;
      m_track_TLP_y[spot] = 0;
      m_track_TLP_z[spot] = 0;
      m_track_TLP_t[spot] = 0;
      m_track_TLP_true_x[spot] = 0;
      m_track_TLP_true_y[spot] = 0;
      m_track_TLP_true_z[spot] = 0;
      m_track_TLP_true_t[spot] = 0;
    }
  }

}
if(m_doCEMC){
    m_Ntowers_CEMC=0;
  for(int spot=0; spot<m_maxNtowersCentral;spot++){

  m_tower_E_CEMC[spot]=0;
  m_tower_iEta_CEMC[spot]=0;
  m_tower_iPhi_CEMC[spot]=0;
  m_tower_trueID_CEMC[spot]=0;
}
  if(m_doClusters){
    m_Nclusters_CEMC=0;
  for(int spot=0; spot<m_maxNclusters;spot++){
    m_clusenergy_CEMC[spot]=0;
    m_cluseta_CEMC[spot]=0;
    m_clustheta_CEMC[spot]=0;
    m_clusphi_CEMC[spot]=0;
    m_cluspt_CEMC[spot]=0;
    m_cluspy_CEMC[spot]=0;
    m_cluspx_CEMC[spot]=0;
    m_cluspz_CEMC[spot]=0;
    m_E_4x4_CEMC[spot]=0;
  }

}
}
if(m_doHCALIN){
    m_Ntowers_HCALIN=0;
  for(int spot=0; spot<m_maxNtowersCentral;spot++){

  m_tower_E_HCALIN[spot]=0;
  m_tower_iEta_HCALIN[spot]=0;
  m_tower_iPhi_HCALIN[spot]=0;
  m_tower_trueID_HCALIN[spot]=0;
}
  if(m_doClusters){
    m_Nclusters_HCALIN=0;
  for(int spot=0; spot<m_maxNclusters;spot++){
    m_clusenergy_HCALIN[spot]=0;
    m_cluseta_HCALIN[spot]=0;
    m_clustheta_HCALIN[spot]=0;
    m_clusphi_HCALIN[spot]=0;
    m_cluspt_HCALIN[spot]=0;
    m_cluspy_HCALIN[spot]=0;
    m_cluspx_HCALIN[spot]=0;
    m_cluspz_HCALIN[spot]=0;

  }

}
}
if(m_doHCALOUT){
    m_Ntowers_HCALOUT=0;
  for(int spot=0; spot<m_maxNtowersCentral;spot++){

  m_tower_E_HCALOUT[spot]=0;
  m_tower_iEta_HCALOUT[spot]=0;
  m_tower_iPhi_HCALOUT[spot]=0;
  m_tower_trueID_HCALOUT[spot]=0;
}
  if(m_doClusters){
    m_Nclusters_HCALOUT=0;
  for(int spot=0; spot<m_maxNclusters;spot++){
    m_clusenergy_HCALOUT[spot]=0;
    m_cluseta_HCALOUT[spot]=0;
    m_clustheta_HCALOUT[spot]=0;
    m_clusphi_HCALOUT[spot]=0;
    m_cluspt_HCALOUT[spot]=0;
    m_cluspy_HCALOUT[spot]=0;
    m_cluspx_HCALOUT[spot]=0;
    m_cluspz_HCALOUT[spot]=0;

  }

}
}
if(m_doVertex){
  m_vertex_x = -1000;
  m_vertex_y = -1000;
  m_vertex_z = -1000;
  m_vertex_NCont = 0;
  m_vertex_true_x = -1000;
  m_vertex_true_y = -1000;
  m_vertex_true_z = -1000;
  if (Verbosity() > 0){ cout << "\t... vertex variables reset" << endl;}
}
if(m_doHits){
  m_nHitsLayers = 0;
  for (int ihit = 0; ihit < m_maxNHits; ihit++)
  {
    m_hits_layerID[ihit] = 0;
    m_hits_trueID[ihit] = 0;
    m_hits_x[ihit] = 0;
    m_hits_y[ihit] = 0;
    m_hits_z[ihit] = 0;
    m_hits_t[ihit] = 0;
  }
  if (Verbosity() > 0){ cout << "\t... hit variables reset" << endl;}
}
}
void TannerEval::resetGeometryArrays(){
  for (int igeo = 0; igeo < m_maxNtowers_GEO; igeo++)
    {
      m_GEO_towers_iEta[igeo] = -10000;
      m_GEO_towers_iPhi[igeo] = -10000;
      m_GEO_towers_Eta[igeo] = -10000;
      m_GEO_towers_Phi[igeo] = -10000;
      m_GEO_towers_x[igeo] = -10000;
      m_GEO_towers_y[igeo] = -10000;
      m_GEO_towers_z[igeo] = -10000;
    }
     m_GEO_Ntowers =-1;
     m_GEO_ID=0;
}

//
// void TannerEval::WaterSeeds(){
//
//     m_outfile = new TFile();
//     m_outfile_geometry = new TFile();
//
//     //HEPMC PARTICLE VARIBLES
//     m_nHEPMCparts = -99;//counter
//
//     m_partid1 = -99;
//     m_partid2 = -99;
//     m_x1 = -99;
//     m_x2 = -99;
//     m_mpi = -99;
//     m_process_id = -99;
//     m_truthenergy = -99;
//     m_trutheta = -99;
//     m_truthphi = -99;
//     m_truthpx = -99;
//     m_truthpy = -99;
//     m_truthpz = -99;
//     m_truthpt = -99;
//     m_truthp = -99;
//     m_truthpid = -99;
//
//
//
//     //PH4G PARTICLE VARIBLES
//     m_nMCparts = -99; //Counter
//
//     m_truthpx_PHG4 = -99;
//     m_truthpy_PHG4 = -99;
//     m_truthpz_PHG4 = -99;
//     m_truthpt_PHG4 = -99;
//     m_truthp_PHG4 = -99;
//     m_truthenergy_PHG4 = -99;
//     m_trutheta_PHG4 = -99;
//     m_truthphi_PHG4 = -99;
//     m_truthpid_PHG4 = -99;
//
//
//
//
//     // Track TTree variables
//     m_nTracks = -99; //counter
//
//     m_tr_px = -99;
//     m_tr_py = -99;
//     m_tr_pz = -99;
//     m_tr_p = -99;
//     m_tr_pt = -99;
//     m_tr_phi = -99;
//     m_tr_eta = -99;
//     m_charge = -99;
//     m_chisq = -99;
//     m_ndf = -99;
//     m_dca = -99;
//     m_tr_x = -99;
//     m_tr_y = -99;
//     m_tr_z = -99;
//     m_truth_is_primary = -99;
//     m_truthtrackpx = -99;
//     m_truthtrackpy = -99;
//     m_truthtrackpz = -99;
//     m_truthtrackp = -99;
//     m_truthtracke = -99;
//     m_truthtrackpt = -99;
//     m_truthtrackphi = -99;
//     m_truthtracketa = -99;
//     m_truthtrackpid = -99;
//
//     /// Cluster and tower variables
//
//
//     //CEMC
//     m_Nclusters_CEMC = -99; //cluster counter
//
//     m_clusenergy_CEMC = -99;
//     m_cluseta_CEMC = -99;
//     m_clustheta_CEMC = -99;
//     m_cluspt_CEMC = -99;
//     m_clusphi_CEMC = -99;
//     m_cluspx_CEMC = -99;
//     m_cluspy_CEMC = -99;
//     m_cluspz_CEMC = -99;
//     m_E_4x4_CEMC = -99;
//     m_clus_trueID_CEMC = -99;
//
//
//     m_Ntowers_CEMC = -99; // tower counter
//
//     m_tower_E_CEMC = -99;
//     m_tower_iEta_CEMC = -99;
//     m_tower_iPhi_CEMC = -99;
//     m_tower_trueID_CEMC = -99;
//
//     //HCALIN
//
//     m_Nclusters_HCALIN = -99; //cluster counter
//
//     m_clusenergy_HCALIN = -99;
//     m_cluseta_HCALIN = -99;
//     m_clustheta_HCALIN = -99;
//     m_cluspt_HCALIN = -99;
//     m_clusphi_HCALIN = -99;
//     m_cluspx_HCALIN = -99;
//     m_cluspy_HCALIN = -99;
//     m_cluspz_HCALIN = -99;
//     m_E_4x4_HCALIN = -99;
//     m_clus_trueID_HCALIN= -99;
//
//
//     m_Ntowers_HCALIN = -99; // tower counter
//
//     m_tower_E_HCALIN = -99;
//     m_tower_iEta_HCALIN = -99;
//     m_tower_iPhi_HCALIN = -99;
//     m_tower_trueID_HCALIN = -99;
//
//
//     //HCALOUT
//     m_Nclusters_HCALOUT = -99; //cluster counter
//
//     m_clusenergy_HCALOUT = -99;
//     m_cluseta_HCALOUT = -99;
//     m_clustheta_HCALOUT = -99;
//     m_cluspt_HCALOUT = -99;
//     m_clusphi_HCALOUT = -99;
//     m_cluspx_HCALOUT = -99;
//     m_cluspy_HCALOUT = -99;
//     m_cluspz_HCALOUT = -99;
//     m_E_4x4_HCALOUT = -99;
//     m_clus_trueID_HCALOUT= -99;
//
//
//     m_Ntowers_HCALOUT = -99; // tower counter
//
//     m_tower_E_HCALOUT = -99;
//     m_tower_iEta_HCALOUT = -99;
//     m_tower_iPhi_HCALOUT = -99;
//     m_tower_trueID_HCALOUT = -99;
//
//     m_GEO_Ntowers=-99;
//     m_GEO_ID=-99;
//     m_GEO_towers_iEta = 99;
//     m_GEO_towers_iPhi = 99;
//     m_GEO_towers_Eta = -99;
//     m_GEO_towers_Phi = -99;
//     m_GEO_towers_x = -99;
//     m_GEO_towers_y = -99;
//     m_GEO_towers_z = -99;
// }
//
// void TannerEval::GrowTrees(){
//
//
//
//
//   if (m_doGeo){
//
//       m_GeoTree = new TTree("GeoTree", "GeoTree");
//       m_GeoTree->Branch("m_GEO_ID", &m_GEO_ID, "nCaloHits/I");
//       m_GeoTree->Branch("m_GEO_Ntowers", &m_GEO_Ntowers, "m_GEO_Ntowers/I");
//       m_GeoTree->Branch("m_GEO_towers_iEta", m_GEO_towers_iEta,"calo_towers_iEta[m_GEO_Ntowers]/I");
//       m_GeoTree->Branch("m_GEO_towers_iPhi", m_GEO_towers_iPhi,"calo_towers_iPhi[m_GEO_Ntowers]/I");
//       m_GeoTree->Branch("m_GEO_towers_Eta", m_GEO_towers_Eta,"calo_towers_Eta[m_GEO_Ntowers]/D");
//       m_GeoTree->Branch("m_GEO_towers_Phi", m_GEO_towers_Phi, "calo_towers_Phi[m_GEO_Ntowers]/D");
//       m_GeoTree->Branch("m_GEO_towers_x", m_GEO_towers_x,"calo_towers_x[m_GEO_Ntowers]/D");
//       m_GeoTree->Branch("m_GEO_towers_y", m_GEO_towers_y,"calo_towers_y[m_GEO_Ntowers]/D");
//       m_GeoTree->Branch("m_GEO_towers_z", m_GEO_towers_z,"calo_towers_z[m_GEO_Ntowers]/D");
//
//
//     }
//
//     m_EventTree = new TTree("EventTree", "EventTree");
//
//  if (m_doEventInfo){
//       // Event level info. This isn't the most efficient way to store this info, but it's straightforward
//       // within the structure of the class, so the size is small compared to the rest of the output.
//       m_EventTree->Branch("m_cross_section", &m_cross_section, "m_cross_section/D");
//       m_EventTree->Branch("m_event_weight", &m_event_weight, "m_event_weight/D");
//       m_EventTree->Branch("m_n_generator_accepted", &m_n_generator_accepted, "m_n_generator_accepted/D");
//     }
//
//  if (m_doTracks){
//       m_EventTree->Branch("m_nTracks", &m_nTracks ,"m_nTracks/I");
//       m_EventTree->Branch("m_tr_px", m_tr_px, "m_tr_px[m_nTracks]/D");
//       m_EventTree->Branch("m_tr_py", m_tr_py, "m_tr_py[m_nTracks]/D");
//       m_EventTree->Branch("m_tr_pz", m_tr_pz, "m_tr_pz[m_nTracks]/D");
//       m_EventTree->Branch("m_tr_p", m_tr_p, "m_tr_p[m_nTracks]/D");
//       m_EventTree->Branch("m_tr_pt", m_tr_pt, "m_tr_pt[m_nTracks]/D");
//       m_EventTree->Branch("m_tr_phi", m_tr_phi, "m_tr_phi[m_nTracks]/D");
//       m_EventTree->Branch("m_tr_eta", m_tr_eta, "m_tr_eta[m_nTracks]/D");
//       m_EventTree->Branch("m_charge", m_charge, "m_charge[m_nTracks]/I");
//       m_EventTree->Branch("m_chisq", m_chisq, "m_chisq[m_nTracks]/D");
//       m_EventTree->Branch("m_ndf", m_ndf, "m_ndf[m_nTracks]/I");
//       m_EventTree->Branch("m_dca", m_dca, "m_dca[m_nTracks]/D");
//       m_EventTree->Branch("m_tr_x", m_tr_x, "m_tr_x[m_nTracks]/D");
//       m_EventTree->Branch("m_tr_y", m_tr_y, "m_tr_y[m_nTracks]/D");
//       m_EventTree->Branch("m_tr_z", m_tr_z, "m_tr_z[m_nTracks]/D");
//       m_EventTree->Branch("m_truth_is_primary", m_truth_is_primary, "m_truth_is_primary[m_nTracks]/I");
//       m_EventTree->Branch("m_truthtrackpx", m_truthtrackpx, "m_truthtrackpx[m_nTracks]/D");
//       m_EventTree->Branch("m_truthtrackpy", m_truthtrackpy, "m_truthtrackpy[m_nTracks]/D");
//       m_EventTree->Branch("m_truthtrackpz", m_truthtrackpz, "m_truthtrackpz[m_nTracks]/D");
//       m_EventTree->Branch("m_truthtrackp", m_truthtrackp, "m_truthtrackp[m_nTracks]/D");
//       m_EventTree->Branch("m_truthtracke", m_truthtracke, "m_truthtracke[m_nTracks]/D");
//       m_EventTree->Branch("m_truthtrackpt", m_truthtrackpt, "m_truthtrackpt[m_nTracks]/D");
//       m_EventTree->Branch("m_truthtrackphi", m_truthtrackphi, "m_truthtrackphi[m_nTracks]/D");
//       m_EventTree->Branch("m_truthtracketa", m_truthtracketa, "m_truthtracketa[m_nTracks]/D");
//       m_EventTree->Branch("m_truthtrackpid", m_truthtrackpid, "m_truthtrackpid[m_nTracks]/I");
//
//     }
//     if (m_doHCALIN){
//     m_EventTree->Branch("m_Nclusters_HCALIN", &m_Nclusters_HCALIN,"m_Nclusters_HCALIN/I");
//     m_EventTree->Branch("m_clusenergy_HCALIN", m_clusenergy_HCALIN, "m_clusenergy_HCALIN[m_Nclusters_HCALIN]/D");
//     m_EventTree->Branch("m_cluseta_HCALIN", m_cluseta_HCALIN, "m_cluseta_HCALIN[m_Nclusters_HCALIN]/D");
//     m_EventTree->Branch("m_clustheta_HCALIN", m_clustheta_HCALIN, "m_clustheta_HCALIN[m_Nclusters_HCALIN]/D");
//     m_EventTree->Branch("m_cluspt_HCALIN", m_cluspt_HCALIN, "m_cluspt_HCALIN[m_Nclusters_HCALIN]/D");
//     m_EventTree->Branch("m_clusphi_HCALIN", m_clusphi_HCALIN, "m_clusphi_HCALIN[m_Nclusters_HCALIN]/D");
//     m_EventTree->Branch("m_cluspx_HCALIN", m_cluspx_HCALIN, "m_cluspx_HCALIN[m_Nclusters_HCALIN]/D");
//     m_EventTree->Branch("m_cluspy_HCALIN", m_cluspy_HCALIN, "m_cluspy_HCALIN[m_Nclusters_HCALIN]/D");
//     m_EventTree->Branch("m_cluspz_HCALIN", m_cluspz_HCALIN, "m_cluspz_HCALIN[m_Nclusters_HCALIN]/D");
//     m_EventTree->Branch("m_E_4x4_HCALIN", m_E_4x4_HCALIN, "m_E_4x4_HCALIN[m_Nclusters_HCALIN]/D");
//     m_EventTree->Branch("m_clus_trueID_HCALIN", m_clus_trueID_HCALIN, "m_clus_trueID_HCALIN[m_Nclusters_HCALIN]/D");
//
//     m_EventTree->Branch("m_Ntowers_HCALIN",&m_Ntowers_HCALIN,"m_Ntowers_HCALIN/I");
//     m_EventTree->Branch("m_tower_E_HCALIN",m_tower_E_HCALIN,"m_tower_E_HCALIN[m_Ntowers_HCALIN]/D");
//     m_EventTree->Branch("m_tower_iEta_HCALIN",m_tower_iEta_HCALIN,"m_tower_iEta_HCALIN[m_Ntowers_HCALIN]/I");
//     m_EventTree->Branch("m_tower_iPhi_HCALIN",m_tower_iPhi_HCALIN,"m_tower_iPhi_HCALIN[m_Ntowers_HCALIN]/I");
//     m_EventTree->Branch("m_tower_trueID_HCALIN",m_tower_trueID_HCALIN,"m_tower_trueID_HCALIN[m_Ntowers_HCALIN]/I");
//
//
//     }
//     if (m_doCEMC){
//     m_EventTree->Branch("m_Nclusters_CEMC", &m_Nclusters_CEMC,"m_Nclusters_CEMC/I");
//     m_EventTree->Branch("m_clusenergy_CEMC", m_clusenergy_CEMC, "m_clusenergy_CEMC[m_Nclusters_CEMC]/D");
//     m_EventTree->Branch("m_cluseta_CEMC", m_cluseta_CEMC, "m_cluseta_CEMC[m_Nclusters_CEMC]/D");
//     m_EventTree->Branch("m_clustheta_CEMC", m_clustheta_CEMC, "m_clustheta_CEMC[m_Nclusters_CEMC]/D");
//     m_EventTree->Branch("m_cluspt_CEMC", m_cluspt_CEMC, "m_cluspt_CEMC[m_Nclusters_CEMC]/D");
//     m_EventTree->Branch("m_clusphi_CEMC", m_clusphi_CEMC, "m_clusphi_CEMC[m_Nclusters_CEMC]/D");
//     m_EventTree->Branch("m_cluspx_CEMC", m_cluspx_CEMC, "m_cluspx_CEMC[m_Nclusters_CEMC]/D");
//     m_EventTree->Branch("m_cluspy_CEMC", m_cluspy_CEMC, "m_cluspy_CEMC[m_Nclusters_CEMC]/D");
//     m_EventTree->Branch("m_cluspz_CEMC", m_cluspz_CEMC, "m_cluspz_CEMC[m_Nclusters_CEMC]/D");
//     m_EventTree->Branch("m_E_4x4_CEMC", m_E_4x4_CEMC, "m_E_4x4_CEMC[m_Nclusters_CEMC]/D");
//     m_EventTree->Branch("m_clus_trueID_CEMC", m_clus_trueID_CEMC, "m_clus_trueID_CEMC[m_Nclusters_CEMC]/D");
//
//     m_EventTree->Branch("m_Ntowers_CEMC",&m_Ntowers_CEMC,"m_Ntowers_CEMC/I");
//     m_EventTree->Branch("m_tower_E_CEMC",m_tower_E_CEMC,"m_tower_E_CEMC[m_Ntowers_CEMC]/D");
//     m_EventTree->Branch("m_tower_iEta_CEMC",m_tower_iEta_CEMC,"m_tower_iEta_CEMC[m_Ntowers_CEMC]/I");
//     m_EventTree->Branch("m_tower_iPhi_CEMC",m_tower_iPhi_CEMC,"m_tower_iPhi_CEMC[m_Ntowers_CEMC]/I");
//     m_EventTree->Branch("m_tower_trueID_CEMC",m_tower_trueID_CEMC,"m_tower_trueID_CEMC[m_Ntowers_CEMC]/I");
//     }
//     if (m_doHCALOUT){
//     m_EventTree->Branch("m_Nclusters_HCALOUT", *m_Nclusters_HCALOUT,"m_Nclusters_HCALOUT/I");
//     m_EventTree->Branch("m_clusenergy_HCALOUT", m_clusenergy_HCALOUT, "m_clusenergy_HCALOUT[m_Nclusters_HCALOUT]/D");
//     m_EventTree->Branch("m_cluseta_HCALOUT", m_cluseta_HCALOUT, "m_cluseta_HCALOUT[m_Nclusters_HCALOUT]/D");
//     m_EventTree->Branch("m_clustheta_HCALOUT", m_clustheta_HCALOUT, "m_clustheta_HCALOUT[m_Nclusters_HCALOUT]/D");
//     m_EventTree->Branch("m_cluspt_HCALOUT", m_cluspt_HCALOUT, "m_cluspt_HCALOUT[m_Nclusters_HCALOUT]/D");
//     m_EventTree->Branch("m_clusphi_HCALOUT", m_clusphi_HCALOUT, "m_clusphi_HCALOUT[m_Nclusters_HCALOUT]/D");
//     m_EventTree->Branch("m_cluspx_HCALOUT", m_cluspx_HCALOUT, "m_cluspx_HCALOUT[m_Nclusters_HCALOUT]/D");
//     m_EventTree->Branch("m_cluspy_HCALOUT", m_cluspy_HCALOUT, "m_cluspy_HCALOUT[m_Nclusters_HCALOUT]/D");
//     m_EventTree->Branch("m_cluspz_HCALOUT", m_cluspz_HCALOUT, "m_cluspz_HCALOUT[m_Nclusters_HCALOUT]/D");
//     m_EventTree->Branch("m_E_4x4_HCALOUT", m_E_4x4_HCALOUT, "m_E_4x4_HCALOUT[m_Nclusters_HCALOUT]/D");
//     m_EventTree->Branch("m_clus_trueID_HCALOUT", m_clus_trueID_HCALOUT, "m_clus_trueID_HCALOUT[m_Nclusters_HCALOUT]/D");
//
//     m_EventTree->Branch("m_Ntowers_HCALOUT",&m_Ntowers_HCALOUT,"m_Ntowers_HCALOUT/I");
//     m_EventTree->Branch("m_tower_E_HCALOUT",m_tower_E_HCALOUT,"m_tower_E_HCALOUT[m_Ntowers_HCALOUT]/D");
//     m_EventTree->Branch("m_tower_iEta_HCALOUT",m_tower_iEta_HCALOUT,"m_tower_iEta_HCALOUT[m_Ntowers_HCALOUT]/I");
//     m_EventTree->Branch("m_tower_iPhi_HCALOUT",m_tower_iPhi_HCALOUT,"m_tower_iPhi_HCALOUT[m_Ntowers_HCALOUT]/I");
//     m_EventTree->Branch("m_tower_trueID_HCALOUT",m_tower_trueID_HCALOUT,"m_tower_trueID_HCALOUT[m_Ntowers_HCALOUT]/I");
//     }
//     if (m_doHEPMCpart){
//        m_EventTree->Branch("m_nHEPMCparts", &m_nHEPMCparts,"m_nHEPMCparts/I");
//        m_EventTree->Branch("m_process_id", &m_process_id, "m_process_id/I");
//        m_EventTree->Branch("m_partid1", m_partid1, "m_partid1[m_nHEPMCparts]/I");
//        m_EventTree->Branch("m_partid2", m_partid2, "m_partid2[m_nHEPMCparts]/I");
//        m_EventTree->Branch("m_x1", m_x1, "m_x1/D");
//        m_EventTree->Branch("m_x2", m_x2, "m_x2/D");
//        m_EventTree->Branch("m_mpi", m_mpi, "m_mpi[m_nHEPMCparts]/I");
//        m_EventTree->Branch("m_truthenergy", m_truthenergy, "m_truthenergy[m_nHEPMCparts]/D");
//        m_EventTree->Branch("m_trutheta", m_trutheta, "m_trutheta[m_nHEPMCparts]/D");
//        m_EventTree->Branch("m_truthphi", m_truthphi, "m_truthphi[m_nHEPMCparts]/D");
//        m_EventTree->Branch("m_truthpx", m_truthpx, "m_truthpx[m_nHEPMCparts]/D");
//        m_EventTree->Branch("m_truthpy", m_truthpy, "m_truthpy[m_nHEPMCparts]/D");
//        m_EventTree->Branch("m_truthpz", m_truthpz, "m_truthpz[m_nHEPMCparts]/D");
//        m_EventTree->Branch("m_truthpt", m_truthpt, "m_truthpt[m_nHEPMCparts]/D");
//        m_EventTree->Branch("m_truthpid", m_truthpid, "m_truthpid[m_nHEPMCparts]/I");
//
//     }
//     if (m_doPHG4part){
//       m_EventTree->Branch("m_nMCparts", &m_nMCparts,"m_nMCparts/I");
//       m_EventTree->Branch("m_truthenergy", m_truthenergy, "m_truthenergy[m_nMCparts]/D");
//       m_EventTree->Branch("m_truthp", m_truthp, "m_truthp[m_nMCparts]/D");
//       m_EventTree->Branch("m_truthpx", m_truthpx, "m_truthpx[m_nMCparts]/D");
//       m_EventTree->Branch("m_truthpy", m_truthpy, "m_truthpy[m_nMCparts]/D");
//       m_EventTree->Branch("m_truthpz", m_truthpz, "m_truthpz[m_nMCparts]/D");
//       m_EventTree->Branch("m_truthpt", m_truthpt, "m_truthpt[m_nMCparts]/D");
//       m_EventTree->Branch("m_truthphi", m_truthphi, "m_truthphi[m_nMCparts]/D");
//       m_EventTree->Branch("m_trutheta", m_trutheta, "m_trutheta[m_nMCparts]/D");
//       m_EventTree->Branch("m_truthpid", m_truthpid, "m_truthpid[m_nMCparts]/I");
//
//     }
//
//
//     }
