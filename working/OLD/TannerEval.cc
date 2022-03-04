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

TannerEval::TannerEval(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_outname(filename)
  , m_doEventInfo(false)
  , m_doGeo(false)
  , m_doTracks(false)
  , m_doHCALIN(false)
  , m_doHCALOUT(false)
  , m_doCEMC(false)
  , m_doHEPMCpart(false)
  , m_doPHG4part(false)
{
    WaterSeeds();
    GrowTrees();

}

  //Module Destructor
TannerEval::~TannerEval()
{
    delete m_EventTree;
    delete m_GeoTree;
    delete m_TrackTree;
    delete m_CEMCTree;
    delete m_HCALINTree;
    delete m_HCALOUTTree;
    delete m_HEPMCTree;
    delete m_PHG4Tree;
}



int TannerEval::Init(PHCompositeNode *topNode){
    if (Verbosity() > 0 ) { cout << "Initializing TannerEval" <<endl;}
    //Create output TFLIES
    m_outfile = new TFile(m_outname.c_str(), "RECREATE");
    if(m_doGeo){  m_outfile_geometry = new TFile("geometry.root","RECREATE"); }
    m_iEvent = 0;


    return Fun4AllReturnCodes::EVENT_OK;

}

  int TannerEval::process_event(PHCompositeNode *topNode){


      if (Verbosity() > 0 ) cout << "Processing event via TannerEval" <<endl;

      //Get Event information
      if (m_doEventInfo){ getEventInfo(topNode); }

      //Resonstructed tracks and truth values
      if(m_doTracks){ getTracks(topNode); }

      // Processing for each calo type
      if (m_doCEMC){
        //check if CEMC eval stack set_exist
        if(!m_CEMC_EvalStack){
          m_CEMC_EvalStack = new CaloEvalStack(topNode,"CEMC");
          m_CEMC_EvalStack->set_verbosity(Verbosity()+1);
          }
        else { m_CEMC_EvalStack->next_event(topNode); }

        cout << "loaded EM Cal evalstack" << endl;

        getCEMC(topNode);

      }

      if (m_doHCALIN){
        if(!m_HCALIN_EvalStack){
          m_HCALIN_EvalStack = new CaloEvalStack(topNode,"HCALIN");
          m_HCALIN_EvalStack->set_verbosity(Verbosity()+1);
        }
        else {m_HCALIN_EvalStack->next_event(topNode); }

        cout << "loaded HCALIN evalstack" << endl;
        getHCALIN(topNode);
        }

      if (m_doHCALOUT){

        if(!m_HCALOUT_EvalStack){
          m_HCALOUT_EvalStack = new CaloEvalStack(topNode,"HCALOUT");
          m_HCALOUT_EvalStack->set_verbosity(Verbosity()+1);
        }

        else { m_HCALOUT_EvalStack->next_event(topNode);}

        cout << "loaded OUT HCal evalstack" << endl;
        getHCALOUT(topNode);

      }

      //counting particles
      if (m_doHEPMCpart){ getHEPMC(topNode);}

      if (m_doPHG4part){getPHG4Truth(topNode);}

      //just in case
      m_iEvent++;
      cout << "Event " << m_iEvent << " complete" <<endl;

      return Fun4AllReturnCodes::EVENT_OK;

  }

  int TannerEval::End(PHCompositeNode *topNode){
    if (Verbosity() > 0) { cout << "Ending TTANEVAL package" << endl; }

    m_outfile->cd();
    if(m_doEventInfo) {  m_EventTree->Write();}
    if(m_doGeo){m_GeoTree->Write();}
    if(m_doTracks) { m_TrackTree->Write(); }
    if(m_doCEMC) { m_CEMCTree->Write(); }
    if(m_doHCALIN){ m_HCALINTree->Write();}
    if(m_doHCALOUT){ m_HCALOUTTree->Write();}
    if(m_doHEPMCpart) {m_HEPMCTree->Write();}
    if(m_doPHG4part) {m_PHG4Tree->Write();}
    m_outfile->Write();
    m_outfile->Close();
    delete m_outfile;

    if(m_doGeo){
      m_outfile_geometry->cd();
      m_GeoTree->Write();
      m_outfile_geometry->Write();
      m_outfile_geometry->Close();
      delete m_outfile_geometry;
    }

    if (Verbosity() > 0)  {
      cout << "========================= " << Name() << "::End() ============================" << endl;
      cout << " " << m_iEvent << " events of output written to: " << m_outname << endl;
      cout << "===========================================================================" << endl;
    }

   //if(m_HCALIN_EvalStack) delete m_HCALIN_EvalStack;
   //if(m_HCALOUT_EvalStack) delete m_HCALOUT_EvalStack;
   //if(m_CEMC_EvalStack) delete m_CEMC_EvalStack;
   //if(m_SVTX_EvalStack) delete m_SVTX_EvalStack;

  cout<< "=============================GOODBYE===================================="<<endl;
  return 0;

  }

void TannerEval::getEventInfo(PHCompositeNode *topNode){

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
      return;
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

    m_EventTree->Fill();
    m_EventTree->Write();
}

void TannerEval::getCEMC(PHCompositeNode *topNode){
        //Do towers for CEMC from method in EventEvaluator.h
    // if (Verbosity() > 0) { cout << "Grabbing CEMC event towers" << endl; }
    // m_Ntowers_CEMC = 0;
    //
    //
    // CaloRawTowerEval* towerevalCEMC = m_CEMC_EvalStack->get_rawtower_eval();
    // string towernodeCEMC = "TOWER_CALIB_CEMC";
    // RawTowerContainer *towersCEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeCEMC.c_str());
    //
    // if (towersCEMC){
    //   m_ALLcalo_Ntowers = 0;
    //   if (Verbosity() > 0) { cout << "saving EMC towers" << endl; }
    //   string towergeomnodeCEMC = "TOWERGEOM_CEMC";
    //   RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeCEMC.c_str());
    //
    //   if (towergeom){
    //   if(m_doGeo){
    //     cout<<"MADE IT TO GEO CEMC LOOP" <<endl;
    //     RawTowerGeomContainer::ConstRange all_towers = towergeom->get_tower_geometries();
    //     for (RawTowerGeomContainer::ConstIterator it = all_towers.first; it != all_towers.second; ++it){
    //
    //       m_ALLcalo_ID = kCEMC;
    //       m_ALLcalo_towers_iEta = it->second->get_bineta();
    //       m_ALLcalo_towers_iPhi = it->second->get_binphi();
    //       m_ALLcalo_towers_Eta = it->second->get_eta();
    //       m_ALLcalo_towers_Phi = it->second->get_phi();
    //       m_ALLcalo_towers_x = it->second->get_center_x();
    //       m_ALLcalo_towers_y = it->second->get_center_y();
    //       m_ALLcalo_towers_z = it->second->get_center_z();
    //       m_ALLcalo_Ntowers++;
    //       m_GeoTree->Fill();
    //     }
    //
    //     cout<<"FILLING "<<endl;
    //
    //     cout<<"CLEARING"<<endl;
    //
    //   }
    //   cout <<"DID CEMC GEO" <<endl;
    //   RawTowerContainer::ConstRange begin_end = towersCEMC->getTowers();
    //   RawTowerContainer::ConstIterator rtiter;
    //   for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter){
    //     RawTower* tower = rtiter->second;
    //     if (tower){
    //       // min energy cut
    //       //if (tower->get_energy() < _reco_e_threshold) continue;
    //
    //       m_tower_iEta_CEMC = tower->get_bineta();
    //       m_tower_iPhi_CEMC = tower->get_binphi();
    //       m_tower_E_CEMC = tower->get_energy();
    //
    //       PHG4Particle* primary = towerevalCEMC->max_truth_primary_particle_by_energy(tower);
    //       if (primary){ m_tower_trueID_CEMC = primary->get_track_id();}
    //         // gflavor = primary->get_pid();
    //         // efromtruth = towerevalCEMC->get_energy_contribution(tower, primary);
    //
    //       else{ m_tower_trueID_CEMC = -10; }
    //       m_CEMCTree->Fill();
    //       m_CEMCTree->Write();
    //       m_Ntowers_CEMC++;
    //     }
    //   }
    //   }
    //   else{
    //      cout << PHWHERE << " ERROR: Can't find " << towergeomnodeCEMC << endl;
    //   // return;
    //    }
    //   if (Verbosity() > 0) { cout << "saved\t" << m_Ntowers_CEMC << "\tCEMC towers" << endl; }
    //    }
    // else {  if (Verbosity() > 0) { cout << PHWHERE << " ERROR: Can't find " << towernodeCEMC << endl; }
    // // return;
    // }

    /// Get the raw cluster container
    /// Note: other cluster containers exist as well. Check out the node tree when
    /// you run a simulation
    RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
    CaloRawClusterEval *clusterevalCEMC = m_CEMC_EvalStack->get_rawcluster_eval();
    m_Nclusters_CEMC=0;
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
    if (vtx == nullptr)
      return;

    /// Trigger emulator
    CaloTriggerInfo *trigger = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");

    /// Can obtain some trigger information if desired
    if(trigger)
      {
        m_E_4x4_CEMC = trigger->get_best_EMCal_4x4_E();
      }

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
      m_clusenergy_CEMC = E_vec_cluster.mag();
      m_cluseta_CEMC = E_vec_cluster.pseudoRapidity();
      m_clustheta_CEMC = E_vec_cluster.getTheta();
      m_cluspt_CEMC = E_vec_cluster.perp();
      m_clusphi_CEMC = E_vec_cluster.getPhi();

      //If pt cut was wanted
      // if (m_cluspt_CEMC < m_mincluspt)
      //   continue;
      m_cluspx_CEMC =  (E_vec_cluster.perp())* cos(E_vec_cluster.getPhi());
      m_cluspy_CEMC = (E_vec_cluster.perp()) * sin(E_vec_cluster.getPhi());
      m_cluspz_CEMC = sqrt(((E_vec_cluster.mag()) *(E_vec_cluster.mag())) - (((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))) - (((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))));

      PHG4Particle* primary = clusterevalCEMC->max_truth_primary_particle_by_energy(cluster);
      if (primary){  m_clus_trueID_CEMC = primary->get_track_id(); }
      else { m_clus_trueID_CEMC = -10; }
      //fill the cluster tree with all CEMC clusters
      m_CEMCTree->Fill();
      m_CEMCTree->Write();
      m_Nclusters_CEMC++;

    }

}

void TannerEval::getHCALIN(PHCompositeNode *topNode){
        //Do towers for HCALIN from method in EventEvaluator.h
  //   if (Verbosity() > 0) { cout << "Grabbing HCALIN event towers" << endl; }
  //
  //   CaloRawTowerEval *towerevalHCALIN = m_HCALIN_EvalStack->get_rawtower_eval();
  //   m_Ntowers_HCALIN = 0;
  //   string towernodeHCALIN = "TOWER_CALIB_HCALIN";
  //   RawTowerContainer *towersHCALIN = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALIN.c_str());
  //
  //   if (towersHCALIN){
  //     m_ALLcalo_Ntowers=0;
  //     if (Verbosity() > 0) { cout << "saving HCALIN towers" << endl; }
  //     string towergeomnodeHCALIN = "TOWERGEOM_HCALIN";
  //     RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALIN.c_str());
  //
  //     if (towergeom){
  //     if(m_doGeo ){
  //       RawTowerGeomContainer::ConstRange all_towers = towergeom->get_tower_geometries();
  //       for (RawTowerGeomContainer::ConstIterator it = all_towers.first; it != all_towers.second; ++it){
  //
  //         m_ALLcalo_ID = kHCALIN;
  //         m_ALLcalo_towers_iEta = it->second->get_bineta();
  //         m_ALLcalo_towers_iPhi = it->second->get_binphi();
  //         m_ALLcalo_towers_Eta = it->second->get_eta();
  //         m_ALLcalo_towers_Phi = it->second->get_phi();
  //         m_ALLcalo_towers_x = it->second->get_center_x();
  //         m_ALLcalo_towers_y = it->second->get_center_y();
  //         m_ALLcalo_towers_z = it->second->get_center_z();
  //         m_ALLcalo_Ntowers++;
  //         m_GeoTree->Fill();
  //       }
  //
  //     }
  //     RawTowerContainer::ConstRange begin_end = towersHCALIN->getTowers();
  //     RawTowerContainer::ConstIterator rtiter;
  //     for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter){
  //       RawTower* tower = rtiter->second;
  //       if (tower){
  //         // min energy cut
  //         //if (tower->get_energy() < _reco_e_threshold) continue;
  //
  //         m_tower_iEta_HCALIN = tower->get_bineta();
  //         m_tower_iPhi_HCALIN = tower->get_binphi();
  //         m_tower_E_HCALIN = tower->get_energy();
  //
  //         PHG4Particle* primary = towerevalHCALIN->max_truth_primary_particle_by_energy(tower);
  //         if (primary){ m_tower_trueID_HCALIN = primary->get_track_id();}
  //           // gflavor = primary->get_pid();
  //           // efromtruth = towerevalHCALIN->get_energy_contribution(tower, primary);
  //
  //         else{ m_tower_trueID_HCALIN = -10; }
  //         m_HCALINTree->Fill();
  //         m_HCALINTree->Write();
  //         m_Ntowers_HCALIN++;
  //       }
  //     }
  //   }
  //   else{
  //     if (Verbosity() > 0) { cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALIN << endl; }
  //     // return;
  //   }
  //   if (Verbosity() > 0) { cout << "saved\t" << m_Ntowers_HCALIN << "\tHCALIN towers" << endl; }
  // }
  //   else {
  //   if (Verbosity() > 0) { cout << PHWHERE << " ERROR: Can't find " << towernodeHCALIN << endl; }
  //   // return;
  // }

    /// Get the raw cluster container
    /// Note: other cluster containers exist as well. Check out the node tree when
    /// you run a simulation
    RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");
    CaloRawClusterEval *clusterevalHCALIN = m_HCALIN_EvalStack->get_rawcluster_eval();
    m_Nclusters_HCALIN=0;
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
      cout << "TannerEval::getHCALINlusters - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
      assert(vertexmap);  // force quit

      return;
    }

    if (vertexmap->empty())
    {
      cout << "TannerEval::getHCALINlusters - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
      return;
    }

    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx == nullptr)
      return;

    /// Trigger emulator
    CaloTriggerInfo *trigger = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");

    /// Can obtain some trigger information if desired
    if(trigger)
      {
      //  m_E_4x4_HCALIN = trigger->get_best_HCALIN_4x4_E();
      }

    RawClusterContainer::ConstRange begin_end = clusters->getClusters();
    RawClusterContainer::ConstIterator clusIter;

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
      m_clusenergy_HCALIN = E_vec_cluster.mag();
      m_cluseta_HCALIN = E_vec_cluster.pseudoRapidity();
      m_clustheta_HCALIN = E_vec_cluster.getTheta();
      m_cluspt_HCALIN = E_vec_cluster.perp();
      m_clusphi_HCALIN = E_vec_cluster.getPhi();

      //If pt cut was wanted
      // if (m_cluspt_HCALIN < m_mincluspt)
      //   continue;
      m_cluspx_HCALIN =  (E_vec_cluster.perp())* cos(E_vec_cluster.getPhi());
      m_cluspy_HCALIN = (E_vec_cluster.perp()) * sin(E_vec_cluster.getPhi());
      m_cluspz_HCALIN = sqrt(((E_vec_cluster.mag()) *(E_vec_cluster.mag())) - (((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))) - (((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))));

      PHG4Particle* primary = clusterevalHCALIN->max_truth_primary_particle_by_energy(cluster);
      if (primary){  m_clus_trueID_HCALIN = primary->get_track_id(); }
      else { m_clus_trueID_HCALIN= -10; }
      //fill the cluster tree with all HCALIN clusters
      m_HCALINTree->Fill();
      m_HCALINTree->Write();
      m_Nclusters_HCALIN++;

    }

  }
void TannerEval::getHCALOUT(PHCompositeNode *topNode){
        //Do towers for HCALOUT from method in EventEvaluator.h
  //   if (Verbosity() > 0) { cout << "Grabbing HCALOUT event towers" << endl; }
  //
  //   CaloRawTowerEval *towerevalHCALOUT = m_HCALOUT_EvalStack->get_rawtower_eval();
  //   m_Ntowers_HCALOUT = 0;
  //   string towernodeHCALOUT = "TOWER_CALIB_HCALOUT";
  //   RawTowerContainer *towersHCALOUT = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALOUT.c_str());
  //
  //   if (towersHCALOUT){
  //     m_ALLcalo_Ntowers=0;
  //     if (Verbosity() > 0) { cout << "saving HCALOUT towers" << endl; }
  //     string towergeomnodeHCALOUT = "TOWERGEOM_HCALOUT";
  //     RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALOUT.c_str());
  //
  //     if (towergeom){
  //     if(m_doGeo){
  //       RawTowerGeomContainer::ConstRange all_towers = towergeom->get_tower_geometries();
  //       for (RawTowerGeomContainer::ConstIterator it = all_towers.first; it != all_towers.second; ++it){
  //
  //         m_ALLcalo_ID = kHCALOUT;
  //         m_ALLcalo_towers_iEta = it->second->get_bineta();
  //         m_ALLcalo_towers_iPhi = it->second->get_binphi();
  //         m_ALLcalo_towers_Eta = it->second->get_eta();
  //         m_ALLcalo_towers_Phi = it->second->get_phi();
  //         m_ALLcalo_towers_x = it->second->get_center_x();
  //         m_ALLcalo_towers_y = it->second->get_center_y();
  //         m_ALLcalo_towers_z = it->second->get_center_z();
  //         m_ALLcalo_Ntowers++;
  //         m_GeoTree->Fill();
  //       }
  //
  //     }
  //     RawTowerContainer::ConstRange begin_end = towersHCALOUT->getTowers();
  //     RawTowerContainer::ConstIterator rtiter;
  //     for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter){
  //       RawTower* tower = rtiter->second;
  //       if (tower){
  //         // min energy cut
  //         //if (tower->get_energy() < _reco_e_threshold) continue;
  //
  //         m_tower_iEta_HCALOUT = tower->get_bineta();
  //         m_tower_iPhi_HCALOUT = tower->get_binphi();
  //         m_tower_E_HCALOUT = tower->get_energy();
  //
  //         PHG4Particle* primary = towerevalHCALOUT->max_truth_primary_particle_by_energy(tower);
  //         if (primary){ m_tower_trueID_HCALOUT = primary->get_track_id();}
  //           // gflavor = primary->get_pid();
  //           // efromtruth = towerevalHCALOUT->get_energy_contribution(tower, primary);
  //
  //         else{ m_tower_trueID_HCALOUT= -10; }
  //           m_HCALOUTTree->Fill();
  //         m_Ntowers_HCALOUT++;
  //       }
  //     }
  //   }
  //   else{
  //     if (Verbosity() > 0) { cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALOUT << endl; }
  //     // return;
  //   }
  //   if (Verbosity() > 0) { cout << "saved\t" << m_Ntowers_HCALOUT << "\tHCALOUT towers" << endl; }
  // }
  //   else {
  //   if (Verbosity() > 0) { cout << PHWHERE << " ERROR: Can't find " << towernodeHCALOUT << endl; }
  //   // return;
  // }

    /// Get the raw cluster container
    /// Note: other cluster containers exist as well. Check out the node tree when
    /// you run a simulation
    RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");
    CaloRawClusterEval *clusterevalHCALOUT = m_HCALIN_EvalStack->get_rawcluster_eval();
    m_Nclusters_HCALOUT=0;
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
      cout << "TannerEval::getHCALOUTlusters - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
      assert(vertexmap);  // force quit

      return;
    }

    if (vertexmap->empty())
    {
      cout << "TannerEval::getHCALOUTlusters - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
      return;
    }

    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx == nullptr)
      return;

    /// Trigger emulator
    CaloTriggerInfo *trigger = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");

    /// Can obtain some trigger information if desired
    if(trigger)
      {
      //  m_E_4x4_HCALOUT = trigger->get_best_HCALOUT_4x4_E();
      }

    RawClusterContainer::ConstRange begin_end = clusters->getClusters();
    RawClusterContainer::ConstIterator clusIter;

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
      m_clusenergy_HCALOUT = E_vec_cluster.mag();
      m_cluseta_HCALOUT = E_vec_cluster.pseudoRapidity();
      m_clustheta_HCALOUT = E_vec_cluster.getTheta();
      m_cluspt_HCALOUT = E_vec_cluster.perp();
      m_clusphi_HCALOUT = E_vec_cluster.getPhi();

      //If pt cut was wanted
      // if (m_cluspt_HCALOUT < m_mincluspt)
      //   continue;
      m_cluspx_HCALOUT =  (E_vec_cluster.perp())* cos(E_vec_cluster.getPhi());
      m_cluspy_HCALOUT = (E_vec_cluster.perp()) * sin(E_vec_cluster.getPhi());
      m_cluspz_HCALOUT = sqrt(((E_vec_cluster.mag()) *(E_vec_cluster.mag())) - (((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* cos(E_vec_cluster.getPhi()))) - (((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))*((E_vec_cluster.perp())* sin(E_vec_cluster.getPhi()))));

      PHG4Particle* primary = clusterevalHCALOUT->max_truth_primary_particle_by_energy(cluster);
      if (primary){  m_clus_trueID_HCALOUT = primary->get_track_id(); }
      else { m_clus_trueID_HCALOUT = -10; }
      //fill the cluster tree with all HCALOUT clusters
      m_HCALOUTTree->Fill();
      m_HCALOUTTree->Fill();
      m_Nclusters_HCALOUT++;


    }


  }
  /**
       * This method gets all of the HEPMC truth particles from the node tree
       * and stores them in a ROOT TTree. The HEPMC truth particles are what,
       * for example, directly comes out of PYTHIA and thus gives you all of
       * the associated parton information
   */
void TannerEval::getHEPMC(PHCompositeNode *topNode){
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
    m_NHEPMCpartsinevent = 0;
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
        m_partid1 = pdfinfo->id1();
        m_partid2 = pdfinfo->id2();
        m_x1 = pdfinfo->x1();
        m_x2 = pdfinfo->x2();

        /// Are there multiple partonic intercations in a p+p event
        m_mpi = truthevent->mpi();

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
          m_truthenergy = (*iter)->momentum().e();
          m_truthpid = (*iter)->pdg_id();

          m_trutheta = (*iter)->momentum().pseudoRapidity();
          m_truthphi = (*iter)->momentum().phi();
          m_truthpx = (*iter)->momentum().px();
          m_truthpy = (*iter)->momentum().py();
          m_truthpz = (*iter)->momentum().pz();
          m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

          /// Fill the truth tree
          m_NHEPMCpartsinevent++;
          m_HEPMCTree->Fill();
          m_HEPMCTree->Write();
        }
      }
    }

}
  /**
   * This function collects the truth PHG4 stable particles that
   * are produced from PYTHIA, or some other event generator. These
   * are the stable particles, e.g. there are not any (for example)
   * partons here.
   */
void TannerEval::getPHG4Truth(PHCompositeNode *topNode){
    /// G4 truth particle node
    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    if (!truthinfo)
    {
      cout << PHWHERE
           << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
           << endl;
      return;
    }

    /// Get the primary particle range
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

    /// Loop over the G4 truth (stable) particles
    m_NPHG4partsinevent = 0;
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
         iter != range.second;
         ++iter)
    {
      /// Get this truth particle
      const PHG4Particle *truth = iter->second;

      /// Get this particles momentum, etc.
      m_truthpx_PHG4 = truth->get_px();
      m_truthpy_PHG4 = truth->get_py();
      m_truthpz_PHG4 = truth->get_pz();
      m_truthp_PHG4 = sqrt(m_truthpx_PHG4 * m_truthpx_PHG4 + m_truthpy_PHG4 * m_truthpy_PHG4 + m_truthpz_PHG4 * m_truthpz_PHG4);

      m_truthenergy_PHG4 = truth->get_e();
      m_truthpt_PHG4 = sqrt(m_truthpx_PHG4 * m_truthpx_PHG4 + m_truthpy_PHG4 * m_truthpy_PHG4);

      m_truthphi_PHG4 = atan(m_truthpy / m_truthpx);

      m_trutheta_PHG4 = atanh(m_truthpz_PHG4 / m_truthenergy_PHG4);
      /// Check for nans
      if (m_trutheta_PHG4 != m_trutheta_PHG4)
        m_trutheta_PHG4 = -99;
      m_truthpid_PHG4 = truth->get_pid();

      /// Fill the g4 truth tree
      m_NPHG4partsinevent++;
    m_PHG4Tree->Fill();
    m_PHG4Tree->Write();
    }

  }

void TannerEval::getTracks(PHCompositeNode *topNode){
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
    if(!m_SVTX_EvalStack)
      {
        m_SVTX_EvalStack = new SvtxEvalStack(topNode);
        m_SVTX_EvalStack->set_verbosity(Verbosity());
      }

    m_SVTX_EvalStack->next_event(topNode);

    /// Get the track evaluator
    SvtxTrackEval *trackeval = m_SVTX_EvalStack->get_track_eval();

    /// Get the range for primary tracks
    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    if (Verbosity() > 1)
    {
      cout << "Get the SVTX tracks" << endl;
    }
    m_Ntracksinevent = 0;
    for (SvtxTrackMap::Iter iter = trackmap->begin();
         iter != trackmap->end();
         ++iter)
    {
      SvtxTrack *track = iter->second;

      /// Get the reconstructed track info
      m_tr_px  = track->get_px();
      m_tr_py  = track->get_py();
      m_tr_pz  = track->get_pz();
      m_tr_p  = sqrt(m_tr_px  * m_tr_px  + m_tr_py  * m_tr_py  + m_tr_pz  * m_tr_pz );

      m_tr_pt  = sqrt(m_tr_px  * m_tr_px  + m_tr_py  * m_tr_py );

      // Make some cuts on the track to clean up sample
      if (m_tr_pt  < 0.5)
        continue;

      m_tr_phi  = track->get_phi();
      m_tr_eta  = track->get_eta();

      m_charge  = track->get_charge();
      m_chisq  = track->get_chisq();
      m_ndf  = track->get_ndf();
      m_dca  = track->get_dca();
      m_tr_x  = track->get_x();
      m_tr_y  = track->get_y();
      m_tr_z  = track->get_z();

      /// Get truth track info that matches this reconstructed track
      PHG4Particle *truthtrack = trackeval->max_truth_particle_by_nclusters(track);
      m_truth_is_primary = truthinfo->is_primary(truthtrack);

      m_truthtrackpx  = truthtrack->get_px();
      m_truthtrackpy  = truthtrack->get_py();
      m_truthtrackpz  = truthtrack->get_pz();
      m_truthtrackp  = sqrt(m_truthtrackpx * m_truthtrackpx + m_truthtrackpy * m_truthtrackpy + m_truthtrackpz * m_truthtrackpz);

      m_truthtracke  = truthtrack->get_e();

      m_truthtrackpt  = sqrt(m_truthtrackpx  * m_truthtrackpx  + m_truthtrackpy  * m_truthtrackpy );
      m_truthtrackphi  = atan(m_truthtrackpy  / m_truthtrackpx );
      m_truthtracketa  = atanh(m_truthtrackpz  / m_truthtrackp );
      m_truthtrackpid  = truthtrack->get_pid();
      m_TrackTree->Fill();
      m_TrackTree->Write();
      m_Ntracksinevent++;

    }

  }




void TannerEval::WaterSeeds(){

    m_outfile = new TFile();
  //  m_outfile_geometry = new TFile();

    //HEPMC PARTICLE VARIBLES
    m_NHEPMCpartsinevent = -99;//counter

    m_partid1 = -99;
    m_partid2 = -99;
    m_x1 = -99;
    m_x2 = -99;
    m_mpi = -99;
    m_process_id = -99;
    m_truthenergy = -99;
    m_trutheta = -99;
    m_truthphi = -99;
    m_truthpx = -99;
    m_truthpy = -99;
    m_truthpz = -99;
    m_truthpt = -99;
    m_truthp = -99;
    m_truthpid = -99;



    //PH4G PARTICLE VARIBLES
    m_NPHG4partsinevent = -99; //Counter

    m_truthpx_PHG4 = -99;
    m_truthpy_PHG4 = -99;
    m_truthpz_PHG4 = -99;
    m_truthpt_PHG4 = -99;
    m_truthp_PHG4 = -99;
    m_truthenergy_PHG4 = -99;
    m_trutheta_PHG4 = -99;
    m_truthphi_PHG4 = -99;
    m_truthpid_PHG4 = -99;




    // Track TTree variables
    m_Ntracksinevent = -99; //counter

    m_tr_px = -99;
    m_tr_py = -99;
    m_tr_pz = -99;
    m_tr_p = -99;
    m_tr_pt = -99;
    m_tr_phi = -99;
    m_tr_eta = -99;
    m_charge = -99;
    m_chisq = -99;
    m_ndf = -99;
    m_dca = -99;
    m_tr_x = -99;
    m_tr_y = -99;
    m_tr_z = -99;
    m_truth_is_primary = -99;
    m_truthtrackpx = -99;
    m_truthtrackpy = -99;
    m_truthtrackpz = -99;
    m_truthtrackp = -99;
    m_truthtracke = -99;
    m_truthtrackpt = -99;
    m_truthtrackphi = -99;
    m_truthtracketa = -99;
    m_truthtrackpid = -99;

    /// Cluster and tower variables


    //CEMC
    m_Nclusters_CEMC = -99; //cluster counter

    m_clusenergy_CEMC = -99;
    m_cluseta_CEMC = -99;
    m_clustheta_CEMC = -99;
    m_cluspt_CEMC = -99;
    m_clusphi_CEMC = -99;
    m_cluspx_CEMC = -99;
    m_cluspy_CEMC = -99;
    m_cluspz_CEMC = -99;
    m_E_4x4_CEMC = -99;
    m_clus_trueID_CEMC = -99;


    m_Ntowers_CEMC = -99; // tower counter

    m_tower_E_CEMC = -99;
    m_tower_iEta_CEMC = -99;
    m_tower_iPhi_CEMC = -99;
    m_tower_trueID_CEMC = -99;

    //HCALIN

    m_Nclusters_HCALIN = -99; //cluster counter

    m_clusenergy_HCALIN = -99;
    m_cluseta_HCALIN = -99;
    m_clustheta_HCALIN = -99;
    m_cluspt_HCALIN = -99;
    m_clusphi_HCALIN = -99;
    m_cluspx_HCALIN = -99;
    m_cluspy_HCALIN = -99;
    m_cluspz_HCALIN = -99;
    m_E_4x4_HCALIN = -99;
    m_clus_trueID_HCALIN= -99;


    m_Ntowers_HCALIN = -99; // tower counter

    m_tower_E_HCALIN = -99;
    m_tower_iEta_HCALIN = -99;
    m_tower_iPhi_HCALIN = -99;
    m_tower_trueID_HCALIN = -99;


    //HCALOUT
    m_Nclusters_HCALOUT = -99; //cluster counter

    m_clusenergy_HCALOUT = -99;
    m_cluseta_HCALOUT = -99;
    m_clustheta_HCALOUT = -99;
    m_cluspt_HCALOUT = -99;
    m_clusphi_HCALOUT = -99;
    m_cluspx_HCALOUT = -99;
    m_cluspy_HCALOUT = -99;
    m_cluspz_HCALOUT = -99;
    m_E_4x4_HCALOUT = -99;
    m_clus_trueID_HCALOUT= -99;


    m_Ntowers_HCALOUT = -99; // tower counter

    m_tower_E_HCALOUT = -99;
    m_tower_iEta_HCALOUT = -99;
    m_tower_iPhi_HCALOUT = -99;
    m_tower_trueID_HCALOUT = -99;

    m_ALLcalo_Ntowers=-99;
    m_ALLcalo_ID=-99;
    m_ALLcalo_towers_iEta = 99;
    m_ALLcalo_towers_iPhi = 99;
    m_ALLcalo_towers_Eta = -99;
    m_ALLcalo_towers_Phi = -99;
    m_ALLcalo_towers_x = -99;
    m_ALLcalo_towers_y = -99;
    m_ALLcalo_towers_z = -99;
}

void TannerEval::GrowTrees(){



 if (m_doEventInfo){
      // Event level info. This isn't the most efficient way to store this info, but it's straightforward
      // within the structure of the class, so the size is small compared to the rest of the output.
      m_EventTree = new TTree("EventTree", "EventTree");
      m_EventTree->Branch("m_cross_section", &m_cross_section, "m_cross_section/D");
      m_EventTree->Branch("m_event_weight", &m_event_weight, "m_event_weight/D");
      m_EventTree->Branch("m_n_generator_accepted", &m_n_generator_accepted, "m_n_generator_accepted/D");
    }
  if (m_doGeo){

      m_GeoTree = new TTree("GeoTree", "GeoTree");
      m_GeoTree->Branch("m_ALLcalo_ID", &m_ALLcalo_ID, "nCaloHits/I");
      m_GeoTree->Branch("m_ALLcalo_Ntowers", &m_ALLcalo_Ntowers, "m_ALLcalo_Ntowers/I");
      m_GeoTree->Branch("m_ALLcalo_towers_iEta", &m_ALLcalo_towers_iEta,"calo_towers_iEta/I");
      m_GeoTree->Branch("m_ALLcalo_towers_iPhi", &m_ALLcalo_towers_iPhi,"calo_towers_iPhi/I");
      m_GeoTree->Branch("m_ALLcalo_towers_Eta", &m_ALLcalo_towers_Eta,"calo_towers_Eta/D");
      m_GeoTree->Branch("m_ALLcalo_towers_Phi", &m_ALLcalo_towers_Phi, "calo_towers_Phi/D");
      m_GeoTree->Branch("m_ALLcalo_towers_x", &m_ALLcalo_towers_x,"calo_towers_x/D");
      m_GeoTree->Branch("m_ALLcalo_towers_y", &m_ALLcalo_towers_y,"calo_towers_y/D");
      m_GeoTree->Branch("m_ALLcalo_towers_z", &m_ALLcalo_towers_z,"calo_towers_z/D");


    }
    if (m_doTracks){
      m_TrackTree = new TTree("TrackTree", "TrackTree");
      m_TrackTree->Branch("m_Ntracksinevent", &m_Ntracksinevent ,"m_Ntracksinevent/I");
      m_TrackTree->Branch("m_tr_px", &m_tr_px, "m_tr_px/D");
      m_TrackTree->Branch("m_tr_py", &m_tr_py, "m_tr_py/D");
      m_TrackTree->Branch("m_tr_pz", &m_tr_pz, "m_tr_pz/D");
      m_TrackTree->Branch("m_tr_p", &m_tr_p, "m_tr_p/D");
      m_TrackTree->Branch("m_tr_pt", &m_tr_pt, "m_tr_pt/D");
      m_TrackTree->Branch("m_tr_phi", &m_tr_phi, "m_tr_phi/D");
      m_TrackTree->Branch("m_tr_eta", &m_tr_eta, "m_tr_eta/D");
      m_TrackTree->Branch("m_charge", &m_charge, "m_charge/I");
      m_TrackTree->Branch("m_chisq", &m_chisq, "m_chisq/D");
      m_TrackTree->Branch("m_ndf", &m_ndf, "m_ndf/I");
      m_TrackTree->Branch("m_dca", &m_dca, "m_dca/D");
      m_TrackTree->Branch("m_tr_x", &m_tr_x, "m_tr_x/D");
      m_TrackTree->Branch("m_tr_y", &m_tr_y, "m_tr_y/D");
      m_TrackTree->Branch("m_tr_z", &m_tr_z, "m_tr_z/D");
      m_TrackTree->Branch("m_truth_is_primary", &m_truth_is_primary, "m_truth_is_primary/I");
      m_TrackTree->Branch("m_truthtrackpx", &m_truthtrackpx, "m_truthtrackpx/D");
      m_TrackTree->Branch("m_truthtrackpy", &m_truthtrackpy, "m_truthtrackpy/D");
      m_TrackTree->Branch("m_truthtrackpz", &m_truthtrackpz, "m_truthtrackpz/D");
      m_TrackTree->Branch("m_truthtrackp", &m_truthtrackp, "m_truthtrackp/D");
      m_TrackTree->Branch("m_truthtracke", &m_truthtracke, "m_truthtracke/D");
      m_TrackTree->Branch("m_truthtrackpt", &m_truthtrackpt, "m_truthtrackpt/D");
      m_TrackTree->Branch("m_truthtrackphi", &m_truthtrackphi, "m_truthtrackphi/D");
      m_TrackTree->Branch("m_truthtracketa", &m_truthtracketa, "m_truthtracketa/D");
      m_TrackTree->Branch("m_truthtrackpid", &m_truthtrackpid, "m_truthtrackpid/I");

    }
    if (m_doHCALIN){
      m_HCALINTree = new TTree("HCALINTree", "HCALINTree");
      m_HCALINTree->Branch("m_Nclusters_HCALIN", &m_Nclusters_HCALIN,"m_Nclusters_HCALIN/I");
      m_HCALINTree->Branch("m_clusenergy_HCALIN", &m_clusenergy_HCALIN, "m_clusenergy_HCALIN/D");
      m_HCALINTree->Branch("m_cluseta_HCALIN", &m_cluseta_HCALIN, "m_cluseta_HCALIN/D");
      m_HCALINTree->Branch("m_clustheta_HCALIN", &m_clustheta_HCALIN, "m_clustheta_HCALIN/D");
      m_HCALINTree->Branch("m_cluspt_HCALIN", &m_cluspt_HCALIN, "m_cluspt_HCALIN/D");
      m_HCALINTree->Branch("m_clusphi_HCALIN", &m_clusphi_HCALIN, "m_clusphi_HCALIN/D");
      m_HCALINTree->Branch("m_cluspx_HCALIN", &m_cluspx_HCALIN, "m_cluspx_HCALIN/D");
      m_HCALINTree->Branch("m_cluspy_HCALIN", &m_cluspy_HCALIN, "m_cluspy_HCALIN/D");
      m_HCALINTree->Branch("m_cluspz_HCALIN", &m_cluspz_HCALIN, "m_cluspz_HCALIN/D");
      m_HCALINTree->Branch("m_E_4x4_HCALIN", &m_E_4x4_HCALIN, "m_E_4x4_HCALIN/D");
      m_HCALINTree->Branch("m_clus_trueID_HCALIN", &m_clus_trueID_HCALIN, "m_clus_trueID_HCALIN/D");

      m_HCALINTree->Branch("m_Ntowers_HCALIN",&m_Ntowers_HCALIN,"m_Ntowers_HCALIN/I");
      m_HCALINTree->Branch("m_tower_E_HCALIN",&m_tower_E_HCALIN,"m_tower_E_HCALIN/D");
      m_HCALINTree->Branch("m_tower_iEta_HCALIN",&m_tower_iEta_HCALIN,"m_tower_iEta_HCALIN/I");
      m_HCALINTree->Branch("m_tower_iPhi_HCALIN",&m_tower_iPhi_HCALIN,"m_tower_iPhi_HCALIN/I");
      m_HCALINTree->Branch("m_tower_trueID_HCALIN",&m_tower_trueID_HCALIN,"m_tower_trueID_HCALIN/I");


    }
    if (m_doCEMC){
      m_CEMCTree = new TTree("CEMCTree", "CEMCTree");
      m_CEMCTree->Branch("m_Nclusters_CEMC", &m_Nclusters_CEMC,"m_Nclusters_CEMC/I");
      m_CEMCTree->Branch("m_clusenergy_CEMC", &m_clusenergy_CEMC, "m_clusenergy_CEMC/D");
      m_CEMCTree->Branch("m_cluseta_CEMC", &m_cluseta_CEMC, "m_cluseta_CEMC/D");
      m_CEMCTree->Branch("m_clustheta_CEMC", &m_clustheta_CEMC, "m_clustheta_CEMC/D");
      m_CEMCTree->Branch("m_cluspt_CEMC", &m_cluspt_CEMC, "m_cluspt_CEMC/D");
      m_CEMCTree->Branch("m_clusphi_CEMC", &m_clusphi_CEMC, "m_clusphi_CEMC/D");
      m_CEMCTree->Branch("m_cluspx_CEMC", &m_cluspx_CEMC, "m_cluspx_CEMC/D");
      m_CEMCTree->Branch("m_cluspy_CEMC", &m_cluspy_CEMC, "m_cluspy_CEMC/D");
      m_CEMCTree->Branch("m_cluspz_CEMC", &m_cluspz_CEMC, "m_cluspz_CEMC/D");
      m_CEMCTree->Branch("m_E_4x4_CEMC", &m_E_4x4_CEMC, "m_E_4x4_CEMC/D");
      m_CEMCTree->Branch("m_clus_trueID_CEMC", &m_clus_trueID_CEMC, "m_clus_trueID_CEMC/D");

      m_CEMCTree->Branch("m_Ntowers_CEMC",&m_Ntowers_CEMC,"m_Ntowers_CEMC/I");
      m_CEMCTree->Branch("m_tower_E_CEMC",&m_tower_E_CEMC,"m_tower_E_CEMC/D");
      m_CEMCTree->Branch("m_tower_iEta_CEMC",&m_tower_iEta_CEMC,"m_tower_iEta_CEMC/I");
      m_CEMCTree->Branch("m_tower_iPhi_CEMC",&m_tower_iPhi_CEMC,"m_tower_iPhi_CEMC/I");
      m_CEMCTree->Branch("m_tower_trueID_CEMC",&m_tower_trueID_CEMC,"m_tower_trueID_CEMC/I");

    }
    if (m_doHCALOUT){
      m_HCALOUTTree = new TTree("HCALOUTTree", "HCALOUTTree");
      m_HCALOUTTree->Branch("m_Nclusters_HCALOUT", &m_Nclusters_HCALOUT,"m_Nclusters_HCALOUT/I");
      m_HCALOUTTree->Branch("m_clusenergy_HCALOUT", &m_clusenergy_HCALOUT, "m_clusenergy_HCALOUT/D");
      m_HCALOUTTree->Branch("m_cluseta_HCALOUT", &m_cluseta_HCALOUT, "m_cluseta_HCALOUT/D");
      m_HCALOUTTree->Branch("m_clustheta_HCALOUT", &m_clustheta_HCALOUT, "m_clustheta_HCALOUT/D");
      m_HCALOUTTree->Branch("m_cluspt_HCALOUT", &m_cluspt_HCALOUT, "m_cluspt_HCALOUT/D");
      m_HCALOUTTree->Branch("m_clusphi_HCALOUT", &m_clusphi_HCALOUT, "m_clusphi_HCALOUT/D");
      m_HCALOUTTree->Branch("m_cluspx_HCALOUT", &m_cluspx_HCALOUT, "m_cluspx_HCALOUT/D");
      m_HCALOUTTree->Branch("m_cluspy_HCALOUT", &m_cluspy_HCALOUT, "m_cluspy_HCALOUT/D");
      m_HCALOUTTree->Branch("m_cluspz_HCALOUT", &m_cluspz_HCALOUT, "m_cluspz_HCALOUT/D");
      m_HCALOUTTree->Branch("m_E_4x4_HCALOUT", &m_E_4x4_HCALOUT, "m_E_4x4_HCALOUT/D");
      m_HCALOUTTree->Branch("m_clus_trueID_HCALOUT", &m_clus_trueID_HCALOUT, "m_clus_trueID_HCALOUT/D");

      m_HCALOUTTree->Branch("m_Ntowers_HCALOUT",&m_Ntowers_HCALOUT,"m_Ntowers_HCALOUT/I");
      m_HCALOUTTree->Branch("m_tower_E_HCALOUT",&m_tower_E_HCALOUT,"m_tower_E_HCALOUT/D");
      m_HCALOUTTree->Branch("m_tower_iEta_HCALOUT",&m_tower_iEta_HCALOUT,"m_tower_iEta_HCALOUT/I");
      m_HCALOUTTree->Branch("m_tower_iPhi_HCALOUT",&m_tower_iPhi_HCALOUT,"m_tower_iPhi_HCALOUT/I");
      m_HCALOUTTree->Branch("m_tower_trueID_HCALOUT",&m_tower_trueID_HCALOUT,"m_tower_trueID_HCALOUT/I");

    }
    if (m_doHEPMCpart){
      m_HEPMCTree = new TTree("HEPMCTree", "HEPMCTree");
      m_HEPMCTree->Branch("m_NHEPMCpartsinevent", &m_NHEPMCpartsinevent,"m_NHEPMCpartsinevent/I");
      m_HEPMCTree->Branch("m_process_id", &m_process_id, "m_process_id/I");
      m_HEPMCTree->Branch("m_partid1", &m_partid1, "m_partid1/I");
      m_HEPMCTree->Branch("m_partid2", &m_partid2, "m_partid2/I");
      m_HEPMCTree->Branch("m_x1", &m_x1, "m_x1/D");
      m_HEPMCTree->Branch("m_x2", &m_x2, "m_x2/D");
      m_HEPMCTree->Branch("m_mpi", &m_mpi, "m_mpi/I");
      m_HEPMCTree->Branch("m_truthenergy", &m_truthenergy, "m_truthenergy/D");
      m_HEPMCTree->Branch("m_trutheta", &m_trutheta, "m_trutheta/D");
      m_HEPMCTree->Branch("m_truthphi", &m_truthphi, "m_truthphi/D");
      m_HEPMCTree->Branch("m_truthpx", &m_truthpx, "m_truthpx/D");
      m_HEPMCTree->Branch("m_truthpy", &m_truthpy, "m_truthpy/D");
      m_HEPMCTree->Branch("m_truthpz", &m_truthpz, "m_truthpz/D");
      m_HEPMCTree->Branch("m_truthpt", &m_truthpt, "m_truthpt/D");
      m_HEPMCTree->Branch("m_truthpid", &m_truthpid, "m_truthpid/I");

    }
    if (m_doPHG4part){
      m_PHG4Tree = new TTree("PHG4Tree","PHG4Tree");
      m_PHG4Tree->Branch("m_NPHG4partsinevent", &m_NPHG4partsinevent,"m_NPHG4partsinevent/I");
      m_PHG4Tree->Branch("m_truthenergy", &m_truthenergy, "m_truthenergy/D");
      m_PHG4Tree->Branch("m_truthp", &m_truthp, "m_truthp/D");
      m_PHG4Tree->Branch("m_truthpx", &m_truthpx, "m_truthpx/D");
      m_PHG4Tree->Branch("m_truthpy", &m_truthpy, "m_truthpy/D");
      m_PHG4Tree->Branch("m_truthpz", &m_truthpz, "m_truthpz/D");
      m_PHG4Tree->Branch("m_truthpt", &m_truthpt, "m_truthpt/D");
      m_PHG4Tree->Branch("m_truthphi", &m_truthphi, "m_truthphi/D");
      m_PHG4Tree->Branch("m_trutheta", &m_trutheta, "m_trutheta/D");
      m_PHG4Tree->Branch("m_truthpid", &m_truthpid, "m_truthpid/I");

    }
  }
