//# include "MEM/MEMAlgo/interface/PyConfig.h"

# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include <sstream>
# include <algorithm> 

# include "MEM/MEMAlgo/interface/Constants.h"
# include "MEM/MEMAlgo/interface/MGIntegration.h"
# include "MEM/MEMAlgo/interface/LHAPDF.h"
# include "MEM/MEMAlgo/interface/Processes.h"



# define evalGauss( x, a0, a1, a2 ) ( (a0) * exp( -0.5 * ((x)-(a1)) * ((x)-(a1)) / ((a2)*(a2)) ) ) 

using namespace std;

extern "C" gsl_rng *gslRNGInit( ) {
  gsl_rng_env_setup ();
  const gsl_rng_type *T = gsl_rng_rand48;
  gsl_rng *r = gsl_rng_alloc (T);
  return r;
}        

extern "C" void gslIntegrate( 
        gsl_monte_function *fct_descr, 
        IntegrationMsg_t* integration,
	int nbrOfDim,
	int nbrOfPoints,
        gsl_rng *rng, 
        gsl_monte_vegas_state *state,
        bool newGrid,
        double *res, double *err, double *chisq) {
  
  if ( newGrid )
    state->stage = 0;
  // else {
  // Don't change the status  assumed a stage >=1
  // }

  // Test boundaries
  double *lower = integration->lowerValues_;
  double *upper = integration->upperValues_;
  int i;
  bool boundariesOk = true;
  for (i=0; i< nbrOfDim; i++) {
    boundariesOk = boundariesOk && (lower[i] < upper[i]);
    if ( !boundariesOk) {
        std::cout << "lower[" << i << "] : " << lower[i] << " - upper[" << i << "] : " << upper[i] << std::endl;
    }
  }
  
  if ( boundariesOk ) {
    std::cout << "gslIntegrate, boundaries OK" << std::endl;
    gsl_monte_vegas_integrate ( fct_descr, integration->lowerValues_, integration->upperValues_, 
            nbrOfDim, nbrOfPoints, rng, state,
            res, err); 
    *chisq = state->chisq;
    //std::cout << "apres gsl_monte_vegas_integrate, res = " << *res << std::endl; // TEMPORAIRE
  } else {
    std::cout << "gslIntegrate, boundaries KO" << std::endl;
    *res   = 0.0;
    *err   = 0.0;
    *chisq = 1.0;
  }
}

void MGIntegration::setDefault() {
  // Set default
   
  // Computation 
  jacobianKind_ = UnknownJacobian;
  verbose_ = 0;
  
  // LHC 
//  double SqrtS_ = 13000; // LHC beam energy   

  // PDF  
  Q_ = Physics::mtop + 0.5*Physics::mHiggs ; // PDF calibration (Gev)

  
}

MGIntegration::MGIntegration() {
  
  // Set default   
  setDefault();
 }

//MGIntegration::MGIntegration( const char *configname, RunConfig &runConfig, int processID ) {
//MGIntegration::MGIntegration( const char *configname, int processID ) {
MGIntegration::MGIntegration( RunConfig runConfig, int processID ) {

  setDefault();
  // Run ttZ integration
  runttZ_integration_ = runConfig.runttZ_integration_;

  // Run ttZ_Zll integration
  runttZ_Zll_integration_ = runConfig.runttZ_Zll_integration_;

  // Run ttW integration
  runttW_integration_ = runConfig.runttW_integration_;

  // Run ttjets integration
  runttjets_integration_ = runConfig.runttjets_integration_;
  
  // Run ttbar_SL integration
  runttbar_SL_integration_ = runConfig.runttbar_SL_integration_;
  
  // Run ttbar_DL integration
  runttbar_DL_integration_ = runConfig.runttbar_DL_integration_;
  
  // Run ttbar_DL_fakelep integration
  runttbar_DL_fakelep_integration_ = runConfig.runttbar_DL_fakelep_integration_;
  
  // Run missing jet integration
  run_missing_jet_integration_ = runConfig.run_missing_jet_integration_;
  
  // Force missing jet integration
  force_missing_jet_integration_ = runConfig.force_missing_jet_integration_;
  
  // Force missing jet integration if no MEM sel perm passing cuts without missing jet
  force_missing_jet_integration_ifnoperm_ = runConfig.force_missing_jet_integration_ifnoperm_;
  
  // Random generator
  flagSameRNG_ = runConfig.flagSameRNG_;
  
  // Computation / run
  flagTFLepTau_ = true; flagTFHadTau_ = true; flagTFFake_ = true;
  flagTFMET_ = true; flagTFJet1_ = true; flagTFJet2_ = true;
  flagTFBJet_leptop_ = true;   flagTFBJet_hadtop_ = true;  flagTFTop_ = true;
  flagTF_fakelep_ = true; flagTF_fakeleptau_ = true;
  flagWME_ = true; flagJac_ = true;
  //
  flagTFLepTau_ = runConfig.flagTFLepTau_; 
  flagTFHadTau_ = runConfig.flagTFHadTau_; 
  flagTFFake_ = runConfig.flagTFFake_;
  flagTFMET_ = runConfig.flagTFMET_; 
  flagTFJet1_ = runConfig.flagTFJet1_; 
  flagTFJet2_ = runConfig.flagTFJet2_;
  flagTFBJet_leptop_ = runConfig.flagTFBJet_leptop_;   
  flagTFBJet_hadtop_ = runConfig.flagTFBJet_hadtop_;  
  flagTFTop_ = runConfig.flagTFTop_;
  flagTF_fakelep_ = runConfig.flagTF_fakelep_; 
  flagTF_fakeleptau_ = runConfig.flagTF_fakeleptau_;
  flagWME_ = runConfig.flagJac_; 
  flagJac_ = runConfig.flagWME_;

  // MET TF
  include_hadrecoil_ = runConfig.include_hadrecoil_;

  // ME version
  MEversion_ = runConfig.MEversion_;

  // Verbose mode  
  verbose_ = runConfig.verbose_;

  // Integral value file

  // MC integration
  nbrOfDim_ttH_ = 0;
  nbrOfDim_ttZ_ = 0;
  nbrOfDim_ttW_ = 0;
  nbrOfDim_ttjets_ = 0;
  nbrOfDim_ttbar_SL_ = 0;
  nbrOfDim_ttbar_DL_ = 0;
  nbrOfDim_ttZ_Zll_ = 0;

  nbrOfDim_ttH_miss_ = 0;
  nbrOfDim_ttZ_miss_ = 0;
  nbrOfDim_ttjets_miss_ = 0;
  nbrOfDim_ttbar_SL_miss_ = 0;
  nbrOfDim_ttZ_Zll_miss_ = 0;

  nbrOfDim_ttH_ = runConfig.nbrOfDimttH_ ;
  nbrOfDim_ttZ_ = runConfig.nbrOfDimttZ_ ;
  nbrOfDim_ttW_ = runConfig.nbrOfDimttW_ ;
  nbrOfDim_ttjets_ = runConfig.nbrOfDimttjets_ ;
  nbrOfDim_ttbar_SL_ = runConfig.nbrOfDimttbar_SL_ ;
  nbrOfDim_ttbar_DL_ = runConfig.nbrOfDimttbar_DL_ ;
  nbrOfDim_ttZ_Zll_ = runConfig.nbrOfDimttZ_Zll_ ;

  if(run_missing_jet_integration_){     
      nbrOfDim_ttH_miss_ = runConfig.nbrOfDimttH_miss_ ;
      nbrOfDim_ttZ_miss_ = runConfig.nbrOfDimttZ_miss_ ;
      nbrOfDim_ttjets_miss_ = runConfig.nbrOfDimttjets_miss_ ;
      nbrOfDim_ttbar_SL_miss_ = runConfig.nbrOfDimttbar_SL_miss_ ;
      nbrOfDim_ttZ_Zll_miss_ = runConfig.nbrOfDimttZ_Zll_miss_ ;
  }

  // GG TODO ???
  // method :
  // nbrOfDim_ = configAsLong("nbrOfDim" )
  
nbrOfPoints_ttH_ = runConfig.nbrOfPoints_ttH_ ;
nbrOfPoints_ttZ_ = runConfig.nbrOfPoints_ttZ_ ;
nbrOfPoints_ttW_ = runConfig.nbrOfPoints_ttW_ ;
nbrOfPoints_ttjets_ = runConfig.nbrOfPoints_ttjets_ ;
nbrOfPoints_ttbar_SL_ = runConfig.nbrOfPoints_ttbar_SL_ ;
nbrOfPoints_ttbar_DL_ = runConfig.nbrOfPoints_ttbar_DL_;
nbrOfPoints_ttZ_Zll_ = runConfig.nbrOfPoints_ttZ_Zll_ ;

  if(run_missing_jet_integration_){
    nbrOfPoints_ttH_miss_ = runConfig.nbrOfPoints_ttH_miss_ ;
    nbrOfPoints_ttZ_miss_ = runConfig.nbrOfPoints_ttZ_miss_ ;
    nbrOfPoints_ttjets_miss_ = runConfig.nbrOfPoints_ttjets_miss_ ;
    nbrOfPoints_ttbar_SL_miss_ = runConfig.nbrOfPoints_ttbar_SL_miss_ ;
    nbrOfPoints_ttZ_Zll_miss_ = runConfig.nbrOfPoints_ttZ_Zll_miss_ ;
    }

    //
  nbrOfPermut_per_jet_ = runConfig.nbrOfPermut_per_jet_;

  use_pT_TFJet_ = runConfig.use_pT_TFJet_;
  
  use_top_compatibility_check_ = runConfig.use_top_compatibility_check_;

  //Integration boundaries
  CI_TFJet_ = runConfig.CI_TFJet_;

  eta_acceptance_ = runConfig.eta_acceptance_ ;  
  jet_radius_ = runConfig.jet_radius_ ;  
  dR_veto_jet_lep_ = runConfig.dR_veto_jet_lep_ ;  
  rel_iso_lep_ = runConfig.rel_iso_lep_ ;  
  pT_cut_ = runConfig.pT_cut_ ;  

 if(run_missing_jet_integration_){
    phi_missing_jet_Boundaries_[0] = runConfig.phi_missing_jet_Boundaries_[0];  
    phi_missing_jet_Boundaries_[1] = runConfig.phi_missing_jet_Boundaries_[1];  
    cosTheta_missing_jet_Boundaries_[0] = runConfig.cosTheta_missing_jet_Boundaries_[0];  
    cosTheta_missing_jet_Boundaries_[1] = runConfig.cosTheta_missing_jet_Boundaries_[1];
 }

  phiNu_tlep_Boundaries_[0] = runConfig.phiNu_tlep_Boundaries_[0];
  phiNu_tlep_Boundaries_[1] = runConfig.phiNu_tlep_Boundaries_[1];

  cosThetaNu_tlep_Boundaries_[0] = runConfig.cosThetaNu_tlep_Boundaries_[0];
  cosThetaNu_tlep_Boundaries_[1] = runConfig.cosThetaNu_tlep_Boundaries_[1];  

  if(runttW_integration_){
  cosThetaNu_W_ttW_Boundaries_[0] = runConfig.cosThetaNu_W_ttW_Boundaries_[0] ;  
  cosThetaNu_W_ttW_Boundaries_[1] = runConfig.cosThetaNu_W_ttW_Boundaries_[1] ;  
  phiNu_W_ttW_Boundaries_[0] = runConfig.phiNu_W_ttW_Boundaries_[0] ;  
  phiNu_W_ttW_Boundaries_[1] = runConfig.phiNu_W_ttW_Boundaries_[1] ;  
  }

  if(runttbar_DL_integration_ || runttbar_DL_fakelep_integration_ || runttW_integration_){
  cosThetaNu_W_ttW_Boundaries_[0] = runConfig.cosThetaNu_ttau_ttbar_DL_ttW_Boundaries_[0] ;  
  cosThetaNu_W_ttW_Boundaries_[1] = runConfig.cosThetaNu_ttau_ttbar_DL_ttW_Boundaries_[1] ;  
  phiNu_W_ttW_Boundaries_[0] = runConfig.phiNu_ttau_ttbar_DL_ttW_Boundaries_[0] ;  
  phiNu_W_ttW_Boundaries_[1] = runConfig.phiNu_ttau_ttbar_DL_ttW_Boundaries_[1] ;  
  }

  force_nonzero_integral_ = runConfig.force_nonzero_integral_;

  sqrtS_ = 13000;
  sqrtS_ = runConfig.sqrtS_ ;

  Q_        = Physics::mtop + 0.5*Physics::mHiggs ;
  Q_ = runConfig.Q_;
 
//  Py_XDECREF(pVariable);
//  closePyConfig( pModule );
  
 
}

MGIntegration::~MGIntegration() {

}

const char* MGIntegration::getJacobianKindAsString( ) {
  const char *str;
  switch( jacobianKind_) {
    case NoJacobian:
      str = "NoJacobian";
      break;
    case MassTauTauBarJacobian:
      str = "MassTauTauBarJacobian";
      break; 
    default :
      str = "MassTauTauBarJacobian";
  }
  return str;
}

void MGIntegration::logConfig( ofstream &stream) { // pas forcement utile
//  int i;
  
  stream << "Configuration" << endl;
  stream << "  //" << endl;
  
  stream << "  Computation/Run" << endl;  
  stream << "    usePDF      : " << usePDF_ << endl; 
  stream << "    jacobianKind: " 
         << " " << getJacobianKindAsString() << endl;  
  stream << "    verbose     : " << verbose_ << endl; 
  stream << "  //" << endl;
  
  stream << "  MC Integration" << endl; 
  stream << "    NbrOfDim ttH  : " <<  nbrOfDim_ttH_ << endl; 
  stream << "    NbrOfDim ttZ   : " <<  nbrOfDim_ttZ_ << endl; 
  stream << "    NbrOfDim ttW   : " <<  nbrOfDim_ttW_ << endl;
  stream << "    NbrOfDim ttjets   : " <<  nbrOfDim_ttjets_ << endl;
  stream << "    NbrOfDim ttbar SL   : " <<  nbrOfDim_ttbar_SL_ << endl;
  stream << "    NbrOfDim ttbar DL   : " <<  nbrOfDim_ttbar_DL_ << endl;
  stream << "    NbrOfDim ttZ Zll   : " <<  nbrOfDim_ttZ_Zll_ << endl; 

  stream << "    NbrOfPoints ttH  : " <<  nbrOfPoints_ttH_ << endl; 
  stream << "    NbrOfPoints ttZ  : " <<  nbrOfPoints_ttZ_ << endl; 
  stream << "    NbrOfPoints ttW  : " <<  nbrOfPoints_ttW_ << endl; 
  stream << "    NbrOfPoints ttjets  : " <<  nbrOfPoints_ttjets_ << endl; 
  stream << "    NbrOfPoints ttbar SL  : " <<  nbrOfPoints_ttbar_SL_ << endl; 
  stream << "    NbrOfPoints ttbar DL  : " <<  nbrOfPoints_ttbar_DL_ << endl; 
  stream << "    NbrOfPoints ttZ Zll : " <<  nbrOfPoints_ttZ_Zll_ << endl; 

  stream << "    NbrOfDim ttH miss  : " <<  nbrOfDim_ttH_miss_ << endl; 
  stream << "    NbrOfDim ttZ miss   : " <<  nbrOfDim_ttZ_miss_ << endl; 
  stream << "    NbrOfDim ttjets miss   : " <<  nbrOfDim_ttjets_miss_ << endl;
  stream << "    NbrOfDim ttbar SL miss   : " <<  nbrOfDim_ttbar_SL_miss_ << endl;
  stream << "    NbrOfDim ttZ Zll miss   : " <<  nbrOfDim_ttZ_Zll_miss_ << endl; 

  stream << "    NbrOfPoints ttH miss  : " <<  nbrOfPoints_ttH_miss_ << endl; 
  stream << "    NbrOfPoints ttZ miss  : " <<  nbrOfPoints_ttZ_miss_ << endl; 
  stream << "    NbrOfPoints ttjets miss  : " <<  nbrOfPoints_ttjets_miss_ << endl; 
  stream << "    NbrOfPoints ttbar SL miss  : " <<  nbrOfPoints_ttbar_SL_miss_ << endl;
  stream << "    NbrOfPoints ttZ Zll miss  : " <<  nbrOfPoints_ttZ_Zll_miss_ << endl; 

  
  stream << "  //" << endl;
  
  stream << "  Physics" << endl;  
  stream << "    SqrtS  (Gev) : " <<  sqrtS_ << endl;   
  stream << "    Q      (Gev) : " <<  Q_ << endl;   
  stream << "    mTau   (Gev) : " <<  Physics::mTau << endl;
  stream << endl;
}


///////////////////////////////////////////////////
////////Helper functions for boundaries ///////////
///////////////////////////////////////////////////


double MGIntegration::u_plus( double mTauTau_2, double cosThetaTauLepTauHad ) {
  double u_plus = mTauTau_2 * TMath::Sqrt( ( mTauTau_2 / ( 4 * Physics::mTau2 ) ) - 1.0 ) 
                / TMath::Sqrt( 1.0 - cosThetaTauLepTauHad * cosThetaTauLepTauHad ); 
  return u_plus;
}


// GG ??? skip if no boundaries

bool MGIntegration::setCosThetaTauLepTauHadBoundaries( int leptonType,
              const TLorentzVector &Lep4P, const TLorentzVector &HadSys4P,
              double PTauLepBoundaries[2], double cosTauLepTauHadBoundaries[2],
              double mTauTauHZ2) {
  
  //To have the boundaries at zero if the determination fails
  PTauLepBoundaries[0] = 0;
  PTauLepBoundaries[1] = 0;
  cosTauLepTauHadBoundaries[0] = 0;
  cosTauLepTauHadBoundaries[1] = 0;  


  TVector3 el  = Lep4P.Vect().Unit();
  TVector3 epi = HadSys4P.Vect().Unit();
  double cosThetaTauLepPi = TMath::Cos( el.Angle(epi) );
//  double sinThetaTauLepPi = TMath::Sqrt(1.0 - cosThetaTauLepPi*cosThetaTauLepPi);
  
  // GG To start with a positive value of gamma2
  double  cosThetaTauLepTauHad =   cosThetaTauLepPi;	        
  
  //
  // Setting the integration boundaries in cosThetaTauLepTauHad and P_TauLep
  //
  
  // PTauLep upper limit
  double u_plus_ = u_plus( mTauTauHZ2, cosThetaTauLepTauHad );
  double ptaulepUpperLimit = std::min( t_plusLepton( Lep4P, leptonType ), u_plus_ );
  //
  double gamma2Max        = 0.;
  double ptaulepGamma2Max = 0.;
  double ptaulepLeft      = std::abs( t_minusLepton( Lep4P, leptonType ) );
  double ptaulepRight     = 0.;
 
  
  for( double ptaulep = ptaulepLeft ; ptaulep <= ptaulepUpperLimit ; ptaulep+= 1.0 ){
    double gamma2  = getGamma2( ptaulep, cosThetaTauLepTauHad, cosThetaTauLepPi, HadSys4P, mTauTauHZ2 );

    if( (gamma2 > 0.0) ) {
      ptaulepRight = ptaulep;
      if(gamma2 > gamma2Max){
        gamma2Max  = gamma2;
        ptaulepGamma2Max = ptaulep;
      }
    }
  }
  

  if( ptaulepRight == 0.0 || (ptaulepRight==ptaulepLeft) ) return false;


  double costhetatautauRight = -10.0;
  double costhetatautauLeft =  10.0;
  //
  for( double costhetatautau = -1; costhetatautau <= 1 ; costhetatautau += 0.0001 ){
 
    //
    double gamma2  = getGamma2( ptaulepGamma2Max, costhetatautau, cosThetaTauLepPi, HadSys4P, mTauTauHZ2 );

    if(gamma2 > 0.0 ) {
      if( costhetatautau < costhetatautauLeft )
        costhetatautauLeft = costhetatautau;     
      if (costhetatautau > costhetatautauRight)
        costhetatautauRight = costhetatautau;           
    }
  }
  

  if ( ( (costhetatautauLeft == 10.0) && (costhetatautauRight == -10.0) ) || (costhetatautauLeft==costhetatautauRight) ) return false;
  
  PTauLepBoundaries[0] = ptaulepLeft;
  PTauLepBoundaries[1] = ptaulepRight;
  cosTauLepTauHadBoundaries[0] = costhetatautauLeft;
  cosTauLepTauHadBoundaries[1] = costhetatautauRight;  
         
  return true;      
} 

//bool MGIntegration::setmTauTau_CosThetaTauLepTauHadBoundaries(){
void MGIntegration::setmTauTau_CosThetaTauLepTauHadBoundaries(){

  for(int perm=0; perm<nbrOfPermut_per_jet_; perm++){

    initVersors(perm);

    double mTauTau = Physics::mZ;
    double PTauLepBoundaries[2];
    double cosTauLepTauHadBoundaries[2];
    bool mass_ok = setCosThetaTauLepTauHadBoundaries( lepton_Tau_Type_,
						      evLep_Tau_4P_, evHadSys_Tau_4P_,
						      PTauLepBoundaries, cosTauLepTauHadBoundaries, 
						      mTauTau*mTauTau);
    
    if(force_nonzero_integral_){
      
      while( !mass_ok && mTauTau<500 ){
	
	mTauTau++;
	
	mass_ok = setCosThetaTauLepTauHadBoundaries( lepton_Tau_Type_,
						     evLep_Tau_4P_, evHadSys_Tau_4P_,
						     PTauLepBoundaries, cosTauLepTauHadBoundaries, 
						     mTauTau*mTauTau);
	
      }
      
    }


    if(integration_type_ == integration_type_wo_miss && !force_missing_jet_integration_){

      mTauTau_ttZ_[perm] = mTauTau ;
      PTauLep_ttZ_Lower_[perm] = PTauLepBoundaries[0];
      PTauLep_ttZ_Upper_[perm] = PTauLepBoundaries[1];
      cosTheta_diTau_ttZ_Lower_[perm] = cosTauLepTauHadBoundaries[0];
      cosTheta_diTau_ttZ_Upper_[perm] = cosTauLepTauHadBoundaries[1];
      include_perm_ttZ_[perm] = mass_ok;

    }

    else if(integration_type_ == integration_type_w_miss || force_missing_jet_integration_){

      for(int i_jet=0; i_jet<n_lightJets_; i_jet++){
	int i = nbrOfPermut_per_jet_*i_jet + perm;

	mTauTau_ttZ_[i] = mTauTau ;
	PTauLep_ttZ_Lower_[i] = PTauLepBoundaries[0];
	PTauLep_ttZ_Upper_[i] = PTauLepBoundaries[1];
	cosTheta_diTau_ttZ_Lower_[i] = cosTauLepTauHadBoundaries[0];
	cosTheta_diTau_ttZ_Upper_[i] = cosTauLepTauHadBoundaries[1];
	include_perm_ttZ_[i] = mass_ok;

      }
    }

    mTauTau = Physics::mHiggs;
  
    mass_ok = setCosThetaTauLepTauHadBoundaries( lepton_Tau_Type_,
						 evLep_Tau_4P_, evHadSys_Tau_4P_,
						 PTauLepBoundaries, cosTauLepTauHadBoundaries, 
						 mTauTau*mTauTau);
    
    if(force_nonzero_integral_){

      while( !mass_ok && mTauTau<500 ){
	
	mTauTau++;
	
	mass_ok = setCosThetaTauLepTauHadBoundaries( lepton_Tau_Type_,
						     evLep_Tau_4P_, evHadSys_Tau_4P_,
						     PTauLepBoundaries, cosTauLepTauHadBoundaries, 
						     mTauTau*mTauTau);
	
      }
      
    }
    
    if(integration_type_ == integration_type_wo_miss && !force_missing_jet_integration_){

      mTauTau_ttH_[perm] = mTauTau ;
      PTauLep_ttH_Lower_[perm] = PTauLepBoundaries[0];
      PTauLep_ttH_Upper_[perm] = PTauLepBoundaries[1];
      cosTheta_diTau_ttH_Lower_[perm] = cosTauLepTauHadBoundaries[0];
      cosTheta_diTau_ttH_Upper_[perm] = cosTauLepTauHadBoundaries[1];
      include_perm_ttH_[perm] = mass_ok;

    }

    else if(integration_type_ == integration_type_w_miss || force_missing_jet_integration_){

      for(int i_jet=0; i_jet<n_lightJets_; i_jet++){
	int i = nbrOfPermut_per_jet_*i_jet + perm;

	mTauTau_ttH_[i] = mTauTau ;
	PTauLep_ttH_Lower_[i] = PTauLepBoundaries[0];
	PTauLep_ttH_Upper_[i] = PTauLepBoundaries[1];
	cosTheta_diTau_ttH_Lower_[i] = cosTauLepTauHadBoundaries[0];
	cosTheta_diTau_ttH_Upper_[i] = cosTauLepTauHadBoundaries[1];
	include_perm_ttH_[i] = mass_ok;

      }
      
    }
    
  }


  if(nbrOfPermut_<nbrOfPermutMax){
      for(int i=nbrOfPermut_; i<nbrOfPermutMax; i++){
	mTauTau_ttZ_[i] = NAN ;
	PTauLep_ttZ_Lower_[i] = NAN;
	PTauLep_ttZ_Upper_[i] = NAN;
	cosTheta_diTau_ttZ_Lower_[i] = NAN;
	cosTheta_diTau_ttZ_Upper_[i] = NAN;

	mTauTau_ttH_[i] = NAN ;
	PTauLep_ttH_Lower_[i] = NAN;
	PTauLep_ttH_Upper_[i] = NAN;
	cosTheta_diTau_ttH_Lower_[i] = NAN;
	cosTheta_diTau_ttH_Upper_[i] = NAN;
      }
  }
    

}

//--- t_plus and s_plus
double MGIntegration::t_plus( const TLorentzVector &vec4P ){

  double E_Pi	= vec4P.E();
  double P_Pi	= vec4P.P();

  double mPi_2 = (E_Pi * E_Pi) - (P_Pi * P_Pi);
  double mTau_2 = Physics::mTau2;

  double tPlus = ( (mTau_2 + mPi_2) * P_Pi + (mTau_2 - mPi_2) * E_Pi ) 
                  / (2 * mPi_2);

  return tPlus;
}

//--- t_minus and s_minus
double MGIntegration::t_minus( const TLorentzVector &vec4P ){
    
  double E_Pi	= vec4P.E();
  double P_Pi	= vec4P.P();

  double mPi_2 = (E_Pi * E_Pi) - (P_Pi * P_Pi);
  double mTau_2 = Physics::mTau2;

  double tMinus = ( (mTau_2 + mPi_2) * P_Pi - (mTau_2 - mPi_2) * E_Pi ) 
                  / (2 * mPi_2);

  return tMinus;
}

void MGIntegration::setTauHadMomentumBoundaries( ) {

  if(evHadSys_Tau_4P_.M()>Physics::mTau){
    //Can happen if tau 4P mis-measured
    //Tau Had TF would anyway be zero

    for( int i=0; i<nbrOfPermut_; i++) { 
      include_perm_ttH_[i] = 0;
      include_perm_ttZ_[i] = 0;
      include_perm_ttW_[i] = 0;
      include_perm_ttbar_DL_[i] = 0;
      include_perm_ttbar_DL_fakelep_tlep_[i] = 0;
      include_perm_ttbar_DL_fakelep_ttau_[i] = 0;
    }

  }

  PTauHad_ttbar_DL_ttW_Boundaries_[0] = t_minus( evHadSys_Tau_4P_ );
  PTauHad_ttbar_DL_ttW_Boundaries_[1] = t_plus( evHadSys_Tau_4P_ );
  return;

}

//--- t_minus and s_minus
double MGIntegration::t_minusLepton( const TLorentzVector &vec4P, int leptonType ){

  double P_Lep	= vec4P.P();
  double m_Lep_2 = 0;
  if( abs(leptonType) == 11) 
    m_Lep_2 = Physics::mElectron*Physics::mElectron ;
  else if( abs(leptonType) == 13)
    m_Lep_2 = Physics::mMuon*Physics::mMuon;
  double m_Tau_2 = Physics::mTau2;
  double E_Lep	= TMath::Sqrt( P_Lep * P_Lep + m_Lep_2 );
  
  double tMinus = ( ( m_Tau_2 + m_Lep_2 ) * P_Lep - E_Lep * ( m_Tau_2 - m_Lep_2 ) ) / ( 2 * m_Lep_2 ) ;
  return tMinus;
}

//--- t_minus 
double MGIntegration::t_plusLepton( const TLorentzVector &vec4P, int leptonType){

  double P_Lep	= vec4P.P();
  double m_Lep_2 = 0;
  if( abs(leptonType) == 11) 
    m_Lep_2 = Physics::mElectron*Physics::mElectron ;
  else if( abs(leptonType) == 13)
    m_Lep_2 = Physics::mMuon*Physics::mMuon; 
  double m_Tau_2 = Physics::mTau2;
  double E_Lep	= TMath::Sqrt( P_Lep * P_Lep + m_Lep_2 );
  
  double tPlus = ( m_Tau_2 / m_Lep_2 ) * E_Lep ;
  return tPlus;
}

pair<double,double> MGIntegration::getE_Quark_Boundaries(TLorentzVector jet4P, TString flavor) {

  double quantile = CI_TFJet_;

  double Elow1 = 0.;
  double Ehigh1 = 0.;
  double Elow2 = 0.;
  double Ehigh2 = 0.;

  double mu = 0. ;
  double sigma = 0. ;

  double chi2cut = TMath::ChisquareQuantile(quantile ,1);

  vector<double> jetTFresolution;
  vector<double> jetTFmean;

  double Egen = jet4P.E();
  if(use_pT_TFJet_)
    Egen=jet4P.Pt();
 
  jetTFresolution = getJetTFresolution(Egen,jet4P.Eta(),flavor);
  jetTFmean = getJetTFmean(Egen,jet4P.Eta(),flavor);

  double Erec = Egen;

  mu = jetTFmean[0];
  sigma = jetTFresolution[0];

  double GeVStep = sigma/10. ;
  while( pow((Erec-mu)/sigma,2) < chi2cut && Egen>0. && Egen<8000.){
    Egen -= GeVStep;
    jetTFresolution = getJetTFresolution(Egen,jet4P.Eta(),flavor);
    jetTFmean = getJetTFmean(Egen,jet4P.Eta(),flavor);
    mu=jetTFmean[0];
    sigma = jetTFresolution[0];
  }
  Elow1 = Egen;

  Egen = Erec;
  jetTFresolution = getJetTFresolution(Egen,jet4P.Eta(),flavor);
  jetTFmean = getJetTFmean(Egen,jet4P.Eta(),flavor);
  mu=jetTFmean[0];
  sigma = jetTFresolution[0];

  while( pow((Erec-mu)/sigma,2) < chi2cut && Egen>0. && Egen<8000.){
    Egen += GeVStep;
    jetTFresolution = getJetTFresolution(Egen,jet4P.Eta(),flavor);
    jetTFmean = getJetTFmean(Egen,jet4P.Eta(),flavor);
    mu=jetTFmean[0];
    sigma = jetTFresolution[0];
  }
  Ehigh1 = Egen;

  Egen = Erec;
  jetTFresolution = getJetTFresolution(Egen,jet4P.Eta(),flavor);
  jetTFmean = getJetTFmean(Egen,jet4P.Eta(),flavor);
  mu=jetTFmean[1];
  sigma = jetTFresolution[1];

  GeVStep = sigma/10. ;
  while( pow((Erec-mu)/sigma,2) < chi2cut && Egen>0. && Egen<8000.){
    Egen -= GeVStep;
    jetTFresolution = getJetTFresolution(Egen,jet4P.Eta(),flavor);
    jetTFmean = getJetTFmean(Egen,jet4P.Eta(),flavor);
    mu=jetTFmean[1];
    sigma = jetTFresolution[1];
  }
  Elow2 = Egen;

  Egen = Erec;
  jetTFresolution = getJetTFresolution(Egen,jet4P.Eta(),flavor);
  jetTFmean = getJetTFmean(Egen,jet4P.Eta(),flavor);
  mu=jetTFmean[1];
  sigma = jetTFresolution[1];

  while( pow((Erec-mu)/sigma,2) < chi2cut && Egen>0. && Egen<8000.){
    Egen += GeVStep;
    jetTFresolution = getJetTFresolution(Egen,jet4P.Eta(),flavor);
    jetTFmean = getJetTFmean(Egen,jet4P.Eta(),flavor);
    mu=jetTFmean[1];
    sigma = jetTFresolution[1];

  }
  Ehigh2 = Egen;

  double Emin=min(Elow1,Elow2);
  if(Emin<0)
    Emin=0;
  double Emax=max(Ehigh1,Ehigh2);

  pair<double,double> E_Bound=make_pair(Emin,Emax);
  if(use_pT_TFJet_){
    E_Bound.first*=cosh(jet4P.Eta());
    E_Bound.second*=cosh(jet4P.Eta());
  }

  return E_Bound;


}

pair<double,double> MGIntegration::getE_fakelep_Boundaries(TLorentzVector lep4P) {


  double quantile = CI_TFJet_;

  double Elow = 0.;
  double Ehigh = 0.;

  double mu = 0. ;
  double sigma = 0. ;

  double chi2cut = TMath::ChisquareQuantile(quantile ,1);

  double Erec = lep4P.Pt();
  double Egen = Erec;

  mu = getFakeLepTFmean(Egen);
  sigma = getFakeLepTFresolution(Egen);


  double GeVStep = sigma/10. ;
  while( pow((Erec-mu)/sigma,2) < chi2cut && Egen>0. && Egen<8000.){
    Egen -= GeVStep;
    mu = getFakeLepTFmean(Egen);
    sigma = getFakeLepTFresolution(Egen);
  }
  Elow = Egen;

  Egen = Erec;

  mu = getFakeLepTFmean(Egen);
  sigma = getFakeLepTFresolution(Egen);
  while( pow((Erec-mu)/sigma,2) < chi2cut && Egen>0. && Egen<8000.){
    Egen += GeVStep;
    mu = getFakeLepTFmean(Egen);
    sigma = getFakeLepTFresolution(Egen);
  }
  Ehigh = Egen;


  double Emin=Elow;
  if(Emin<0)
    Emin=0;
  double Emax=Ehigh;

  pair<double,double> E_Bound=make_pair(Emin,Emax);
  E_Bound.first*=cosh(lep4P.Eta());
  E_Bound.second*=cosh(lep4P.Eta());

  return E_Bound;


}

void MGIntegration::setE_lightQuark_Boundaries() {

  if(integration_type_ == integration_type_wo_miss && !force_missing_jet_integration_){

    pair<double,double> E_Bound=getE_Quark_Boundaries(evJet1_4P_,"light");
    for(int i=0; i<nbrOfPermut_; i++){
      EQuark1_Lower_[i] = E_Bound.first;
      EQuark1_Upper_[i] = E_Bound.second;
    }
    for(int i=nbrOfPermut_; i<nbrOfPermutMax; i++){
      EQuark1_Lower_[i] = NAN;
      EQuark1_Upper_[i] = NAN;
    }

  }

  else if(integration_type_ == integration_type_w_miss || force_missing_jet_integration_){

    for(int i=0; i<nbrOfPermut_; i++){
      initVersors_miss(i);
      pair<double,double> E_Bound=getE_Quark_Boundaries(evJet1_4P_,"light");
      EQuark1_Lower_[i] = E_Bound.first;
      EQuark1_Upper_[i] = E_Bound.second;
    }
    if(nbrOfPermut_<nbrOfPermutMax){
      for(int i=nbrOfPermut_; i<nbrOfPermutMax; i++){
	EQuark1_Lower_[i] = NAN;
	EQuark1_Upper_[i] = NAN;
      }
    }

  }

}

void MGIntegration::checkCompatibility_TopHad(){

  //Check if M(bjj) compatible with Mtop within uncertainty

  pair<double,double> boundsJet1 = getE_Quark_Boundaries(evJet1_4P_, "light");
  pair<double,double> boundsJet2 = getE_Quark_Boundaries(evJet2_4P_, "light");

  TLorentzVector jet1Low, jet2Low;
  jet1Low.SetPtEtaPhiE ( boundsJet1.first/cosh(evJet1_4P_.Eta()), evJet1_4P_.Eta(), evJet1_4P_.Phi(), boundsJet1.first );
  jet2Low.SetPtEtaPhiE ( boundsJet2.first/cosh(evJet2_4P_.Eta()), evJet2_4P_.Eta(), evJet2_4P_.Phi(), boundsJet2.first );

  TLorentzVector jet1High, jet2High;
  jet1High.SetPtEtaPhiE ( boundsJet1.second/cosh(evJet1_4P_.Eta()), evJet1_4P_.Eta(), evJet1_4P_.Phi(), boundsJet1.second );
  jet2High.SetPtEtaPhiE ( boundsJet2.second/cosh(evJet2_4P_.Eta()), evJet2_4P_.Eta(), evJet2_4P_.Phi(), boundsJet2.second );

  double Mt = Physics::mtop;

  for(int perm=0; perm<nbrOfPermut_; perm++){
    
    initVersors(perm);
    
    pair<double,double> boundsBJet = getE_Quark_Boundaries(evBJet_hadtop_4P_, "b");

    TLorentzVector BJetLow, BJetHigh;

    BJetLow.SetPtEtaPhiE ( boundsBJet.first/cosh(evBJet_hadtop_4P_.Eta()), evBJet_hadtop_4P_.Eta(), evBJet_hadtop_4P_.Phi(), boundsBJet.first );
    BJetHigh.SetPtEtaPhiE ( boundsBJet.second/cosh(evBJet_hadtop_4P_.Eta()), evBJet_hadtop_4P_.Eta(), evBJet_hadtop_4P_.Phi(), boundsBJet.second );

    double massLow = (jet1Low+jet2Low+BJetLow).M();
    double massHigh = (jet1High+jet2High+BJetHigh).M();

    if(!(massLow<=Mt && Mt<=massHigh)){
      include_perm_ttH_[perm] = 0;
      include_perm_ttZ_[perm] = 0;
      include_perm_ttjets_[perm] = 0;
      include_perm_ttbar_SL_[perm] = 0;
      include_perm_ttZ_Zll_[perm] = 0;
    }

  }


}

void MGIntegration::checkCompatibility_TopHad_missing_jet(){

  //Check if M(bj)<=Mtop within uncertainty

  double Mt = Physics::mtop;

  for(int perm=0; perm<nbrOfPermut_; perm++){
    
    initVersors_miss(perm);

    /*std::cout << "evJet1_4P_ : " << std::endl; //TEMPORAIRE
    evJet1_4P_.Print(); //TEMPORAIRE
    std::cout << "evBJet_hadtop_4P_ : " << std::endl; //TEMPORAIRE
    evBJet_hadtop_4P_.Print(); //TEMPORAIRE */
    
    pair<double,double> boundsJet1 = getE_Quark_Boundaries(evJet1_4P_, "light");

    TLorentzVector jet1Low;
    jet1Low.SetPtEtaPhiE ( boundsJet1.first/cosh(evJet1_4P_.Eta()), evJet1_4P_.Eta(), evJet1_4P_.Phi(), boundsJet1.first );

    pair<double,double> boundsBJet = getE_Quark_Boundaries(evBJet_hadtop_4P_, "b");

    TLorentzVector BJetLow;
    BJetLow.SetPtEtaPhiE ( boundsBJet.first/cosh(evBJet_hadtop_4P_.Eta()), evBJet_hadtop_4P_.Eta(), evBJet_hadtop_4P_.Phi(), boundsBJet.first );

    double MbjLow=(BJetLow+jet1Low).M();
    //std::cout << "MbjLow : " << MbjLow << std::endl;  //TEMPORAIRE
    
    if(Mt<MbjLow){
      include_perm_ttH_[perm] = 0;
      include_perm_ttZ_[perm] = 0;
      include_perm_ttjets_[perm] = 0;
      include_perm_ttbar_SL_[perm] = 0;
      include_perm_ttZ_Zll_[perm] = 0;
    }

  }

}

void MGIntegration::checkCompatibility_TopLep(){

  //Check if M(b+l)<=Mtop within uncertainty
  double Mt = Physics::mtop;

  for(int perm=0; perm<nbrOfPermut_; perm++){
    
    if(integration_type_==integration_type_wo_miss)
      initVersors(perm);
    else
      initVersors_miss(perm);

    pair<double,double> boundsBJet = getE_Quark_Boundaries(evBJet_leptop_4P_, "b");

    TLorentzVector BJetLow;
    BJetLow.SetPtEtaPhiE ( boundsBJet.first/cosh(evBJet_leptop_4P_.Eta()), evBJet_leptop_4P_.Eta(), evBJet_leptop_4P_.Phi(), boundsBJet.first );

    double MblLow=(BJetLow+evLep_top_4P_).M();

    if(Mt<MblLow){
      include_perm_ttH_[perm] = 0;
      include_perm_ttZ_[perm] = 0;
      include_perm_ttW_[perm] = 0;
      include_perm_ttjets_[perm] = 0;
      include_perm_ttbar_SL_[perm] = 0;
      include_perm_ttbar_DL_[perm] = 0;
      include_perm_ttZ_Zll_[perm] = 0;
    }

  }

}

void MGIntegration::checkCompatibility_TopLep_fakelep(){

  //Check if M(l+l)<=Mtop within uncertainty
  double Mt = Physics::mtop;

  for(int perm=0; perm<nbrOfPermut_; perm++){
    
    if(integration_type_==integration_type_wo_miss)
      initVersors(perm);
    else
      initVersors_miss(perm);

    pair<double,double> boundsFakeLep = getE_fakelep_Boundaries(evLep_Tau_4P_);

    TLorentzVector FakeLepLow;
    FakeLepLow.SetPtEtaPhiE ( boundsFakeLep.first/cosh(evLep_Tau_4P_.Eta()), evLep_Tau_4P_.Eta(), evLep_Tau_4P_.Phi(), boundsFakeLep.first );

    double MllLow=(FakeLepLow+evLep_top_4P_).M();

    if(Mt<MllLow){
      include_perm_ttbar_DL_fakelep_tlep_[perm] = 0;
    }

  }

}

void MGIntegration::checkCompatibility_TopTau(){

  //Check if M(b+tau)<=Mtop within uncertainty
  double Mt = Physics::mtop;

  for(int perm=0; perm<nbrOfPermut_; perm++){
    
    if(integration_type_==integration_type_wo_miss)
      initVersors(perm);
    else
      initVersors_miss(perm);

    pair<double,double> boundsBJet = getE_Quark_Boundaries(evBJet_hadtop_4P_, "b");

    TLorentzVector BJetLow;
    BJetLow.SetPtEtaPhiE ( boundsBJet.first/cosh(evBJet_hadtop_4P_.Eta()), evBJet_hadtop_4P_.Eta(), evBJet_hadtop_4P_.Phi(), boundsBJet.first );

    double MblLow=(BJetLow+evHadSys_Tau_4P_).M();

    if(Mt<MblLow){
      include_perm_ttW_[perm] = 0;
      include_perm_ttbar_DL_[perm] = 0;
    }

  }

}

void MGIntegration::checkCompatibility_TopTau_fakelep(){

  //Check if M(l+tau)<=Mtop within uncertainty
  double Mt = Physics::mtop;

  for(int perm=0; perm<nbrOfPermut_; perm++){
    
    if(integration_type_==integration_type_wo_miss)
      initVersors(perm);
    else
      initVersors_miss(perm);

    pair<double,double> boundsFakeLep = getE_fakelep_Boundaries(evLep_Tau_4P_);

    TLorentzVector FakeLepLow;
    FakeLepLow.SetPtEtaPhiE ( boundsFakeLep.first/cosh(evLep_Tau_4P_.Eta()), evLep_Tau_4P_.Eta(), evLep_Tau_4P_.Phi(), boundsFakeLep.first );


    double MltauLow=(FakeLepLow+evHadSys_Tau_4P_).M();

    if(Mt<MltauLow){
      include_perm_ttbar_DL_fakelep_ttau_[perm] = 0;
    }

  }

}

void MGIntegration::checkCompatibility_Zll(){

  //Check if lepton reconstructed with Z mass constraint is compatible with reco tau

  for(int perm=0; perm<nbrOfPermut_; perm++){
    
    if(integration_type_==integration_type_wo_miss)
      initVersors(perm);
    else
      initVersors_miss(perm);

    TLorentzVector lep_Z_faketau_4P;

    double MZ = Physics::mZ;
    double CosTheta_ll = TMath::Cos(evLep_Tau_4P_.Angle(evHadSys_Tau_4P_.Vect()));
    double E_gen = MZ*MZ / (2*evLep_Tau_4P_.E()*(1-CosTheta_ll));  
    double pt_gen=E_gen/TMath::CosH(evHadSys_Tau_4P_.Eta());
    
    double pt_reco=evHadSys_Tau_4P_.Pt();

    double mu = getFakeTauLepTFmean(pt_gen,lepton_Tau_Type_);
    double sigma = getFakeTauLepTFresolution(pt_gen,lepton_Tau_Type_);

    double chi = (pt_reco-mu)/sigma;
    
    if(abs(chi)>5) include_perm_ttZ_Zll_[perm] = 0;

  }


}


///////////////////////////////////////////////////
////////Helper functions for di-tau system ////////
///////////////////////////////////////////////////


double MGIntegration::getPTauHad( double _mtautau2, double PTauLep, double cosThetaTauLepTauHad ) {
  double mTau2 = Physics::mTau2;
  double PTauLep2 = PTauLep * PTauLep;
  double M_2 = 0.5 * _mtautau2 - mTau2;
  double E_TauLep = TMath::Sqrt( PTauLep * PTauLep + mTau2 );
  double sinThetaTauLepTauHad = TMath::Sqrt( 1.0 - cosThetaTauLepTauHad * cosThetaTauLepTauHad );
  double sin_2 = sinThetaTauLepTauHad * sinThetaTauLepTauHad;
  double PTauHad = ( M_2 * PTauLep * cosThetaTauLepTauHad 
                      +  E_TauLep * TMath::Sqrt( M_2 * M_2 - mTau2 * (mTau2 + PTauLep2 * sin_2 ) )  
                    ) / (mTau2 + PTauLep2  * sin_2 );
  return PTauHad;
}

double MGIntegration::getGamma2( double PTauLep, double cosThetaTauLepTauHad, double cosThetaTauLepPi, 
        const TLorentzVector &HadSys4P, double mTauTauHZ2 ) {
  
  double M_2 = 0.5 * mTauTauHZ2 - Physics::mTau2; 
  //
  double sinThetaTauLepTauHad2 = 1.0 - cosThetaTauLepTauHad*cosThetaTauLepTauHad;
  //
  double PTauLep2 = PTauLep * PTauLep;
  double ETauLep = TMath::Sqrt( PTauLep2 + Physics::mTau2 );

  double PTauHad = ( M_2 * PTauLep * cosThetaTauLepTauHad 
                  +  ETauLep * TMath::Sqrt( M_2 * M_2 - Physics::mTau2 * (Physics::mTau2 + PTauLep2 * sinThetaTauLepTauHad2 ) )  
                ) / (Physics::mTau2 + PTauLep2  * sinThetaTauLepTauHad2 );   

//  double ETauHad = TMath::Sqrt( Physics::mTau2 + PTauHad * PTauHad );
  double cosThetaTauHadPi = getCosThetaTauHadPi( PTauHad, HadSys4P );
  
  double alpha, beta, gamma, gamma2;
  const char *err;
  getAlphaBetaGamma( cosThetaTauLepTauHad, cosThetaTauLepPi, cosThetaTauHadPi,
                     alpha, beta, gamma, gamma2, &err);
  return gamma2;
}

//--- cosThetaTauHadPi 
// ??? GG verify arg of the fct
double MGIntegration::getCosThetaTauHadPi( double P_TauHad, const TLorentzVector &HadSys4P) {
  double mTau_2 = Physics::mTau2;
  double E_Tau 	= TMath::Sqrt( P_TauHad*P_TauHad + mTau_2 );

  double E_Pi	= HadSys4P.E();
  double P_Pi	= HadSys4P.P();
  double mPi_2  = ( E_Pi*E_Pi - P_Pi*P_Pi );
 
  double cosThetaTauHadPi = ( 2*E_Tau*E_Pi - (mTau_2 + mPi_2) )
                       / ( 2*P_TauHad*P_Pi );

  return cosThetaTauHadPi;
}

void MGIntegration::getAlphaBetaGamma(double cosThetaTauLepTauHad, 
        double cosThetaTauLepPi, double cosThetaTauHadPi,
        double &alpha, double &beta, double &gamma, double &gamma2, const char **error){
// GG ???
  double sinThetaTauLepPi2 = 1.0 - cosThetaTauLepPi*cosThetaTauLepPi;
    
  alpha = ( ( cosThetaTauLepTauHad - cosThetaTauLepPi * cosThetaTauHadPi     ) / ( sinThetaTauLepPi2 ) );  
  beta  = ( ( cosThetaTauHadPi     - cosThetaTauLepPi * cosThetaTauLepTauHad ) / ( sinThetaTauLepPi2 ) );
 
  gamma2 = ( 1 - ( alpha * alpha ) - ( beta * beta ) - 2 * alpha * beta * cosThetaTauLepPi );  
  gamma =  TMath::Sqrt( gamma2 );
    if( TMath::IsNaN( gamma ) ) { if (*error == 0) *error = "Gamma";}  
}


///////////////////////////////////////////////////
////////    Helper functions for tops      ////////
///////////////////////////////////////////////////


///Can also be used for nu from leptonic top with CosTheta_qq->CosTheta_lnu and Eq->El
double MGIntegration::getEqbar_Enu(double CosTheta_qq, double Eq){

  double MW = Physics::mW;
  double Eqbar= MW*MW / (2*Eq*(1-CosTheta_qq));
  return Eqbar;

}

double MGIntegration::getEb(TLorentzVector W, TLorentzVector bjet){

  double Mt = Physics::mtop;
  double MW = Physics::mW;
  double mb = Physics::mb;

  double deltaM2=0.5*(Mt*Mt-MW*MW-mb*mb);
  double EW=W.E();
  TVector3 eb=(bjet.Vect()).Unit();
  double Web=W.Vect()*eb;

  double Eb_plus= ( EW*deltaM2 + fabs(Web)*TMath::Sqrt(deltaM2*deltaM2 - mb*mb*(EW*EW - Web*Web)) ) / (EW*EW - Web*Web);
  double Eb_minus= ( EW*deltaM2 - fabs(Web)*TMath::Sqrt(deltaM2*deltaM2 - mb*mb*(EW*EW - Web*Web)) ) / (EW*EW - Web*Web);

  if(Web>0){
    if(Eb_plus > deltaM2/EW)
      return Eb_plus;
    else if(Eb_minus > deltaM2/EW)
      return Eb_minus;
  }
  
  else{
    if(Eb_plus < deltaM2/EW)
      return Eb_plus;
    else if(Eb_minus < deltaM2/EW)
      return Eb_minus;
  }

  return 0;

}


///////////////////////////////////////////////////
////////   Helper functions for ttjets     ////////
///////////////////////////////////////////////////


double MGIntegration::getPFake_upperBoundary_CI(double quantile,
						const TLorentzVector &evHadSys4P){

  //Same as for fake TF
  double mu=-0.331309;
  double sigma = 0.179943 ;
  double chi2cut = TMath::ChisquareQuantile(quantile ,1);
  
  double Ehigh = evHadSys4P.E() / (1+mu-sigma*sqrt(chi2cut));
  
  return Ehigh;
  
}

void MGIntegration::setPOutFakeBoundaries( ) {

  POutFakeBoundaries_ttjets_[0] = 0.2 * evHadSys_Tau_4P_.E();
  POutFakeBoundaries_ttjets_[1] = getPFake_upperBoundary_CI(CI_TFJet_,evHadSys_Tau_4P_);
  
}


///////////////////////////////////////////////////
///////////    Transfer functions       ///////////
///////////////////////////////////////////////////


//--- Jet TFs 

vector<double> MGIntegration::getJetTFresolution( double Egen, double eta, TString flavor) const {

  double cst1[3] = {0,0,0}, cst2[3] = {0,0,0};

  if(flavor=="light"){

    if (use_pT_TFJet_) {
      if ( fabs(eta) < 1.5){
	cst1[0] = 0;
	cst1[1] = 1.53;
	cst1[2] = 0.13;
	
	cst2[0] = 0;
	cst2[1] = 0.78;
	cst2[2] = 0.05;      
      } else if( fabs(eta)>1.5 && fabs(eta)<3.0){
	cst1[0] = 4.26;
	cst1[1] = 0.58;
	cst1[2] = 0.0;
	
	cst2[0] = 2.56;
	cst2[1] = 1.99;
	cst2[2] = 0.0;
      } else {
	cst1[0] = 0.0;
	cst1[1] = 0.58;
	cst1[2] = 0.12;
	
	cst2[0] = 0.0;
	cst2[1] = 0.6;
	cst2[2] = 0.05;
      }
    } else {
      if(fabs(eta)<1.5){
	cst1[0] = 0.06;
	cst1[1] = 0.99;
	cst1[2] = 0.02;
	
	cst2[0] = 10.50;
	cst2[1] = 2.40;
	cst2[2] = 0.0;
      } else if( fabs(eta)>1.5 && fabs(eta)<3.0){
	cst1[0] = 19.65;
	cst1[1] = 3.56;
	cst1[2] = 0.0;
	
	cst2[0] = 0.0;
	cst2[1] = 1.22;
	cst2[2] = 0.08;
      } else{
	cst1[0] = 60.65;
	cst1[1] = 0.0;
	cst1[2] = 0.1;
	
	cst2[0] = 0.01;
	cst2[1] = 1.70;
	cst2[2] = 0.04;
      }
    }

  }


  else if(flavor=="b"){

    if(use_pT_TFJet_) {
      if (fabs(eta)<1.5) {
	cst1[0] = 0.0;
	cst1[1] = 1.07;
	cst1[2] = 0.05;
	
	cst2[0] = 0.0;
	cst2[1] = 0.0;
	cst2[2] = 0.23;  
      } else if(fabs(eta)>1.5 && fabs(eta)<3.0){
	cst1[0] = 0.;
	cst1[1] = 1.14;
	cst1[2] = 0.;
	
	cst2[0] = 3.86;
	cst2[1] = 0.47;
	cst2[2] = 0.21;  
      } else {
	cst1[0] = 0.0;
	cst1[1] = 0.0;
	cst1[2] = 0.0;
	
	cst2[0] = 0.0;
	cst2[1] = 0.0;
	cst2[2] = 0.0;  
      }
    } else {
      if (fabs(eta)<1.5) {
	cst1[0] = 0.0;
	cst1[1] = 0.89;
	cst1[2] = 0.16;
	
	cst2[0] = 0.0;
	cst2[1] = 0.91;
	cst2[2] = 0.04;  
      } else if (fabs(eta)>1.5 && fabs(eta)<3.0) {
	cst1[0] = 0.0;
	cst1[1] = 0.0;
	cst1[2] = 0.19;
	
	cst2[0] = 0.0;
	cst2[1] = 0.0;
	cst2[2] = 0.12;  
      } else{
	cst1[0] = 0.0;
	cst1[1] = 0.0;
	cst1[2] = 0.0;
	
	cst2[0] = 0.0;
	cst2[1] = 0.0;
	cst2[2] = 0.0;  
      } 
    }
           
  }



  vector<double> res;
  res.push_back( sqrt( cst1[0]*cst1[0] + ( cst1[1]*cst1[1]*fabs(Egen) ) + (cst1[2]*Egen)*(cst1[2]*Egen) ) );
  res.push_back( sqrt( cst2[0]*cst2[0] + ( cst2[1]*cst2[1]*fabs(Egen) ) + (cst2[2]*Egen)*(cst2[2]*Egen) ) );

  return res;


}

vector<double> MGIntegration::getJetTFmean( double Egen, double eta, TString flavor) const {

  double cst1[2] = {0,0}, cst2[2] = {0,0};

  if(flavor=="light"){
    
    if(use_pT_TFJet_){
      
      if(fabs(eta)<1.5){   
	cst1[0] = -2.33;
	cst1[1] = 0.99;
	
	cst2[0] = 9.39;
	cst2[1] = 0.94;
      }
      
      else if(fabs(eta)>1.5 && fabs(eta)<3.0){
	cst1[0] = 10.79;
	cst1[1] = 0.94;
	
	cst2[0] = 51.65;
	cst2[1] = 0.82;
      } else{
	cst1[0] = -1.83;
	cst1[1] = 0.64;
	
	cst2[0] = 3.52;
	cst2[1] = 0.62;
      }
      
    } else{
      if(fabs(eta)<1.5){   
	cst1[0] = 1.95;
	cst1[1] = 0.94;
	
	cst2[0] = -4.12;
	cst2[1] = 1.01;
      } else if(fabs(eta)>1.5 && fabs(eta)<3.0){
	cst1[0] = 58.65;
	cst1[1] = 0.77;
	
	cst2[0] = -4.46;
	cst2[1] = 0.97;
      } else {
	cst1[0] = 122.31;
	cst1[1] = 0.58;
	
	cst2[0] = 56.86;
	cst2[1] = 0.59;
      }
    }
    
  }

  else if(flavor=="b"){

    if (use_pT_TFJet_) {
      if(fabs(eta)<1.5) {
	cst1[0] = -6.24;
	cst1[1] = 0.99;
	
	cst2[0] = -8.12;
	cst2[1] = 0.97;      
      } else if (fabs(eta)>1.5 && fabs(eta)<2.5) {
	cst1[0] = -6.24;
	cst1[1] = 0.98;
	
	cst2[0] = -10.31;
	cst2[1] = 0.94;    
      } else {
	cst1[0] = 0.0;
      cst1[1] = 0.0;
      
      cst2[0] = 0.0;
      cst2[1] = 0.0;
      }
    } else {
      if (fabs(eta)<1.5) {
	cst1[0] = -0.73;
	cst1[1] = 0.91;
	
	cst2[0] = -5.91;
	cst2[1] = 1.0;    
      } else if (fabs(eta)>1.5 && fabs(eta)<2.5) {
	cst1[0] = 26.42;
	cst1[1] = 0.77;
	
	cst2[0] = -28.78;
	cst2[1] = 1.03;   
      } else {
	cst1[0] = 0.0;
	cst1[1] = 0.0;
	
	cst2[0] = 0.0;
	cst2[1] = 0.0;    
      }
    }
    
    
  }


  vector<double> mean;
  mean.push_back( cst1[0] + cst1[1]*Egen );
  mean.push_back( cst2[0] + cst2[1]*Egen );
  return mean;

}

double MGIntegration::getJetTFfraction( double eta, TString flavor) const {

  double f=0.;

  if(flavor=="light"){

    if(use_pT_TFJet_){
      if(fabs(eta)<1.5)  
	f=0.81;
      else if(fabs(eta)>1.5 && fabs(eta)<3.0)
	f=0.8;
      else
	f=0.63;
    }
    else{
    if(fabs(eta)<1.5)  
      f=0.66;
    else if(fabs(eta)>1.5 && fabs(eta)<3.0)
      f=0.52;
    else
      f=0.5;
    }
    
  }

  
  else if(flavor=="b"){

    if(use_pT_TFJet_){
      if(fabs(eta)<1.5)  
	f=0.72;
      else if(fabs(eta)>1.5 && fabs(eta)<3.0)
	f=0.71;
    }
    else{
      if(fabs(eta)<1.5)  
	f=0.67;
      else if(fabs(eta)>1.5 && fabs(eta)<3.0)
      f=0.6;
    }

  }


  return f;

}

double MGIntegration::getJetTF( const TLorentzVector &quark4P, const TLorentzVector &evJet4P, TString flavor) const {

  double f=getJetTFfraction(quark4P.Eta(),flavor);

  double Egen=quark4P.E();
  double Erec=evJet4P.E();
  if(use_pT_TFJet_){
    Egen=quark4P.Pt();
    Erec=evJet4P.Pt();
  }

  vector<double> res  = getJetTFresolution(Egen, quark4P.Eta(), flavor);
  vector<double> mean = getJetTFmean(Egen, quark4P.Eta(), flavor);  

  double sigmaV0 = res[0]; 
  double sigmaV1 = res[1]; 
  double mu0=mean[0]; 
  double mu1=mean[1];
  double a0 = 1/(sigmaV0*TMath::Sqrt(2*TMath::Pi()) );
  double a1 = 1/(sigmaV1*TMath::Sqrt(2*TMath::Pi()) );

  double TFgaus = f * evalGauss(Erec,a0,mu0,sigmaV0) + (1-f) * evalGauss(Erec,a1,mu1,sigmaV1);

  return TFgaus;  
}

double MGIntegration::getFakeLepTFresolution(double pT_gen) const {

  double cst[3];
  cst[0] = 6.26;
  cst[1] = 0.;
  cst[2] = 0.14;
  
  double res = TMath::Sqrt( pow(cst[0],2) + pow(cst[1]*TMath::Sqrt(pT_gen),2) + pow(cst[2]*pT_gen,2) );
  
  return res;

}

double MGIntegration::getFakeLepTFmean(double pT_gen) const {

  double cst[2];
  cst[0] = -10.74;
  cst[1] = 0.58;

  double mean = cst[0]+cst[1]*pT_gen;

  return mean;

}

double MGIntegration::getFakeLepTF( const TLorentzVector &quark4P, const TLorentzVector &evLep4P) const {

  double pT_gen = quark4P.Pt();
  double pT_reco = evLep4P.Pt();

  double mu = getFakeLepTFmean(pT_gen);
  double sigma = getFakeLepTFresolution(pT_gen);
  double a = 1/(sigma*TMath::Sqrt(2*TMath::Pi()) );

  double TF = evalGauss(pT_reco,a,mu,sigma);

  return TF;

}

double MGIntegration::getFakeTauLepTFresolution(double pT_gen,int leptype) const {

  double cst[2] = {0,0};

  if(abs(leptype)==11){
    cst[0] = 4.69;
    cst[1] = 0.04;
  }
  else if(abs(leptype)==13){
    cst[0] = 4.46;
    cst[1] = 0.02;
  }

  double res = cst[0]+cst[1]*pT_gen;

  return res;

}

double MGIntegration::getFakeTauLepTFmean(double pT_gen,int leptype) const {

  double cst[2] = {0,0};

  if(abs(leptype)==11){
    cst[0] = 0.84;
    cst[1] = 0.96;
  }
  else if(abs(leptype)==13){
    cst[0] = 1.17;
    cst[1] = 0.97;
  }

  double mean = cst[0]+cst[1]*pT_gen;

  return mean;

}

double MGIntegration::getFakeTauLepTFn(double pT_gen,int leptype) const {

  double cst[2] = {0,0};

  if(abs(leptype)==11){
    cst[0] = 0.;
    cst[1] = 0.;
  }
  else if(abs(leptype)==13){
    cst[0] = 2.29;
    cst[1] = -0.0151;
  }

  double n = cst[0]+cst[1]*pT_gen;
  if(n<0.1) n=0.1;

  return n;

}

double MGIntegration::getFakeTauLepTF( const TLorentzVector &lep4P, const TLorentzVector &evTau4P, int leptype) const {

  double pT_gen = lep4P.Pt();
  double pT_reco = evTau4P.Pt();

  double mu = getFakeTauLepTFmean(pT_gen,leptype);
  double sigma = getFakeTauLepTFresolution(pT_gen,leptype);
  double n = getFakeTauLepTFn(pT_gen,leptype);

  double TF = 1.;

  if(abs(leptype)==11){

    double I = 1/(TMath::Pi()*sigma);
    
    TF = I / (1 + pow((pT_reco-mu)/sigma , 2));

  }


  else if (abs(leptype)==13){

    double chi = (pT_reco-mu)/sigma;
    double alpha = 1.2;
    double A = pow(n/alpha,n) * exp(-alpha*alpha/2.);
    double B = n/alpha - alpha;
    //double C = n/alpha * 1./(n-1.) * exp(-alpha*alpha/2.);
    double C = exp(-alpha*alpha/2.)/(n-1.) *(pow(n/alpha,1-n)-pow(B+5.,1-n));
    double D = TMath::Sqrt(TMath::Pi()/2.) * (1 + TMath::Erf(alpha/TMath::Sqrt(2.)));
    double N = 1./(sigma*(C+D));

    if(abs(chi)>5)
      TF = 0;
    else if ( chi < alpha )
      TF = N * exp( - 0.5*chi*chi );
    else
      TF = N * A * pow(B + chi,-n);       

  }

  return TF;

}

double MGIntegration::getMissJetAcceptance( const TLorentzVector &quark4P) const {

  double A=0;
  double dR_lep1=quark4P.DeltaR(evLep1_4P_);
  double dR_lep2=quark4P.DeltaR(evLep2_4P_);
  double dR_tau=quark4P.DeltaR(evHadSys_Tau_4P_);
  double dR_jet1=quark4P.DeltaR(evJet1_4P_);
  double dR_BJet1=quark4P.DeltaR(evBJet1_4P_);
  double dR_BJet2=quark4P.DeltaR(evBJet2_4P_);
  double dR_jet_min=min(dR_jet1,min(dR_BJet1,dR_BJet2));
  
  //Jet out of eta acceptance
  if(abs(quark4P.Eta())>eta_acceptance_)
    A=1;

  //Jet would have made fail the lepton isolation
  else if(dR_lep1<dR_veto_jet_lep_ && (quark4P.E()/evLep1_4P_.E()) > rel_iso_lep_)
    A=0;

  else if(dR_lep2<dR_veto_jet_lep_ && (quark4P.E()/evLep2_4P_.E()) > rel_iso_lep_)
    A=0;

  else if(dR_tau<dR_veto_jet_lep_)
    A=0;

  //Jet would be merged with other jet
  else if(dR_jet_min<jet_radius_)
    A=1;

  //Jet did not pass pT cut
  else{

    double f = getJetTFfraction(quark4P.Eta(),"light");

    double Egen=quark4P.E();
    if(use_pT_TFJet_)
      Egen=quark4P.Pt();

    vector<double> res = getJetTFresolution(Egen,quark4P.Eta(),"light");
    vector<double> mean = getJetTFmean(Egen,quark4P.Eta(),"light");

    double sigmaV0 = res[0]; 
    double sigmaV1 = res[1]; 
    double mu0=mean[0]; 
    double mu1=mean[1];

    double Eup=pT_cut_*cosh(quark4P.Eta());
    if(use_pT_TFJet_)
      Eup=pT_cut_;

    A = 1 / TMath::Sqrt(TMath::Pi()) * 
      ( f * ( TMath::Erf( (Eup-mu0)/(TMath::Sqrt(2)*sigmaV0) ) + TMath::Erf( (mu0)/(TMath::Sqrt(2)*sigmaV0) ) )
	+ (1-f) * ( TMath::Erf( (Eup-mu1)/(TMath::Sqrt(2)*sigmaV1) ) + TMath::Erf( (mu1)/(TMath::Sqrt(2)*sigmaV1) ) ) );

  }
    

  return A;

}

//--- TauLeptonicTF Without angular dependance
double MGIntegration::getTauLeptonicTF( double P_TauLep, const TLorentzVector &Lep4P) const {

  // Energy part
  double E_Lep = Lep4P.E();
  //
  double mTau_2 = Physics::mTau2;
  double E_Tau  = TMath::Sqrt( P_TauLep * P_TauLep + mTau_2);

  double z = E_Lep / E_Tau;
  double TauLepTF = (1-z) * (5 + 5*z - 4 * z*z) / E_Tau ;
  
  return TauLepTF;
}

//--- TauHadronicTF
double MGIntegration::getTauHadronicTF( double P_TauHad, const TLorentzVector &HadSys4P ) const {
  double TauHadTF;
  double mATau = Physics::mTau;
  double E_ATau   = TMath::Sqrt( P_TauHad * P_TauHad + mATau * mATau );
  double z = HadSys4P.E() / E_ATau ; 

  if ( z < ( HadSys4P.M() / mATau ) * ( HadSys4P.M() / mATau ) )
  	TauHadTF = 0;
  else {
    TauHadTF = ( (1. / E_ATau ) * ( 1 / ( 1 - ( HadSys4P.M() / mATau ) * ( HadSys4P.M() / mATau ) ) ) );    
  }

  return TauHadTF;
}

//--- Fake TauHadronic TF
double MGIntegration::getFakeTF( const TLorentzVector &part4P, const TLorentzVector &evHadSys4P ) const {

  double xcur=(evHadSys4P.E()-part4P.E())/part4P.E();
  
  double alpha=2.208933;
  double n=0.647207;
  double mu=-0.331309;
  double sigma=0.179943;

  double A = pow((n/alpha),n) * exp(-alpha*alpha/2);
  double B = n/alpha - alpha;
  double C = n/alpha * 1./(n-1) * exp(-alpha*alpha/2);
  double D = sqrt(PI/2.) * ( 1 + TMath::Erf(alpha/sqrt(2)) );
  double N = 1 / ( sigma * (C+D) );

  double TauHadTF = 1;

  if ((xcur-mu)/sigma < alpha)
    TauHadTF = N*exp(-(xcur-mu)*(xcur-mu)/ (2*sigma*sigma));
  else
    TauHadTF = N*A*pow((B + (xcur-mu)/sigma),-n);

  return TauHadTF;

}

//From Feynman diagram of top decay
//Can be used for leptonic or hadronic top decay
double MGIntegration::getTopTF(const TLorentzVector &t, const TLorentzVector &b, const TLorentzVector &l, const TLorentzVector &nu){

  return (b*nu)*(t*l);

}


///////////////////////////////////////////////////
///////////        Jacobian             ///////////
///////////////////////////////////////////////////

//--- Jacobian term from tau

double MGIntegration::getJacobian_Tau( const TLorentzVector &tau) const {
    
  double jacobian = tau.P()*tau.P()/tau.E();
  return jacobian;

}


//--- Jacobian term from the di-tau system

double MGIntegration::getJacobian_diTau( const TLorentzVector &TauLep4P, const TLorentzVector &TauHad4P,
        double gamma, double cosTheta_TauLepTauHad,  double sinTheta_TauLepPi) const {
  double E_TauLep = TauLep4P.E();
  double P_TauLep = TauLep4P.P();
  double E_TauHad = TauHad4P.E();
  double P_TauHad = TauHad4P.P();
    
  double jacobian = 1. / std::abs( 
                          ( E_TauLep * P_TauHad / E_TauHad  - P_TauLep * cosTheta_TauLepTauHad ) 
                          * sinTheta_TauLepPi * gamma  
                         );
  return jacobian;
}



double MGIntegration::getJacobian_tophad(TLorentzVector b, TLorentzVector W, double Eqbar, double CosTheta_qq){

  double b_Pmag=b.P();
  double Eb=b.E();
  TVector3 eb=(b.Vect()).Unit();
  double Web=W.Vect()*eb;
  double EW=W.E();
  
  double Jac= b_Pmag*Eqbar/(1-CosTheta_qq) * 1 / fabs(Eb/b_Pmag * Web - EW) ;
  return Jac;

}


double MGIntegration::getJacobian_toplep(TLorentzVector b, TLorentzVector W, double El, double Enu){

  double b_Pmag=b.P();
  double Eb=b.E();
  TVector3 eb=(b.Vect()).Unit();
  double Web=W.Vect()*eb;
  double EW=W.E();
  
  double Jac= b_Pmag*Enu*Enu/El * 1 / fabs(Eb/b_Pmag * Web - EW) ;
  return Jac;

}



double MGIntegration::getJacobian_W( double Enu, double Elep) const {

  return Enu*Enu/Elep;

}


///////////////////////////////////////////////////
///////////   Integrand utilities       ///////////
///////////////////////////////////////////////////

void MGIntegration::setEventParameters( const IntegrationMsg_t &data, bool force_missing_jet_integration_forevent) {


  std::cout << "MGIntegration::setEventParameters" << std::endl;
/*  std::cout << "data.evLep1_4P_ : " << data.evLep1_4P_[0] << ", " << data.evLep1_4P_[1] << ", " << data.evLep1_4P_[2] << ", " << data.evLep1_4P_[3] << ", " << std::endl;
  std::cout << "data.evLep1_4P_ : " << data.lepton1_Type_ << std::endl;
  std::cout << "data.evLep2_4P_ : " << data.evLep2_4P_[0] << ", " << data.evLep2_4P_[1] << ", " << data.evLep2_4P_[2] << ", " << data.evLep2_4P_[3] << ", " << std::endl;
  std::cout << "data.evLep2_4P_ : " << data.lepton2_Type_ << std::endl;*/
  setLep_4Ps( data.evLep1_4P_ , data.lepton1_Type_,
	      data.evLep2_4P_ , data.lepton2_Type_);

/*  std::cout << "data.evHadSys_Tau_4P_ : " << data.evHadSys_Tau_4P_[0] << ", " << data.evHadSys_Tau_4P_[1] << ", " << data.evHadSys_Tau_4P_[2] << ", " << data.evHadSys_Tau_4P_[3] << ", " << std::endl;
  std::cout << "data.HadtauDecayMode_ : " << data.HadtauDecayMode_ << std::endl;*/
  setHadronicTau4P( data.evHadSys_Tau_4P_ , data.HadtauDecayMode_);
 
  event_type_ = data.event_type_;

  if(!force_missing_jet_integration_forevent)
    integration_type_ = data.integration_type_;

  setlightJet_4Ps( data.evJet1_4P_ , data.evJet2_4P_);
/*    std::cout << "nbrOfPermut_per_jet_ setlightJet_4Ps : " << nbrOfPermut_per_jet_ << std::endl; // TEMPORAIRE
    std::cout << "nbrOfPermut_         setlightJet_4Ps : " << nbrOfPermut_ << std::endl; // TEMPORAIRE */
  setlightJets_list_4P( data.evJets_4P_, data.n_lightJets_ );
/*    std::cout << "nbrOfPermut_per_jet_ setlightJets_list_4P : " << nbrOfPermut_per_jet_ << std::endl; // TEMPORAIRE
    std::cout << "nbrOfPermut_         setlightJets_list_4P : " << nbrOfPermut_ << std::endl; // TEMPORAIRE //
    std::cout << "*-*-*-" << std::endl; // TEMPORAIRE */

  setBJet_4Ps( data.evBJet1_4P_ , data.evBJet2_4P_);

  setMET_4P( data.evRecoMET4P_);
  setMETCov( data.evV_ );

  for( int i=0; i<nbrOfPermut_per_jet_; i++){
    include_perm_ttbar_DL_[i] = 1;
    include_perm_ttbar_DL_fakelep_tlep_[i] = 1;
    include_perm_ttbar_DL_fakelep_ttau_[i] = 1;
    include_perm_ttW_[i] = 1;
  }
  for( int i=nbrOfPermut_per_jet_; i<nbrOfPermutMax; i++){
    include_perm_ttbar_DL_[i] = 0; 
    include_perm_ttbar_DL_fakelep_tlep_[i] = 0; 
    include_perm_ttbar_DL_fakelep_ttau_[i] = 0; 
    include_perm_ttW_[i] = 0;
  }

  for( int i=0; i<nbrOfPermut_; i++) { 
    include_perm_ttH_[i] = 1;
    include_perm_ttZ_[i] = 1;
    include_perm_ttjets_[i] = 1;
    include_perm_ttbar_SL_[i] = 1;
    include_perm_ttZ_Zll_[i] = 1;
  }

  if(nbrOfPermut_<nbrOfPermutMax){
    for( int i=nbrOfPermut_; i<nbrOfPermutMax; i++) { 
      include_perm_ttH_[i] = 0;
      include_perm_ttZ_[i] = 0;
      include_perm_ttjets_[i] = 0;
      include_perm_ttbar_SL_[i] = 0;
      include_perm_ttZ_Zll_[i] = 0;
    }
  }

  setE_lightQuark_Boundaries();
  setPOutFakeBoundaries();
  
  setmTauTau_CosThetaTauLepTauHadBoundaries();
  setTauHadMomentumBoundaries();

  if(use_top_compatibility_check_){

    if(integration_type_ == integration_type_wo_miss && !force_missing_jet_integration_)
      checkCompatibility_TopHad();
    else if(integration_type_ == integration_type_w_miss || force_missing_jet_integration_)
      checkCompatibility_TopHad_missing_jet();

    checkCompatibility_TopLep();
    checkCompatibility_TopTau();
    checkCompatibility_TopLep_fakelep();
    checkCompatibility_TopTau_fakelep();
  }

  checkCompatibility_Zll();


  // Check if at least one perm included
  if(integration_type_ == integration_type_wo_miss && force_missing_jet_integration_ifnoperm_){

    bool atleast_one_perm_incl = false;
    for( int i=0; i<nbrOfPermut_; i++)
      atleast_one_perm_incl |= include_perm_ttH_[i];
    if(!atleast_one_perm_incl){
      integration_type_ = integration_type_w_miss;
      setEventParameters(data,true);
    }
  
  
  }


  std::cout << "call values initialization" << std::endl;
  // MPI process source
  MPIInfo_  = data.MPIInfo_;
  eventID_ = data.eventID_;
  nRun_   = data.nRun_;
  nLumi_  = data.nLumi_;
  nEvent_ = data.nEvent_;

  for( int i=0; i<nbrOfPermutMax; i++) { 

    integralttH_[i] = NAN;
    integralttZ_[i] = NAN;
    integralttW_[i] = NAN;
    integralttjets_[i] = NAN;
    integralttbar_SL_[i] = NAN;
    integralttbar_DL_[i] = NAN;
    integralttbar_DL_fakelep_tlep_[i] = NAN;
    integralttbar_DL_fakelep_ttau_[i] = NAN;
    integralttZ_Zll_[i] = NAN;


    stderrttH_[i] = NAN;
    stderrttZ_[i] = NAN;
    stderrttW_[i] = NAN;
    stderrttjets_[i] = NAN;
    stderrttbar_SL_[i] = NAN;
    stderrttbar_DL_[i] = NAN;
    stderrttbar_DL_fakelep_tlep_[i] = NAN;
    stderrttbar_DL_fakelep_ttau_[i] = NAN;
    stderrttZ_Zll_[i] = NAN;

    chiSquarettH_[i] = NAN;
    chiSquarettZ_[i] = NAN;
    chiSquarettW_[i] = NAN;
    chiSquarettjets_[i] = NAN;
    chiSquarettbar_SL_[i] = NAN;
    chiSquarettbar_DL_[i] = NAN;
    chiSquarettbar_DL_fakelep_tlep_[i] = NAN;
    chiSquarettbar_DL_fakelep_ttau_[i] = NAN;
    chiSquarettZ_Zll_[i] = NAN;

    integrationEfficiencyttH_[i] = -1;
    integrationEfficiencyttZ_[i] = -1;
    integrationEfficiencyttW_[i] = -1;
    integrationEfficiencyttjets_[i] = -1;
    integrationEfficiencyttbar_SL_[i] = -1;
    integrationEfficiencyttbar_DL_[i] = -1;
    integrationEfficiencyttbar_DL_fakelep_tlep_[i] = -1;
    integrationEfficiencyttbar_DL_fakelep_ttau_[i] = -1;
    integrationEfficiencyttZ_Zll_[i] = -1;

    totalDrawsttH_[i] = -1;
    totalDrawsttZ_[i] = -1;
    totalDrawsttW_[i] = -1;
    totalDrawsttjets_[i] = -1;
    totalDrawsttbar_SL_[i] = -1;
    totalDrawsttbar_DL_[i] = -1;
    totalDrawsttbar_DL_fakelep_tlep_[i] = -1;
    totalDrawsttbar_DL_fakelep_ttau_[i] = -1;
    totalDrawsttZ_Zll_[i] = -1;

    compTimettH_[i] = NAN;
    compTimettZ_[i] = NAN;
    compTimettW_[i] = NAN;
    compTimettjets_[i] = NAN;    
    compTimettbar_SL_[i] = NAN;
    compTimettbar_DL_[i] = NAN;
    compTimettbar_DL_fakelep_tlep_[i] = NAN;
    compTimettbar_DL_fakelep_ttau_[i] = NAN;
    compTimettZ_Zll_[i] = NAN;

  }
  
    std::cout << "nbrOfPermut_         sortie : " << nbrOfPermut_ << std::endl; // TEMPORAIRE
  std::cout << "end of setEventParameters" << std::endl;
  
}

void MGIntegration::copyBoundaries( IntegrationMsg_t *data ) {
  
  data->event_type_ = event_type_;
  data->integration_type_ = integration_type_;
  std::cout << "MGIntegration::copyBoundaries" << std::endl;
  std::cout << "event_type_       : " << event_type_ << std::endl;
  //std::cout << "integration_type_ : " << integration_type_ << std::endl; // TEMPORAIRE

  if(integration_type_ == integration_type_wo_miss){
    data->nbrOfDim_ttH_ = nbrOfDim_ttH_;
    data->nbrOfDim_ttZ_ = nbrOfDim_ttZ_;
    data->nbrOfDim_ttW_ = nbrOfDim_ttW_;
    data->nbrOfDim_ttjets_ = nbrOfDim_ttjets_;
    data->nbrOfDim_ttbar_SL_ = nbrOfDim_ttbar_SL_;
    data->nbrOfDim_ttZ_Zll_ = nbrOfDim_ttZ_Zll_;

    data->nbrOfPoints_ttH_ = nbrOfPoints_ttH_;
    data->nbrOfPoints_ttZ_ = nbrOfPoints_ttZ_;
    data->nbrOfPoints_ttW_ = nbrOfPoints_ttW_;
    data->nbrOfPoints_ttjets_ = nbrOfPoints_ttjets_;
    data->nbrOfPoints_ttbar_SL_ = nbrOfPoints_ttbar_SL_;
    data->nbrOfPoints_ttZ_Zll_ = nbrOfPoints_ttZ_Zll_;

  }  

  else if(integration_type_ == integration_type_w_miss){
    data->nbrOfDim_ttH_miss_ = nbrOfDim_ttH_miss_;
    data->nbrOfDim_ttZ_miss_ = nbrOfDim_ttZ_miss_;
    data->nbrOfDim_ttjets_miss_ = nbrOfDim_ttjets_miss_;
    data->nbrOfDim_ttbar_SL_miss_ = nbrOfDim_ttbar_SL_miss_;
    data->nbrOfDim_ttZ_Zll_miss_ = nbrOfDim_ttZ_Zll_miss_;

    data->nbrOfPoints_ttH_miss_ = nbrOfPoints_ttH_miss_;
    data->nbrOfPoints_ttZ_miss_ = nbrOfPoints_ttZ_miss_;
    data->nbrOfPoints_ttjets_miss_ = nbrOfPoints_ttjets_miss_;
    data->nbrOfPoints_ttbar_SL_miss_ = nbrOfPoints_ttbar_SL_miss_;
    data->nbrOfPoints_ttZ_Zll_miss_ = nbrOfPoints_ttZ_Zll_miss_;

  }    

  data->nbrOfDim_ttbar_DL_ = nbrOfDim_ttbar_DL_;
  data->nbrOfPoints_ttbar_DL_ = nbrOfPoints_ttbar_DL_;
  std::cout << "\tnbrOfDim_ttbar_DL_ : " << nbrOfDim_ttbar_DL_ << std::endl;
  std::cout << "\tnbrOfPoints_ttbar_DL_ : " << nbrOfPoints_ttbar_DL_ << std::endl;

  data->nbrOfPermut_ = nbrOfPermut_;
  std::cout << "\tcopyBoundaries nbrOfPermut_ : " << nbrOfPermut_ << std::endl;

  // Boundaries
  for( int i=0; i< DimensionMax; i++) {
    data->lowerValues_[i] = lowerValues_[i];
    data->upperValues_[i] = upperValues_[i];
  }
  
  for( int i=0; i<nbrOfPermutMax; i++) { 

    data->EQuark1_Lower_[i] = EQuark1_Lower_[i];
    data->EQuark1_Upper_[i] = EQuark1_Upper_[i];

    data->include_perm_ttH_[i] = include_perm_ttH_[i];
    data->mTauTau_ttH_[i] = mTauTau_ttH_[i];
    data->PTauLep_ttH_Lower_[i] = PTauLep_ttH_Lower_[i];
    data->PTauLep_ttH_Upper_[i] = PTauLep_ttH_Upper_[i];
    data->cosTheta_diTau_ttH_Lower_[i] = cosTheta_diTau_ttH_Lower_[i];
    data->cosTheta_diTau_ttH_Upper_[i] = cosTheta_diTau_ttH_Upper_[i];

    data->include_perm_ttZ_[i] = include_perm_ttZ_[i];
    data->mTauTau_ttZ_[i] = mTauTau_ttZ_[i];
    data->PTauLep_ttZ_Lower_[i] = PTauLep_ttZ_Lower_[i];
    data->PTauLep_ttZ_Upper_[i] = PTauLep_ttZ_Upper_[i];
    data->cosTheta_diTau_ttZ_Lower_[i] = cosTheta_diTau_ttZ_Lower_[i];
    data->cosTheta_diTau_ttZ_Upper_[i] = cosTheta_diTau_ttZ_Upper_[i];

    data->include_perm_ttW_[i] = include_perm_ttW_[i];

    data->include_perm_ttjets_[i] = include_perm_ttjets_[i];

    data->include_perm_ttbar_SL_[i] = include_perm_ttbar_SL_[i];

    data->include_perm_ttbar_DL_[i] = include_perm_ttbar_DL_[i];
    data->include_perm_ttbar_DL_fakelep_tlep_[i] = include_perm_ttbar_DL_fakelep_tlep_[i];
    data->include_perm_ttbar_DL_fakelep_ttau_[i] = include_perm_ttbar_DL_fakelep_ttau_[i];

    data->include_perm_ttZ_Zll_[i] = include_perm_ttZ_Zll_[i];


  }

  for( int i=0; i<2; i++) { 
    data->phi_missing_jet_Boundaries_[i] = phi_missing_jet_Boundaries_[i];
    data->cosTheta_missing_jet_Boundaries_[i] = cosTheta_missing_jet_Boundaries_[i];
    data->phiNu_tlep_Boundaries_[i] = phiNu_tlep_Boundaries_[i];
    data->cosThetaNu_tlep_Boundaries_[i] = cosThetaNu_tlep_Boundaries_[i];
    data->POutFakeBoundaries_ttjets_[i] = POutFakeBoundaries_ttjets_[i];

    data->phiNu_ttau_ttbar_DL_ttW_Boundaries_[i] = phiNu_ttau_ttbar_DL_ttW_Boundaries_[i];
    data->cosThetaNu_ttau_ttbar_DL_ttW_Boundaries_[i] = cosThetaNu_ttau_ttbar_DL_ttW_Boundaries_[i];
    data->PTauHad_ttbar_DL_ttW_Boundaries_[i] = PTauHad_ttbar_DL_ttW_Boundaries_[i];
    data->cosThetaNu_W_ttW_Boundaries_[i] = cosThetaNu_W_ttW_Boundaries_[i];
    data->phiNu_W_ttW_Boundaries_[i] = phiNu_W_ttW_Boundaries_[i];

  }

    
}

double MGIntegration::evalttH(const double* x ) {


  double mb=Physics::mb;

 
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@                   Take the input variable form VEGAS                 @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double  P_TauLep	  = x[PTauLep_id];        //module of the momentum of the tau
  double  cosThetaTauLepTauHad = x[cosThetaTauLepTauHad_id]; //cos of the TauLep/TauHad angle
  double  EQuark1	  = x[EQuark1_id];     //module of the momentum of the final state Quark1
  double  cosThetaNu = x[cosThetaNu_tlep_id]; // cosTheta of the nu from the leptonic top
  double  phiNu = x[phiNu_tlep_id]; // phi of the nu from the leptonic top  

  double cosTheta_missing_jet = 0;
  double phi_missing_jet      = 0;
  if(integration_type_ == integration_type_w_miss){
    cosTheta_missing_jet = x[cosTheta_missing_jet_id];  // cosTheta of the quark associated to the missing jet
    phi_missing_jet      = x[phi_missing_jet_id];   // phi of the quark associated to the missing jet
  }

  const char* error = 0;
  // To remove ???
//  char str[256];
/*    std::cout << "Vegas point :" 
	       << " " << " cosTheta_missing_jet = "    << cosTheta_missing_jet
	       << " " << " phi_missing_jet = "    << phi_missing_jet
	       << " " << " P_TauLep = "       << P_TauLep
	       << " " << " cosThetaTLepTHad = " << cosThetaTauLepTauHad
	       << " " << " EQuark1 = "    << EQuark1
	       << " " << " cosThetaNu = "    << cosThetaNu
	       << " " << " phiNu = "    << phiNu
	       << " " << endl;*/

  if (verbose_ >= IntegrandLevel){
    (*stream_) << "Vegas point :" 
	       << " " << " cosTheta_missing_jet = "    << cosTheta_missing_jet
	       << " " << " phi_missing_jet = "    << phi_missing_jet
	       << " " << " P_TauLep = "       << P_TauLep
	       << " " << " cosThetaTLepTHad = " << cosThetaTauLepTauHad
	       << " " << " EQuark1 = "    << EQuark1
	       << " " << " cosThetaNu = "    << cosThetaNu
	       << " " << " phiNu = "    << phiNu
	       << " " << endl;
//    std::cout << " IntegrandLevel <= verbose_" << std::endl;
  }
  else {
//      std::cout << " IntegrandLevel > verbose_" << std::endl;
  }


  //get the 3D momentum vector
  TVector3 Lep_Tau_3P    = evLep_Tau_4P_.Vect();
  TVector3 HadSys_Tau_3P = evHadSys_Tau_4P_.Vect();  
  
  //Normalization of the 3D momentum vector:
  TVector3  e_lep = Lep_Tau_3P.Unit();
  TVector3  e_pi  = HadSys_Tau_3P.Unit();
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the Leptonic Tau 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TVector3 taulep = P_TauLep * e_lep;
  
  double E_TauLep = TMath::Sqrt( P_TauLep*P_TauLep + Physics::mTau2 );
  TLorentzVector TauLep4P(taulep, E_TauLep);
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the Hadronic Tau 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
  // Frame Initialization         
  // 2nd ref. sys. (e_tau, e_pi, e_yy)
  TVector3 e_tau;
  TVector3 e_yy;
  //Normalization of the tau 3D momentum
  e_tau  = taulep.Unit();
  e_yy   = e_tau.Cross(e_pi);
  e_yy   = e_yy.Unit();
  
  // Evaluation of the Hadronic Tau 4-Vector 
  double P_TauHad = getPTauHad( m_TauTau_2_, P_TauLep, cosThetaTauLepTauHad );
  double cosThetaTauLepPi = TMath::Cos( taulep.Angle(e_pi) );
  double sinThetaTauLepPi = TMath::Sqrt( 1.0 - cosThetaTauLepPi*cosThetaTauLepPi);
  double cosThetaTauHadPi = getCosThetaTauHadPi( P_TauHad, evHadSys_Tau_4P_ );

  //
  double alpha, beta, gamma, gamma2;
  
  getAlphaBetaGamma( cosThetaTauLepTauHad, 
        cosThetaTauLepPi, cosThetaTauHadPi,
        alpha, beta, gamma, gamma2, &error);
 
  //
  TVector3 tauhad_plus = P_TauHad * (alpha*e_tau + beta*e_pi + gamma*e_yy);
  TVector3 tauhad_minus = P_TauHad * (alpha*e_tau + beta*e_pi - gamma*e_yy);
  //  
  double E_TauHad = TMath::Sqrt( P_TauHad*P_TauHad + Physics::mTau2 );
  //
  TLorentzVector TauHad_plus_4P(tauhad_plus, E_TauHad); 
  TLorentzVector TauHad_minus_4P(tauhad_minus, E_TauHad); 
  if (verbose_ >= IntegrandLevel){
    (*stream_) << "Building Hadronic Tau :" << endl;
    (*stream_) << " .... alpha, beta, gamma " 
                << " " << alpha
                << " " << beta
                << " " << gamma
               << endl;
  }
  if (verbose_ >= IntegrandLevel){
    (*stream_) << "Build Leptonic/Hadronic Tau :" << endl;
    (*stream_) << " .... leptonic Tau " 
            <<  TauLep4P.Px() << ", "
            <<  TauLep4P.Py() << ", "
            <<  TauLep4P.Pz() << ", "
            <<  TauLep4P.E()     
            << endl;
   (*stream_) << " .... hadronic Tau plus" 
            <<  TauHad_plus_4P.Px() << ", "
            <<  TauHad_plus_4P.Py() << ", "
            <<  TauHad_plus_4P.Pz() << ", "
            <<  TauHad_plus_4P.E()    
            << endl;
   (*stream_) << " .... hadronic Tau minus" 
            <<  TauHad_minus_4P.Px() << ", "
            <<  TauHad_minus_4P.Py() << ", "
            <<  TauHad_minus_4P.Pz() << ", "
            <<  TauHad_minus_4P.E()    
            << endl;
  }
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                          Constrain on Leptonic Tau                   @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if ( P_TauLep >  u_plus(m_TauTau_2_, cosThetaTauLepTauHad) ) {
    if (error == 0) error = "outLepTau";
  }  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                          Constrain on Hadronic Tau                   @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  if ( (P_TauHad <  evSMinus_)  || (P_TauHad >  evSPlus_) ) {
    if (error == 0) error = "outHadTau";
  }  




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the leptonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_lep_reco_4P;
  TLorentzVector W_lep_4P;
  TLorentzVector nu_top_4P;
  TLorentzVector t_lep_4P;

  if(error==0){
    
    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu),phiNu);
    double CosThetalnu_top=TMath::Cos(evLep_top_4P_.Angle(enu));
    double Enu=getEqbar_Enu(CosThetalnu_top,evLep_top_4P_.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_top_4P=nu_top_4P_dummy;
    W_lep_4P=nu_top_4P+evLep_top_4P_;

    
    double Eb_lep=getEb(W_lep_4P, evBJet_leptop_4P_);
    double b_lep_Pmag=TMath::Sqrt(Eb_lep*Eb_lep - mb*mb);
    double pT_bl=b_lep_Pmag/TMath::CosH(evBJet_leptop_4P_.Eta());
    if(Eb_lep>0){
      b_lep_reco_4P.SetPtEtaPhiE(pT_bl,evBJet_leptop_4P_.Eta(),evBJet_leptop_4P_.Phi(),Eb_lep);
      t_lep_4P=W_lep_4P+b_lep_reco_4P;
    }
    else{
      error = "leptop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build leptonic top :" << endl;
      (*stream_) << " .... b from lep top " 
              <<  b_lep_reco_4P.Px() << ", "
              <<  b_lep_reco_4P.Py() << ", "
              <<  b_lep_reco_4P.Pz() << ", "
              <<  b_lep_reco_4P.E()     
              << endl;
     (*stream_) << " .... lep top" 
              <<  t_lep_4P.Px() << ", "
              <<  t_lep_4P.Py() << ", "
              <<  t_lep_4P.Pz() << ", "
              <<  t_lep_4P.E()    
              << endl;
    }

  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the hadronic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector q1_reco_4P;
  TLorentzVector q2_reco_4P;
  TLorentzVector b_had_reco_4P;
  TLorentzVector W_had_4P;
  TLorentzVector t_had_4P;

  if(error==0){

    double pt_q1=EQuark1/TMath::CosH(evJet1_4P_.Eta());
    q1_reco_4P.SetPtEtaPhiE(pt_q1,evJet1_4P_.Eta(),evJet1_4P_.Phi(),EQuark1);
    
//    std::cout << "integration_type_ : " << integration_type_ << std::endl;// TEMPORAIRE
    
    if(integration_type_ == integration_type_wo_miss){
      double CosTheta_qq=TMath::Cos(evJet1_4P_.Angle(evJet2_4P_.Vect()));
      double EQuark2=getEqbar_Enu(CosTheta_qq,EQuark1);
      double pt_q2=EQuark2/TMath::CosH(evJet2_4P_.Eta());
      q2_reco_4P.SetPtEtaPhiE(pt_q2,evJet2_4P_.Eta(),evJet2_4P_.Phi(),EQuark2);
    }
    else if(integration_type_ == integration_type_w_miss){
      TVector3 e_q2;
      e_q2.SetMagThetaPhi(1,TMath::ACos(cosTheta_missing_jet),phi_missing_jet);
      double CosTheta_qq=TMath::Cos(evJet1_4P_.Angle(e_q2));
      double EQuark2=getEqbar_Enu(CosTheta_qq,EQuark1);
      double pt_q2=EQuark2/TMath::CosH(e_q2.Eta());
      q2_reco_4P.SetPtEtaPhiE(pt_q2,e_q2.Eta(),e_q2.Phi(),EQuark2);
    }

    W_had_4P=q1_reco_4P+q2_reco_4P;
    /*std::cout << "q1_reco_4P : "  ; //TEMPORAIRE
    q1_reco_4P.Print();
    std::cout << "q2_reco_4P : "  ; //TEMPORAIRE
    q2_reco_4P.Print();
    std::cout << "W_had_4P : "  ; //TEMPORAIRE
    W_had_4P.Print();*/
    
    double Eb_had=getEb(W_had_4P, evBJet_hadtop_4P_);
    double b_had_Pmag=TMath::Sqrt(Eb_had*Eb_had - mb*mb);
    double pT_bh=b_had_Pmag/TMath::CosH(evBJet_hadtop_4P_.Eta());
    if(Eb_had>0){
      b_had_reco_4P.SetPtEtaPhiE(pT_bh,evBJet_hadtop_4P_.Eta(),evBJet_hadtop_4P_.Phi(),Eb_had);
      t_had_4P=W_had_4P+b_had_reco_4P;
    }
    else{
      error = "hadtop";
    }
    /*std::cout << "pT_bh : " << pT_bh << std::endl; //TEMPORAIRE
    std::cout << "evBJet_hadtop_4P_.Eta() : " << evBJet_hadtop_4P_.Eta() << std::endl; //TEMPORAIRE
    std::cout << "evBJet_hadtop_4P_.Phi() : " << evBJet_hadtop_4P_.Phi() << std::endl; //TEMPORAIRE
    std::cout << "Eb_had : " << Eb_had << std::endl; //TEMPORAIRE */
    
    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build hadronic top :" << endl;
      (*stream_) << " .... b from had top " 
              <<  b_had_reco_4P.Px() << ", "
              <<  b_had_reco_4P.Py() << ", "
              <<  b_had_reco_4P.Pz() << ", "
              <<  b_had_reco_4P.E()     
              << endl;
      (*stream_) << " .... q1 from had top " 
              <<  q1_reco_4P.Px() << ", "
              <<  q1_reco_4P.Py() << ", "
              <<  q1_reco_4P.Pz() << ", "
              <<  q1_reco_4P.E()     
              << endl;
      (*stream_) << " .... q2 from had top " 
              <<  q2_reco_4P.Px() << ", "
              <<  q2_reco_4P.Py() << ", "
              <<  q2_reco_4P.Pz() << ", "
              <<  q2_reco_4P.E()     
              << endl;
      (*stream_) << " .... had top " 
		 <<  t_had_4P.Px() << ", "
		 <<  t_had_4P.Py() << ", "
		 <<  t_had_4P.Pz() << ", "
		 <<  t_had_4P.E()    
		 << endl;
    }


  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                Recoil                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector tot_plus_4P;
  TLorentzVector tot_minus_4P;
  TLorentzVector tot4P;
  TLorentzVector rho4P;
  TLorentzVector TauHad4P;

  if(error==0){
    
    tot_plus_4P = TauLep4P + TauHad_plus_4P + t_lep_4P + t_had_4P ;
    tot_minus_4P = TauLep4P + TauHad_minus_4P + t_lep_4P + t_had_4P ;
    tot4P = (tot_plus_4P.Pt() < tot_minus_4P.Pt()) ? tot_plus_4P : tot_minus_4P;
    TauHad4P = (tot_plus_4P.Pt() < tot_minus_4P.Pt()) ? TauHad_plus_4P : TauHad_minus_4P;
    if(integration_type_ == integration_type_wo_miss)
      rho4P = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_Tau_4P_ + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_ + evJet1_4P_ + evJet2_4P_ ;
    else if(integration_type_ == integration_type_w_miss)
      rho4P = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_Tau_4P_ + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_ + evJet1_4P_ ;


  }

  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        JET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

  double T_Jet1 = 1.0;
  double T_Jet2 = 1.0;
  double T_BJet_leptop = 1.0;
  double T_BJet_hadtop = 1.0;

  if ( error == 0 ) {

    if (flagTFJet1_) {

      if(integration_type_ == integration_type_wo_miss)
	T_Jet1 = getJetTF( q1_reco_4P, evJet1_4P_, "light");

      
      else if(integration_type_ == integration_type_w_miss){
	if(q1_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector q1q2_reco_4P = q1_reco_4P;
	  q1q2_reco_4P.SetE(q1_reco_4P.E() + q2_reco_4P.E());
	  T_Jet1 = getJetTF( q1q2_reco_4P, evJet1_4P_, "light");
	}
	else
	  T_Jet1 = getJetTF( q1_reco_4P, evJet1_4P_, "light");
      }
    }

    //
    if (flagTFJet2_) {
      if(integration_type_ == integration_type_wo_miss)
	T_Jet2 = getJetTF( q2_reco_4P, evJet2_4P_,"light");	
      else if(integration_type_ == integration_type_w_miss)
	T_Jet2 = getMissJetAcceptance( q2_reco_4P);
    } 
    //
    if (flagTFBJet_leptop_) {
      if(integration_type_ == integration_type_wo_miss)
	T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");
      
      else if(integration_type_ == integration_type_w_miss){
	if(b_lep_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector b_lep_q2_reco_4P = b_lep_reco_4P;
	  b_lep_q2_reco_4P.SetE(b_lep_reco_4P.E() + q2_reco_4P.E());
	  T_BJet_leptop = getJetTF( b_lep_q2_reco_4P, evBJet_leptop_4P_,"b");
	}
	else
	  T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");
      }
    }
    //
    if (flagTFBJet_hadtop_) {
        if(integration_type_ == integration_type_wo_miss)
            T_BJet_hadtop = getJetTF( b_had_reco_4P, evBJet_hadtop_4P_,"b");
     
        else if(integration_type_ == integration_type_w_miss){
            if(b_had_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
                // Two quarks are merged in a single jet: add two energies for TF
                TLorentzVector b_had_q2_reco_4P = b_had_reco_4P;
                b_had_q2_reco_4P.SetE(b_had_reco_4P.E() + q2_reco_4P.E());
                T_BJet_hadtop = getJetTF( b_had_q2_reco_4P, evBJet_leptop_4P_,"b");
            }
            else
                T_BJet_hadtop = getJetTF( b_had_reco_4P, evBJet_hadtop_4P_,"b");
        }
    } 
    /*std::cout << "b_had_reco_4P : " << std::endl; //TEMPORAIRE
    b_had_reco_4P.Print();                      //TEMPORAIRE
    std::cout << "evBJet_hadtop_4P_ : " << std::endl;  //TEMPORAIRE
    evBJet_hadtop_4P_.Print();  */                      //TEMPORAIRE


    // Log
/*    if (outgoingJetsOut_) {
      (*outgoingJetsOut_) << setiosflags(ios::fixed) << scientific;
      (*outgoingJetsOut_) 
         << " " << setw(21) << evJet1_4P_.Pt()
         << " " << setw(21) << evJet1_4P_.Eta()
         << " " << setw(21) << evJet1_4P_.Phi()
         << " " << setw(21) << evJet1_4P_.P()
         << " " << setw(21) << evJet1_4P_.E()      
         << " " << setw(21) << q1_reco_4P.Pt()
         << " " << setw(21) << q1_reco_4P.Eta()
         << " " << setw(21) << q1_reco_4P.Phi()
         << " " << setw(21) << q1_reco_4P.P()
         << " " << setw(21) << q1_reco_4P.E()
         << " " << setw(21) << T_Jet1
         << " " << setw(21) << evJet2_4P_.Pt()
         << " " << setw(21) << evJet2_4P_.Eta()
         << " " << setw(21) << evJet2_4P_.Phi()
         << " " << setw(21) << evJet2_4P_.P()
         << " " << setw(21) << evJet2_4P_.E()      
         << " " << setw(21) << q2_reco_4P.Pt()
         << " " << setw(21) << q2_reco_4P.Eta()
         << " " << setw(21) << q2_reco_4P.Phi()
         << " " << setw(21) << q2_reco_4P.P()
         << " " << setw(21) << q2_reco_4P.E()
         << " " << setw(21) << T_Jet2
	 << " " << setw(21) << evBJet_leptop_4P_.Pt()
         << " " << setw(21) << evBJet_leptop_4P_.Eta()
         << " " << setw(21) << evBJet_leptop_4P_.Phi()
         << " " << setw(21) << evBJet_leptop_4P_.P()
         << " " << setw(21) << evBJet_leptop_4P_.E()      
         << " " << setw(21) << b_lep_reco_4P.Pt()
         << " " << setw(21) << b_lep_reco_4P.Eta()
         << " " << setw(21) << b_lep_reco_4P.Phi()
         << " " << setw(21) << b_lep_reco_4P.P()
         << " " << setw(21) << b_lep_reco_4P.E()
         << " " << setw(21) << T_BJet_leptop
	 << " " << setw(21) << evBJet_hadtop_4P_.Pt()
         << " " << setw(21) << evBJet_hadtop_4P_.Eta()
         << " " << setw(21) << evBJet_hadtop_4P_.Phi()
         << " " << setw(21) << evBJet_hadtop_4P_.P()
         << " " << setw(21) << evBJet_hadtop_4P_.E()      
         << " " << setw(21) << b_had_reco_4P.Pt()
         << " " << setw(21) << b_had_reco_4P.Eta()
         << " " << setw(21) << b_had_reco_4P.Phi()
         << " " << setw(21) << b_had_reco_4P.P()
         << " " << setw(21) << b_had_reco_4P.E()
         << " " << setw(21) << T_BJet_hadtop
         << endl;
    } */
  } 


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        MET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_MET=1.0;
  if ( error == 0 ) {
    if (flagTFMET_) {

      TVector3 RhoT;
      TVector3 PT;
      
      if(include_hadrecoil_){
	
	RhoT.SetXYZ((rho4P).Px(), (rho4P).Py(), 0);
	PT.SetXYZ((tot4P).Px(), (tot4P).Py(), 0);

      }

      else{
	
	TLorentzVector tot4P_nohad =  TauLep4P + TauHad4P + W_lep_4P ;
	TLorentzVector rho4P_nohad = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_Tau_4P_ + evLep_top_4P_ ;
	RhoT.SetXYZ((rho4P_nohad).Px(), (rho4P_nohad).Py(), 0);
	PT.SetXYZ((tot4P_nohad).Px(), (tot4P_nohad).Py(), 0);


      }

      /*std::cout << "evV_[0] : " << evV_[0] << " - "; //TEMPORAIRE
      std::cout << "evV_[1] : " << evV_[1] << " - "; //TEMPORAIRE
      std::cout << "evV_[2] : " << evV_[2] << " - "; //TEMPORAIRE
      std::cout << "evV_[3] : " << evV_[3] << std::endl; //TEMPORAIRE */
      
      TVector3 Rec = RhoT-PT;
      // Order : (0,0), (0,1), (1,0), (1,1) ???
      double sqrtDetV = TMath::Sqrt( std::abs( evV_[0] * evV_[3] - evV_[1] * evV_[2]));
      double f = Rec.X()*(evV_[0]*Rec.X() + evV_[1]*Rec.Y()) 
               + Rec.Y()*(evV_[2]*Rec.X() + evV_[3]*Rec.Y());
      T_MET =  TMath::Exp( -0.5*f ) * sqrtDetV / ( 2 * TMath::Pi() ) ;
    
    /*std::cout << "Rec : " ;  //TEMPORAIRE
    Rec.Print();                        //TEMPORAIRE
    std::cout << "RhoT : " ;  //TEMPORAIRE
    RhoT.Print();                        //TEMPORAIRE
    std::cout << "PT : " ;  //TEMPORAIRE
    PT.Print();                        //TEMPORAIRE
    std::cout << "f : " << f ;  //TEMPORAIRE
    std::cout << " - T_MET : " << T_MET << std::endl;  //TEMPORAIRE */
    

    }
  } 

 
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                   LEPTONIC TAU Transfer function                     @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_lepTau = 1.0;
  if ( error == 0 ) {
    if ( flagTFLepTau_ )
      T_lepTau = getTauLeptonicTF(P_TauLep, evLep_Tau_4P_);
  }
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                    HADRONIC TAU Transfer function                    @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_hadTau = 1.0;
  if ( error == 0 ) {
    if ( flagTFHadTau_ )
      T_hadTau = getTauHadronicTF(P_TauHad, evHadSys_Tau_4P_);  
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                       TOP Transfer functions                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_lepTop = 1.0;
  double T_hadTop = 1.0;

  if ( error == 0 ) {
    if ( flagTFTop_ ){
      T_lepTop = getTopTF(t_lep_4P,b_lep_reco_4P,evLep_top_4P_,nu_top_4P);

      //We don't know which q is the quark or the antiquark
      T_hadTop = getTopTF(t_had_4P,b_had_reco_4P,q1_reco_4P,q2_reco_4P)
	+ getTopTF(t_had_4P,b_had_reco_4P,q2_reco_4P,q1_reco_4P);

    }
  }

 

  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Jacobian Term                             @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
 
  double Jac = 1.0;
  
  if ( error == 0 ) {  
    if( flagJac_ ) {

      double Jac_taulep = getJacobian_Tau( TauLep4P );
      double Jac_tauhad = getJacobian_Tau( TauHad4P );
      double Jac_diTau = getJacobian_diTau( TauLep4P, TauHad4P, gamma, cosThetaTauLepTauHad, sinThetaTauLepPi );
      double cosTheta_qq=TMath::Cos(q1_reco_4P.Angle(q2_reco_4P.Vect()));
      double Jac_tophad = getJacobian_tophad( b_had_reco_4P, W_had_4P, q2_reco_4P.E(), cosTheta_qq);
      double Jac_toplep = getJacobian_toplep( b_lep_reco_4P, W_lep_4P, evLep_top_4P_.E(), nu_top_4P.E());
      Jac = Jac_taulep * Jac_tauhad * Jac_diTau * Jac_tophad * Jac_toplep;

      if( TMath::IsNaN( Jac ) ) { if (error == 0) error = "Jacobian";}  
    }
  }




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                 BOOST                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double totPt = 0.0;
  TLorentzVector tot4P_boosted;
  TLorentzVector TauLep4P_boosted;
  TLorentzVector TauHad4P_boosted;
  TLorentzVector t_had_4P_boosted;
  TLorentzVector t_lep_4P_boosted;

  if ( error == 0 ) {
    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Boosting :" << endl;
      (*stream_) << "  initial total momentum tot4P " 
              <<  tot4P.Px() << ", "
              <<  tot4P.Py() << ", "
              <<  tot4P.Pz() << ", "
              <<  tot4P.E() 
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
    totPt = tot4P.Pt();
    
    // Set transverse boost
    TVector3 boost3D = tot4P.BoostVector();
    boost3D.SetZ(0);
    
    //
    // Application of the boost to the final state 4-momenta
    // WARRNING : only following 4-momenta are now booosted
    tot4P_boosted=tot4P;
    tot4P_boosted.Boost   ( -boost3D );

    // ??? implulsion > energy 
    if( TMath::IsNaN(tot4P_boosted.E()) ) { if (error == 0) error = "Boost";}  

    TauLep4P_boosted=TauLep4P;
    TauHad4P_boosted=TauHad4P;
    t_had_4P_boosted=t_had_4P;
    t_lep_4P_boosted=t_lep_4P;

    TauLep4P_boosted.Boost( -boost3D );
    TauHad4P_boosted.Boost( -boost3D );
    t_had_4P_boosted.Boost( -boost3D );
    t_lep_4P_boosted.Boost( -boost3D );

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "  Boosted total momentum tot4P " 
              <<  tot4P_boosted.Px() << ", "
              <<  tot4P_boosted.Py() << ", "
              <<  tot4P_boosted.Pz() << ", "
              <<  tot4P_boosted.E() 
              << endl;
      (*stream_) << "  Boost magnitude " 
              <<  totPt
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        THE BJORKEN FRACTION                          @@@
  //@@@                         Matrix Elements                              @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double wME = 1.0;
  double x_a;
  double x_b;
  x_a = (tot4P_boosted.E() + tot4P_boosted.Pz()) / sqrtS_ ;
  x_b = (tot4P_boosted.E() - tot4P_boosted.Pz()) / sqrtS_ ;
  
  if ( (x_a > 1.0) || (x_a < 0.0) || (x_b > 1.0 ) ||  (x_b < 0.0 ) ) 
    error = "BFrac";
  
  if ( error == 0 ) {

    TLorentzVector top_4P;
    TLorentzVector Atop_4P;
    TLorentzVector Tau_4P;
    TLorentzVector ATau_4P;

    if(lepton_Tau_charge_>0){
      ATau_4P=TauLep4P_boosted;
      Tau_4P=TauHad4P_boosted;
    }
    else{
      Tau_4P=TauLep4P_boosted;
      ATau_4P=TauHad4P_boosted;
    }

    if(lepton_top_charge_>0){
      top_4P=t_lep_4P_boosted;
      Atop_4P=t_had_4P_boosted;
    }
    else{
      top_4P=t_had_4P_boosted;
      Atop_4P=t_lep_4P_boosted;
    }


    // Parton momentum fraction
    //
    // MadGraph quadri-vector construction
    double inPart_a[4], inPart_b[4]; 
    double outTau[4], outATau[4];
    double outtop[4], outAtop[4];
    double inE;
    // Incoming quarks
    inE =  0.5 * x_a * sqrtS_;
    inPart_a[0] = inE; inPart_a[1] = 0.; inPart_a[2] =  0; inPart_a[3] = inE; 
    inE =  0.5 * x_b * sqrtS_;
    inPart_b[0] = inE; inPart_b[1] = 0.; inPart_b[2] =  0; inPart_b[3] =-inE; 
    // Tau / ATau
    outTau[0] = Tau_4P.E(); outTau[1] = Tau_4P.Px(); outTau[2] = Tau_4P.Py(); outTau[3] = Tau_4P.Pz();
    outATau[0] = ATau_4P.E(); outATau[1] = ATau_4P.Px(); outATau[2] = ATau_4P.Py(); outATau[3] = ATau_4P.Pz();
    // top / Atop
    outtop[0] = top_4P.E(); outtop[1] = top_4P.Px(); outtop[2] = top_4P.Py(); outtop[3] = top_4P.Pz();
    outAtop[0] = Atop_4P.E(); outAtop[1] = Atop_4P.Px(); outAtop[2] = Atop_4P.Py(); outAtop[3] = Atop_4P.Pz();

    // Quadri-vector order for MadGraph
    //   p[0] incoming parton0
    //   p[1] incoming parton1
    //   p[2] outgoing top
    //   p[3] outgoing Atop
    //   p[4] outgoing tau-
    //   p[5] outgoing tau+
    vector<double*> p(6,(double*) 0);
    p[0] = inPart_a ; p[1] = inPart_b; 
    p[2] = outtop; p[3] = outAtop; 
    p[4] = outTau ; p[5] = outATau   ; 

    // Matrix Element
    if (flagWME_) {
      if (signalME_){
	//if(MEversion_==1)
	wME = get_ttHWeightedME( x_a, x_b, p);
      }
      else{
	//if(MEversion_==1)
	//wME = get_ttZWeightedME ( x_a, x_b, p );
	wME = get_ttZZonlyWeightedME ( x_a, x_b, p );
      }
      if( TMath::IsNaN(wME) ) { if (error == 0) error = "WeightME";}
    }
  }
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Weight evaluation                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    double Eval =  T_Jet1 * T_Jet2 * T_BJet_hadtop * T_BJet_leptop * T_MET * T_lepTau * T_hadTau * T_lepTop* T_hadTop * Jac * wME;
/*        std::cout << "Integrand value :" ;
        std::cout<< " T_Jet1 = " << T_Jet1 ;
        std::cout<< " T_Jet2 = " << T_Jet2 ;
        std::cout<< " T_BJet_hadtop = " << T_BJet_hadtop ;
        std::cout<< " T_BJet_leptop = " << T_BJet_leptop ;
        std::cout<< " T_MET = " << T_MET ;
        std::cout<< " T_lepTau = " << T_lepTau ;
        std::cout<< " T_hadTau = " << T_hadTau ;
        std::cout<< " T_lepTop = " << T_lepTop ;
        std::cout<< " T_hadTop = " << T_hadTop ;
        std::cout<< " Jac = " << Jac ;
        std::cout<< " wME = " << wME  ;
        std::cout<< " Eval = " << Eval << std::endl;*/
        //std::cout<< " error = " << error << std::endl;/**/
    if (error) {
        Eval = 0;
    }

  if (verbose_ >= IntegrandLevel){
    // if( Eval > 10.0e-40) {
    (*stream_) << "Integrand value :" 
       << " T_Jet1 = " << T_Jet1
       << " T_Jet2 = " << T_Jet2
       << " T_BJet_hadtop = " << T_BJet_hadtop
       << " T_BJet_leptop = " << T_BJet_leptop
       << " T_MET = " << T_MET
       << " T_lepTau = " << T_lepTau
       << " T_hadTau = " << T_hadTau
       << " T_lepTop = " << T_lepTop
       << " T_hadTop = " << T_hadTop
       << " Jac = " << Jac
       << " wME = " << wME 
       << " Eval = "    << Eval
       << " Error= " << error
       << endl;
    // }
  }
/*  const char *err_str = "Ok";
  if (error) 
    err_str = error;
  else{*/
  if ( !error) {
    if (signalME_){
      integr_EfficiencyttH_++;
    }
    else{ 
      integr_EfficiencyttZ_++;
    }
  }

  if (signalME_) 
    tot_DrawsttH_++;
  else 
    tot_DrawsttZ_++;

//  std::cout << "Eval : " << Eval << std::endl;
  return Eval;
}

double MGIntegration::evalttW(const double* x ) {


  double mb=Physics::mb;

 
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@                   Take the input variable form VEGAS                 @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double  cosThetaNu_tlep = x[cosThetaNu_tlep_ttbar_DL_ttW_id]; // cosTheta of the nu from the leptonic top
  double  phiNu_tlep = x[phiNu_tlep_ttbar_DL_ttW_id];           // phi of the nu from the leptonic top
  double  cosThetaNu_ttau = x[cosThetaNu_ttau_ttbar_DL_ttW_id]; // cosTheta of the nu from the tauonic top
  double  phiNu_ttau = x[phiNu_ttau_ttbar_DL_ttW_id];           // phi of the nu from the tauonic top
  double P_TauHad = x[PTauHad_ttbar_DL_ttW_id];
  double  cosThetaNu_W	  = x[cosThetaNu_W_ttW_id];    //cosTheta of the nu from W
  double phiNu_W          = x[phiNu_W_ttW_id];         //phi of the nu from W  



  const char* error = 0;
  // To remove ???
//  char str[256];
  

 if (verbose_ >= IntegrandLevel){
    (*stream_) << "Vegas point :" 
	       << " " << " cosThetaNu_tlep = "    << cosThetaNu_tlep
	       << " " << " phiNu_tlep = "    << phiNu_tlep
	       << " " << " cosThetaNu_ttau = "    << cosThetaNu_ttau
	       << " " << " phiNu_ttau = "    << phiNu_ttau
	       << " " << " P_TauHad = "    << P_TauHad
	       << " " << " cosThetaNu_W = "       << cosThetaNu_W
	       << " " << " phiNu_W = " << phiNu_W
	       << " " << endl;
  }


  //get the 3D momentum vector
  TVector3 Lep_W_3P    = evLep_Tau_4P_.Vect();
  TVector3 e_lep = Lep_W_3P.Unit();


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                    Evaluation of the W 4-Vector                      @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  TVector3 enu_W;
  enu_W.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_W),phiNu_W);
  double CosThetalnu_W=TMath::Cos(evLep_Tau_4P_.Angle(enu_W));
  double Enu_W=getEqbar_Enu(CosThetalnu_W,evLep_Tau_4P_.E());
  TLorentzVector nu_W_4P(Enu_W*enu_W,Enu_W);
  TLorentzVector W_4P=nu_W_4P+evLep_Tau_4P_;

  if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build W :" << endl;
      (*stream_) << " .... nu " 
              <<  nu_W_4P.Px() << ", "
              <<  nu_W_4P.Py() << ", "
              <<  nu_W_4P.Pz() << ", "
              <<  nu_W_4P.E()     
              << endl;
     (*stream_) << " .... W" 
              <<  W_4P.Px() << ", "
              <<  W_4P.Py() << ", "
              <<  W_4P.Pz() << ", "
              <<  W_4P.E()    
              << endl;
  }




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the leptonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_lep_reco_4P;
  TLorentzVector W_lep_4P;
  TLorentzVector nu_lep_4P;
  TLorentzVector t_lep_4P;

  if(error==0){
    
    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_tlep),phiNu_tlep);
    double CosThetalnu_top=TMath::Cos(evLep_top_4P_.Angle(enu));
    double Enu=getEqbar_Enu(CosThetalnu_top,evLep_top_4P_.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_lep_4P=nu_top_4P_dummy;
    W_lep_4P=nu_lep_4P+evLep_top_4P_;
    
    double Eb_lep=getEb(W_lep_4P, evBJet_leptop_4P_);
    double b_lep_Pmag=TMath::Sqrt(Eb_lep*Eb_lep - mb*mb);
    double pT_bl=b_lep_Pmag/TMath::CosH(evBJet_leptop_4P_.Eta());
    if(Eb_lep>0){
      b_lep_reco_4P.SetPtEtaPhiE(pT_bl,evBJet_leptop_4P_.Eta(),evBJet_leptop_4P_.Phi(),Eb_lep);
      t_lep_4P=W_lep_4P+b_lep_reco_4P;
    }
    else{
      error = "leptop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build leptonic top :" << endl;
      (*stream_) << " .... b from lep top " 
              <<  b_lep_reco_4P.Px() << ", "
              <<  b_lep_reco_4P.Py() << ", "
              <<  b_lep_reco_4P.Pz() << ", "
              <<  b_lep_reco_4P.E()     
              << endl;
     (*stream_) << " .... lep top" 
              <<  t_lep_4P.Px() << ", "
              <<  t_lep_4P.Py() << ", "
              <<  t_lep_4P.Pz() << ", "
              <<  t_lep_4P.E()    
              << endl;
    }

  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the tauonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_tau_reco_4P;
  TLorentzVector W_tau_4P;
  TLorentzVector nu_tau_4P;
  TLorentzVector tau_4P;
  TLorentzVector t_tau_4P;

  if(error==0){

    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_ttau),phiNu_ttau);
    double pT_tau = P_TauHad / TMath::CosH(evHadSys_Tau_4P_.Eta());
    double m_tau = Physics::mTau;
    tau_4P.SetPtEtaPhiM(pT_tau,evHadSys_Tau_4P_.Eta(),evHadSys_Tau_4P_.Phi(),m_tau);

    double CosThetataunu_top=TMath::Cos(tau_4P.Angle(enu));
    double Enu=getEqbar_Enu(CosThetataunu_top,tau_4P.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_tau_4P=nu_top_4P_dummy;
    W_tau_4P=nu_tau_4P+tau_4P;

    double Eb_tau=getEb(W_tau_4P, evBJet_hadtop_4P_);
    double b_tau_Pmag=TMath::Sqrt(Eb_tau*Eb_tau - mb*mb);
    double pT_btau=b_tau_Pmag/TMath::CosH(evBJet_hadtop_4P_.Eta());
    if(Eb_tau>0){
      b_tau_reco_4P.SetPtEtaPhiE(pT_btau,evBJet_hadtop_4P_.Eta(),evBJet_hadtop_4P_.Phi(),Eb_tau);
      t_tau_4P=W_tau_4P+b_tau_reco_4P;
    }
    else{
      error = "tautop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build tauonic top :" << endl;
      (*stream_) << " .... b from tau top " 
              <<  b_tau_reco_4P.Px() << ", "
              <<  b_tau_reco_4P.Py() << ", "
              <<  b_tau_reco_4P.Pz() << ", "
              <<  b_tau_reco_4P.E()     
              << endl;
      (*stream_) << " .... tau" 
              <<  tau_4P.Px() << ", "
              <<  tau_4P.Py() << ", "
              <<  tau_4P.Pz() << ", "
              <<  tau_4P.E()    
              << endl;
     (*stream_) << " .... tau top" 
              <<  t_tau_4P.Px() << ", "
              <<  t_tau_4P.Py() << ", "
              <<  t_tau_4P.Pz() << ", "
              <<  t_tau_4P.E()    
              << endl;
    }


  }



  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                Recoil                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector tot4P;
  TLorentzVector rho4P;

  if(error==0){
    
    tot4P = W_4P + t_lep_4P + t_tau_4P ;
    rho4P = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_Tau_4P_ + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_;

  }

  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        JET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

  double T_BJet_leptop = 1.0;
  double T_BJet_tautop = 1.0;

  if ( error == 0 ) {

    if (flagTFBJet_leptop_)
      T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");   
    if (flagTFBJet_hadtop_)
      T_BJet_tautop = getJetTF( b_tau_reco_4P, evBJet_hadtop_4P_,"b");

    // Log
/*    if (outgoingJetsOut_ttbar_DL_) {
      (*outgoingJetsOut_ttbar_DL_) << setiosflags(ios::fixed) << scientific;
      (*outgoingJetsOut_ttbar_DL_) 
	 << " " << setw(21) << evBJet_leptop_4P_.Pt()
         << " " << setw(21) << evBJet_leptop_4P_.Eta()
         << " " << setw(21) << evBJet_leptop_4P_.Phi()
         << " " << setw(21) << evBJet_leptop_4P_.P()
         << " " << setw(21) << evBJet_leptop_4P_.E()      
         << " " << setw(21) << b_lep_reco_4P.Pt()
         << " " << setw(21) << b_lep_reco_4P.Eta()
         << " " << setw(21) << b_lep_reco_4P.Phi()
         << " " << setw(21) << b_lep_reco_4P.P()
         << " " << setw(21) << b_lep_reco_4P.E()
         << " " << setw(21) << T_BJet_leptop
	 << " " << setw(21) << evBJet_hadtop_4P_.Pt()
         << " " << setw(21) << evBJet_hadtop_4P_.Eta()
         << " " << setw(21) << evBJet_hadtop_4P_.Phi()
         << " " << setw(21) << evBJet_hadtop_4P_.P()
         << " " << setw(21) << evBJet_hadtop_4P_.E()      
         << " " << setw(21) << b_tau_reco_4P.Pt()
         << " " << setw(21) << b_tau_reco_4P.Eta()
         << " " << setw(21) << b_tau_reco_4P.Phi()
         << " " << setw(21) << b_tau_reco_4P.P()
         << " " << setw(21) << b_tau_reco_4P.E()
         << " " << setw(21) << T_BJet_tautop
         << endl;
    } */
  }




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        MET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_MET=1.0;
  if ( error == 0 ) {
    if (flagTFMET_) {

      TVector3 RhoT;
      TVector3 PT;
      
      if(include_hadrecoil_){
	
	RhoT.SetXYZ((rho4P).Px(), (rho4P).Py(), 0);
	PT.SetXYZ((tot4P).Px(), (tot4P).Py(), 0);

      }

      else{
	
	TLorentzVector tot4P_nohad = W_4P + W_lep_4P + W_tau_4P;
	TLorentzVector rho4P_nohad = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_Tau_4P_ + evLep_top_4P_ ;
	RhoT.SetXYZ((rho4P_nohad).Px(), (rho4P_nohad).Py(), 0);
	PT.SetXYZ((tot4P_nohad).Px(), (tot4P_nohad).Py(), 0);


      }


      TVector3 Rec = RhoT-PT;
      // Order : (0,0), (0,1), (1,0), (1,1) ???
      double sqrtDetV = TMath::Sqrt( std::abs( evV_[0] * evV_[3] - evV_[1] * evV_[2]));
      double f = Rec.X()*(evV_[0]*Rec.X() + evV_[1]*Rec.Y()) 
               + Rec.Y()*(evV_[2]*Rec.X() + evV_[3]*Rec.Y());
      T_MET =  TMath::Exp( -0.5*f ) * sqrtDetV / ( 2 * TMath::Pi() ) ;
    

    }
  } 



  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                       TOP Transfer functions                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_lepTop = 1.0;
  double T_tauTop = 1.0;

  if ( error == 0 ) {
    if ( flagTFTop_ ){
      T_lepTop = getTopTF(t_lep_4P,b_lep_reco_4P,evLep_top_4P_,nu_lep_4P);
      T_tauTop = getTopTF(t_tau_4P,b_tau_reco_4P,tau_4P,nu_tau_4P);
    }
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                    HADRONIC TAU Transfer function                    @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_hadTau = 1.0;
  if ( error == 0 ) {
    if ( flagTFHadTau_ )
      T_hadTau = getTauHadronicTF(P_TauHad, evHadSys_Tau_4P_);  
  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Jacobian Term                             @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
 
  double Jac = 1.0;
  
  if ( error == 0 ) {  
    if( flagJac_ ) {

      double Jac_W=getJacobian_W( nu_W_4P.E(), evLep_Tau_4P_.E() );
      double Jac_tau = getJacobian_Tau( tau_4P ) * tau_4P.E() ; //Energy already taken into account in toptau Jacpobian 
      double Jac_toptau = getJacobian_toplep( b_tau_reco_4P, W_tau_4P, tau_4P.E(), nu_tau_4P.E());
      double Jac_toplep = getJacobian_toplep( b_lep_reco_4P, W_lep_4P, evLep_top_4P_.E(), nu_lep_4P.E());
      Jac = Jac_tau * Jac_toptau * Jac_toplep * Jac_W;

      if( TMath::IsNaN( Jac ) ) { if (error == 0) error = "Jacobian";}  
    }
  }




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                 BOOST                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double totPt = 0.0;
  TLorentzVector tot4P_boosted;
  TLorentzVector W_4P_boosted;
  TLorentzVector t_tau_4P_boosted;
  TLorentzVector t_lep_4P_boosted;

  if ( error == 0 ) {
    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Boosting :" << endl;
      (*stream_) << "  initial total momentum tot4P " 
              <<  tot4P.Px() << ", "
              <<  tot4P.Py() << ", "
              <<  tot4P.Pz() << ", "
              <<  tot4P.E() 
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
    totPt = tot4P.Pt();
    
    // Set transverse boost
    TVector3 boost3D = tot4P.BoostVector();
    boost3D.SetZ(0);
    
    //
    // Application of the boost to the final state 4-momenta
    // WARRNING : only following 4-momenta are now booosted
    tot4P_boosted=tot4P;
    tot4P_boosted.Boost   ( -boost3D );

    // ??? implulsion > energy 
    if( TMath::IsNaN(tot4P_boosted.E()) ) { if (error == 0) error = "Boost";}  

    W_4P_boosted=W_4P;
    t_tau_4P_boosted=t_tau_4P;
    t_lep_4P_boosted=t_lep_4P;

    W_4P_boosted.Boost( -boost3D );
    t_tau_4P_boosted.Boost( -boost3D );
    t_lep_4P_boosted.Boost( -boost3D );

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "  Boosted total momentum tot4P " 
              <<  tot4P_boosted.Px() << ", "
              <<  tot4P_boosted.Py() << ", "
              <<  tot4P_boosted.Pz() << ", "
              <<  tot4P_boosted.E() 
              << endl;
      (*stream_) << "  Boost magnitude " 
              <<  totPt
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        THE BJORKEN FRACTION                          @@@
  //@@@                         Matrix Elements                              @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double wME = 1.0;
  double x_a;
  double x_b;
  x_a = (tot4P_boosted.E() + tot4P_boosted.Pz()) / sqrtS_ ;
  x_b = (tot4P_boosted.E() - tot4P_boosted.Pz()) / sqrtS_ ;
  
  if ( (x_a > 1.0) || (x_a < 0.0) || (x_b > 1.0 ) ||  (x_b < 0.0 ) ) 
    error = "BFrac";
  
  if ( error == 0 ) {

    TLorentzVector top_4P;
    TLorentzVector Atop_4P;

    if(lepton_top_charge_>0){
      top_4P=t_lep_4P_boosted;
      Atop_4P=t_tau_4P_boosted;
    }
    else{
      top_4P=t_tau_4P_boosted;
      Atop_4P=t_lep_4P_boosted;
    }


    // Parton momentum fraction
    //
    // MadGraph quadri-vector construction
    double inPart_a[4], inPart_b[4]; 
    double outW[4];
    double outtop[4], outAtop[4];
    double inE;
    // Incoming quarks
    inE =  0.5 * x_a * sqrtS_;
    inPart_a[0] = inE; inPart_a[1] = 0.; inPart_a[2] =  0; inPart_a[3] = inE; 
    inE =  0.5 * x_b * sqrtS_;
    inPart_b[0] = inE; inPart_b[1] = 0.; inPart_b[2] =  0; inPart_b[3] =-inE; 
    // W / Fake
    outW[0] = W_4P.E(); outW[1] = W_4P.Px(); outW[2] = W_4P.Py(); outW[3] = W_4P.Pz();
    // top / Atop
    outtop[0] = top_4P.E(); outtop[1] = top_4P.Px(); outtop[2] = top_4P.Py(); outtop[3] = top_4P.Pz();
    outAtop[0] = Atop_4P.E(); outAtop[1] = Atop_4P.Px(); outAtop[2] = Atop_4P.Py(); outAtop[3] = Atop_4P.Pz();

    // Quadri-vector order for MadGraph
    //   p[0] incoming parton0
    //   p[1] incoming parton1
    //   p[2] outgoing top
    //   p[3] outgoing Atop
    //   p[4] outgoing W
    vector<double*> p(6,(double*) 0);
    p[0] = inPart_a ; p[1] = inPart_b; 
    p[2] = outtop; p[3] = outAtop; 
    p[4] = outW   ; 

    // Matrix Element
    if (flagWME_) {
      if (lepton_Tau_charge_<0){	
	wME = get_ttWmWeightedME( x_a, x_b, p);
      }
      else{
	wME = get_ttWpWeightedME ( x_a, x_b, p );
      }
      if( TMath::IsNaN(wME) ) { if (error == 0) error = "WeightME";}
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Weight evaluation                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  double Eval =  T_hadTau * T_BJet_tautop * T_BJet_leptop * T_MET * T_lepTop * T_tauTop * Jac * wME;

  if (error) Eval = 0;

  if (verbose_ >= IntegrandLevel){
    // if( Eval > 10.0e-40) {
    (*stream_) << "Integrand value :" 
	       << " T_hadTau = " << T_hadTau
	       << " T_BJet_tautop = " << T_BJet_tautop
	       << " T_BJet_leptop = " << T_BJet_leptop
	       << " T_MET = " << T_MET 
	       << " T_lepTop = " << T_lepTop
	       << " T_tauTop = " << T_tauTop
	       << " Jac = " << Jac
	       << " wME = " << wME 
	       << " Eval = "    << Eval
	       << " Error= " << error
	       << endl;
    // }
  }


/*  const char *err_str = "Ok";
  if (error) 
    err_str = error;
  else{*/
  if ( !error) {
    integr_EfficiencyttW_++;
  }


  tot_DrawsttW_++;

/*      
  if (valuesOut_ttW_) {
    writeFctValues_ttW(eventID_, lepton_Tau_charge_,
		       cosThetaNu_tlep, phiNu_tlep,
		       cosThetaNu_ttau, phiNu_ttau,
		       P_TauHad,
		       cosThetaNu_W, phiNu_W,
		       totPt,
		       rho4P, tot4P,
		       T_BJet_tautop, T_BJet_leptop,
		       T_MET, T_hadTau,
		       T_lepTop, T_tauTop,
		       Jac,  wME, Eval, err_str);
  }

  if (valuesOutRoot_ttW_) {
    valuesOutRoot_ttW_->fill(eventID_, lepton_Tau_charge_,	       
			     cosThetaNu_tlep, phiNu_tlep,
			     cosThetaNu_ttau, phiNu_ttau,
			     P_TauHad,
			     cosThetaNu_W, phiNu_W,
			     totPt,
			     rho4P, tot4P,
			     T_BJet_tautop, T_BJet_leptop,
			     T_MET, T_hadTau,
			     T_lepTop, T_tauTop,
			     Jac,  wME, Eval, err_str);

  }

*/

  return Eval;
}

/*

void MGIntegration::writeFctValues_ttW(int64_t eventID, bool lep_sign, 
				       double cosThetaNu_tlep, double phiNu_tlep,
				       double cosThetaNu_ttau, double phiNu_ttau,
				       double PTauHad,
				       double cosThetaNu_W, double phiNu_W,
				       double boost, 
				       TLorentzVector rho4P, TLorentzVector boost4P, 
				       double TFBjet_tautop, double TFBjet_leptop,
				       double TFMET, double TFHadTau,
				       double TFLepTop, double TFTauTop,
				       double Jac, double wME, 
				       double eval,				
				       const char *errStr) const {
  if (valuesOut_ttW_) {
    (*valuesOut_ttW_) << setiosflags(ios::fixed) << scientific;
    (*valuesOut_ttW_) 
      << " " << setw(9) << eventID
      << " " << setw(9) << lep_sign
      << " " << setw(15) << cosThetaNu_tlep
      << " " << setw(15) << phiNu_tlep    
      << " " << setw(15) << cosThetaNu_ttau
      << " " << setw(15) << phiNu_ttau
      << " " << setw(15) << PTauHad
      << " " << setw(15) << cosThetaNu_W
      << " " << setw(15) << phiNu_W
      << " " << setw(15) << boost
      << " " << setw(15) << rho4P.Px()
      << " " << setw(15) << rho4P.Py()
      << " " << setw(15) << rho4P.Pz()
      << " " << setw(15) << rho4P.E()
      << " " << setw(15) << boost4P.Px()
      << " " << setw(15) << boost4P.Py()
      << " " << setw(15) << boost4P.Pz()
      << " " << setw(15) << boost4P.E()
      << " " << setw(15) << TFBjet_tautop
      << " " << setw(15) << TFBjet_leptop
      << " " << setw(15) << TFMET
      << " " << setw(15) << TFHadTau
      << " " << setw(15) << TFLepTop
      << " " << setw(15) << TFTauTop
      << " " << setw(15) << Jac
      << " " << setw(15) << wME 
      << " " << setw(15) << eval
      << " " << setw(15) << errStr
      << endl; 


  }  
}

*/

double MGIntegration::evalttjets(const double* x ) {


  double mb=Physics::mb;

 
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@                   Take the input variable form VEGAS                 @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double  EQuark1	  = x[EQuark1_ttjets_id];         //module of the momentum of the final state Quark1
  double  cosThetaNu_tlep = x[cosThetaNu_tlep_ttjets_id]; // cosTheta of the nu from the leptonic top
  double  phiNu_tlep = x[phiNu_tlep_ttjets_id];           // phi of the nu from the leptonic top
  double  P_outFake	  = x[POutFake_ttjets_id];        //module of the momentum of the final state parton wich fakes a tauh

  double cosTheta_missing_jet = 0;
  double phi_missing_jet      = 0;
  if(integration_type_ == integration_type_w_miss){
    cosTheta_missing_jet = x[cosTheta_missing_jet_ttjets_id];  // cosTheta of the quark associated to the missing jet
    phi_missing_jet      = x[phi_missing_jet_ttjets_id];   // phi of the quark associated to the missing jet
  }


  const char* error = 0;
  // To remove ???
//  char str[256];
  

  if (verbose_ >= IntegrandLevel){
    (*stream_) << "Vegas point :" 
	       << " " << " cosTheta_missing_jet = "    << cosTheta_missing_jet
	       << " " << " phi_missing_jet = "    << phi_missing_jet
	       << " " << " EQuark1 = "    << EQuark1
	       << " " << " cosThetaNu_tlep = "    << cosThetaNu_tlep
	       << " " << " phiNu_tlep = "    << phiNu_tlep
	       << " " << " P_outFake = "    << P_outFake
	       << " " << endl;
  }



  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the parton fake 4-Vector                 @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector outFake_4P;
  double PtFake   =   P_outFake * TMath::Sin( evHadSys_Tau_4P_.Theta() );
  outFake_4P.SetPtEtaPhiE(PtFake, evHadSys_Tau_4P_.Eta(), evHadSys_Tau_4P_.Phi(), P_outFake);  


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the leptonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_lep_reco_4P;
  TLorentzVector W_lep_4P;
  TLorentzVector nu_top_4P;
  TLorentzVector t_lep_4P;

  if(error==0){
    
    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_tlep),phiNu_tlep);
    double CosThetalnu_top=TMath::Cos(evLep_top_4P_.Angle(enu));
    double Enu=getEqbar_Enu(CosThetalnu_top,evLep_top_4P_.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_top_4P=nu_top_4P_dummy;
    W_lep_4P=nu_top_4P+evLep_top_4P_;
    
    double Eb_lep=getEb(W_lep_4P, evBJet_leptop_4P_);
    double b_lep_Pmag=TMath::Sqrt(Eb_lep*Eb_lep - mb*mb);
    double pT_bl=b_lep_Pmag/TMath::CosH(evBJet_leptop_4P_.Eta());
    if(Eb_lep>0){
      b_lep_reco_4P.SetPtEtaPhiE(pT_bl,evBJet_leptop_4P_.Eta(),evBJet_leptop_4P_.Phi(),Eb_lep);
      t_lep_4P=W_lep_4P+b_lep_reco_4P;
    }
    else{
      error = "leptop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build leptonic top :" << endl;
      (*stream_) << " .... b from lep top " 
              <<  b_lep_reco_4P.Px() << ", "
              <<  b_lep_reco_4P.Py() << ", "
              <<  b_lep_reco_4P.Pz() << ", "
              <<  b_lep_reco_4P.E()     
              << endl;
     (*stream_) << " .... lep top" 
              <<  t_lep_4P.Px() << ", "
              <<  t_lep_4P.Py() << ", "
              <<  t_lep_4P.Pz() << ", "
              <<  t_lep_4P.E()    
              << endl;
    }

  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the hadronic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector q1_reco_4P;
  TLorentzVector q2_reco_4P;
  TLorentzVector b_had_reco_4P;
  TLorentzVector W_had_4P;
  TLorentzVector t_had_4P;

  if(error==0){

    double pt_q1=EQuark1/TMath::CosH(evJet1_4P_.Eta());
    q1_reco_4P.SetPtEtaPhiE(pt_q1,evJet1_4P_.Eta(),evJet1_4P_.Phi(),EQuark1);

    if(integration_type_ == integration_type_wo_miss){
      double CosTheta_qq=TMath::Cos(evJet1_4P_.Angle(evJet2_4P_.Vect()));
      double EQuark2=getEqbar_Enu(CosTheta_qq,EQuark1);
      double pt_q2=EQuark2/TMath::CosH(evJet2_4P_.Eta());
      q2_reco_4P.SetPtEtaPhiE(pt_q2,evJet2_4P_.Eta(),evJet2_4P_.Phi(),EQuark2);
    }
    else if(integration_type_ == integration_type_w_miss){
      TVector3 e_q2;
      e_q2.SetMagThetaPhi(1,TMath::ACos(cosTheta_missing_jet),phi_missing_jet);
      double CosTheta_qq=TMath::Cos(evJet1_4P_.Angle(e_q2));
      double EQuark2=getEqbar_Enu(CosTheta_qq,EQuark1);
      double pt_q2=EQuark2/TMath::CosH(e_q2.Eta());
      q2_reco_4P.SetPtEtaPhiE(pt_q2,e_q2.Eta(),e_q2.Phi(),EQuark2);
    }
    
    W_had_4P=q1_reco_4P+q2_reco_4P;
    double Eb_had=getEb(W_had_4P, evBJet_hadtop_4P_);
    double b_had_Pmag=TMath::Sqrt(Eb_had*Eb_had - mb*mb);
    double pT_bh=b_had_Pmag/TMath::CosH(evBJet_hadtop_4P_.Eta());
    if(Eb_had>0){
      b_had_reco_4P.SetPtEtaPhiE(pT_bh,evBJet_hadtop_4P_.Eta(),evBJet_hadtop_4P_.Phi(),Eb_had);
      t_had_4P=W_had_4P+b_had_reco_4P;
    }
    else{
      error = "hadtop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build hadronic top :" << endl;
      (*stream_) << " .... b from had top " 
              <<  b_had_reco_4P.Px() << ", "
              <<  b_had_reco_4P.Py() << ", "
              <<  b_had_reco_4P.Pz() << ", "
              <<  b_had_reco_4P.E()     
              << endl;
      (*stream_) << " .... q1 from had top " 
              <<  q1_reco_4P.Px() << ", "
              <<  q1_reco_4P.Py() << ", "
              <<  q1_reco_4P.Pz() << ", "
              <<  q1_reco_4P.E()     
              << endl;
      (*stream_) << " .... q2 from had top " 
              <<  q2_reco_4P.Px() << ", "
              <<  q2_reco_4P.Py() << ", "
              <<  q2_reco_4P.Pz() << ", "
              <<  q2_reco_4P.E()     
              << endl;
      (*stream_) << " .... had top " 
		 <<  t_had_4P.Px() << ", "
		 <<  t_had_4P.Py() << ", "
		 <<  t_had_4P.Pz() << ", "
		 <<  t_had_4P.E()    
		 << endl;
    }


  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                Recoil                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector tot4P;
  TLorentzVector rho4P;

  if(error==0){
    
    tot4P = outFake_4P + t_lep_4P + t_had_4P ;
    if(integration_type_ == integration_type_wo_miss)
      rho4P = evRecoMET4P_ + evHadSys_Tau_4P_  + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_ + evJet1_4P_ + evJet2_4P_ ;
    else if(integration_type_ == integration_type_wo_miss)
      rho4P = evRecoMET4P_ + evHadSys_Tau_4P_  + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_ + evJet1_4P_ ;

  }

  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        JET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

  double T_Jet1 = 1.0;
  double T_Jet2 = 1.0;
  double T_BJet_leptop = 1.0;
  double T_BJet_hadtop = 1.0;
  double T_Fake = 1.0;

  if ( error == 0 ) {
    if (flagTFJet1_) {
      if(integration_type_ == integration_type_wo_miss)
	T_Jet1 = getJetTF( q1_reco_4P, evJet1_4P_,"light");
      
      else if(integration_type_ == integration_type_w_miss){
	if(q1_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector q1q2_reco_4P = q1_reco_4P;
	  q1q2_reco_4P.SetE(q1_reco_4P.E() + q2_reco_4P.E());
	  T_Jet1 = getJetTF( q1q2_reco_4P, evJet1_4P_,"light");
	}
	else
	  T_Jet1 = getJetTF( q1_reco_4P, evJet1_4P_,"light");
      }
    }

    //
    if (flagTFJet2_) {
      if(integration_type_ == integration_type_wo_miss)
	T_Jet2 = getJetTF( q2_reco_4P, evJet2_4P_,"light");
      else if(integration_type_ == integration_type_w_miss)
	T_Jet2 = getMissJetAcceptance( q2_reco_4P);
    } 
    //

    if (flagTFBJet_leptop_) {
      if(integration_type_ == integration_type_wo_miss)
	T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");
      
      else if(integration_type_ == integration_type_w_miss){
	if(b_lep_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector b_lep_q2_reco_4P = b_lep_reco_4P;
	  b_lep_q2_reco_4P.SetE(b_lep_reco_4P.E() + q2_reco_4P.E());
	  T_BJet_leptop = getJetTF( b_lep_q2_reco_4P, evBJet_leptop_4P_,"b");
	}
	else
	  T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");
      }
    }

    //
    if (flagTFBJet_hadtop_) {
      if(integration_type_ == integration_type_wo_miss)
      	T_BJet_hadtop = getJetTF( b_had_reco_4P, evBJet_hadtop_4P_,"b");

      else if(integration_type_ == integration_type_w_miss){
	if(b_had_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector b_had_q2_reco_4P = b_had_reco_4P;
	  b_had_q2_reco_4P.SetE(b_had_reco_4P.E() + q2_reco_4P.E());
	  T_BJet_hadtop = getJetTF( b_had_q2_reco_4P, evBJet_leptop_4P_,"b");
	}
	else
	  T_BJet_hadtop = getJetTF( b_had_reco_4P, evBJet_hadtop_4P_,"b");
      }
    }

    //
    if (flagTFFake_) {
      T_Fake = getFakeTF( outFake_4P, evHadSys_Tau_4P_);
    } 

    // Log
/*    if (outgoingJetsOut_ttjets_) {
      (*outgoingJetsOut_ttjets_) << setiosflags(ios::fixed) << scientific;
      (*outgoingJetsOut_ttjets_) 
         << " " << setw(21) << evJet1_4P_.Pt()
         << " " << setw(21) << evJet1_4P_.Eta()
         << " " << setw(21) << evJet1_4P_.Phi()
         << " " << setw(21) << evJet1_4P_.P()
         << " " << setw(21) << evJet1_4P_.E()      
         << " " << setw(21) << q1_reco_4P.Pt()
         << " " << setw(21) << q1_reco_4P.Eta()
         << " " << setw(21) << q1_reco_4P.Phi()
         << " " << setw(21) << q1_reco_4P.P()
         << " " << setw(21) << q1_reco_4P.E()
         << " " << setw(21) << T_Jet1
         << " " << setw(21) << evJet2_4P_.Pt()
         << " " << setw(21) << evJet2_4P_.Eta()
         << " " << setw(21) << evJet2_4P_.Phi()
         << " " << setw(21) << evJet2_4P_.P()
         << " " << setw(21) << evJet2_4P_.E()      
         << " " << setw(21) << q2_reco_4P.Pt()
         << " " << setw(21) << q2_reco_4P.Eta()
         << " " << setw(21) << q2_reco_4P.Phi()
         << " " << setw(21) << q2_reco_4P.P()
         << " " << setw(21) << q2_reco_4P.E()
         << " " << setw(21) << T_Jet2
	 << " " << setw(21) << evBJet_leptop_4P_.Pt()
         << " " << setw(21) << evBJet_leptop_4P_.Eta()
         << " " << setw(21) << evBJet_leptop_4P_.Phi()
         << " " << setw(21) << evBJet_leptop_4P_.P()
         << " " << setw(21) << evBJet_leptop_4P_.E()      
         << " " << setw(21) << b_lep_reco_4P.Pt()
         << " " << setw(21) << b_lep_reco_4P.Eta()
         << " " << setw(21) << b_lep_reco_4P.Phi()
         << " " << setw(21) << b_lep_reco_4P.P()
         << " " << setw(21) << b_lep_reco_4P.E()
         << " " << setw(21) << T_BJet_leptop
	 << " " << setw(21) << evBJet_hadtop_4P_.Pt()
         << " " << setw(21) << evBJet_hadtop_4P_.Eta()
         << " " << setw(21) << evBJet_hadtop_4P_.Phi()
         << " " << setw(21) << evBJet_hadtop_4P_.P()
         << " " << setw(21) << evBJet_hadtop_4P_.E()      
         << " " << setw(21) << b_had_reco_4P.Pt()
         << " " << setw(21) << b_had_reco_4P.Eta()
         << " " << setw(21) << b_had_reco_4P.Phi()
         << " " << setw(21) << b_had_reco_4P.P()
         << " " << setw(21) << b_had_reco_4P.E()
         << " " << setw(21) << T_BJet_hadtop
         << endl;
    } */
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        MET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_MET=1.0;
  if ( error == 0 ) {
    if (flagTFMET_) {

      TVector3 RhoT;
      TVector3 PT;
      
      if(include_hadrecoil_){
	
	RhoT.SetXYZ((rho4P).Px(), (rho4P).Py(), 0);
	PT.SetXYZ((tot4P).Px(), (tot4P).Py(), 0);

      }

      else{
	
	TLorentzVector tot4P_nohad = outFake_4P + W_lep_4P ;
	TLorentzVector rho4P_nohad = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_top_4P_ ;
	RhoT.SetXYZ((rho4P_nohad).Px(), (rho4P_nohad).Py(), 0);
	PT.SetXYZ((tot4P_nohad).Px(), (tot4P_nohad).Py(), 0);


      }


      TVector3 Rec = RhoT-PT;
      // Order : (0,0), (0,1), (1,0), (1,1) ???
      double sqrtDetV = TMath::Sqrt( std::abs( evV_[0] * evV_[3] - evV_[1] * evV_[2]));
      double f = Rec.X()*(evV_[0]*Rec.X() + evV_[1]*Rec.Y()) 
               + Rec.Y()*(evV_[2]*Rec.X() + evV_[3]*Rec.Y());
      T_MET =  TMath::Exp( -0.5*f ) * sqrtDetV / ( 2 * TMath::Pi() ) ;
    

    }
  } 



  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                       TOP Transfer functions                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_lepTop = 1.0;
  double T_hadTop = 1.0;

  if ( error == 0 ) {
    if ( flagTFTop_ ){
      T_lepTop = getTopTF(t_lep_4P,b_lep_reco_4P,evLep_top_4P_,nu_top_4P);

      //We don't know which q is the quark or the antiquark
      T_hadTop = getTopTF(t_had_4P,b_had_reco_4P,q1_reco_4P,q2_reco_4P)
	+ getTopTF(t_had_4P,b_had_reco_4P,q2_reco_4P,q1_reco_4P);

    }
  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Jacobian Term                             @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
 
  double Jac = 1.0;
  
  if ( error == 0 ) {  
    if( flagJac_ ) {

      double cosTheta_qq=TMath::Cos(q1_reco_4P.Angle(q2_reco_4P.Vect()));
      double Jac_tophad = getJacobian_tophad( b_had_reco_4P, W_had_4P, q2_reco_4P.E(), cosTheta_qq);
      double Jac_toplep = getJacobian_toplep( b_lep_reco_4P, W_lep_4P, evLep_top_4P_.E(), nu_top_4P.E());
      Jac = Jac_tophad * Jac_toplep;

      if( TMath::IsNaN( Jac ) ) { if (error == 0) error = "Jacobian";}  
    }
  }




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                 BOOST                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double totPt = 0.0;
  TLorentzVector tot4P_boosted;
  TLorentzVector outFake_4P_boosted;
  TLorentzVector t_had_4P_boosted;
  TLorentzVector t_lep_4P_boosted;

  if ( error == 0 ) {
    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Boosting :" << endl;
      (*stream_) << "  initial total momentum tot4P " 
              <<  tot4P.Px() << ", "
              <<  tot4P.Py() << ", "
              <<  tot4P.Pz() << ", "
              <<  tot4P.E() 
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
    totPt = tot4P.Pt();
    
    // Set transverse boost
    TVector3 boost3D = tot4P.BoostVector();
    boost3D.SetZ(0);
    
    //
    // Application of the boost to the final state 4-momenta
    // WARRNING : only following 4-momenta are now booosted
    tot4P_boosted=tot4P;
    tot4P_boosted.Boost   ( -boost3D );

    // ??? implulsion > energy 
    if( TMath::IsNaN(tot4P_boosted.E()) ) { if (error == 0) error = "Boost";}  

    outFake_4P_boosted=outFake_4P;   
    t_had_4P_boosted=t_had_4P;
    t_lep_4P_boosted=t_lep_4P;

    outFake_4P_boosted.Boost( -boost3D );
    t_had_4P_boosted.Boost( -boost3D );
    t_lep_4P_boosted.Boost( -boost3D );

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "  Boosted total momentum tot4P " 
              <<  tot4P_boosted.Px() << ", "
              <<  tot4P_boosted.Py() << ", "
              <<  tot4P_boosted.Pz() << ", "
              <<  tot4P_boosted.E() 
              << endl;
      (*stream_) << "  Boost magnitude " 
              <<  totPt
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        THE BJORKEN FRACTION                          @@@
  //@@@                         Matrix Elements                              @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double wME = 1.0;
  double x_a;
  double x_b;
  x_a = (tot4P_boosted.E() + tot4P_boosted.Pz()) / sqrtS_ ;
  x_b = (tot4P_boosted.E() - tot4P_boosted.Pz()) / sqrtS_ ;
  
  if ( (x_a > 1.0) || (x_a < 0.0) || (x_b > 1.0 ) ||  (x_b < 0.0 ) ) 
    error = "BFrac";
  
  if ( error == 0 ) {

    TLorentzVector top_4P;
    TLorentzVector Atop_4P;

    if(lepton_top_charge_>0){
      top_4P=t_lep_4P_boosted;
      Atop_4P=t_had_4P_boosted;
    }
    else{
      top_4P=t_had_4P_boosted;
      Atop_4P=t_lep_4P_boosted;
    }


    // Parton momentum fraction
    //
    // MadGraph quadri-vector construction
    double inPart_a[4], inPart_b[4]; 
    double outFake[4];
    double outtop[4], outAtop[4];
    double inE;
    // Incoming quarks
    inE =  0.5 * x_a * sqrtS_;
    inPart_a[0] = inE; inPart_a[1] = 0.; inPart_a[2] =  0; inPart_a[3] = inE; 
    inE =  0.5 * x_b * sqrtS_;
    inPart_b[0] = inE; inPart_b[1] = 0.; inPart_b[2] =  0; inPart_b[3] =-inE; 
    // Fake
    outFake[0] = outFake_4P.E(); outFake[1] = outFake_4P.Px(); outFake[2] = outFake_4P.Py(); outFake[3] = outFake_4P.Pz();
    // top / Atop
    outtop[0] = top_4P.E(); outtop[1] = top_4P.Px(); outtop[2] = top_4P.Py(); outtop[3] = top_4P.Pz();
    outAtop[0] = Atop_4P.E(); outAtop[1] = Atop_4P.Px(); outAtop[2] = Atop_4P.Py(); outAtop[3] = Atop_4P.Pz();

    // Quadri-vector order for MadGraph
    //   p[0] incoming parton0
    //   p[1] incoming parton1
    //   p[2] outgoing top
    //   p[3] outgoing Atop
    //   p[4] outgoing fake
    vector<double*> p(5,(double*) 0);
    p[0] = inPart_a ; p[1] = inPart_b; 
    p[2] = outtop; p[3] = outAtop; 
    p[4] = outFake ;

    // Matrix Element
    if (flagWME_) {
      wME = get_ttjetsWeightedME( x_a, x_b, p);
      if( TMath::IsNaN(wME) ) { if (error == 0) error = "WeightME";}
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Weight evaluation                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  double Eval =  T_Jet1 * T_Jet2 * T_BJet_hadtop * T_BJet_leptop * T_MET * T_Fake * T_lepTop * T_hadTop * Jac * wME;
  if (error) Eval = 0;

  if (verbose_ >= IntegrandLevel){
    // if( Eval > 10.0e-40) {
    (*stream_) << "Integrand value :" 
       << " T_Jet1 = " << T_Jet1
       << " T_Jet2 = " << T_Jet2
       << " T_BJet_hadtop = " << T_BJet_hadtop
       << " T_BJet_leptop = " << T_BJet_leptop
       << " T_MET = " << T_MET 
       << " T_fake = " << T_Fake
       << " T_lepTop = " << T_lepTop
       << " T_hadTop = " << T_hadTop
       << " Jac = " << Jac
       << " wME = " << wME 
       << " Eval = "    << Eval
       << " Error= " << error
       << endl;
    // }
  }
/*  const char *err_str = "Ok";
  if (error) 
    err_str = error;
  else{*/
  if ( !error) {
    integr_Efficiencyttjets_++;
  }


  tot_Drawsttjets_++;

/*      
  if (valuesOut_ttjets_) {
    writeFctValues_ttjets(eventID_,
			  cosTheta_missing_jet, phi_missing_jet,
			  P_outFake,
			  EQuark1,
			  cosThetaNu_tlep, phiNu_tlep,
			  totPt,
			  rho4P, tot4P,
			  T_Jet1, T_Jet2,
			  T_BJet_hadtop, T_BJet_leptop,
			  T_MET, T_Fake,
			  T_lepTop, T_hadTop,
			  Jac,  wME, Eval, err_str);
  }

  if (valuesOutRoot_ttjets_) {
    valuesOutRoot_ttjets_->fill(eventID_,
				cosTheta_missing_jet, phi_missing_jet,
				P_outFake,
				EQuark1,
				cosThetaNu_tlep, phiNu_tlep,
				totPt,
				rho4P, tot4P,
				T_Jet1, T_Jet2,
				T_BJet_hadtop, T_BJet_leptop,
				T_MET, T_Fake,
				T_lepTop, T_hadTop,
				Jac,  wME, Eval, err_str);
    
  }

*/

  return Eval;
}

/*

void MGIntegration::writeFctValues_ttjets(int64_t eventID,
					  double cosTheta_miss_jet, double phi_miss_jet,
					  double POutFake,
					  double EQuark1,
					  double cosThetaNu_tlep, double phiNu_tlep,
					  double boost, 
					  TLorentzVector rho4P, TLorentzVector boost4P,
					  double TFJet1, double TFJet2,
					  double TFBjet_hadtop, double TFBjet_leptop,
					  double TFMET, 
					  double TFFake,
					  double TFLepTop, double TFHadTop,
					  double Jac, double wME, 
					  double eval,
					  const char *errStr) const {
  if (valuesOut_ttjets_) {
    (*valuesOut_ttjets_) << setiosflags(ios::fixed) << scientific;
    (*valuesOut_ttjets_) 
       << " " << setw(9) << eventID
       << " " << setw(15) << cosTheta_miss_jet
       << " " << setw(15) << phi_miss_jet
       << " " << setw(15) << POutFake      
       << " " << setw(15) << EQuark1
       << " " << setw(15) << cosThetaNu_tlep
       << " " << setw(15) << phiNu_tlep           
       << " " << setw(15) << rho4P.Px()
       << " " << setw(15) << rho4P.Py()
       << " " << setw(15) << rho4P.Pz()
       << " " << setw(15) << rho4P.E()
       << " " << setw(15) << boost4P.Px()
       << " " << setw(15) << boost4P.Py()
       << " " << setw(15) << boost4P.Pz()
       << " " << setw(15) << boost4P.E()
       << " " << setw(15) << boost
       << " " << setw(15) << TFJet1
       << " " << setw(15) << TFJet2
       << " " << setw(15) << TFBjet_hadtop
       << " " << setw(15) << TFBjet_leptop
       << " " << setw(15) << TFMET
       << " " << setw(15) << TFFake
       << " " << setw(15) << TFLepTop
       << " " << setw(15) << TFHadTop
       << " " << setw(15) << Jac
       << " " << setw(15) << wME 
       << " " << setw(15) << eval
       << " " << setw(15) << errStr
       << endl;
  }  
}

*/

double MGIntegration::evalttbar_SL(const double* x ) {


  double mb=Physics::mb;

 
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@                   Take the input variable form VEGAS                 @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double  EQuark1	  = x[EQuark1_ttbar_SL_ttZ_Zll_id];         //module of the momentum of the final state Quark1
  double  cosThetaNu_tlep = x[cosThetaNu_tlep_ttbar_SL_ttZ_Zll_id]; // cosTheta of the nu from the leptonic top
  double  phiNu_tlep = x[phiNu_tlep_ttbar_SL_ttZ_Zll_id];           // phi of the nu from the leptonic top

  double cosTheta_missing_jet = 0;
  double phi_missing_jet      = 0;
  if(integration_type_ == integration_type_w_miss){
    cosTheta_missing_jet = x[cosTheta_missing_jet_ttbar_SL_ttZ_Zll_id];  // cosTheta of the quark associated to the missing jet
    phi_missing_jet      = x[phi_missing_jet_ttbar_SL_ttZ_Zll_id];   // phi of the quark associated to the missing jet
  }


  const char* error = 0;
  // To remove ???
//  char str[256];
  

  if (verbose_ >= IntegrandLevel){
    (*stream_) << "Vegas point :" 
	       << " " << " cosTheta_missing_jet = "    << cosTheta_missing_jet
	       << " " << " phi_missing_jet = "    << phi_missing_jet
	       << " " << " EQuark1 = "    << EQuark1
	       << " " << " cosThetaNu_tlep = "    << cosThetaNu_tlep
	       << " " << " phiNu_tlep = "    << phiNu_tlep
	       << " " << endl;
  }



  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the leptonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_lep_reco_4P;
  TLorentzVector W_lep_4P;
  TLorentzVector nu_top_4P;
  TLorentzVector t_lep_4P;

  if(error==0){
    
    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_tlep),phiNu_tlep);
    double CosThetalnu_top=TMath::Cos(evLep_top_4P_.Angle(enu));
    double Enu=getEqbar_Enu(CosThetalnu_top,evLep_top_4P_.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_top_4P=nu_top_4P_dummy;
    W_lep_4P=nu_top_4P+evLep_top_4P_;
    
    double Eb_lep=getEb(W_lep_4P, evBJet_leptop_4P_);
    double b_lep_Pmag=TMath::Sqrt(Eb_lep*Eb_lep - mb*mb);
    double pT_bl=b_lep_Pmag/TMath::CosH(evBJet_leptop_4P_.Eta());
    if(Eb_lep>0){
      b_lep_reco_4P.SetPtEtaPhiE(pT_bl,evBJet_leptop_4P_.Eta(),evBJet_leptop_4P_.Phi(),Eb_lep);
      t_lep_4P=W_lep_4P+b_lep_reco_4P;
    }
    else{
      error = "leptop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build leptonic top :" << endl;
      (*stream_) << " .... b from lep top " 
              <<  b_lep_reco_4P.Px() << ", "
              <<  b_lep_reco_4P.Py() << ", "
              <<  b_lep_reco_4P.Pz() << ", "
              <<  b_lep_reco_4P.E()     
              << endl;
     (*stream_) << " .... lep top" 
              <<  t_lep_4P.Px() << ", "
              <<  t_lep_4P.Py() << ", "
              <<  t_lep_4P.Pz() << ", "
              <<  t_lep_4P.E()    
              << endl;
    }

  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the hadronic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector q1_reco_4P;
  TLorentzVector q2_reco_4P;
  TLorentzVector b_had_reco_4P;
  TLorentzVector W_had_4P;
  TLorentzVector t_had_4P;

  if(error==0){

    double pt_q1=EQuark1/TMath::CosH(evJet1_4P_.Eta());
    q1_reco_4P.SetPtEtaPhiE(pt_q1,evJet1_4P_.Eta(),evJet1_4P_.Phi(),EQuark1);

    if(integration_type_ == integration_type_wo_miss){
      double CosTheta_qq=TMath::Cos(evJet1_4P_.Angle(evJet2_4P_.Vect()));
      double EQuark2=getEqbar_Enu(CosTheta_qq,EQuark1);
      double pt_q2=EQuark2/TMath::CosH(evJet2_4P_.Eta());
      q2_reco_4P.SetPtEtaPhiE(pt_q2,evJet2_4P_.Eta(),evJet2_4P_.Phi(),EQuark2);
    }
    else if(integration_type_ == integration_type_w_miss){
      TVector3 e_q2;
      e_q2.SetMagThetaPhi(1,TMath::ACos(cosTheta_missing_jet),phi_missing_jet);
      double CosTheta_qq=TMath::Cos(evJet1_4P_.Angle(e_q2));
      double EQuark2=getEqbar_Enu(CosTheta_qq,EQuark1);
      double pt_q2=EQuark2/TMath::CosH(e_q2.Eta());
      q2_reco_4P.SetPtEtaPhiE(pt_q2,e_q2.Eta(),e_q2.Phi(),EQuark2);
    }
    
    W_had_4P=q1_reco_4P+q2_reco_4P;
    double Eb_had=getEb(W_had_4P, evBJet_hadtop_4P_);
    double b_had_Pmag=TMath::Sqrt(Eb_had*Eb_had - mb*mb);
    double pT_bh=b_had_Pmag/TMath::CosH(evBJet_hadtop_4P_.Eta());
    if(Eb_had>0){
      b_had_reco_4P.SetPtEtaPhiE(pT_bh,evBJet_hadtop_4P_.Eta(),evBJet_hadtop_4P_.Phi(),Eb_had);
      t_had_4P=W_had_4P+b_had_reco_4P;
    }
    else{
      error = "hadtop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build hadronic top :" << endl;
      (*stream_) << " .... b from had top " 
              <<  b_had_reco_4P.Px() << ", "
              <<  b_had_reco_4P.Py() << ", "
              <<  b_had_reco_4P.Pz() << ", "
              <<  b_had_reco_4P.E()     
              << endl;
      (*stream_) << " .... q1 from had top " 
              <<  q1_reco_4P.Px() << ", "
              <<  q1_reco_4P.Py() << ", "
              <<  q1_reco_4P.Pz() << ", "
              <<  q1_reco_4P.E()     
              << endl;
      (*stream_) << " .... q2 from had top " 
              <<  q2_reco_4P.Px() << ", "
              <<  q2_reco_4P.Py() << ", "
              <<  q2_reco_4P.Pz() << ", "
              <<  q2_reco_4P.E()     
              << endl;
      (*stream_) << " .... had top " 
		 <<  t_had_4P.Px() << ", "
		 <<  t_had_4P.Py() << ", "
		 <<  t_had_4P.Pz() << ", "
		 <<  t_had_4P.E()    
		 << endl;
    }


  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                Recoil                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector tot4P;
  TLorentzVector rho4P;

  if(error==0){
    
    tot4P = t_lep_4P + t_had_4P ;
    if(integration_type_ == integration_type_wo_miss)
      rho4P = evRecoMET4P_ + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_ + evJet1_4P_ + evJet2_4P_ ;
    else if(integration_type_ == integration_type_wo_miss)
      rho4P = evRecoMET4P_ + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_ + evJet1_4P_ ;

  }

  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        JET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

  double T_Jet1 = 1.0;
  double T_Jet2 = 1.0;
  double T_BJet_leptop = 1.0;
  double T_BJet_hadtop = 1.0;

  if ( error == 0 ) {
    if (flagTFJet1_) {
      if(integration_type_ == integration_type_wo_miss)
	T_Jet1 = getJetTF( q1_reco_4P, evJet1_4P_,"light");
      
      else if(integration_type_ == integration_type_w_miss){
	if(q1_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector q1q2_reco_4P = q1_reco_4P;
	  q1q2_reco_4P.SetE(q1_reco_4P.E() + q2_reco_4P.E());
	  T_Jet1 = getJetTF( q1q2_reco_4P, evJet1_4P_,"light");
	}
	else
	  T_Jet1 = getJetTF( q1_reco_4P, evJet1_4P_,"light");
      }
    }

    //
    if (flagTFJet2_) {
      if(integration_type_ == integration_type_wo_miss)
	T_Jet2 = getJetTF( q2_reco_4P, evJet2_4P_,"light");
      else if(integration_type_ == integration_type_w_miss)
	T_Jet2 = getMissJetAcceptance( q2_reco_4P);
    } 
    //

    if (flagTFBJet_leptop_) {
      if(integration_type_ == integration_type_wo_miss)
	T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");
      
      else if(integration_type_ == integration_type_w_miss){
	if(b_lep_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector b_lep_q2_reco_4P = b_lep_reco_4P;
	  b_lep_q2_reco_4P.SetE(b_lep_reco_4P.E() + q2_reco_4P.E());
	  T_BJet_leptop = getJetTF( b_lep_q2_reco_4P, evBJet_leptop_4P_,"b");
	}
	else
	  T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");

      }
    }

    //
    if (flagTFBJet_hadtop_) {
      if(integration_type_ == integration_type_wo_miss)
	T_BJet_hadtop = getJetTF( b_had_reco_4P, evBJet_hadtop_4P_,"b");
      
      else if(integration_type_ == integration_type_w_miss){
	if(b_had_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector b_had_q2_reco_4P = b_had_reco_4P;
	  b_had_q2_reco_4P.SetE(b_had_reco_4P.E() + q2_reco_4P.E());
	  T_BJet_hadtop = getJetTF( b_had_q2_reco_4P, evBJet_leptop_4P_,"b");
	}
	else
	  T_BJet_hadtop = getJetTF( b_had_reco_4P, evBJet_hadtop_4P_,"b");
      }
    }


    // Log
 /*   if (outgoingJetsOut_ttbar_SL_) {
      (*outgoingJetsOut_ttbar_SL_) << setiosflags(ios::fixed) << scientific;
      (*outgoingJetsOut_ttbar_SL_) 
         << " " << setw(21) << evJet1_4P_.Pt()
         << " " << setw(21) << evJet1_4P_.Eta()
         << " " << setw(21) << evJet1_4P_.Phi()
         << " " << setw(21) << evJet1_4P_.P()
         << " " << setw(21) << evJet1_4P_.E()      
         << " " << setw(21) << q1_reco_4P.Pt()
         << " " << setw(21) << q1_reco_4P.Eta()
         << " " << setw(21) << q1_reco_4P.Phi()
         << " " << setw(21) << q1_reco_4P.P()
         << " " << setw(21) << q1_reco_4P.E()
         << " " << setw(21) << T_Jet1
         << " " << setw(21) << evJet2_4P_.Pt()
         << " " << setw(21) << evJet2_4P_.Eta()
         << " " << setw(21) << evJet2_4P_.Phi()
         << " " << setw(21) << evJet2_4P_.P()
         << " " << setw(21) << evJet2_4P_.E()      
         << " " << setw(21) << q2_reco_4P.Pt()
         << " " << setw(21) << q2_reco_4P.Eta()
         << " " << setw(21) << q2_reco_4P.Phi()
         << " " << setw(21) << q2_reco_4P.P()
         << " " << setw(21) << q2_reco_4P.E()
         << " " << setw(21) << T_Jet2
	 << " " << setw(21) << evBJet_leptop_4P_.Pt()
         << " " << setw(21) << evBJet_leptop_4P_.Eta()
         << " " << setw(21) << evBJet_leptop_4P_.Phi()
         << " " << setw(21) << evBJet_leptop_4P_.P()
         << " " << setw(21) << evBJet_leptop_4P_.E()      
         << " " << setw(21) << b_lep_reco_4P.Pt()
         << " " << setw(21) << b_lep_reco_4P.Eta()
         << " " << setw(21) << b_lep_reco_4P.Phi()
         << " " << setw(21) << b_lep_reco_4P.P()
         << " " << setw(21) << b_lep_reco_4P.E()
         << " " << setw(21) << T_BJet_leptop
	 << " " << setw(21) << evBJet_hadtop_4P_.Pt()
         << " " << setw(21) << evBJet_hadtop_4P_.Eta()
         << " " << setw(21) << evBJet_hadtop_4P_.Phi()
         << " " << setw(21) << evBJet_hadtop_4P_.P()
         << " " << setw(21) << evBJet_hadtop_4P_.E()      
         << " " << setw(21) << b_had_reco_4P.Pt()
         << " " << setw(21) << b_had_reco_4P.Eta()
         << " " << setw(21) << b_had_reco_4P.Phi()
         << " " << setw(21) << b_had_reco_4P.P()
         << " " << setw(21) << b_had_reco_4P.E()
         << " " << setw(21) << T_BJet_hadtop
         << endl;
    } */
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        MET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_MET=1.0;
  if ( error == 0 ) {
    if (flagTFMET_) {

      TVector3 RhoT;
      TVector3 PT;
      
      if(include_hadrecoil_){
	
	RhoT.SetXYZ((rho4P).Px(), (rho4P).Py(), 0);
	PT.SetXYZ((tot4P).Px(), (tot4P).Py(), 0);

      }

      else{
	
	TLorentzVector tot4P_nohad = W_lep_4P ;
	TLorentzVector rho4P_nohad = evRecoMET4P_ + evLep_top_4P_ ;
	RhoT.SetXYZ((rho4P_nohad).Px(), (rho4P_nohad).Py(), 0);
	PT.SetXYZ((tot4P_nohad).Px(), (tot4P_nohad).Py(), 0);


      }


      TVector3 Rec = RhoT-PT;
      // Order : (0,0), (0,1), (1,0), (1,1) ???
      double sqrtDetV = TMath::Sqrt( std::abs( evV_[0] * evV_[3] - evV_[1] * evV_[2]));
      double f = Rec.X()*(evV_[0]*Rec.X() + evV_[1]*Rec.Y()) 
               + Rec.Y()*(evV_[2]*Rec.X() + evV_[3]*Rec.Y());
      T_MET =  TMath::Exp( -0.5*f ) * sqrtDetV / ( 2 * TMath::Pi() ) ;
    

    }
  } 



  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                       TOP Transfer functions                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_lepTop = 1.0;
  double T_hadTop = 1.0;

  if ( error == 0 ) {
    if ( flagTFTop_ ){
      T_lepTop = getTopTF(t_lep_4P,b_lep_reco_4P,evLep_top_4P_,nu_top_4P);

      //We don't know which q is the quark or the antiquark
      T_hadTop = getTopTF(t_had_4P,b_had_reco_4P,q1_reco_4P,q2_reco_4P)
	+ getTopTF(t_had_4P,b_had_reco_4P,q2_reco_4P,q1_reco_4P);

    }
  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Jacobian Term                             @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
 
  double Jac = 1.0;
  
  if ( error == 0 ) {  
    if( flagJac_ ) {

      double cosTheta_qq=TMath::Cos(q1_reco_4P.Angle(q2_reco_4P.Vect()));
      double Jac_tophad = getJacobian_tophad( b_had_reco_4P, W_had_4P, q2_reco_4P.E(), cosTheta_qq);
      double Jac_toplep = getJacobian_toplep( b_lep_reco_4P, W_lep_4P, evLep_top_4P_.E(), nu_top_4P.E());
      Jac = Jac_tophad * Jac_toplep;

      if( TMath::IsNaN( Jac ) ) { if (error == 0) error = "Jacobian";}  
    }
  }




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                 BOOST                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double totPt = 0.0;
  TLorentzVector tot4P_boosted;
  TLorentzVector t_had_4P_boosted;
  TLorentzVector t_lep_4P_boosted;

  if ( error == 0 ) {
    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Boosting :" << endl;
      (*stream_) << "  initial total momentum tot4P " 
              <<  tot4P.Px() << ", "
              <<  tot4P.Py() << ", "
              <<  tot4P.Pz() << ", "
              <<  tot4P.E() 
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
    totPt = tot4P.Pt();
    
    // Set transverse boost
    TVector3 boost3D = tot4P.BoostVector();
    boost3D.SetZ(0);
    
    //
    // Application of the boost to the final state 4-momenta
    // WARRNING : only following 4-momenta are now booosted
    tot4P_boosted=tot4P;
    tot4P_boosted.Boost   ( -boost3D );

    // ??? implulsion > energy 
    if( TMath::IsNaN(tot4P_boosted.E()) ) { if (error == 0) error = "Boost";}  

    t_had_4P_boosted=t_had_4P;
    t_lep_4P_boosted=t_lep_4P;

    t_had_4P_boosted.Boost( -boost3D );
    t_lep_4P_boosted.Boost( -boost3D );

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "  Boosted total momentum tot4P " 
              <<  tot4P_boosted.Px() << ", "
              <<  tot4P_boosted.Py() << ", "
              <<  tot4P_boosted.Pz() << ", "
              <<  tot4P_boosted.E() 
              << endl;
      (*stream_) << "  Boost magnitude " 
              <<  totPt
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        THE BJORKEN FRACTION                          @@@
  //@@@                         Matrix Elements                              @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double wME = 1.0;
  double x_a;
  double x_b;
  x_a = (tot4P_boosted.E() + tot4P_boosted.Pz()) / sqrtS_ ;
  x_b = (tot4P_boosted.E() - tot4P_boosted.Pz()) / sqrtS_ ;
  
  if ( (x_a > 1.0) || (x_a < 0.0) || (x_b > 1.0 ) ||  (x_b < 0.0 ) ) 
    error = "BFrac";
  
  if ( error == 0 ) {

    TLorentzVector top_4P;
    TLorentzVector Atop_4P;

    if(lepton_top_charge_>0){
      top_4P=t_lep_4P_boosted;
      Atop_4P=t_had_4P_boosted;
    }
    else{
      top_4P=t_had_4P_boosted;
      Atop_4P=t_lep_4P_boosted;
    }


    // Parton momentum fraction
    //
    // MadGraph quadri-vector construction
    double inPart_a[4], inPart_b[4]; 
    double outtop[4], outAtop[4];
    double inE;
    // Incoming quarks
    inE =  0.5 * x_a * sqrtS_;
    inPart_a[0] = inE; inPart_a[1] = 0.; inPart_a[2] =  0; inPart_a[3] = inE; 
    inE =  0.5 * x_b * sqrtS_;
    inPart_b[0] = inE; inPart_b[1] = 0.; inPart_b[2] =  0; inPart_b[3] =-inE; 
    // top / Atop
    outtop[0] = top_4P.E(); outtop[1] = top_4P.Px(); outtop[2] = top_4P.Py(); outtop[3] = top_4P.Pz();
    outAtop[0] = Atop_4P.E(); outAtop[1] = Atop_4P.Px(); outAtop[2] = Atop_4P.Py(); outAtop[3] = Atop_4P.Pz();

    // Quadri-vector order for MadGraph
    //   p[0] incoming parton0
    //   p[1] incoming parton1
    //   p[2] outgoing top
    //   p[3] outgoing Atop
    vector<double*> p(4,(double*) 0);
    p[0] = inPart_a ; p[1] = inPart_b; 
    p[2] = outtop; p[3] = outAtop; 

    // Matrix Element
    if (flagWME_) {
      wME = get_ttbarWeightedME( x_a, x_b, p);
      if( TMath::IsNaN(wME) ) { if (error == 0) error = "WeightME";}
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Weight evaluation                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  double Eval =  T_Jet1 * T_Jet2 * T_BJet_hadtop * T_BJet_leptop * T_MET * T_lepTop * T_hadTop * Jac * wME;
  if (error) Eval = 0;

  if (verbose_ >= IntegrandLevel){
    // if( Eval > 10.0e-40) {
    (*stream_) << "Integrand value :" 
       << " T_Jet1 = " << T_Jet1
       << " T_Jet2 = " << T_Jet2
       << " T_BJet_hadtop = " << T_BJet_hadtop
       << " T_BJet_leptop = " << T_BJet_leptop
       << " T_MET = " << T_MET 
       << " T_lepTop = " << T_lepTop
       << " T_hadTop = " << T_hadTop
       << " Jac = " << Jac
       << " wME = " << wME 
       << " Eval = "    << Eval
       << " Error= " << error
       << endl;
    // }
  }
/*  const char *err_str = "Ok";
  if (error) 
    err_str = error;
  else{*/
  if ( !error) {
    integr_Efficiencyttbar_SL_++;
  }


  tot_Drawsttbar_SL_++;

/*      
  if (valuesOut_ttbar_SL_) {
    writeFctValues_ttbar_SL(eventID_,
			  cosTheta_missing_jet, phi_missing_jet,
			  EQuark1,
			  cosThetaNu_tlep, phiNu_tlep,
			  totPt,
			  rho4P, tot4P,
			  T_Jet1, T_Jet2,
			  T_BJet_hadtop, T_BJet_leptop,
			  T_MET,
			  T_lepTop, T_hadTop,
			  Jac,  wME, Eval, err_str);
  }

  if (valuesOutRoot_ttbar_SL_) {
    valuesOutRoot_ttbar_SL_->fill(eventID_,
				cosTheta_missing_jet, phi_missing_jet,
				EQuark1,
				cosThetaNu_tlep, phiNu_tlep,
				totPt,
				rho4P, tot4P,
				T_Jet1, T_Jet2,
				T_BJet_hadtop, T_BJet_leptop,
				T_MET,
				T_lepTop, T_hadTop,
				Jac,  wME, Eval, err_str);
    
  }
*/


  return Eval;
}

/*

void MGIntegration::writeFctValues_ttbar_SL(int64_t eventID,
					  double cosTheta_miss_jet, double phi_miss_jet,
					  double EQuark1,
					  double cosThetaNu_tlep, double phiNu_tlep,
					  double boost, 
					  TLorentzVector rho4P, TLorentzVector boost4P,
					  double TFJet1, double TFJet2,
					  double TFBjet_hadtop, double TFBjet_leptop,
					  double TFMET, 
					  double TFLepTop, double TFHadTop,
					  double Jac, double wME, 
					  double eval,
					  const char *errStr) const {
  if (valuesOut_ttbar_SL_) {
    (*valuesOut_ttbar_SL_) << setiosflags(ios::fixed) << scientific;
    (*valuesOut_ttbar_SL_) 
       << " " << setw(9) << eventID
       << " " << setw(15) << cosTheta_miss_jet
       << " " << setw(15) << phi_miss_jet
       << " " << setw(15) << EQuark1
       << " " << setw(15) << cosThetaNu_tlep
       << " " << setw(15) << phiNu_tlep           
       << " " << setw(15) << rho4P.Px()
       << " " << setw(15) << rho4P.Py()
       << " " << setw(15) << rho4P.Pz()
       << " " << setw(15) << rho4P.E()
       << " " << setw(15) << boost4P.Px()
       << " " << setw(15) << boost4P.Py()
       << " " << setw(15) << boost4P.Pz()
       << " " << setw(15) << boost4P.E()
       << " " << setw(15) << boost
       << " " << setw(15) << TFJet1
       << " " << setw(15) << TFJet2
       << " " << setw(15) << TFBjet_hadtop
       << " " << setw(15) << TFBjet_leptop
       << " " << setw(15) << TFMET
       << " " << setw(15) << TFLepTop
       << " " << setw(15) << TFHadTop
       << " " << setw(15) << Jac
       << " " << setw(15) << wME 
       << " " << setw(15) << eval
       << " " << setw(15) << errStr
       << endl;
  }  
}
*/

double MGIntegration::evalttbar_DL(const double* x ) {


  double mb=Physics::mb;

 
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@                   Take the input variable form VEGAS                 @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double  cosThetaNu_tlep = x[cosThetaNu_tlep_ttbar_DL_ttW_id]; // cosTheta of the nu from the leptonic top
  double  phiNu_tlep = x[phiNu_tlep_ttbar_DL_ttW_id];           // phi of the nu from the leptonic top
  double  cosThetaNu_ttau = x[cosThetaNu_ttau_ttbar_DL_ttW_id]; // cosTheta of the nu from the tauonic top
  double  phiNu_ttau = x[phiNu_ttau_ttbar_DL_ttW_id];           // phi of the nu from the tauonic top
  double P_TauHad = x[PTauHad_ttbar_DL_ttW_id];


  const char* error = 0;
  // To remove ???
//  char str[256];
  

  if (verbose_ >= IntegrandLevel){
    (*stream_) << "Vegas point :" 
	       << " " << " cosThetaNu_tlep = "    << cosThetaNu_tlep
	       << " " << " phiNu_tlep = "    << phiNu_tlep
	       << " " << " cosThetaNu_ttau = "    << cosThetaNu_ttau
	       << " " << " phiNu_ttau = "    << phiNu_ttau
	       << " " << " P_TauHad = "    << P_TauHad
	       << " " << endl;
  }



  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the leptonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_lep_reco_4P;
  TLorentzVector W_lep_4P;
  TLorentzVector nu_lep_4P;
  TLorentzVector t_lep_4P;

  if(error==0){
    
    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_tlep),phiNu_tlep);
    double CosThetalnu_top=TMath::Cos(evLep_top_4P_.Angle(enu));
    double Enu=getEqbar_Enu(CosThetalnu_top,evLep_top_4P_.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_lep_4P=nu_top_4P_dummy;
    W_lep_4P=nu_lep_4P+evLep_top_4P_;
    
    double Eb_lep=getEb(W_lep_4P, evBJet_leptop_4P_);
    double b_lep_Pmag=TMath::Sqrt(Eb_lep*Eb_lep - mb*mb);
    double pT_bl=b_lep_Pmag/TMath::CosH(evBJet_leptop_4P_.Eta());
    if(Eb_lep>0){
      b_lep_reco_4P.SetPtEtaPhiE(pT_bl,evBJet_leptop_4P_.Eta(),evBJet_leptop_4P_.Phi(),Eb_lep);
      t_lep_4P=W_lep_4P+b_lep_reco_4P;
    }
    else{
      error = "leptop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build leptonic top :" << endl;
      (*stream_) << " .... b from lep top " 
              <<  b_lep_reco_4P.Px() << ", "
              <<  b_lep_reco_4P.Py() << ", "
              <<  b_lep_reco_4P.Pz() << ", "
              <<  b_lep_reco_4P.E()     
              << endl;
     (*stream_) << " .... lep top" 
              <<  t_lep_4P.Px() << ", "
              <<  t_lep_4P.Py() << ", "
              <<  t_lep_4P.Pz() << ", "
              <<  t_lep_4P.E()    
              << endl;
    }

  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the tauonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_tau_reco_4P;
  TLorentzVector W_tau_4P;
  TLorentzVector nu_tau_4P;
  TLorentzVector tau_4P;
  TLorentzVector t_tau_4P;

  if(error==0){

    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_ttau),phiNu_ttau);
    double pT_tau = P_TauHad / TMath::CosH(evHadSys_Tau_4P_.Eta());
    double m_tau = Physics::mTau;
    tau_4P.SetPtEtaPhiM(pT_tau,evHadSys_Tau_4P_.Eta(),evHadSys_Tau_4P_.Phi(),m_tau);

    double CosThetataunu_top=TMath::Cos(tau_4P.Angle(enu));
    double Enu=getEqbar_Enu(CosThetataunu_top,tau_4P.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_tau_4P=nu_top_4P_dummy;
    W_tau_4P=nu_tau_4P+tau_4P;
    
    double Eb_tau=getEb(W_tau_4P, evBJet_hadtop_4P_);
    double b_tau_Pmag=TMath::Sqrt(Eb_tau*Eb_tau - mb*mb);
    double pT_btau=b_tau_Pmag/TMath::CosH(evBJet_hadtop_4P_.Eta());
    if(Eb_tau>0){
      b_tau_reco_4P.SetPtEtaPhiE(pT_btau,evBJet_hadtop_4P_.Eta(),evBJet_hadtop_4P_.Phi(),Eb_tau);
      t_tau_4P=W_tau_4P+b_tau_reco_4P;
    }
    else{
      error = "tautop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build tauonic top :" << endl;
      (*stream_) << " .... b from tau top " 
              <<  b_tau_reco_4P.Px() << ", "
              <<  b_tau_reco_4P.Py() << ", "
              <<  b_tau_reco_4P.Pz() << ", "
              <<  b_tau_reco_4P.E()     
              << endl;
      (*stream_) << " .... tau" 
              <<  tau_4P.Px() << ", "
              <<  tau_4P.Py() << ", "
              <<  tau_4P.Pz() << ", "
              <<  tau_4P.E()    
              << endl;
     (*stream_) << " .... tau top" 
              <<  t_tau_4P.Px() << ", "
              <<  t_tau_4P.Py() << ", "
              <<  t_tau_4P.Pz() << ", "
              <<  t_tau_4P.E()    
              << endl;
    }


  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                Recoil                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector tot4P;
  TLorentzVector rho4P;

  if(error==0){
    
    tot4P = t_lep_4P + t_tau_4P ;
    rho4P = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_;

  }

  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        JET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

  double T_BJet_leptop = 1.0;
  double T_BJet_tautop = 1.0;

  if ( error == 0 ) {

    if (flagTFBJet_leptop_)
      T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");      

    //
    if (flagTFBJet_hadtop_)
      T_BJet_tautop = getJetTF( b_tau_reco_4P, evBJet_hadtop_4P_,"b");


    // Log
/*    if (outgoingJetsOut_ttbar_DL_) {
      (*outgoingJetsOut_ttbar_DL_) << setiosflags(ios::fixed) << scientific;
      (*outgoingJetsOut_ttbar_DL_) 
	 << " " << setw(21) << evBJet_leptop_4P_.Pt()
         << " " << setw(21) << evBJet_leptop_4P_.Eta()
         << " " << setw(21) << evBJet_leptop_4P_.Phi()
         << " " << setw(21) << evBJet_leptop_4P_.P()
         << " " << setw(21) << evBJet_leptop_4P_.E()      
         << " " << setw(21) << b_lep_reco_4P.Pt()
         << " " << setw(21) << b_lep_reco_4P.Eta()
         << " " << setw(21) << b_lep_reco_4P.Phi()
         << " " << setw(21) << b_lep_reco_4P.P()
         << " " << setw(21) << b_lep_reco_4P.E()
         << " " << setw(21) << T_BJet_leptop
	 << " " << setw(21) << evBJet_hadtop_4P_.Pt()
         << " " << setw(21) << evBJet_hadtop_4P_.Eta()
         << " " << setw(21) << evBJet_hadtop_4P_.Phi()
         << " " << setw(21) << evBJet_hadtop_4P_.P()
         << " " << setw(21) << evBJet_hadtop_4P_.E()      
         << " " << setw(21) << b_tau_reco_4P.Pt()
         << " " << setw(21) << b_tau_reco_4P.Eta()
         << " " << setw(21) << b_tau_reco_4P.Phi()
         << " " << setw(21) << b_tau_reco_4P.P()
         << " " << setw(21) << b_tau_reco_4P.E()
         << " " << setw(21) << T_BJet_tautop
         << endl;
    }*/ 
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        MET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_MET=1.0;
  if ( error == 0 ) {
    if (flagTFMET_) {

      TVector3 RhoT;
      TVector3 PT;
      
      if(include_hadrecoil_){
	
	RhoT.SetXYZ((rho4P).Px(), (rho4P).Py(), 0);
	PT.SetXYZ((tot4P).Px(), (tot4P).Py(), 0);

      }

      else{
	
	TLorentzVector tot4P_nohad = W_lep_4P + W_tau_4P;
	TLorentzVector rho4P_nohad = evRecoMET4P_ + evLep_top_4P_ + evHadSys_Tau_4P_;
	RhoT.SetXYZ((rho4P_nohad).Px(), (rho4P_nohad).Py(), 0);
	PT.SetXYZ((tot4P_nohad).Px(), (tot4P_nohad).Py(), 0);


      }


      TVector3 Rec = RhoT-PT;
      // Order : (0,0), (0,1), (1,0), (1,1) ???
      double sqrtDetV = TMath::Sqrt( std::abs( evV_[0] * evV_[3] - evV_[1] * evV_[2]));
      double f = Rec.X()*(evV_[0]*Rec.X() + evV_[1]*Rec.Y()) 
               + Rec.Y()*(evV_[2]*Rec.X() + evV_[3]*Rec.Y());
      T_MET =  TMath::Exp( -0.5*f ) * sqrtDetV / ( 2 * TMath::Pi() ) ;
    

    }
  } 


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                    HADRONIC TAU Transfer function                    @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_hadTau = 1.0;
  if ( error == 0 ) {
    if ( flagTFHadTau_ )
      T_hadTau = getTauHadronicTF(P_TauHad, evHadSys_Tau_4P_);  
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                       TOP Transfer functions                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_lepTop = 1.0;
  double T_tauTop = 1.0;

  if ( error == 0 ) {
    if ( flagTFTop_ ){
      T_lepTop = getTopTF(t_lep_4P,b_lep_reco_4P,evLep_top_4P_,nu_lep_4P);
      T_tauTop = getTopTF(t_tau_4P,b_tau_reco_4P,tau_4P,nu_tau_4P);
    }
  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Jacobian Term                             @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
 
  double Jac = 1.0;
  
  if ( error == 0 ) {  
    if( flagJac_ ) {

      double Jac_tau = getJacobian_Tau( tau_4P ) * tau_4P.E() ; //Energy already taken into account in toptau Jacpobian 
      double Jac_toptau = getJacobian_toplep( b_tau_reco_4P, W_tau_4P, tau_4P.E(), nu_tau_4P.E());
      double Jac_toplep = getJacobian_toplep( b_lep_reco_4P, W_lep_4P, evLep_top_4P_.E(), nu_lep_4P.E());
      Jac = Jac_tau * Jac_toptau * Jac_toplep;

      if( TMath::IsNaN( Jac ) ) { if (error == 0) error = "Jacobian";}  
    }
  }




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                 BOOST                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double totPt = 0.0;
  TLorentzVector tot4P_boosted;
  TLorentzVector t_tau_4P_boosted;
  TLorentzVector t_lep_4P_boosted;

  if ( error == 0 ) {
    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Boosting :" << endl;
      (*stream_) << "  initial total momentum tot4P " 
              <<  tot4P.Px() << ", "
              <<  tot4P.Py() << ", "
              <<  tot4P.Pz() << ", "
              <<  tot4P.E() 
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
    totPt = tot4P.Pt();
    
    // Set transverse boost
    TVector3 boost3D = tot4P.BoostVector();
    boost3D.SetZ(0);
    
    //
    // Application of the boost to the final state 4-momenta
    // WARRNING : only following 4-momenta are now booosted
    tot4P_boosted=tot4P;
    tot4P_boosted.Boost   ( -boost3D );

    // ??? implulsion > energy 
    if( TMath::IsNaN(tot4P_boosted.E()) ) { if (error == 0) error = "Boost";}  

    t_tau_4P_boosted=t_tau_4P;
    t_lep_4P_boosted=t_lep_4P;

    t_tau_4P_boosted.Boost( -boost3D );
    t_lep_4P_boosted.Boost( -boost3D );

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "  Boosted total momentum tot4P " 
              <<  tot4P_boosted.Px() << ", "
              <<  tot4P_boosted.Py() << ", "
              <<  tot4P_boosted.Pz() << ", "
              <<  tot4P_boosted.E() 
              << endl;
      (*stream_) << "  Boost magnitude " 
              <<  totPt
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        THE BJORKEN FRACTION                          @@@
  //@@@                         Matrix Elements                              @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double wME = 1.0;
  double x_a;
  double x_b;
  x_a = (tot4P_boosted.E() + tot4P_boosted.Pz()) / sqrtS_ ;
  x_b = (tot4P_boosted.E() - tot4P_boosted.Pz()) / sqrtS_ ;
  
  if ( (x_a > 1.0) || (x_a < 0.0) || (x_b > 1.0 ) ||  (x_b < 0.0 ) ) 
    error = "BFrac";
  
  if ( error == 0 ) {

    TLorentzVector top_4P;
    TLorentzVector Atop_4P;

    if(lepton_top_charge_>0){
      top_4P=t_lep_4P_boosted;
      Atop_4P=t_tau_4P_boosted;
    }
    else{
      top_4P=t_tau_4P_boosted;
      Atop_4P=t_lep_4P_boosted;
    }


    // Parton momentum fraction
    //
    // MadGraph quadri-vector construction
    double inPart_a[4], inPart_b[4]; 
    double outtop[4], outAtop[4];
    double inE;
    // Incoming quarks
    inE =  0.5 * x_a * sqrtS_;
    inPart_a[0] = inE; inPart_a[1] = 0.; inPart_a[2] =  0; inPart_a[3] = inE; 
    inE =  0.5 * x_b * sqrtS_;
    inPart_b[0] = inE; inPart_b[1] = 0.; inPart_b[2] =  0; inPart_b[3] =-inE; 
    // top / Atop
    outtop[0] = top_4P.E(); outtop[1] = top_4P.Px(); outtop[2] = top_4P.Py(); outtop[3] = top_4P.Pz();
    outAtop[0] = Atop_4P.E(); outAtop[1] = Atop_4P.Px(); outAtop[2] = Atop_4P.Py(); outAtop[3] = Atop_4P.Pz();

    // Quadri-vector order for MadGraph
    //   p[0] incoming parton0
    //   p[1] incoming parton1
    //   p[2] outgoing top
    //   p[3] outgoing Atop
    //   p[4] outgoing fake
    vector<double*> p(5,(double*) 0);
    p[0] = inPart_a ; p[1] = inPart_b; 
    p[2] = outtop; p[3] = outAtop; 

    // Matrix Element
    if (flagWME_) {
      wME = get_ttbarWeightedME( x_a, x_b, p);
      if( TMath::IsNaN(wME) ) { if (error == 0) error = "WeightME";}
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Weight evaluation                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  double Eval =  T_hadTau * T_BJet_tautop * T_BJet_leptop * T_MET * T_lepTop * T_tauTop * Jac * wME;
  if (error) Eval = 0;

  if (verbose_ >= IntegrandLevel){
    // if( Eval > 10.0e-40) {
    (*stream_) << "Integrand value :" 
       << " T_hadTau = " << T_hadTau
       << " T_BJet_tautop = " << T_BJet_tautop
       << " T_BJet_leptop = " << T_BJet_leptop
       << " T_MET = " << T_MET 
       << " T_lepTop = " << T_lepTop
       << " T_tauTop = " << T_tauTop
       << " Jac = " << Jac
       << " wME = " << wME 
       << " Eval = "    << Eval
       << " Error= " << error
       << endl;
    // }
  }
/*  const char *err_str = "Ok";
  if (error) 
    err_str = error;
  else{*/
  if ( !error) {
    integr_Efficiencyttbar_DL_++;
  }


  tot_Drawsttbar_DL_++;

/*      
  if (valuesOut_ttbar_DL_) {
    writeFctValues_ttbar_DL(eventID_,			  
			    cosThetaNu_tlep, phiNu_tlep,
			    cosThetaNu_ttau, phiNu_ttau,
			    P_TauHad,
			    totPt,
			    rho4P, tot4P,
			    T_BJet_tautop, T_BJet_leptop,
			    T_MET, T_hadTau,
			    T_lepTop, T_tauTop,
			    Jac,  wME, Eval, err_str);
  }

  if (valuesOutRoot_ttbar_DL_) {
    valuesOutRoot_ttbar_DL_->fill(eventID_,		       
				  cosThetaNu_tlep, phiNu_tlep,
				  cosThetaNu_ttau, phiNu_ttau,
				  P_TauHad,
				  totPt,
				  rho4P, tot4P,
				  T_BJet_tautop, T_BJet_leptop,
				  T_MET, T_hadTau,
				  T_lepTop, T_tauTop,
				  Jac,  wME, Eval, err_str);
    
  }
*/

  return Eval;
}

/*
void MGIntegration::writeFctValues_ttbar_DL(int64_t eventID,
					    double cosThetaNu_tlep, double phiNu_tlep,
					    double cosThetaNu_ttau, double phiNu_ttau,
					    double PTauHad,
					     double boost, 
					    TLorentzVector rho4P, TLorentzVector boost4P, 
					    double TFBjet_tautop, double TFBjet_leptop,
					    double TFMET, double TFHadTau,
					    double TFLepTop, double TFTauTop,
					    double Jac, double wME, 
					    double eval,
					    const char *errStr) const {
  if (valuesOut_ttbar_DL_) {
    (*valuesOut_ttbar_DL_) << setiosflags(ios::fixed) << scientific;
    (*valuesOut_ttbar_DL_) 
      << " " << setw(9) << eventID
      << " " << setw(15) << cosThetaNu_tlep
      << " " << setw(15) << phiNu_tlep    
      << " " << setw(15) << cosThetaNu_ttau
      << " " << setw(15) << phiNu_ttau
      << " " << setw(15) << PTauHad
      << " " << setw(15) << rho4P.Px()
      << " " << setw(15) << rho4P.Py()
      << " " << setw(15) << rho4P.Pz()
      << " " << setw(15) << rho4P.E()
      << " " << setw(15) << boost4P.Px()
      << " " << setw(15) << boost4P.Py()
      << " " << setw(15) << boost4P.Pz()
      << " " << setw(15) << boost4P.E()
      << " " << setw(15) << boost
      << " " << setw(15) << TFBjet_tautop
      << " " << setw(15) << TFBjet_leptop
      << " " << setw(15) << TFMET
      << " " << setw(15) << TFHadTau
      << " " << setw(15) << TFLepTop
      << " " << setw(15) << TFTauTop
      << " " << setw(15) << Jac
      << " " << setw(15) << wME 
      << " " << setw(15) << eval
      << " " << setw(15) << errStr
      << endl;
  }  
}
*/

double MGIntegration::evalttbar_DL_fakelep_tlep(const double* x ) {


  double mb=Physics::mb;

 
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@                   Take the input variable form VEGAS                 @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double  cosThetaNu_tlep = x[cosThetaNu_tlep_ttbar_DL_ttW_id]; // cosTheta of the nu from the leptonic top
  double  phiNu_tlep = x[phiNu_tlep_ttbar_DL_ttW_id];           // phi of the nu from the leptonic top
  double  cosThetaNu_ttau = x[cosThetaNu_ttau_ttbar_DL_ttW_id]; // cosTheta of the nu from the tauonic top
  double  phiNu_ttau = x[phiNu_ttau_ttbar_DL_ttW_id];           // phi of the nu from the tauonic top
  double P_TauHad = x[PTauHad_ttbar_DL_ttW_id];


  const char* error = 0;
  // To remove ???
//  char str[256];
  

  if (verbose_ >= IntegrandLevel){
    (*stream_) << "Vegas point :" 
	       << " " << " cosThetaNu_tlep = "    << cosThetaNu_tlep
	       << " " << " phiNu_tlep = "    << phiNu_tlep
	       << " " << " cosThetaNu_ttau = "    << cosThetaNu_ttau
	       << " " << " phiNu_ttau = "    << phiNu_ttau
	       << " " << " P_TauHad = "    << P_TauHad
	       << " " << endl;
  }



  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the leptonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_lep_reco_4P;
  TLorentzVector W_lep_4P;
  TLorentzVector nu_lep_4P;
  TLorentzVector t_lep_4P;

  if(error==0){
    
    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_tlep),phiNu_tlep);
    double CosThetalnu_top=TMath::Cos(evLep_top_4P_.Angle(enu));
    double Enu=getEqbar_Enu(CosThetalnu_top,evLep_top_4P_.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_lep_4P=nu_top_4P_dummy;
    W_lep_4P=nu_lep_4P+evLep_top_4P_;
    
    // b quark is reconstructed as fake lepton
    double Eb_lep=getEb(W_lep_4P, evLep_Tau_4P_);
    double b_lep_Pmag=TMath::Sqrt(Eb_lep*Eb_lep - mb*mb);
    double pT_bl=b_lep_Pmag/TMath::CosH(evLep_Tau_4P_.Eta());
    if(Eb_lep>0){
      b_lep_reco_4P.SetPtEtaPhiE(pT_bl,evLep_Tau_4P_.Eta(),evLep_Tau_4P_.Phi(),Eb_lep);
      t_lep_4P=W_lep_4P+b_lep_reco_4P;
    }
    else{
      error = "leptop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build leptonic top :" << endl;
      (*stream_) << " .... b from lep top " 
              <<  b_lep_reco_4P.Px() << ", "
              <<  b_lep_reco_4P.Py() << ", "
              <<  b_lep_reco_4P.Pz() << ", "
              <<  b_lep_reco_4P.E()     
              << endl;
     (*stream_) << " .... lep top" 
              <<  t_lep_4P.Px() << ", "
              <<  t_lep_4P.Py() << ", "
              <<  t_lep_4P.Pz() << ", "
              <<  t_lep_4P.E()    
              << endl;
    }

  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the tauonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_tau_reco_4P;
  TLorentzVector W_tau_4P;
  TLorentzVector nu_tau_4P;
  TLorentzVector tau_4P;
  TLorentzVector t_tau_4P;

  if(error==0){

    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_ttau),phiNu_ttau);
    double pT_tau = P_TauHad / TMath::CosH(evHadSys_Tau_4P_.Eta());
    double m_tau = Physics::mTau;
    tau_4P.SetPtEtaPhiM(pT_tau,evHadSys_Tau_4P_.Eta(),evHadSys_Tau_4P_.Phi(),m_tau);

    double CosThetataunu_top=TMath::Cos(tau_4P.Angle(enu));
    double Enu=getEqbar_Enu(CosThetataunu_top,tau_4P.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_tau_4P=nu_top_4P_dummy;
    W_tau_4P=nu_tau_4P+tau_4P;
    
    double Eb_tau=getEb(W_tau_4P, evBJet_hadtop_4P_);
    double b_tau_Pmag=TMath::Sqrt(Eb_tau*Eb_tau - mb*mb);
    double pT_btau=b_tau_Pmag/TMath::CosH(evBJet_hadtop_4P_.Eta());
    if(Eb_tau>0){
      b_tau_reco_4P.SetPtEtaPhiE(pT_btau,evBJet_hadtop_4P_.Eta(),evBJet_hadtop_4P_.Phi(),Eb_tau);
      t_tau_4P=W_tau_4P+b_tau_reco_4P;
    }
    else{
      error = "tautop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build tauonic top :" << endl;
      (*stream_) << " .... b from tau top " 
              <<  b_tau_reco_4P.Px() << ", "
              <<  b_tau_reco_4P.Py() << ", "
              <<  b_tau_reco_4P.Pz() << ", "
              <<  b_tau_reco_4P.E()     
              << endl;
      (*stream_) << " .... tau" 
              <<  tau_4P.Px() << ", "
              <<  tau_4P.Py() << ", "
              <<  tau_4P.Pz() << ", "
              <<  tau_4P.E()    
              << endl;
     (*stream_) << " .... tau top" 
              <<  t_tau_4P.Px() << ", "
              <<  t_tau_4P.Py() << ", "
              <<  t_tau_4P.Pz() << ", "
              <<  t_tau_4P.E()    
              << endl;
    }


  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                Recoil                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector tot4P;
  TLorentzVector rho4P;

  if(error==0){
    
    tot4P = t_lep_4P + t_tau_4P ;
    rho4P = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_top_4P_ + evLep_Tau_4P_ + evBJet_hadtop_4P_;

  }

  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        JET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

  double T_fakelep_leptop = 1.0;
  double T_BJet_tautop = 1.0;

  if ( error == 0 ) {

    if (flagTF_fakelep_)
      T_fakelep_leptop = getFakeLepTF( b_lep_reco_4P, evLep_Tau_4P_);      

    //
    if (flagTFBJet_hadtop_)
      T_BJet_tautop = getJetTF( b_tau_reco_4P, evBJet_hadtop_4P_,"b");

    // Log
 /*   if (outgoingJetsOut_ttbar_DL_) {
      (*outgoingJetsOut_ttbar_DL_) << setiosflags(ios::fixed) << scientific;
      (*outgoingJetsOut_ttbar_DL_) 
	 << " " << setw(21) << evLep_Tau_4P_.Pt()
         << " " << setw(21) << evLep_Tau_4P_.Eta()
         << " " << setw(21) << evLep_Tau_4P_.Phi()
         << " " << setw(21) << evLep_Tau_4P_.P()
         << " " << setw(21) << evLep_Tau_4P_.E()      
         << " " << setw(21) << b_lep_reco_4P.Pt()
         << " " << setw(21) << b_lep_reco_4P.Eta()
         << " " << setw(21) << b_lep_reco_4P.Phi()
         << " " << setw(21) << b_lep_reco_4P.P()
         << " " << setw(21) << b_lep_reco_4P.E()
         << " " << setw(21) << T_fakelep_leptop
	 << " " << setw(21) << evBJet_hadtop_4P_.Pt()
         << " " << setw(21) << evBJet_hadtop_4P_.Eta()
         << " " << setw(21) << evBJet_hadtop_4P_.Phi()
         << " " << setw(21) << evBJet_hadtop_4P_.P()
         << " " << setw(21) << evBJet_hadtop_4P_.E()      
         << " " << setw(21) << b_tau_reco_4P.Pt()
         << " " << setw(21) << b_tau_reco_4P.Eta()
         << " " << setw(21) << b_tau_reco_4P.Phi()
         << " " << setw(21) << b_tau_reco_4P.P()
         << " " << setw(21) << b_tau_reco_4P.E()
         << " " << setw(21) << T_BJet_tautop
         << endl;
    } */
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        MET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_MET=1.0;
  if ( error == 0 ) {
    if (flagTFMET_) {

      TVector3 RhoT;
      TVector3 PT;
      
      if(include_hadrecoil_){
	
	RhoT.SetXYZ((rho4P).Px(), (rho4P).Py(), 0);
	PT.SetXYZ((tot4P).Px(), (tot4P).Py(), 0);

      }

      else{
	
	TLorentzVector tot4P_nohad = W_lep_4P + W_tau_4P;
	TLorentzVector rho4P_nohad = evRecoMET4P_ + evLep_top_4P_ + evHadSys_Tau_4P_;
	RhoT.SetXYZ((rho4P_nohad).Px(), (rho4P_nohad).Py(), 0);
	PT.SetXYZ((tot4P_nohad).Px(), (tot4P_nohad).Py(), 0);


      }


      TVector3 Rec = RhoT-PT;
      // Order : (0,0), (0,1), (1,0), (1,1) ???
      double sqrtDetV = TMath::Sqrt( std::abs( evV_[0] * evV_[3] - evV_[1] * evV_[2]));
      double f = Rec.X()*(evV_[0]*Rec.X() + evV_[1]*Rec.Y()) 
               + Rec.Y()*(evV_[2]*Rec.X() + evV_[3]*Rec.Y());
      T_MET =  TMath::Exp( -0.5*f ) * sqrtDetV / ( 2 * TMath::Pi() ) ;
    

    }
  } 


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                    HADRONIC TAU Transfer function                    @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_hadTau = 1.0;
  if ( error == 0 ) {
    if ( flagTFHadTau_ )
      T_hadTau = getTauHadronicTF(P_TauHad, evHadSys_Tau_4P_);  
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                       TOP Transfer functions                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_lepTop = 1.0;
  double T_tauTop = 1.0;

  if ( error == 0 ) {
    if ( flagTFTop_ ){
      T_lepTop = getTopTF(t_lep_4P,b_lep_reco_4P,evLep_top_4P_,nu_lep_4P);
      T_tauTop = getTopTF(t_tau_4P,b_tau_reco_4P,tau_4P,nu_tau_4P);
    }
  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Jacobian Term                             @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
 
  double Jac = 1.0;
  
  if ( error == 0 ) {  
    if( flagJac_ ) {

      double Jac_tau = getJacobian_Tau( tau_4P ) * tau_4P.E() ; //Energy already taken into account in toptau Jacpobian 
      double Jac_toptau = getJacobian_toplep( b_tau_reco_4P, W_tau_4P, tau_4P.E(), nu_tau_4P.E());
      double Jac_toplep = getJacobian_toplep( b_lep_reco_4P, W_lep_4P, evLep_top_4P_.E(), nu_lep_4P.E());
      Jac = Jac_tau * Jac_toptau * Jac_toplep;

      if( TMath::IsNaN( Jac ) ) { if (error == 0) error = "Jacobian";}  
    }
  }




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                 BOOST                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double totPt = 0.0;
  TLorentzVector tot4P_boosted;
  TLorentzVector t_tau_4P_boosted;
  TLorentzVector t_lep_4P_boosted;

  if ( error == 0 ) {
    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Boosting :" << endl;
      (*stream_) << "  initial total momentum tot4P " 
              <<  tot4P.Px() << ", "
              <<  tot4P.Py() << ", "
              <<  tot4P.Pz() << ", "
              <<  tot4P.E() 
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
    totPt = tot4P.Pt();
    
    // Set transverse boost
    TVector3 boost3D = tot4P.BoostVector();
    boost3D.SetZ(0);
    
    //
    // Application of the boost to the final state 4-momenta
    // WARRNING : only following 4-momenta are now booosted
    tot4P_boosted=tot4P;
    tot4P_boosted.Boost   ( -boost3D );

    // ??? implulsion > energy 
    if( TMath::IsNaN(tot4P_boosted.E()) ) { if (error == 0) error = "Boost";}  

    t_tau_4P_boosted=t_tau_4P;
    t_lep_4P_boosted=t_lep_4P;

    t_tau_4P_boosted.Boost( -boost3D );
    t_lep_4P_boosted.Boost( -boost3D );

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "  Boosted total momentum tot4P " 
              <<  tot4P_boosted.Px() << ", "
              <<  tot4P_boosted.Py() << ", "
              <<  tot4P_boosted.Pz() << ", "
              <<  tot4P_boosted.E() 
              << endl;
      (*stream_) << "  Boost magnitude " 
              <<  totPt
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        THE BJORKEN FRACTION                          @@@
  //@@@                         Matrix Elements                              @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double wME = 1.0;
  double x_a;
  double x_b;
  x_a = (tot4P_boosted.E() + tot4P_boosted.Pz()) / sqrtS_ ;
  x_b = (tot4P_boosted.E() - tot4P_boosted.Pz()) / sqrtS_ ;
  
  if ( (x_a > 1.0) || (x_a < 0.0) || (x_b > 1.0 ) ||  (x_b < 0.0 ) ) 
    error = "BFrac";
  
  if ( error == 0 ) {

    TLorentzVector top_4P;
    TLorentzVector Atop_4P;

    if(lepton_top_charge_>0){
      top_4P=t_lep_4P_boosted;
      Atop_4P=t_tau_4P_boosted;
    }
    else{
      top_4P=t_tau_4P_boosted;
      Atop_4P=t_lep_4P_boosted;
    }


    // Parton momentum fraction
    //
    // MadGraph quadri-vector construction
    double inPart_a[4], inPart_b[4]; 
    double outtop[4], outAtop[4];
    double inE;
    // Incoming quarks
    inE =  0.5 * x_a * sqrtS_;
    inPart_a[0] = inE; inPart_a[1] = 0.; inPart_a[2] =  0; inPart_a[3] = inE; 
    inE =  0.5 * x_b * sqrtS_;
    inPart_b[0] = inE; inPart_b[1] = 0.; inPart_b[2] =  0; inPart_b[3] =-inE; 
    // top / Atop
    outtop[0] = top_4P.E(); outtop[1] = top_4P.Px(); outtop[2] = top_4P.Py(); outtop[3] = top_4P.Pz();
    outAtop[0] = Atop_4P.E(); outAtop[1] = Atop_4P.Px(); outAtop[2] = Atop_4P.Py(); outAtop[3] = Atop_4P.Pz();

    // Quadri-vector order for MadGraph
    //   p[0] incoming parton0
    //   p[1] incoming parton1
    //   p[2] outgoing top
    //   p[3] outgoing Atop
    //   p[4] outgoing fake
    vector<double*> p(5,(double*) 0);
    p[0] = inPart_a ; p[1] = inPart_b; 
    p[2] = outtop; p[3] = outAtop; 

    // Matrix Element
    if (flagWME_) {
      wME = get_ttbarWeightedME( x_a, x_b, p);
      if( TMath::IsNaN(wME) ) { if (error == 0) error = "WeightME";}
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Weight evaluation                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  double Eval =  T_hadTau * T_BJet_tautop * T_fakelep_leptop * T_MET * T_lepTop * T_tauTop * Jac * wME;
  if (error) Eval = 0;

  if (verbose_ >= IntegrandLevel){
    // if( Eval > 10.0e-40) {
    (*stream_) << "Integrand value :" 
       << " T_hadTau = " << T_hadTau
       << " T_BJet_tautop = " << T_BJet_tautop
       << " T_fakelep_leptop = " << T_fakelep_leptop
       << " T_MET = " << T_MET 
       << " T_lepTop = " << T_lepTop
       << " T_tauTop = " << T_tauTop
       << " Jac = " << Jac
       << " wME = " << wME 
       << " Eval = "    << Eval
       << " Error= " << error
       << endl;
    // }
  }
/*  const char *err_str = "Ok";
  if (error) 
    err_str = error;
  else{*/
  if ( !error) {
    integr_Efficiencyttbar_DL_fakelep_tlep_++;
  }


  tot_Drawsttbar_DL_fakelep_tlep_++;


  return Eval;
}

double MGIntegration::evalttbar_DL_fakelep_ttau(const double* x ) {


  double mb=Physics::mb;

 
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@                   Take the input variable form VEGAS                 @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double  cosThetaNu_tlep = x[cosThetaNu_tlep_ttbar_DL_ttW_id]; // cosTheta of the nu from the leptonic top
  double  phiNu_tlep = x[phiNu_tlep_ttbar_DL_ttW_id];           // phi of the nu from the leptonic top
  double  cosThetaNu_ttau = x[cosThetaNu_ttau_ttbar_DL_ttW_id]; // cosTheta of the nu from the tauonic top
  double  phiNu_ttau = x[phiNu_ttau_ttbar_DL_ttW_id];           // phi of the nu from the tauonic top
  double P_TauHad = x[PTauHad_ttbar_DL_ttW_id];


  const char* error = 0;
  // To remove ???
//  char str[256];
  

  if (verbose_ >= IntegrandLevel){
    (*stream_) << "Vegas point :" 
	       << " " << " cosThetaNu_tlep = "    << cosThetaNu_tlep
	       << " " << " phiNu_tlep = "    << phiNu_tlep
	       << " " << " cosThetaNu_ttau = "    << cosThetaNu_ttau
	       << " " << " phiNu_ttau = "    << phiNu_ttau
	       << " " << " P_TauHad = "    << P_TauHad
	       << " " << endl;
  }



  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the leptonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_lep_reco_4P;
  TLorentzVector W_lep_4P;
  TLorentzVector nu_lep_4P;
  TLorentzVector t_lep_4P;

  if(error==0){
    
    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_tlep),phiNu_tlep);
    double CosThetalnu_top=TMath::Cos(evLep_top_4P_.Angle(enu));
    double Enu=getEqbar_Enu(CosThetalnu_top,evLep_top_4P_.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_lep_4P=nu_top_4P_dummy;
    W_lep_4P=nu_lep_4P+evLep_top_4P_;
    
    double Eb_lep=getEb(W_lep_4P, evBJet_leptop_4P_);
    double b_lep_Pmag=TMath::Sqrt(Eb_lep*Eb_lep - mb*mb);
    double pT_bl=b_lep_Pmag/TMath::CosH(evBJet_leptop_4P_.Eta());
    if(Eb_lep>0){
      b_lep_reco_4P.SetPtEtaPhiE(pT_bl,evBJet_leptop_4P_.Eta(),evBJet_leptop_4P_.Phi(),Eb_lep);
      t_lep_4P=W_lep_4P+b_lep_reco_4P;
    }
    else{
      error = "leptop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build leptonic top :" << endl;
      (*stream_) << " .... b from lep top " 
              <<  b_lep_reco_4P.Px() << ", "
              <<  b_lep_reco_4P.Py() << ", "
              <<  b_lep_reco_4P.Pz() << ", "
              <<  b_lep_reco_4P.E()     
              << endl;
     (*stream_) << " .... lep top" 
              <<  t_lep_4P.Px() << ", "
              <<  t_lep_4P.Py() << ", "
              <<  t_lep_4P.Pz() << ", "
              <<  t_lep_4P.E()    
              << endl;
    }

  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the tauonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_tau_reco_4P;
  TLorentzVector W_tau_4P;
  TLorentzVector nu_tau_4P;
  TLorentzVector tau_4P;
  TLorentzVector t_tau_4P;

  if(error==0){

    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu_ttau),phiNu_ttau);
    double pT_tau = P_TauHad / TMath::CosH(evHadSys_Tau_4P_.Eta());
    double m_tau = Physics::mTau;
    tau_4P.SetPtEtaPhiM(pT_tau,evHadSys_Tau_4P_.Eta(),evHadSys_Tau_4P_.Phi(),m_tau);

    double CosThetataunu_top=TMath::Cos(tau_4P.Angle(enu));
    double Enu=getEqbar_Enu(CosThetataunu_top,tau_4P.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_tau_4P=nu_top_4P_dummy;
    W_tau_4P=nu_tau_4P+tau_4P;
    
    // b quark is reconstructed as fake lepton
    double Eb_tau=getEb(W_tau_4P, evLep_Tau_4P_);
    double b_tau_Pmag=TMath::Sqrt(Eb_tau*Eb_tau - mb*mb);
    double pT_btau=b_tau_Pmag/TMath::CosH(evLep_Tau_4P_.Eta());
    if(Eb_tau>0){
      b_tau_reco_4P.SetPtEtaPhiE(pT_btau,evLep_Tau_4P_.Eta(),evLep_Tau_4P_.Phi(),Eb_tau);
      t_tau_4P=W_tau_4P+b_tau_reco_4P;
    }
    else{
      error = "tautop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build tauonic top :" << endl;
      (*stream_) << " .... b from tau top " 
              <<  b_tau_reco_4P.Px() << ", "
              <<  b_tau_reco_4P.Py() << ", "
              <<  b_tau_reco_4P.Pz() << ", "
              <<  b_tau_reco_4P.E()     
              << endl;
      (*stream_) << " .... tau" 
              <<  tau_4P.Px() << ", "
              <<  tau_4P.Py() << ", "
              <<  tau_4P.Pz() << ", "
              <<  tau_4P.E()    
              << endl;
     (*stream_) << " .... tau top" 
              <<  t_tau_4P.Px() << ", "
              <<  t_tau_4P.Py() << ", "
              <<  t_tau_4P.Pz() << ", "
              <<  t_tau_4P.E()    
              << endl;
    }


  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                Recoil                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector tot4P;
  TLorentzVector rho4P;

  if(error==0){
    
    tot4P = t_lep_4P + t_tau_4P ;
    rho4P = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_top_4P_ + evLep_Tau_4P_ + evBJet_leptop_4P_;

  }

  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        JET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

  double T_BJet_leptop = 1.0;
  double T_fakelep_tautop = 1.0;

  if ( error == 0 ) {

    if (flagTFBJet_hadtop_)
      T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");

    //
    if (flagTF_fakelep_)
      T_fakelep_tautop = getFakeLepTF( b_tau_reco_4P, evLep_Tau_4P_);


    // Log
/*    if (outgoingJetsOut_ttbar_DL_) {
      (*outgoingJetsOut_ttbar_DL_) << setiosflags(ios::fixed) << scientific;
      (*outgoingJetsOut_ttbar_DL_) 
	 << " " << setw(21) << evBJet_leptop_4P_.Pt()
         << " " << setw(21) << evBJet_leptop_4P_.Eta()
         << " " << setw(21) << evBJet_leptop_4P_.Phi()
         << " " << setw(21) << evBJet_leptop_4P_.P()
         << " " << setw(21) << evBJet_leptop_4P_.E()	       
         << " " << setw(21) << b_lep_reco_4P.Pt()
         << " " << setw(21) << b_lep_reco_4P.Eta()
         << " " << setw(21) << b_lep_reco_4P.Phi()
         << " " << setw(21) << b_lep_reco_4P.P()
         << " " << setw(21) << b_lep_reco_4P.E()
         << " " << setw(21) << T_BJet_leptop
	 << " " << setw(21) << evLep_Tau_4P_.Pt()
         << " " << setw(21) << evLep_Tau_4P_.Eta()
         << " " << setw(21) << evLep_Tau_4P_.Phi()
         << " " << setw(21) << evLep_Tau_4P_.P()
         << " " << setw(21) << evLep_Tau_4P_.E()	      
         << " " << setw(21) << b_tau_reco_4P.Pt()
         << " " << setw(21) << b_tau_reco_4P.Eta()
         << " " << setw(21) << b_tau_reco_4P.Phi()
         << " " << setw(21) << b_tau_reco_4P.P()
         << " " << setw(21) << b_tau_reco_4P.E()
         << " " << setw(21) << T_fakelep_tautop
         << endl;
    } */
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        MET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_MET=1.0;
  if ( error == 0 ) {
    if (flagTFMET_) {

      TVector3 RhoT;
      TVector3 PT;
      
      if(include_hadrecoil_){
	
	RhoT.SetXYZ((rho4P).Px(), (rho4P).Py(), 0);
	PT.SetXYZ((tot4P).Px(), (tot4P).Py(), 0);

      }

      else{
	
	TLorentzVector tot4P_nohad = W_lep_4P + W_tau_4P;
	TLorentzVector rho4P_nohad = evRecoMET4P_ + evLep_top_4P_ + evHadSys_Tau_4P_;
	RhoT.SetXYZ((rho4P_nohad).Px(), (rho4P_nohad).Py(), 0);
	PT.SetXYZ((tot4P_nohad).Px(), (tot4P_nohad).Py(), 0);


      }


      TVector3 Rec = RhoT-PT;
      // Order : (0,0), (0,1), (1,0), (1,1) ???
      double sqrtDetV = TMath::Sqrt( std::abs( evV_[0] * evV_[3] - evV_[1] * evV_[2]));
      double f = Rec.X()*(evV_[0]*Rec.X() + evV_[1]*Rec.Y()) 
               + Rec.Y()*(evV_[2]*Rec.X() + evV_[3]*Rec.Y());
      T_MET =  TMath::Exp( -0.5*f ) * sqrtDetV / ( 2 * TMath::Pi() ) ;
    

    }
  } 


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                    HADRONIC TAU Transfer function                    @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_hadTau = 1.0;
  if ( error == 0 ) {
    if ( flagTFHadTau_ )
      T_hadTau = getTauHadronicTF(P_TauHad, evHadSys_Tau_4P_);  
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                       TOP Transfer functions                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_lepTop = 1.0;
  double T_tauTop = 1.0;

  if ( error == 0 ) {
    if ( flagTFTop_ ){
      T_lepTop = getTopTF(t_lep_4P,b_lep_reco_4P,evLep_top_4P_,nu_lep_4P);
      T_tauTop = getTopTF(t_tau_4P,b_tau_reco_4P,tau_4P,nu_tau_4P);
    }
  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Jacobian Term                             @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
 
  double Jac = 1.0;
  
  if ( error == 0 ) {  
    if( flagJac_ ) {

      double Jac_tau = getJacobian_Tau( tau_4P ) * tau_4P.E() ; //Energy already taken into account in toptau Jacpobian 
      double Jac_toptau = getJacobian_toplep( b_tau_reco_4P, W_tau_4P, tau_4P.E(), nu_tau_4P.E());
      double Jac_toplep = getJacobian_toplep( b_lep_reco_4P, W_lep_4P, evLep_top_4P_.E(), nu_lep_4P.E());
      Jac = Jac_tau * Jac_toptau * Jac_toplep;

      if( TMath::IsNaN( Jac ) ) { if (error == 0) error = "Jacobian";}  
    }
  }




  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                 BOOST                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double totPt = 0.0;
  TLorentzVector tot4P_boosted;
  TLorentzVector t_tau_4P_boosted;
  TLorentzVector t_lep_4P_boosted;

  if ( error == 0 ) {
    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Boosting :" << endl;
      (*stream_) << "  initial total momentum tot4P " 
              <<  tot4P.Px() << ", "
              <<  tot4P.Py() << ", "
              <<  tot4P.Pz() << ", "
              <<  tot4P.E() 
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
    totPt = tot4P.Pt();
    
    // Set transverse boost
    TVector3 boost3D = tot4P.BoostVector();
    boost3D.SetZ(0);
    
    //
    // Application of the boost to the final state 4-momenta
    // WARRNING : only following 4-momenta are now booosted
    tot4P_boosted=tot4P;
    tot4P_boosted.Boost   ( -boost3D );

    // ??? implulsion > energy 
    if( TMath::IsNaN(tot4P_boosted.E()) ) { if (error == 0) error = "Boost";}  

    t_tau_4P_boosted=t_tau_4P;
    t_lep_4P_boosted=t_lep_4P;

    t_tau_4P_boosted.Boost( -boost3D );
    t_lep_4P_boosted.Boost( -boost3D );

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "  Boosted total momentum tot4P " 
              <<  tot4P_boosted.Px() << ", "
              <<  tot4P_boosted.Py() << ", "
              <<  tot4P_boosted.Pz() << ", "
              <<  tot4P_boosted.E() 
              << endl;
      (*stream_) << "  Boost magnitude " 
              <<  totPt
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        THE BJORKEN FRACTION                          @@@
  //@@@                         Matrix Elements                              @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double wME = 1.0;
  double x_a;
  double x_b;
  x_a = (tot4P_boosted.E() + tot4P_boosted.Pz()) / sqrtS_ ;
  x_b = (tot4P_boosted.E() - tot4P_boosted.Pz()) / sqrtS_ ;
  
  if ( (x_a > 1.0) || (x_a < 0.0) || (x_b > 1.0 ) ||  (x_b < 0.0 ) ) 
    error = "BFrac";
  
  if ( error == 0 ) {

    TLorentzVector top_4P;
    TLorentzVector Atop_4P;

    if(lepton_top_charge_>0){
      top_4P=t_lep_4P_boosted;
      Atop_4P=t_tau_4P_boosted;
    }
    else{
      top_4P=t_tau_4P_boosted;
      Atop_4P=t_lep_4P_boosted;
    }


    // Parton momentum fraction
    //
    // MadGraph quadri-vector construction
    double inPart_a[4], inPart_b[4]; 
    double outtop[4], outAtop[4];
    double inE;
    // Incoming quarks
    inE =  0.5 * x_a * sqrtS_;
    inPart_a[0] = inE; inPart_a[1] = 0.; inPart_a[2] =  0; inPart_a[3] = inE; 
    inE =  0.5 * x_b * sqrtS_;
    inPart_b[0] = inE; inPart_b[1] = 0.; inPart_b[2] =  0; inPart_b[3] =-inE; 
    // top / Atop
    outtop[0] = top_4P.E(); outtop[1] = top_4P.Px(); outtop[2] = top_4P.Py(); outtop[3] = top_4P.Pz();
    outAtop[0] = Atop_4P.E(); outAtop[1] = Atop_4P.Px(); outAtop[2] = Atop_4P.Py(); outAtop[3] = Atop_4P.Pz();

    // Quadri-vector order for MadGraph
    //   p[0] incoming parton0
    //   p[1] incoming parton1
    //   p[2] outgoing top
    //   p[3] outgoing Atop
    //   p[4] outgoing fake
    vector<double*> p(5,(double*) 0);
    p[0] = inPart_a ; p[1] = inPart_b; 
    p[2] = outtop; p[3] = outAtop; 

    // Matrix Element
    if (flagWME_) {
      wME = get_ttbarWeightedME( x_a, x_b, p);
      if( TMath::IsNaN(wME) ) { if (error == 0) error = "WeightME";}
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Weight evaluation                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  double Eval =  T_hadTau * T_fakelep_tautop * T_BJet_leptop * T_MET * T_lepTop * T_tauTop * Jac * wME;
  if (error) Eval = 0;

  if (verbose_ >= IntegrandLevel){
    // if( Eval > 10.0e-40) {
    (*stream_) << "Integrand value :" 
       << " T_hadTau = " << T_hadTau
       << " T_fakelep_tautop = " << T_fakelep_tautop
       << " T_BJet_leptop = " << T_BJet_leptop
       << " T_MET = " << T_MET 
       << " T_lepTop = " << T_lepTop
       << " T_tauTop = " << T_tauTop
       << " Jac = " << Jac
       << " wME = " << wME 
       << " Eval = "    << Eval
       << " Error= " << error
       << endl;
    // }
  }
/*  const char *err_str = "Ok";
  if (error) 
    err_str = error;
  else{*/
  if ( !error) {
    integr_Efficiencyttbar_DL_fakelep_ttau_++;
  }


  tot_Drawsttbar_DL_fakelep_ttau_++;

  return Eval;
}

double MGIntegration::evalttZ_Zll(const double* x ) {

  double mb=Physics::mb;

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  //@@@                   Take the input variable form VEGAS                 @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double  EQuark1	  = x[EQuark1_ttbar_SL_ttZ_Zll_id];     //module of the momentum of the final state Quark1
  double  cosThetaNu = x[cosThetaNu_tlep_ttbar_SL_ttZ_Zll_id]; // cosTheta of the nu from the leptonic top
  double  phiNu = x[phiNu_tlep_ttbar_SL_ttZ_Zll_id]; // phi of the nu from the leptonic top  

  double cosTheta_missing_jet = 0;
  double phi_missing_jet      = 0;
  if(integration_type_ == integration_type_w_miss){
    cosTheta_missing_jet = x[cosTheta_missing_jet_ttbar_SL_ttZ_Zll_id];  // cosTheta of the quark associated to the missing jet
    phi_missing_jet      = x[phi_missing_jet_ttbar_SL_ttZ_Zll_id];   // phi of the quark associated to the missing jet
  }

  const char* error = 0;
  // To remove ???
//  char str[256];

  if (verbose_ >= IntegrandLevel){
    (*stream_) << "Vegas point :" 
	       << " " << " cosTheta_missing_jet = "    << cosTheta_missing_jet
	       << " " << " phi_missing_jet = "    << phi_missing_jet
	       << " " << " EQuark1 = "    << EQuark1
	       << " " << " cosThetaNu = "    << cosThetaNu
	       << " " << " phiNu = "    << phiNu
	       << " " << endl;
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@              Evaluation of the fake lepton -> tau_h 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector lep_Z_faketau_4P;

  double MZ = Physics::mZ;
  double CosTheta_ll = TMath::Cos(evLep_Tau_4P_.Angle(evHadSys_Tau_4P_.Vect()));
  double E_fake = MZ*MZ / (2*evLep_Tau_4P_.E()*(1-CosTheta_ll));  
  double pt_fake=E_fake/TMath::CosH(evHadSys_Tau_4P_.Eta());

  lep_Z_faketau_4P.SetPtEtaPhiE(pt_fake,evHadSys_Tau_4P_.Eta(),evHadSys_Tau_4P_.Phi(),E_fake);

  TLorentzVector Z_4P = lep_Z_faketau_4P + evLep_Tau_4P_;

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the leptonic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector b_lep_reco_4P;
  TLorentzVector W_lep_4P;
  TLorentzVector nu_top_4P;
  TLorentzVector t_lep_4P;

  if(error==0){
    
    TVector3 enu;
    enu.SetMagThetaPhi(1,TMath::ACos(cosThetaNu),phiNu);
    double CosThetalnu_top=TMath::Cos(evLep_top_4P_.Angle(enu));
    double Enu=getEqbar_Enu(CosThetalnu_top,evLep_top_4P_.E());
    TLorentzVector nu_top_4P_dummy(Enu*enu,Enu);
    nu_top_4P=nu_top_4P_dummy;
    W_lep_4P=nu_top_4P+evLep_top_4P_;

    
    double Eb_lep=getEb(W_lep_4P, evBJet_leptop_4P_);
    double b_lep_Pmag=TMath::Sqrt(Eb_lep*Eb_lep - mb*mb);
    double pT_bl=b_lep_Pmag/TMath::CosH(evBJet_leptop_4P_.Eta());
    if(Eb_lep>0){
      b_lep_reco_4P.SetPtEtaPhiE(pT_bl,evBJet_leptop_4P_.Eta(),evBJet_leptop_4P_.Phi(),Eb_lep);
      t_lep_4P=W_lep_4P+b_lep_reco_4P;
    }
    else{
      error = "leptop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build leptonic top :" << endl;
      (*stream_) << " .... b from lep top " 
              <<  b_lep_reco_4P.Px() << ", "
              <<  b_lep_reco_4P.Py() << ", "
              <<  b_lep_reco_4P.Pz() << ", "
              <<  b_lep_reco_4P.E()     
              << endl;
     (*stream_) << " .... lep top" 
              <<  t_lep_4P.Px() << ", "
              <<  t_lep_4P.Py() << ", "
              <<  t_lep_4P.Pz() << ", "
              <<  t_lep_4P.E()    
              << endl;
    }

  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@               Evaluation of the hadronic top 4-Vector                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector q1_reco_4P;
  TLorentzVector q2_reco_4P;
  TLorentzVector b_had_reco_4P;
  TLorentzVector W_had_4P;
  TLorentzVector t_had_4P;

  if(error==0){

    double pt_q1=EQuark1/TMath::CosH(evJet1_4P_.Eta());
    q1_reco_4P.SetPtEtaPhiE(pt_q1,evJet1_4P_.Eta(),evJet1_4P_.Phi(),EQuark1);
    
    if(integration_type_ == integration_type_wo_miss){
      double CosTheta_qq=TMath::Cos(evJet1_4P_.Angle(evJet2_4P_.Vect()));
      double EQuark2=getEqbar_Enu(CosTheta_qq,EQuark1);
      double pt_q2=EQuark2/TMath::CosH(evJet2_4P_.Eta());
      q2_reco_4P.SetPtEtaPhiE(pt_q2,evJet2_4P_.Eta(),evJet2_4P_.Phi(),EQuark2);
    }
    else if(integration_type_ == integration_type_w_miss){
      TVector3 e_q2;
      e_q2.SetMagThetaPhi(1,TMath::ACos(cosTheta_missing_jet),phi_missing_jet);
      double CosTheta_qq=TMath::Cos(evJet1_4P_.Angle(e_q2));
      double EQuark2=getEqbar_Enu(CosTheta_qq,EQuark1);
      double pt_q2=EQuark2/TMath::CosH(e_q2.Eta());
      q2_reco_4P.SetPtEtaPhiE(pt_q2,e_q2.Eta(),e_q2.Phi(),EQuark2);
    }

    W_had_4P=q1_reco_4P+q2_reco_4P;
    double Eb_had=getEb(W_had_4P, evBJet_hadtop_4P_);
    double b_had_Pmag=TMath::Sqrt(Eb_had*Eb_had - mb*mb);
    double pT_bh=b_had_Pmag/TMath::CosH(evBJet_hadtop_4P_.Eta());
    if(Eb_had>0){
      b_had_reco_4P.SetPtEtaPhiE(pT_bh,evBJet_hadtop_4P_.Eta(),evBJet_hadtop_4P_.Phi(),Eb_had);
      t_had_4P=W_had_4P+b_had_reco_4P;
    }
    else{
      error = "hadtop";
    }

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Build hadronic top :" << endl;
      (*stream_) << " .... b from had top " 
              <<  b_had_reco_4P.Px() << ", "
              <<  b_had_reco_4P.Py() << ", "
              <<  b_had_reco_4P.Pz() << ", "
              <<  b_had_reco_4P.E()     
              << endl;
      (*stream_) << " .... q1 from had top " 
              <<  q1_reco_4P.Px() << ", "
              <<  q1_reco_4P.Py() << ", "
              <<  q1_reco_4P.Pz() << ", "
              <<  q1_reco_4P.E()     
              << endl;
      (*stream_) << " .... q2 from had top " 
              <<  q2_reco_4P.Px() << ", "
              <<  q2_reco_4P.Py() << ", "
              <<  q2_reco_4P.Pz() << ", "
              <<  q2_reco_4P.E()     
              << endl;
      (*stream_) << " .... had top " 
		 <<  t_had_4P.Px() << ", "
		 <<  t_had_4P.Py() << ", "
		 <<  t_had_4P.Pz() << ", "
		 <<  t_had_4P.E()    
		 << endl;
    }


  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                Recoil                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  TLorentzVector tot4P;
  TLorentzVector rho4P;

  if(error==0){
    
    tot4P = Z_4P + t_lep_4P + t_had_4P ;

    if(integration_type_ == integration_type_wo_miss)
      rho4P = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_Tau_4P_ + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_ + evJet1_4P_ + evJet2_4P_ ;
    else if(integration_type_ == integration_type_w_miss)
      rho4P = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_Tau_4P_ + evLep_top_4P_ + evBJet_leptop_4P_ + evBJet_hadtop_4P_ + evJet1_4P_ ;

  }


  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        JET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

  double T_Jet1 = 1.0;
  double T_Jet2 = 1.0;
  double T_BJet_leptop = 1.0;
  double T_BJet_hadtop = 1.0;

  if ( error == 0 ) {

    if (flagTFJet1_) {
      if(integration_type_ == integration_type_wo_miss)
	T_Jet1 = getJetTF( q1_reco_4P, evJet1_4P_,"light");
      
      else if(integration_type_ == integration_type_w_miss){
	if(q1_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector q1q2_reco_4P = q1_reco_4P;
	  q1q2_reco_4P.SetE(q1_reco_4P.E() + q2_reco_4P.E());
	  T_Jet1 = getJetTF( q1q2_reco_4P, evJet1_4P_,"light");
	}
	else
	  T_Jet1 = getJetTF( q1_reco_4P, evJet1_4P_,"light");
      }
    }

    //
    if (flagTFJet2_) {
      if(integration_type_ == integration_type_wo_miss)
	T_Jet2 = getJetTF( q2_reco_4P, evJet2_4P_,"light");
      else if(integration_type_ == integration_type_w_miss)
	T_Jet2 = getMissJetAcceptance( q2_reco_4P);
    } 
    //
    if (flagTFBJet_leptop_) {
      if(integration_type_ == integration_type_wo_miss)
	T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");
      
      else if(integration_type_ == integration_type_w_miss){
	if(b_lep_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector b_lep_q2_reco_4P = b_lep_reco_4P;
	  b_lep_q2_reco_4P.SetE(b_lep_reco_4P.E() + q2_reco_4P.E());
	  T_BJet_leptop = getJetTF( b_lep_q2_reco_4P, evBJet_leptop_4P_,"b");
	}
	else
	  T_BJet_leptop = getJetTF( b_lep_reco_4P, evBJet_leptop_4P_,"b");

      }
    }
    //
    if (flagTFBJet_hadtop_) {
      if(integration_type_ == integration_type_wo_miss)
	T_BJet_hadtop = getJetTF( b_had_reco_4P, evBJet_hadtop_4P_,"b");
      
      else if(integration_type_ == integration_type_w_miss){
	if(b_had_reco_4P.DeltaR(q2_reco_4P)<jet_radius_){
	  // Two quarks are merged in a single jet: add two energies for TF
	  TLorentzVector b_had_q2_reco_4P = b_had_reco_4P;
	  b_had_q2_reco_4P.SetE(b_had_reco_4P.E() + q2_reco_4P.E());
	  T_BJet_hadtop = getJetTF( b_had_q2_reco_4P, evBJet_leptop_4P_,"b");
	}
	else
	  T_BJet_hadtop = getJetTF( b_had_reco_4P, evBJet_hadtop_4P_,"b");
      }
    } 


    // Log
/*    if (outgoingJetsOut_) {
      (*outgoingJetsOut_) << setiosflags(ios::fixed) << scientific;
      (*outgoingJetsOut_) 
         << " " << setw(21) << evJet1_4P_.Pt()
         << " " << setw(21) << evJet1_4P_.Eta()
         << " " << setw(21) << evJet1_4P_.Phi()
         << " " << setw(21) << evJet1_4P_.P()
         << " " << setw(21) << evJet1_4P_.E()      
         << " " << setw(21) << q1_reco_4P.Pt()
         << " " << setw(21) << q1_reco_4P.Eta()
         << " " << setw(21) << q1_reco_4P.Phi()
         << " " << setw(21) << q1_reco_4P.P()
         << " " << setw(21) << q1_reco_4P.E()
         << " " << setw(21) << T_Jet1
         << " " << setw(21) << evJet2_4P_.Pt()
         << " " << setw(21) << evJet2_4P_.Eta()
         << " " << setw(21) << evJet2_4P_.Phi()
         << " " << setw(21) << evJet2_4P_.P()
         << " " << setw(21) << evJet2_4P_.E()      
         << " " << setw(21) << q2_reco_4P.Pt()
         << " " << setw(21) << q2_reco_4P.Eta()
         << " " << setw(21) << q2_reco_4P.Phi()
         << " " << setw(21) << q2_reco_4P.P()
         << " " << setw(21) << q2_reco_4P.E()
         << " " << setw(21) << T_Jet2
	 << " " << setw(21) << evBJet_leptop_4P_.Pt()
         << " " << setw(21) << evBJet_leptop_4P_.Eta()
         << " " << setw(21) << evBJet_leptop_4P_.Phi()
         << " " << setw(21) << evBJet_leptop_4P_.P()
         << " " << setw(21) << evBJet_leptop_4P_.E()      
         << " " << setw(21) << b_lep_reco_4P.Pt()
         << " " << setw(21) << b_lep_reco_4P.Eta()
         << " " << setw(21) << b_lep_reco_4P.Phi()
         << " " << setw(21) << b_lep_reco_4P.P()
         << " " << setw(21) << b_lep_reco_4P.E()
         << " " << setw(21) << T_BJet_leptop
	 << " " << setw(21) << evBJet_hadtop_4P_.Pt()
         << " " << setw(21) << evBJet_hadtop_4P_.Eta()
         << " " << setw(21) << evBJet_hadtop_4P_.Phi()
         << " " << setw(21) << evBJet_hadtop_4P_.P()
         << " " << setw(21) << evBJet_hadtop_4P_.E()      
         << " " << setw(21) << b_had_reco_4P.Pt()
         << " " << setw(21) << b_had_reco_4P.Eta()
         << " " << setw(21) << b_had_reco_4P.Phi()
         << " " << setw(21) << b_had_reco_4P.P()
         << " " << setw(21) << b_had_reco_4P.E()
         << " " << setw(21) << T_BJet_hadtop
         << endl;
    } */
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        MET Transfer function                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_MET=1.0;
  if ( error == 0 ) {
    if (flagTFMET_) {

      TVector3 RhoT;
      TVector3 PT;
      
      if(include_hadrecoil_){
	
	RhoT.SetXYZ((rho4P).Px(), (rho4P).Py(), 0);
	PT.SetXYZ((tot4P).Px(), (tot4P).Py(), 0);

      }

      else{
	
	TLorentzVector tot4P_nohad =  Z_4P + W_lep_4P ;
	TLorentzVector rho4P_nohad = evRecoMET4P_ + evHadSys_Tau_4P_ + evLep_Tau_4P_ + evLep_top_4P_ ;
	RhoT.SetXYZ((rho4P_nohad).Px(), (rho4P_nohad).Py(), 0);
	PT.SetXYZ((tot4P_nohad).Px(), (tot4P_nohad).Py(), 0);


      }


      TVector3 Rec = RhoT-PT;
      // Order : (0,0), (0,1), (1,0), (1,1) ???
      double sqrtDetV = TMath::Sqrt( std::abs( evV_[0] * evV_[3] - evV_[1] * evV_[2]));
      double f = Rec.X()*(evV_[0]*Rec.X() + evV_[1]*Rec.Y()) 
               + Rec.Y()*(evV_[2]*Rec.X() + evV_[3]*Rec.Y());
      T_MET =  TMath::Exp( -0.5*f ) * sqrtDetV / ( 2 * TMath::Pi() ) ;
    

    }
  } 

 
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                   FAKE TAU from lep Transfer function                     @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_fake_leptau = 1.0;
  if ( error == 0 ) {
    if ( flagTF_fakeleptau_ )
      T_fake_leptau = getFakeTauLepTF( lep_Z_faketau_4P, evHadSys_Tau_4P_, lepton_Tau_Type_);
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                       TOP Transfer functions                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double T_lepTop = 1.0;
  double T_hadTop = 1.0;

  if ( error == 0 ) {
    if ( flagTFTop_ ){
      T_lepTop = getTopTF(t_lep_4P,b_lep_reco_4P,evLep_top_4P_,nu_top_4P);

      //We don't know which q is the quark or the antiquark
      T_hadTop = getTopTF(t_had_4P,b_had_reco_4P,q1_reco_4P,q2_reco_4P)
	+ getTopTF(t_had_4P,b_had_reco_4P,q2_reco_4P,q1_reco_4P);

    }
  }

  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Jacobian Term                             @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
 
  double Jac = 1.0;
  
  if ( error == 0 ) {  
    if( flagJac_ ) {

      
      double Jac_Z = getJacobian_W( lep_Z_faketau_4P.E(), evLep_Tau_4P_.E() );
      double cosTheta_qq=TMath::Cos(q1_reco_4P.Angle(q2_reco_4P.Vect()));
      double Jac_tophad = getJacobian_tophad( b_had_reco_4P, W_had_4P, q2_reco_4P.E(), cosTheta_qq);
      double Jac_toplep = getJacobian_toplep( b_lep_reco_4P, W_lep_4P, evLep_top_4P_.E(), nu_top_4P.E());
      Jac = Jac_Z * Jac_tophad * Jac_toplep;

      if( TMath::IsNaN( Jac ) ) { if (error == 0) error = "Jacobian";}  
    }
  }


  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                                 BOOST                                @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double totPt = 0.0;
  TLorentzVector tot4P_boosted;
  TLorentzVector Lep1_4P_boosted;
  TLorentzVector Lep2_4P_boosted;
  TLorentzVector t_had_4P_boosted;
  TLorentzVector t_lep_4P_boosted;

  if ( error == 0 ) {
    if (verbose_ >= IntegrandLevel){
      (*stream_) << "Boosting :" << endl;
      (*stream_) << "  initial total momentum tot4P " 
              <<  tot4P.Px() << ", "
              <<  tot4P.Py() << ", "
              <<  tot4P.Pz() << ", "
              <<  tot4P.E() 
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
    totPt = tot4P.Pt();
    
    // Set transverse boost
    TVector3 boost3D = tot4P.BoostVector();
    boost3D.SetZ(0);
    
    //
    // Application of the boost to the final state 4-momenta
    // WARRNING : only following 4-momenta are now booosted
    tot4P_boosted=tot4P;
    tot4P_boosted.Boost   ( -boost3D );

    // ??? implulsion > energy 
    if( TMath::IsNaN(tot4P_boosted.E()) ) { if (error == 0) error = "Boost";}  

    Lep1_4P_boosted=evLep_Tau_4P_;
    Lep2_4P_boosted=lep_Z_faketau_4P;
    t_had_4P_boosted=t_had_4P;
    t_lep_4P_boosted=t_lep_4P;

    Lep1_4P_boosted.Boost( -boost3D );
    Lep2_4P_boosted.Boost( -boost3D );
    t_had_4P_boosted.Boost( -boost3D );
    t_lep_4P_boosted.Boost( -boost3D );

    if (verbose_ >= IntegrandLevel){
      (*stream_) << "  Boosted total momentum tot4P " 
              <<  tot4P_boosted.Px() << ", "
              <<  tot4P_boosted.Py() << ", "
              <<  tot4P_boosted.Pz() << ", "
              <<  tot4P_boosted.E() 
              << endl;
      (*stream_) << "  Boost magnitude " 
              <<  totPt
              << endl;
      (*stream_) << "  initial total momentum tot4Pt "  << tot4P.Pt() << endl;    
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                        THE BJORKEN FRACTION                          @@@
  //@@@                         Matrix Elements                              @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  double wME = 1.0;
  double x_a;
  double x_b;
  x_a = (tot4P_boosted.E() + tot4P_boosted.Pz()) / sqrtS_ ;
  x_b = (tot4P_boosted.E() - tot4P_boosted.Pz()) / sqrtS_ ;
  
  if ( (x_a > 1.0) || (x_a < 0.0) || (x_b > 1.0 ) ||  (x_b < 0.0 ) ) 
    error = "BFrac";
  
  if ( error == 0 ) {

    TLorentzVector top_4P;
    TLorentzVector Atop_4P;
    TLorentzVector Lep_4P;
    TLorentzVector ALep_4P;

    if(lepton_Tau_charge_>0){
      ALep_4P=Lep1_4P_boosted;
      Lep_4P=Lep2_4P_boosted;
    }
    else{
      Lep_4P=Lep1_4P_boosted;
      ALep_4P=Lep2_4P_boosted;
    }

    if(lepton_top_charge_>0){
      top_4P=t_lep_4P_boosted;
      Atop_4P=t_had_4P_boosted;
    }
    else{
      top_4P=t_had_4P_boosted;
      Atop_4P=t_lep_4P_boosted;
    }


    // Parton momentum fraction
    //
    // MadGraph quadri-vector construction
    double inPart_a[4], inPart_b[4]; 
    double outLep[4], outALep[4];
    double outtop[4], outAtop[4];
    double inE;
    // Incoming quarks
    inE =  0.5 * x_a * sqrtS_;
    inPart_a[0] = inE; inPart_a[1] = 0.; inPart_a[2] =  0; inPart_a[3] = inE; 
    inE =  0.5 * x_b * sqrtS_;
    inPart_b[0] = inE; inPart_b[1] = 0.; inPart_b[2] =  0; inPart_b[3] =-inE; 
    // Lep / ALep
    outLep[0] = Lep_4P.E(); outLep[1] = Lep_4P.Px(); outLep[2] = Lep_4P.Py(); outLep[3] = Lep_4P.Pz();
    outALep[0] = ALep_4P.E(); outALep[1] = ALep_4P.Px(); outALep[2] = ALep_4P.Py(); outALep[3] = ALep_4P.Pz();
    // top / Atop
    outtop[0] = top_4P.E(); outtop[1] = top_4P.Px(); outtop[2] = top_4P.Py(); outtop[3] = top_4P.Pz();
    outAtop[0] = Atop_4P.E(); outAtop[1] = Atop_4P.Px(); outAtop[2] = Atop_4P.Py(); outAtop[3] = Atop_4P.Pz();

    // Quadri-vector order for MadGraph
    //   p[0] incoming parton0
    //   p[1] incoming parton1
    //   p[2] outgoing top
    //   p[3] outgoing Atop
    //   p[4] outgoing lep-
    //   p[5] outgoing lep+
    vector<double*> p(6,(double*) 0);
    p[0] = inPart_a ; p[1] = inPart_b; 
    p[2] = outtop; p[3] = outAtop; 
    p[4] = outLep ; p[5] = outALep   ; 

    // Matrix Element
    if (flagWME_) {
      wME = get_ttZZonlyZllWeightedME( x_a, x_b, p);
      if( TMath::IsNaN(wME) ) { if (error == 0) error = "WeightME";}
    }
  }
  
  
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  //@@@                            Weight evaluation                         @@@
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  double Eval =  T_Jet1 * T_Jet2 * T_BJet_hadtop * T_BJet_leptop * T_MET * T_fake_leptau * T_lepTop* T_hadTop * Jac * wME;
  if (error) Eval = 0;

  if (verbose_ >= IntegrandLevel){
    // if( Eval > 10.0e-40) {
    (*stream_) << "Integrand value :" 
       << " T_Jet1 = " << T_Jet1
       << " T_Jet2 = " << T_Jet2
       << " T_BJet_hadtop = " << T_BJet_hadtop
       << " T_BJet_leptop = " << T_BJet_leptop
       << " T_MET = " << T_MET
       << " T_fake_leptau = " << T_fake_leptau
       << " T_lepTop = " << T_lepTop
       << " T_hadTop = " << T_hadTop
       << " Jac = " << Jac
       << " wME = " << wME 
       << " Eval = "    << Eval
       << " Error= " << error
       << endl;
    // }
  }
/*  const char *err_str = "Ok";
  if (error) 
    err_str = error;
  else{*/
  if ( !error) {
    integr_EfficiencyttZ_Zll_++;
  }

  tot_DrawsttZ_Zll_++;

  return Eval;
}

extern "C" double wrapper_evalttH( double *x, size_t dim, void* param) {
  return static_cast<MGIntegration*>(param)->evalttH(x);
}

extern "C" double wrapper_evalttW( double *x, size_t dim, void* param) {
  return static_cast<MGIntegration*>(param)->evalttW(x);
}

extern "C" double wrapper_evalttjets( double *x, size_t dim, void* param) {
  return static_cast<MGIntegration*>(param)->evalttjets(x);
}

extern "C" double wrapper_evalttbar_SL( double *x, size_t dim, void* param) {
  return static_cast<MGIntegration*>(param)->evalttbar_SL(x);
}

extern "C" double wrapper_evalttbar_DL( double *x, size_t dim, void* param) {
  return static_cast<MGIntegration*>(param)->evalttbar_DL(x);
}

extern "C" double wrapper_evalttbar_DL_fakelep_tlep( double *x, size_t dim, void* param) {
  return static_cast<MGIntegration*>(param)->evalttbar_DL_fakelep_tlep(x);
}

extern "C" double wrapper_evalttbar_DL_fakelep_ttau( double *x, size_t dim, void* param) {
  return static_cast<MGIntegration*>(param)->evalttbar_DL_fakelep_ttau(x);
}

extern "C" double wrapper_evalttZ_Zll( double *x, size_t dim, void* param) {
  return static_cast<MGIntegration*>(param)->evalttZ_Zll(x);
}
