#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TFile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <Pythia8/Pythia.h>
#include <TStopwatch.h>

#include "fastjet/ClusterSequence.hh"
#include <iostream>

#include <TClonesArray.h>
#include "src/JHistos.h"
#include "src/AliJCard.h"

#define DEBUG 0

using namespace fastjet;
using namespace std;
using namespace Pythia8; 
double dotproduct(TLorentzVector jet,TLorentzVector particle);
void IdentifyParticleType(Event &part, int ith, vector<bool>& ispid);
enum { kCharged, kPion, kKaon, kProton, kCPlus, kCMinus, kPhoton, kNPType };//enumerating the option for the switch
// If true the run uses pion mass, and if false the run uses correct mass.
const bool usePionMass = false;

class MyUserInfo : public PseudoJet::UserInfoBase{
	public:
		// default ctor
		MyUserInfo(const int & pdg_id_in,const vector<bool> &pType) :
			_pdg_id(pdg_id_in){ ispid = pType;}

		/// access to the PDG id
		int pdg_id() const { return _pdg_id;}
		void PrintIsPid() const { 
			for(unsigned int i=0;i<ispid.size();i++) {
			 cout << ispid.at(i)<<":";
			}
			cout << endl;
		}
		bool IsType(int i) const { return ispid[i];}

	protected:
		int _pdg_id;         // the associated pdg id
		vector<bool> ispid;
};


int main(int argc, char **argv) {

	if(argc<4){
		cout<<"usage: pythia pythia.config card.input <output.root> [random_seed]"<<endl;exit(1);//usage of the code, requieres 4 inputs
	}
	TStopwatch timer; 
	timer.Start();   

	char* pythiaconfig  = argv[1];
	char* cardfilename  = argv[2];
	double pTHatMin     = atof(argv[3]);
	double pTHatMax     = atof(argv[4]);
	TString outputs = argv[5];
	Int_t random_seed = argc>6 ? atoi(argv[6]) : 0;//placing the inputs into variables


	TFile *fout = new TFile(outputs.Data(),"RECREATE");
	fout->cd();//opening of the output file

	//---------------------
	//Pythia initialization 
	//---------------------
	Pythia pythia;   // Generator.
	//Event& event      = pythia.event;
	ParticleData& pdt = pythia.particleData;

	// Read in commands from external file.
	pythia.readFile(pythiaconfig);   

	// Extract settings to be used in the main program.
	int    nEvent  = pythia.mode("Main:numberOfEvents");
	bool   showCS  = pythia.flag("Main:showChangedSettings");
	bool   showCPD = pythia.flag("Main:showChangedParticleData");
	double energy  = pythia.mode("Beams:eCM");

	pythia.readString(Form("PhaseSpace:pTHatMin ==%f",pTHatMin));
	pythia.readString(Form("PhaseSpace:pTHatMax ==%f",pTHatMax));
	cout<<"Events="<<nEvent <<" RNDM seed "<< random_seed << endl;

	pythia.readString("Random:setSeed = on");
	pythia.readString(Form("Random:seed=%02d",random_seed));

	// Initialize. Beam parameters set in .cmnd file.
	pythia.init();

	// List changed data. 
	if (showCS)  pythia.settings.listChanged();
	if (showCPD) pdt.listChanged();

	//-------------------------------------------------------
	// Histograms and tools
	//-------------------------------------------------------
	AliJCard *fcard = new AliJCard(cardfilename);
	fcard->PrintOut();
	JHistos *fhistos = new JHistos(fcard);
	fhistos->CreateQAHistos(); 
	fhistos->CreateDiJetHistos(); 

    fhistos->fhDiJetInvM->Sumw2(); // This is NOT done in the JHistos.cxx.
    
    // Logarithmic binning.
    int nBINS2=100;
    double logBinsX2[nBINS2+2], limL2=0.1, LimH2=2000;
    double logBW2 = (log(LimH2)-log(limL2))/nBINS2;
    for(int ij=0;ij<=nBINS2;ij++) logBinsX2[ij+1]=limL2*exp(ij*logBW2);
    logBinsX2[0]=0;

    TH1D *hDiJetInvMDeltaPhiCut = new TH1D("hDiJetInvMDeltaPhiCut","Di-jet inv M with DeltaPhi cut.",nBINS2, logBinsX2 );
    hDiJetInvMDeltaPhiCut->Sumw2();
	TH1D *hCrossSectionInfo = new TH1D("hCrossSection","CrossSectionInfo",8,0,8);
	hCrossSectionInfo->Sumw2();
	TH1D *hDeltaPhi = new TH1D("hDeltaPhi","DeltaPhi",100,0,2*TMath::Pi());
	hDeltaPhi->Sumw2();
	TH2D *hLjetPtVsNtrial = new TH2D("hLjetPtVsNtrial","Leading jet pT vs Ntrial",1000,10., 300.,1000,0,1e9);
	hLjetPtVsNtrial->Sumw2();

	//------------------------------------------------------------------
	// Define jet reconstruction
	//------------------------------------------------------------------
	vector<PseudoJet> finalparticles;
	vector<PseudoJet> chparticles;

	double PartMinPtCutForJet = 0.15;// atlas 0.5 cms/alice 0.15
	double etaMaxCutForJet = 0.5;// atlas 1.2 cms 2.5 alice 0.5
	double coneR = 0.4; // atlas 0.6, cms 0.7 alice 0.4
	double etaMaxCutForPart = etaMaxCutForJet+coneR;///
	double MinJetPt = 20.; // Min Jet Pt cut to disregard low pt jets

	JetDefinition jet_def(antikt_algorithm, coneR); 

	//TClonesArray *inputList = new TClonesArray("AliJBaseTrack",1500);
	//--------------------------------------------------------
	//         B e g i n    e v e n t    l o o p.
	//--------------------------------------------------------
	cout<<"Let's start" <<endl; 
	int ieout = nEvent/20;
	if (ieout<1) ieout=1;
	int EventCounter = 0;
	Int_t nTried = 0; 
	Int_t prev_nTried = 0;
	Int_t nTrial = 0;
	Int_t nAccepted = 0;
	Float_t sigmaGen = 0.0;
	Float_t ebeweight = 1.0;

	for(int iEvent = 0; iEvent < nEvent; ++iEvent) {//begin event loop

		if (!pythia.next()) continue;
		finalparticles.clear();
		chparticles.clear();
		int icent = 0;
		nTried = pythia.info.nTried();
		nTrial = nTried - prev_nTried;
		prev_nTried = nTried;
		sigmaGen = pythia.info.sigmaGen();
		ebeweight = 1.0; //no event-by-event weight at all. //sigmaGen/nTrial;
		fhistos->fhiCentr->Fill(icent);//filling the histogram counting events
		hCrossSectionInfo->Fill(7.5,ebeweight);
		if(iEvent % ieout == 0) cout << iEvent << "\t" << int(float(iEvent)/nEvent*100) << "%, nTried:" << nTried << ", nTrial:" << nTrial << ", sigma:" << sigmaGen << endl;

		for (int i = 0; i < pythia.event.size(); ++i) {//loop over all the particles in the event
			// Building input particle list for the jet reconstruction
			// Make sure "applying same cut done in ATLAS paper ?  flavor of the parciles ,eta , min pt cut 
			if (pythia.event[i].isFinal() && TMath::Abs(pythia.event[i].eta()) < etaMaxCutForPart && pythia.event[i].pT()>PartMinPtCutForJet){
				fhistos->fhJetPartPt->Fill(pythia.event[i].pT());			    
				fhistos->fhJetPartEta->Fill(pythia.event[i].eta());			    
				fhistos->fhJetPartPhi->Fill(pythia.event[i].phi());
				// Building the charged particle list to be used for manual constituent construction instead of using fastjet consistuents.
				if( pythia.event[i].isCharged() && pythia.event[i].isHadron() ) { // Only check if it is charged and hadron since the acceptance was checked in the previous if
					fhistos->fhchPartPt->Fill(pythia.event[i].pT());
					fhistos->fhchPartPhi->Fill(pythia.event[i].phi()); 
					fhistos->fhchPartEta->Fill(pythia.event[i].eta());//filling particle eta and pt histograms
                    if (usePionMass) { // Using pion+ mass for every particle.
                        chparticles.push_back(PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), TMath::Sqrt(TMath::Power(pythia.event[i].pz(),2.0) + pythia.event[i].pT2() + TMath::Power(0.13957018,2.0))));
                    } else { // Using real mass for every particle.
                        chparticles.push_back(PseudoJet(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));
                    }
				}
				vector<bool> IsPid(kNPType,false);
				IdentifyParticleType(pythia.event, i, IsPid);
				fastjet::PseudoJet particle(pythia.event[i].px(), pythia.event[i].py() , pythia.event[i].pz(), pythia.event[i].e());
				particle.set_user_index(pythia.event[i].id()); // check IdentifyParticleType, particle id into the index to build up hisgrams based of id selection later
				particle.set_user_info(new MyUserInfo(pythia.event[i].id(), IsPid));
				finalparticles.push_back( particle );
			} // end of finalparticles
		}

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// Run the clustering, Reconstruct jets
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		ClusterSequence cs(chparticles, jet_def);
		vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// Loop over jets and fill various histos 
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		for (unsigned i = 0; i < jets.size(); i++) {
			// jet eta cut
			if(fabs(jets[i].eta())>etaMaxCutForJet) continue; 
			fhistos->fhJetPseudorapidity->Fill(jets[i].eta(), ebeweight);  
			fhistos->fhJetRapidity->Fill(jets[i].rap(), ebeweight);  
			fhistos->fhJetPt->Fill(jets[i].pt(), ebeweight);
		}//end of the jet loop
		if( jets.size()> 0) {
			hLjetPtVsNtrial->Fill(jets[0].pt(),nTrial);
		}
		// DiJet 
		if(jets.size()>1) {
			if(jets[1].pt()>20.0) {
				PseudoJet dijet = jets[0] + jets[1];
				double mjj = dijet.m();
				double ptpair = dijet.pt();
				int mbin = fcard->GetBin(kDiJetType, mjj);
				fhistos->fhDiJetInvM->Fill(mjj, ebeweight);
				fhistos->fhDiJetInvMBins[mbin]->Fill(mjj, ebeweight);
				if(mbin>-1) fhistos->fhDiJetPtPairRaw[mbin]->Fill(ptpair, ebeweight);
				if(mbin>-1) fhistos->fhDiJetPtPair[mbin]->Fill(ptpair, ptpair<1e-8?1./1e-8:ebeweight/ptpair);
				double dPhi = jets[1].delta_phi_to(jets[0]);
				double dPhi2  = dPhi<0?dPhi+TMath::TwoPi():dPhi;
				if(mbin>-1) fhistos->fhDiJetDeltaPhi[mbin]->Fill(dPhi2, ebeweight);
				hDeltaPhi->Fill(dPhi2, ebeweight);

                // If subleading jet is on the opposite hemisphere compared to leading jet.
                if(dPhi2 > TMath::Pi()/2 && dPhi2 < 3/2*TMath::Pi()) {
                   hDiJetInvMDeltaPhiCut->Fill(mjj,ebeweight); 
                }
			}
		}

		EventCounter++;
		if(iEvent == nEvent-1) cout << nEvent << "\t" << "100%, nTried:" << pythia.info.nTried() << ", sigma:" << pythia.info.sigmaGen() << endl ;
	}//event loop

	nTried = pythia.info.nTried();
	nAccepted = pythia.info.nAccepted();
	sigmaGen = pythia.info.sigmaGen();
	//double sigmaErr = pythia.info.sigmaErr();
	hCrossSectionInfo->Fill(0.5,nTried);
	cout << "nTried after loop:" << nTried << endl;// print also inside event loop and in the macro.
	hCrossSectionInfo->Fill(1.5,nAccepted);
	//cout << "nAccepted after loop:" << nAccepted << endl;
	hCrossSectionInfo->Fill(2.5,sigmaGen);
	cout << "sigma after loop:" << sigmaGen << endl;
	hCrossSectionInfo->Fill(3.5,EventCounter);
	hCrossSectionInfo->Fill(4.5,energy);
	hCrossSectionInfo->Fill(5.5,1); // for counting # of merged
	hCrossSectionInfo->Fill(6.5,pythia.info.weightSum()); // for counting # of merged

	fout->Write();
	fcard->WriteCard(fout); // Write card into the file
	fout->Close();
	cout << EventCounter << " events are analyzed successfully."<< endl;
	timer.Print(); 
	return 0;
}


// CHECK FUNCTION
double dotproduct(TLorentzVector jets,TLorentzVector particles){//takes in two Lorentz vectors and transform it into TVector3 and makes the dot product
	double jetpx=jets.Px();
	double jetpy=jets.Py();
	double jetpz=jets.Pz();

	double partpx=particles.Px();
	double partpy=particles.Py();
	double partpz=particles.Pz();

	TVector3 jet=TVector3(jetpx,jetpy,jetpz);
	TVector3 particle=TVector3(partpx,partpy,partpz);

	double dotproduct= particle.Dot(jet);

	//delete jet;
	//delete particle;

	return dotproduct;
}

void IdentifyParticleType(Event &part, int i, vector<bool>& ispid) {
	// input pid from pythia
	// out  enum { kCharged, kPion, kKaon, kProton, kCPlus, kCMinus };//enumerating the option for the switch
	// only for charged particles
	if(part[i].isCharged() && part[i].isHadron()) ispid[kCharged]=true; // if not photon/electro/muon 
	if(part[i].id()==22) ispid[kPhoton]=true;
	if(TMath::Abs(part[i].id())==211 ) ispid[kPion]=true;
	if(part[i].id() == 311 || TMath::Abs(part[i].id()) == 321) ispid[kKaon]=true;
	if(TMath::Abs(part[i].id()) == 2212) ispid[kProton]=true;
	if(ispid[kCharged] && part[i].charge()>0) ispid[kCPlus]=true; 
	if(ispid[kCharged] && part[i].charge()<0) ispid[kCMinus]=true; 
}
