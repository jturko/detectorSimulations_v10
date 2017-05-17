//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file analysis/shared/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
// $Id: EventAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* run, HistoManager* histo)
    :G4UserEventAction(),
      fRunAct(run),fHistoManager(histo),
      fPrintModulo(0)
{
    pLabAngle = -1;
    numberOfHits = 0;
    numberOfSteps = 0;
    fPrintModulo = 1000;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
    const G4int numEvents = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
    evtNb = evt->GetEventID();
    if (evtNb%fPrintModulo == 0)
        //    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
        printf( " ---> Ev.# %5d :: %.2f %% \r", evtNb, 100.*G4double(evtNb)/G4double(numEvents));
    G4cout.flush();
    
    //std::cout << "Event # = " << evtNb << std::endl;

    ClearVariables();

    SetTotScintPhotons(0);
    SetQuartzScintPhotons(0);
    
    pParticleMap.clear();
    pLabAngle = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
    //    FillParticleType() ;
    FillGriffinCryst() ;
    Fill8piCryst() ;
    FillLaBrCryst() ;
    FillAncillaryBgo() ;
    FillSceptar() ;
    FillGridCell() ;

    //G4cout << "numberOfHits = " << numberOfHits << G4endl;
    for (G4int i = 0 ; i < numberOfHits; i++) {
		fHistoManager->FillHitNtuple(hitTrackerI[0][i], hitTrackerI[1][i], hitTrackerI[2][i], hitTrackerI[3][i],  hitTrackerI[4][i], hitTrackerI[5][i], hitTrackerI[6][i], hitTrackerI[7][i], hitTrackerI[8][i], hitTrackerD[0][i]/keV, hitTrackerD[1][i]/mm, hitTrackerD[2][i]/mm, hitTrackerD[3][i]/mm, hitTrackerD[4][i]/second, hitTrackerI[9][i], GetTotScintPhotons(), GetQuartzScintPhotons(), hitTrackerD[5][i]/keV, hitTrackerD[6][i]/keV, hitTrackerD[7][i]/keV, hitTrackerD[8][i]/keV, hitTrackerD[9][i]/keV, hitTrackerD[10][i]/keV, hitTrackerD[11][i]/keV, hitTrackerD[12][i]/keV, hitTrackerD[13][i]/keV, hitTrackerD[14][i]/keV , hitTrackerD[15][i]/keV, hitTrackerD[16][i]);
    }
    for (G4int i = 0 ; i < numberOfSteps; i++) {
		fHistoManager->FillStepNtuple(stepTrackerI[0][i], stepTrackerI[1][i], stepTrackerI[2][i], stepTrackerI[3][i],  stepTrackerI[4][i], stepTrackerI[5][i], stepTrackerI[6][i], stepTrackerI[7][i], stepTrackerI[8][i], stepTrackerD[0][i]/keV, stepTrackerD[1][i]/mm, stepTrackerD[2][i]/mm, stepTrackerD[3][i]/mm, stepTrackerD[4][i]/second, stepTrackerI[9][i], GetTotScintPhotons(), GetQuartzScintPhotons(), stepTrackerD[5][i]/keV, stepTrackerD[6][i]/keV, stepTrackerD[7][i]/keV, stepTrackerD[8][i]/keV, stepTrackerD[9][i]/keV, stepTrackerD[10][i]/keV, stepTrackerD[11][i]/keV, stepTrackerD[12][i]/keV, stepTrackerD[13][i]/keV, stepTrackerD[14][i]/keV , stepTrackerD[15][i]/keV, stepTrackerD[16][i]);
    }

    ClearVariables();

    const G4int numEvents = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
    if(numEvents < 10) G4cout << "totScintPhotons = " << GetTotScintPhotons() << "    quartzScintPhotons = " << GetQuartzScintPhotons() << G4endl;
    
    pParticleMap.clear();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EventAction::ClearVariables()
{
    if(fHistoManager->GetStepTrackerBool()) {
        stepIndex = 0;
        numberOfSteps = 0;
        for (G4int i = 0 ; i < MAXSTEPS; i++) {
            for (G4int j = 0 ; j < NUMSTEPVARS; j++) {
                stepTrackerI[j][i] = 0;
                stepTrackerD[j][i] = 0.0;
            }
        }
    }

    if(fHistoManager->GetHitTrackerBool()) {
        hitIndex = 0;
        numberOfHits = 0;
        pTrackID = -1;
        pParentID = -1;

        for (G4int i = 0 ; i < MAXHITS; i++) {
            pHitMnemonic[i] = "XXX00XX00X";
            for (G4int j = 0 ; j < NUMSTEPVARS; j++) {
                hitTrackerI[j][i] = 0;
                hitTrackerD[j][i] = 0.0;
            }
        }
    }

    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++) {
        particleTypes[i]                  = 0;
    }
    for (G4int i = 0 ; i < MAXNUMDET; i++) {
        EightPiCrystEnergyDet[i]  = 0 ;
        EightPiCrystTrackDet[i]   = 0 ;

        LaBrCrystEnergyDet[i]  = 0 ;
        LaBrCrystTrackDet[i]   = 0 ;

        AncillaryBgoEnergyDet[i]  = 0 ;
        AncillaryBgoTrackDet[i]   = 0 ;

        SceptarEnergyDet[i]  = 0 ;
        SceptarTrackDet[i]   = 0 ;

        GridCellElectronEKinDet[i]  = 0 ;
        GridCellElectronTrackDet[i]   = 0 ;
        GridCellGammaEKinDet[i]  = 0 ;
        GridCellGammaTrackDet[i]   = 0 ;
        GridCellNeutronEKinDet[i]  = 0 ;
        GridCellNeutronTrackDet[i]   = 0 ;

    }
    for (G4int i = 0 ; i < MAXNUMDETGRIFFIN; i++) {
        for (G4int j = 0 ; j < MAXNUMCRYGRIFFIN; j++) {

            GriffinCrystEnergyDet[i][j]                 = 0;
            GriffinCrystTrackDet[i][j]                  = 0;
            GriffinSuppressorBackEnergyDet[i][j]        = 0;
            GriffinSuppressorBackTrackDet[i][j]         = 0;
            GriffinSuppressorLeftExtensionEnergyDet[i][j]   = 0;
            GriffinSuppressorLeftExtensionTrackDet[i][j]    = 0;
            GriffinSuppressorLeftSideEnergyDet[i][j]        = 0;
            GriffinSuppressorLeftSideTrackDet[i][j]         = 0;
            GriffinSuppressorRightExtensionEnergyDet[i][j]   = 0;
            GriffinSuppressorRightExtensionTrackDet[i][j]    = 0;
            GriffinSuppressorRightSideEnergyDet[i][j]        = 0;
            GriffinSuppressorRightSideTrackDet[i][j]         = 0;
        }
    }
    // NOTE: Clear the variables from the new Fill___Cryst functions.

    fHistoManager->clear_Edep();
    fHistoManager->clear_Ekin();
    fHistoManager->clear_Ptype();
}


void EventAction::FillParticleType()
{
    G4int numParticleTypes = 0;
    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++)
    {
        if (particleTypes[i] != 0)
        { // if particle type 'i' has non-zero counts
            for (G4int j = 0 ; j< particleTypes[i]; j++)
            { // loop over the number of time we saw it
                G4cout << "particleTypes[" << i << "] = " << particleTypes[i] << G4endl;
                fHistoManager->FillHisto(astats_particle_type_in_each_step, i);
            }
        }
    }

    // Fill the number of particle types in the event
    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++)
    {
        if (particleTypes[i] != 0)
            numParticleTypes++;
    }
    fHistoManager->FillHisto(astats_particle_type_in_each_event, numParticleTypes);
}

void EventAction::FillGriffinCryst()
{
    G4double  energySum = 0;
    G4double  energySumDet = 0;
    G4bool SuppressorBackFired[MAXNUMDETGRIFFIN] = {0};
    G4bool SuppressorExtensionFired[MAXNUMDETGRIFFIN] = {0};
    G4bool SuppressorSideFired[MAXNUMDETGRIFFIN] = {0};
    G4bool SuppressorFired = false;

    // Fill Griffin Histos
    for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
        energySumDet = 0;
        // Find if any suppressors were fired
        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            if (GriffinSuppressorBackEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorBackFired[i] = true;
            }
            if (GriffinSuppressorLeftExtensionEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorExtensionFired[i] = true;
            }
            if (GriffinSuppressorRightExtensionEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorExtensionFired[i] = true;
            }
            if (GriffinSuppressorLeftSideEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorSideFired[i] = true;
            }
            if (GriffinSuppressorRightSideEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorSideFired[i] = true;
            }
            if ( !SuppressorFired && ( SuppressorBackFired[i] || SuppressorExtensionFired[i] || SuppressorSideFired[i] ) ) {
                SuppressorFired = true;
            }
        }

        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            if(GriffinCrystEnergyDet[i][j] > MINENERGYTHRES) {
                // fill energies in each crystal
                if(WRITEEDEPHISTOS)     fHistoManager->FillHisto((griffin_crystal_unsup_edep_det0_cry0+(MAXNUMDETGRIFFIN*j))+i, GriffinCrystEnergyDet[i][j]);
                if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(griffin_crystal_unsup_edep_cry, GriffinCrystEnergyDet[i][j]);
                if(!SuppressorBackFired[i] && !SuppressorExtensionFired[i] && !SuppressorSideFired[i]) { // Suppressor fired?
                    if(WRITEEDEPHISTOS) fHistoManager->FillHisto((griffin_crystal_sup_edep_det0_cry0+(MAXNUMDETGRIFFIN*j))+i, GriffinCrystEnergyDet[i][j]);
                    if(WRITEEDEPHISTOS) fHistoManager->FillHisto(griffin_crystal_sup_edep_cry, GriffinCrystEnergyDet[i][j]);
                }
                energySumDet += GriffinCrystEnergyDet[i][j];
            }
        }
        if(energySumDet > MINENERGYTHRES) {
            // fill energies in each detector
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(griffin_crystal_unsup_edep_det0+i, energySumDet);
            // fill standard energy and track spectra
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(griffin_crystal_unsup_edep, energySumDet);
            if(!SuppressorBackFired[i] && !SuppressorExtensionFired[i] && !SuppressorSideFired[i]) {
                // fill energies in each detector
                if(WRITEEDEPHISTOS) fHistoManager->FillHisto(griffin_crystal_sup_edep_det0+i, energySumDet);
                // fill standard energy and track spectra
                if(WRITEEDEPHISTOS) fHistoManager->FillHisto(griffin_crystal_sup_edep, energySumDet);
            }
        }
        energySum += energySumDet;
    }

    if(energySum > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(griffin_crystal_unsup_edep_sum, energySum);
        if(!SuppressorFired) {
            if(WRITEEDEPHISTOS) fHistoManager->FillHisto(griffin_crystal_sup_edep_sum, energySum);
        }
    }
}

void EventAction::Fill8piCryst()
{
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(EightPiCrystEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(Eightpi_crystal_edep, EightPiCrystEnergyDet[j]);
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(Eightpi_crystal_edep_det0+j, EightPiCrystEnergyDet[j]);
            energySumDet += EightPiCrystEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(Eightpi_crystal_edep_sum, energySumDet);
    }

}


void EventAction::FillLaBrCryst()
{
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(LaBrCrystEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(labr_crystal_edep, LaBrCrystEnergyDet[j]);
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(labr_crystal_edep_det0+j, LaBrCrystEnergyDet[j]);
            energySumDet += LaBrCrystEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(labr_crystal_edep_sum, energySumDet);
    }
}

void EventAction::FillAncillaryBgo()
{
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(AncillaryBgoEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(ancillary_bgo_crystal_edep, AncillaryBgoEnergyDet[j]);
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(ancillary_bgo_crystal_edep_det0+j, AncillaryBgoEnergyDet[j]);
            energySumDet += AncillaryBgoEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(ancillary_bgo_crystal_edep_sum, energySumDet);
    }
}

void EventAction::FillSceptar()
{
    G4double energySumDet = 0;
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(SceptarEnergyDet[j] > MINENERGYTHRES) {
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(sceptar_edep, SceptarEnergyDet[j]);
            if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(sceptar_edep_det0+j, SceptarEnergyDet[j]);
            energySumDet += SceptarEnergyDet[j];
        }
    }
    if(energySumDet > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     fHistoManager->FillHisto(sceptar_edep_sum, energySumDet);
    }

}

void EventAction::FillGridCell()
{
    for (G4int j=0; j < MAXNUMDET; j++) {
        if(WRITEEKINHISTOS && GridCellElectronEKinDet[j] > MINENERGYTHRES)     fHistoManager->FillHisto(gridcell_electron_ekin_det0+j, GridCellElectronEKinDet[j]);
        if(WRITEEKINHISTOS && GridCellGammaEKinDet[j] > MINENERGYTHRES)     fHistoManager->FillHisto(gridcell_gamma_ekin_det0+j, GridCellGammaEKinDet[j]);
        if(WRITEEKINHISTOS && GridCellNeutronEKinDet[j] > MINENERGYTHRES)     fHistoManager->FillHisto(gridcell_neutron_ekin_det0+j, GridCellNeutronEKinDet[j]);
    }
}

void EventAction::AddHitTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ, G4int numScintPhotons, G4int numQuartzPhotons, G4double eDepD, G4double eDepC, G4double eDepP, G4double eDepA, G4double eDepE, G4double eDepN, G4double eDepOther, G4double eDepBe, G4double eDepB, G4double eDepT, G4double eDepG, G4double lab_angle)
{
    G4bool newhit = true;
    for (G4int i = 0 ; i < numberOfHits; i++) {
        if(pHitMnemonic[i] == mnemonic) {
            // sum the new enery
            hitTrackerD[0][i] = hitTrackerD[0][i] + depEnergy;
            hitTrackerI[4][i] = particleType;
            hitTrackerD[1][i] = posx;
            hitTrackerD[2][i] = posy;
            hitTrackerD[3][i] = posz;
            hitTrackerD[5][i] = hitTrackerD[5][i] + eDepD;
            hitTrackerD[6][i] = hitTrackerD[6][i] + eDepC;
            hitTrackerD[7][i] = hitTrackerD[7][i] + eDepP;
            hitTrackerD[8][i] = hitTrackerD[8][i] + eDepA;
            hitTrackerD[9][i] = hitTrackerD[9][i] + eDepE;
            hitTrackerD[10][i] = hitTrackerD[10][i] + eDepN;
            hitTrackerD[11][i] = hitTrackerD[11][i] + eDepOther;
            hitTrackerD[12][i] = hitTrackerD[12][i] + eDepBe;
            hitTrackerD[13][i] = hitTrackerD[13][i] + eDepB;
            hitTrackerD[14][i] = hitTrackerD[14][i] + eDepT;
            hitTrackerD[15][i] = hitTrackerD[15][i] + eDepG;
            hitTrackerD[16][i] = pLabAngle;
            newhit = false;
            break;
        }
    }
    if (newhit) { // new hit
        pHitMnemonic[hitIndex] = mnemonic;
        pTrackID = trackID;
        pParentID = parentID;
        hitTrackerI[0][hitIndex] = eventNumber;
        hitTrackerI[1][hitIndex] = trackID;
        hitTrackerI[2][hitIndex] = parentID;
        hitTrackerI[3][hitIndex] = stepNumber;
        hitTrackerI[4][hitIndex] = particleType;
        hitTrackerI[5][hitIndex] = processType;
        hitTrackerI[6][hitIndex] = systemID;
        hitTrackerI[7][hitIndex] = cryNumber;
        hitTrackerI[8][hitIndex] = detNumber;
        hitTrackerI[9][hitIndex] = targetZ;
        hitTrackerI[10][hitIndex] = numScintPhotons;
        hitTrackerI[11][hitIndex] = numQuartzPhotons;
        hitTrackerD[0][hitIndex] = depEnergy;
        hitTrackerD[1][hitIndex] = posx;
        hitTrackerD[2][hitIndex] = posy;
        hitTrackerD[3][hitIndex] = posz;
        hitTrackerD[4][hitIndex] = time;
        hitTrackerD[5][hitIndex] = eDepD;
        hitTrackerD[6][hitIndex] = eDepC;
        hitTrackerD[7][hitIndex] = eDepP;
        hitTrackerD[8][hitIndex] = eDepA;
        hitTrackerD[9][hitIndex] = eDepE;
        hitTrackerD[10][hitIndex] = eDepN;
        hitTrackerD[11][hitIndex] = eDepOther;
        hitTrackerD[12][hitIndex] = eDepBe;
        hitTrackerD[13][hitIndex] = eDepB;
        hitTrackerD[14][hitIndex] = eDepT;
        hitTrackerD[15][hitIndex] = eDepG;
        hitTrackerD[16][hitIndex] = pLabAngle;

        hitIndex++;
        numberOfHits = hitIndex;

        if(numberOfHits >= MAXHITS) {
            G4cout << "ERROR! Too many hits!" << G4endl;
        }
    }

    // push back the new element(s) of the vector(s) w/ edep, ekin, and particleType
    fHistoManager->push_back_Edep(depEnergy);
    fHistoManager->push_back_Ekin(0.);
    fHistoManager->push_back_Ptype(particleType);
}

void EventAction::AddHitTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ, G4int numScintPhotons, G4int numQuartzPhotons, G4double eDepD, G4double eDepC, G4double eDepP, G4double eDepA, G4double eDepE, G4double eDepN, G4double eDepOther, G4double eDepBe, G4double eDepB, G4double eDepT, G4double eDepG, G4double lab_angle, G4double kinEnergy) // with ekin
{
    G4bool newhit = true;
    for (G4int i = 0 ; i < numberOfHits; i++) {
        if(pHitMnemonic[i] == mnemonic) {
            // sum the new enery
            hitTrackerD[0][i] = hitTrackerD[0][i] + depEnergy;
            hitTrackerD[1][i] = posx;
            hitTrackerD[2][i] = posy;
            hitTrackerD[3][i] = posz;
            hitTrackerI[4][i] = particleType;
            hitTrackerD[5][i] = hitTrackerD[5][i] + eDepD;
            hitTrackerD[6][i] = hitTrackerD[6][i] + eDepC;
            hitTrackerD[7][i] = hitTrackerD[7][i] + eDepP;
            hitTrackerD[8][i] = hitTrackerD[8][i] + eDepA;
            hitTrackerD[9][i] = hitTrackerD[9][i] + eDepE;
            hitTrackerD[10][i] = hitTrackerD[10][i] + eDepN;
            hitTrackerD[11][i] = hitTrackerD[11][i] + eDepOther;
            hitTrackerD[12][i] = hitTrackerD[12][i] + eDepBe;
            hitTrackerD[13][i] = hitTrackerD[13][i] + eDepB;
            hitTrackerD[14][i] = hitTrackerD[14][i] + eDepT;
            hitTrackerD[15][i] = hitTrackerD[15][i] + eDepG;
            hitTrackerD[16][i] = pLabAngle;
            newhit = false;
            break;
        }
    }
    if (newhit) { // new hit
        pHitMnemonic[hitIndex] = mnemonic;
        pTrackID = trackID;
        pParentID = parentID;
        hitTrackerI[0][hitIndex] = eventNumber;
        hitTrackerI[1][hitIndex] = trackID;
        hitTrackerI[2][hitIndex] = parentID;
        hitTrackerI[3][hitIndex] = stepNumber;
        hitTrackerI[4][hitIndex] = particleType;
        hitTrackerI[5][hitIndex] = processType;
        hitTrackerI[6][hitIndex] = systemID;
        hitTrackerI[7][hitIndex] = cryNumber;
        hitTrackerI[8][hitIndex] = detNumber;
        hitTrackerI[9][hitIndex] = targetZ;
        hitTrackerI[10][hitIndex] = numScintPhotons;
        hitTrackerI[11][hitIndex] = numQuartzPhotons;
        hitTrackerD[0][hitIndex] = depEnergy;
        hitTrackerD[1][hitIndex] = posx;
        hitTrackerD[2][hitIndex] = posy;
        hitTrackerD[3][hitIndex] = posz;
        hitTrackerD[4][hitIndex] = time;
        hitTrackerD[5][hitIndex] = eDepD;
        hitTrackerD[6][hitIndex] = eDepC;
        hitTrackerD[7][hitIndex] = eDepP;
        hitTrackerD[8][hitIndex] = eDepA;
        hitTrackerD[9][hitIndex] = eDepE;
        hitTrackerD[10][hitIndex] = eDepN;
        hitTrackerD[11][hitIndex] = eDepOther;
        hitTrackerD[12][hitIndex] = eDepBe;
        hitTrackerD[13][hitIndex] = eDepB;
        hitTrackerD[14][hitIndex] = eDepT;
        hitTrackerD[15][hitIndex] = eDepG;
        hitTrackerD[16][hitIndex] = pLabAngle;

        hitIndex++;
        numberOfHits = hitIndex;

        if(numberOfHits >= MAXHITS) {
            G4cout << "ERROR! Too many hits!" << G4endl;
        }
    }

    // push back the new element(s) of the vector(s) w/ edep, ekin, and particleType
    fHistoManager->push_back_Edep(depEnergy);
    fHistoManager->push_back_Ekin(kinEnergy);
    fHistoManager->push_back_Ptype(particleType);
}


void EventAction::AddStepTracker(G4String mnemonic, G4int eventNumber, G4int trackID, G4int parentID, G4int stepNumber, G4int particleType, G4int processType, G4int systemID, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, G4int targetZ, G4int numScintPhotons, G4int numQuartzPhotons, G4double eDepD, G4double eDepC, G4double eDepP, G4double eDepA, G4double eDepE, G4double eDepN, G4double eDepOther, G4double eDepBe, G4double eDepB, G4double eDepT, G4double eDepG, G4double lab_angle)
{
    G4bool newstep = true;
    if (newstep) { // new step
        stepTrackerI[0][stepIndex] = eventNumber;
        stepTrackerI[1][stepIndex] = trackID;
        stepTrackerI[2][stepIndex] = parentID;
        stepTrackerI[3][stepIndex] = stepNumber;
        stepTrackerI[4][stepIndex] = particleType;
        stepTrackerI[5][stepIndex] = processType;
        stepTrackerI[6][stepIndex] = systemID;
        stepTrackerI[7][stepIndex] = cryNumber;
        stepTrackerI[8][stepIndex] = detNumber;
        stepTrackerI[9][stepIndex] = targetZ;
        stepTrackerI[10][stepIndex] = numScintPhotons;
        stepTrackerI[11][stepIndex] = numQuartzPhotons;
        stepTrackerD[0][stepIndex] = depEnergy;
        stepTrackerD[1][stepIndex] = posx;
        stepTrackerD[2][stepIndex] = posy;
        stepTrackerD[3][stepIndex] = posz;
        stepTrackerD[4][stepIndex] = time;
        stepTrackerD[5][stepIndex] = eDepD;
        stepTrackerD[6][stepIndex] = eDepC;
        stepTrackerD[7][stepIndex] = eDepP;
        stepTrackerD[8][stepIndex] = eDepA;
        stepTrackerD[9][stepIndex] = eDepE;
        stepTrackerD[10][stepIndex] = eDepN;
        stepTrackerD[11][stepIndex] = eDepOther;
        stepTrackerD[12][stepIndex] = eDepBe;
        stepTrackerD[13][stepIndex] = eDepB;
        stepTrackerD[14][stepIndex] = eDepT;
        stepTrackerD[15][stepIndex] = eDepG;
        stepTrackerD[16][stepIndex] = pLabAngle;

        stepIndex++;

        numberOfSteps = stepIndex;
        if(numberOfSteps >= MAXSTEPS) {
            G4cout << "ERROR! Too many steps!" << G4endl;
        }
    }
}



