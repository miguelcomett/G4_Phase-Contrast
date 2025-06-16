#ifndef PepiSteppingAction_h
#define PepiSteppingAction_h 1
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "PepiDetectorConstruction.hh"
#include "PepiEventAction.hh"

class PepiSteppingAction : public G4UserSteppingAction
{
public:
   PepiSteppingAction(PepiEventAction* eventAction);
   ~PepiSteppingAction();
	
   virtual void UserSteppingAction(const G4Step*);

private:
   PepiEventAction* fEventAction;
};

#endif
