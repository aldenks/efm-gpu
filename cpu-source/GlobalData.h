#ifndef GLOBALDATA_H
#define	GLOBALDATA_H

#include "Network.h"

#define null 0
#define ZERO 10e-5
#define NEG_ZERO -ZERO

//Network
extern Network network;
//Total number of reactions in the network
extern int numReactions;
//Total number of metbaolites in the network
extern int numMetabolites;
//Boolean flag for external metabolites (true if a metabolite is external)
extern vector<bool> externalMetabolites;
//Boolean flag for reversible reactions (true if a reaction is reversible)
extern vector<bool> reversible;
//Number of remaining metabolites
extern int numMetabolitesRemaining;
//Reaction names
extern vector<string> reactions;
//Metabolites names
extern vector<string> metabolites;

#endif	/* GLOBALDATA_H */
