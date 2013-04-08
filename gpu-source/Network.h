#ifndef NETWORK_H
#define	NETWORK_H

#include <cstdio>
#include <string>
#include <vector>

using namespace std;

#define LINE_SIZE 102400
#define TOKEN_SIZE 1024

class Network {
private:

   int readLine(FILE* file, char* line, int length) {
      int count, result;
      char c;
      for (count = 0; count < length; count++) {
         result = fscanf(file, "%c", &c);
         if (result == EOF) {
            break;
         }
         if (c == '\r') {
            count--;
            continue;
         }
         if (c == '\n') {
            break;
         }
         line[count] = c;
      }
      line[count] = 0;
      return count;
   }

   void nextToken(char* str, char* token, int& pointer) {
      int count = 0;
      while (str[pointer] != '\t' && str[pointer] != 0) {
         token[count] = str[pointer];
         count++;
         pointer++;
      }
      token[count] = 0;
      if (str[pointer] == '\t') {
         pointer++;
      }
   }

   bool validateNetwork() {
      unsigned int reactionCount = reactions.size();
      unsigned int metaboliteCount = metabolites.size();
      bool valid = (reactionCount == reversible.size());
      valid = valid && (metaboliteCount == external.size());
      valid = valid && (metaboliteCount == s.size());
      for (unsigned int m = 0; valid && (m < metaboliteCount); m++) {
         valid = (reactionCount == s[m].size());
      }
      return valid;
   }

   void splitReversibleReactions() {
      for (unsigned int r = 0; r < reversible.size(); r++) {
         if (reversible[r]) {
            string rxn(reactions[r].c_str());
            rxn.resize(rxn.size() + 2);
            rxn[rxn.size() - 2] = '_';
            rxn[rxn.size() - 1] = 'R';
            reactions.insert(reactions.begin() + r + 1, rxn);
            reversible.insert(reversible.begin() + r + 1, true);
            for (unsigned int m = 0; m < metabolites.size(); m++) {
               s[m].insert(s[m].begin() + r + 1, (s[m][r] != 0) ? -s[m][r] : 0);
            }
            reversiblePairs.push_back(r);
            r++;
         }
      }
   }
public:
   //Stoichiometric matrix of the network containing metabolite coefficients
   vector<vector<float> > s;
   //Names of metabolites in the network
   vector<string> metabolites;
   //Names of reactions in the network
   vector<string> reactions;
   //A boolean flag for metabolites. If true the metabolite is external.
   vector<bool> external;
   //A boolean flag for reactions. If true the reaction was reversible.
   vector<bool> reversible;
   //A list of reversible reaction pairs. Two consecutive elements represent a reversible reaction
   vector<int> reversiblePairs;

   Network() {
   }

   ~Network() {
   }

   bool readNetworkFile(const char* filename) {
      FILE* file = fopen(filename, "r");
      if (file == NULL) {
         fprintf(stderr, "Unable to open the following file for reading: %s\r\n", filename);
         return false;
      }
      char line [LINE_SIZE];
      char token [TOKEN_SIZE];
      int pointer;
      int lineNumber = 0;
      int tokenNumber;
      while (readLine(file, line, LINE_SIZE) > 0) {
         pointer = 0;
         tokenNumber = 0;
         if (lineNumber >= 2) {
            s.push_back(vector<float>());
         }
         while (line[pointer] != 0) {
            nextToken(line, token, pointer);
            switch (lineNumber) {
                  //Reaction names
               case 0:
                  if (tokenNumber >= 2) {
                     reactions.push_back(string(token));
                  }
                  break;
                  //Reaction reversibility
               case 1:
                  if (tokenNumber >= 2) {
                     reversible.push_back(strcmp(token, "1") == 0);
                  }
                  break;
                  //Metabolites
               default:
                  if (tokenNumber == 0) {
                     metabolites.push_back(string(token));
                  } else if (tokenNumber == 1) {
                     external.push_back(strcmp(token, "1") == 0);
                  } else {
                     s[lineNumber - 2].push_back(atof(token));
                  }
            }
            tokenNumber++;
         }
         lineNumber++;
      }
      fclose(file);
      if (!validateNetwork()) {
         fprintf(stderr, "The following file is not a valid network file : %s\r\n", filename);
         return false;
      }
      splitReversibleReactions();
      return true;
   }
};

#endif	/* NETWORK_H */

