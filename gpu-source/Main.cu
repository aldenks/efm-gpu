#include "Setup.h"
#include "EFMGenerator.h"

void execute(const char* file) {
   if (!network.readNetworkFile(file)) {
      fprintf(stderr, "Error loading network file\n");
      return;
   }
   if (!setup()) {
      return;
   }
   generateEFMs();
   freeResources();
}

int main(int argc, char** argv) {
   if (argc == 2) {
      execute(argv[1]);
   } else {
      printf("Please specify the network file.\n");
   }
   return (EXIT_SUCCESS);
}
