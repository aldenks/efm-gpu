#include "Setup.h"
#include "EFMGenerator.h"
#include <conio.h>

void execute(const char* file) {
   printf("Reading network file : %s\n", file);
   if (!network.readNetworkFile(file)) {
      fprintf(stderr, "Error loading network file\n");
      return;
   }
   printf("Network file read successfully\n");
   //network.print();
   printf("Setting up memory\n");
   if (!setup()) {
      printf("Memory setup not successful\n");
      return;
   }
   printf("Memory setup successful\n");
   printf("Generating EFMs\n");
   generateEFMs();
   printf("Freeing resources\n");
   freeResources();
}

int main(int argc, char** argv) {
   //if (argc == 2) {
   //   execute(argv[1]);
   //} else {
   //   printf("Please specify the network file.\n");
   //}
   execute("E:\\Developer\\GPU\\EFM\\source\\ecoli-irrev-compressed.xls");
   printf("Press any key to exit . . .\n");
   getch();
   return (EXIT_SUCCESS);
}
