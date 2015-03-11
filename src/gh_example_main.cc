#ifndef __CINT__

#include "GHExample.h"
#include <time.h>

using namespace std;

/**
 * @brief the main routine
 * @param argc number of parameters
 * @param argv the parameters as strings
 * @return exit code
 */
int main(int argc, char *argv[])
{

    clock_t start, end;
    start = clock();

    // Create instance of analysis class
    GHExample* analysis = new GHExample;

    // Perform basic configuration
    if(!analysis->BaseConfig(argc, argv, "GoAT", "Physics"))
    {
        system("man ./documents/goat.man");
        return 0;
    }

    // Perform full initialisation
    if(!analysis->Init())
    {
        cout << "ERROR: Init failed!" << endl;
        return 0;
    }

    // Run over files
    analysis->TraverseFiles();

    end = clock();
    cout << endl;
    cout << "Time required for execution: "
    << (Double_t)(end-start)/CLOCKS_PER_SEC
    << " seconds." << "\n\n";

    if(analysis) delete analysis;

    return 0;
}

#endif
