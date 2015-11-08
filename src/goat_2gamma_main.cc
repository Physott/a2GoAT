#ifndef __CINT__

#include "GoAT_2gamma.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
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

    // Initialise strings
    std::string configfile = "";
    std::string serverfile = "";
    std::string dir_in = "";
    std::string dir_out = "";
    std::string file_in = "";
    std::string file_out = "";
    std::string pre_in = "";
    std::string pre_out = "";

    Int_t length;
    std::string flag;

    Bool_t overwrite = kTRUE;

    if(argc == 1)
    {
        system("man ./documents/goat.man");
        return 0;
    }
    else if(argc == 2)
    {
        flag = argv[1];
        flag.erase(0,1);
        if(strcmp(flag.c_str(), "help") == 0) system("man ./documents/goat.man");
        else configfile= argv[1];
    }
    else
    {
        for(int i=1; i<argc; i++)
        {
            flag = argv[i];
            if(flag.find_first_of("-") == 0)
            {
                i++;
                flag.erase(0,1);
                if(strcmp(flag.c_str(), "s") == 0) serverfile = argv[i];
                else if(strcmp(flag.c_str(), "d") == 0) dir_in = argv[i];
                else if(strcmp(flag.c_str(), "D") == 0) dir_out = argv[i];
                else if(strcmp(flag.c_str(), "f") == 0) file_in = argv[i];
                else if(strcmp(flag.c_str(), "F") == 0) file_out = argv[i];
                else if(strcmp(flag.c_str(), "p") == 0) pre_in = argv[i];
                else if(strcmp(flag.c_str(), "P") == 0) pre_out = argv[i];
                else if(strcmp(flag.c_str(), "n") == 0)
                {
                        overwrite = kFALSE;
                    i--;
                }
                else
                {
                    cout << "Unknown flag " << flag << endl;
                    return 0;
                }
            }
            else configfile = argv[i];
        }
    }

    // Check that config file exists:
    ifstream cfile(configfile.c_str());
    if(!cfile)
    {
        cout << "Config file '" << configfile << "' could not be found." << endl;
        return 0;
    }

    // If server file is specified, check that it exists
    if(serverfile.length() > 0)
    {
        // Check that file exists:
        ifstream sfile(serverfile.c_str());
        if(!sfile)
        {
            cout << "Server file '" << serverfile << "' could not be found" << endl;
            return 0;
        }
    }
    // If no server file is specified, allow for checking in the config file
    else serverfile = configfile;

    // Create instance of GoAT_Detailed class
    GTreeManager* goat = new GoAT_2gamma;

    // If unset, scan server or config file for file settings
    if(dir_in.length() == 0)
    {
        flag = goat->ReadConfig("Input-Directory",0,(Char_t*)serverfile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) dir_in = flag;
    }

    if(dir_out.length() == 0)
    {
        flag = goat->ReadConfig("Output-Directory",0,(Char_t*)serverfile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) dir_out = flag;
    }

    if(file_in.length() == 0)
    {
        flag = goat->ReadConfig("Input-File",0,(Char_t*)serverfile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) file_in = flag;
    }

    if(file_out.length() == 0)
    {
        flag = goat->ReadConfig("Output-File",0,(Char_t*)serverfile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) file_out = flag;
    }

    if(pre_in.length() == 0)
    {
        flag = goat->ReadConfig("Input-Prefix",0,(Char_t*)serverfile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) pre_in = flag;
    }

    if(pre_out.length() == 0)
    {
        flag = goat->ReadConfig("Output-Prefix",0,(Char_t*)serverfile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) pre_out = flag;
    }
    // Finished scanning for file settings

    // Fix directories to include final slash if not there
    if(dir_in.find_last_of("/") != (dir_in.length()-1)) dir_in += "/";
    if(dir_out.find_last_of("/") != (dir_out.length()-1)) dir_out += "/";

    // Output user settings (Set to defaults if still unspecified)
    cout << endl << "User inputs" << endl;
    cout << "Config file:      '" << configfile << "' chosen" << endl;
    if(dir_in.length() != 0)  	cout << "Input directory:  '" << dir_in << "' chosen" << endl;
    if(dir_out.length() != 0)  	cout << "Output directory: '" << dir_out << "' chosen" << endl;
    if(file_in.length() != 0)  	cout << "Input file:       '" << file_in << "' chosen" << endl;
    if(file_out.length() != 0) 	cout << "Output file:      '" << file_out << "' chosen" << endl;
    if(pre_in.length() != 0)  	cout << "Input prefix:     '" << pre_in << "' chosen" << endl;
    else { pre_in = "Acqu"; 	cout << "Input prefix:     '" << pre_in << "' chosen by default" << endl; }
    if(pre_out.length() != 0)  	cout << "Output prefix:    '" << pre_out << "' chosen" << endl;
    else { pre_out = "GoAT"; 	cout << "Output prefix:    '" << pre_out << "' chosen by default" << endl; }
    cout << endl;

    // Perform full initialisation
    goat->SetConfigFile(configfile.c_str());
    if(!goat->Init())
    {
        cout << "ERROR: GoAT Init failed!" << endl;
        return 0;
    }

    std::string file;
    std::string prefix;
    std::string suffix;

    Int_t files_found = 0;
    // If input file is specified, use it
    if(file_in.length() > 0)
    {
        cout << "Searching for input file(s)" << endl;
        file = file_in;
        length = file.length();
        // File should at least have '.root' at the end
        if(length >= 5)
        {
            // Add input directory to it
            file_in = dir_in+file_in;
            cout << "Input file  '" << file_in << "' chosen" << endl;

            // If output file is specified, use it
            if(file_out.length() > 0) file_out = dir_out+file_out;
            // If output file is not specified, build it
            else
            {
                // If output directory is not specified, build it
                if(dir_out.length() == 0)
                {
                    prefix = file.substr(0,file.find_last_of("/")+1);
                    dir_out = dir_in+prefix;
                }
                // If input prefix doesn't match, simply prepend output prefix to the file name
                if(file.find(pre_in)>file.length()) suffix = ("_"+file.substr(file.find_last_of("/")+1,length-(file.find_last_of("/")+1)));
                // If input prefix does match, switch prefixes
                else suffix = file.substr(file.find_last_of("/")+1+pre_in.length(),length-(file.find_last_of("/")+1+pre_in.length()));
                // Build output file name
                file_out = dir_out+pre_out+suffix;
            }

            cout << "Output file '" << file_out << "' chosen" << endl << endl;
            if(!goat->StartFile(file_in.c_str(), file_out.c_str())) cout << "ERROR: GoAT failed on file " << file_in << "!" << endl;
            files_found++;
        }
    }
    // Otherwise scan input directory for matching files
    else
    {
        cout << "Searching input directory for files matching input prefix" << endl;
        cout << "Input prefix  '" << pre_in << "' chosen" << endl;
        cout << "Output prefix '" << pre_out << "' chosen" << endl;

        // If output directory is not specified, use the input directory
        if(dir_in.length()  == 0) dir_in = "./";
        if(dir_out.length() == 0) dir_out = dir_in;

        // Create list of files in input directory
        TSystemFile *sys_file;
        TSystemDirectory *sys_dir = new TSystemDirectory("files",dir_in.c_str());
        TList *file_list = sys_dir->GetListOfFiles();
        file_list->Sort();
        TIter file_iter(file_list);

        // Iterate over files
        while((sys_file=(TSystemFile*)file_iter()))
        {
            file = sys_file->GetName();
            length = file.length();
            // File should at least have '.root' at the end
            if(length >= 5)
            {
                //Check that prefixes and suffixes match
                prefix = file.substr(0,pre_in.length());
                suffix = file.substr(length-5,5);
                if(((strcmp(prefix.c_str(),pre_in.c_str()) == 0)) && (strcmp(suffix.c_str(),".root") == 0))
                {
                    // Build input file name
                    file_in = dir_in+file;
                    // Build output file name
                    suffix = file.substr(pre_in.length(),length-pre_in.length());
                    file_out = dir_out+pre_out+suffix;

                    // Check for previously created output file
                    if((gSystem->IsFileInIncludePath(file_out.c_str())) && !overwrite) continue;

                    files_found++;

                    // Run GoAT
                    if(!goat->StartFile(file_in.c_str(), file_out.c_str()))
                        cout << "ERROR: GoAT failed on file " << file_in << "!" << endl;
                }
            }
        }
    }
    if (files_found == 0)
    {
        cout << "ERROR: No AcquRoot files found!" << endl;
        return 0;
    }

    end = clock();
    cout << endl;
    cout << "Time required for execution: "
    << (double)(end-start)/CLOCKS_PER_SEC
    << " seconds." << "\n\n";

    if(goat)    delete goat;

    return 0;
}



#endif
