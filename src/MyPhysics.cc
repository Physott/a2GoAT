#ifndef __CINT__

#include "MyPhysics.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include <time.h>
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

	if(argc == 1)
	{
		cout << "Please provide a config file" << endl;
		return 0;
	}
	else if(argc == 2) configfile = argv[1];
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

    // Create instance of MyPhysics class
    MyPhysics* petap = new MyPhysics;

	// If unset, scan server or config file for file settings
	if(dir_in.length() == 0)
	{
        flag = petap->ReadConfig("Input-Directory",0,(Char_t*)serverfile.c_str());
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) dir_in = flag;
	}
	
	if(dir_out.length() == 0)
	{	
        flag = petap->ReadConfig("Output-Directory",0,(Char_t*)serverfile.c_str());
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) dir_out = flag;
	}
	
	if(file_in.length() == 0)
	{	
        flag = petap->ReadConfig("Input-File",0,(Char_t*)serverfile.c_str());
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) file_in = flag;
	}
	
	if(file_out.length() == 0)
	{	
        flag = petap->ReadConfig("Output-File",0,(Char_t*)serverfile.c_str());
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) file_out = flag;
	}
	
	if(pre_in.length() == 0)
	{	
        flag = petap->ReadConfig("Input-Prefix",0,(Char_t*)serverfile.c_str());
		flag.erase(0,flag.find_first_not_of(" "));
		if(strcmp(flag.c_str(),"nokey") != 0) pre_in = flag;
	}
	
	if(pre_out.length() == 0)
	{	
        flag = petap->ReadConfig("Output-Prefix",0,(Char_t*)serverfile.c_str());
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
	else { pre_in = "GoAT"; 	cout << "Input prefix:     '" << pre_in << "' chosen by default" << endl; }
	if(pre_out.length() != 0)  	cout << "Output prefix:    '" << pre_out << "' chosen" << endl;	
	else { pre_out = "Eta"; 	cout << "Output prefix:    '" << pre_out << "' chosen by default" << endl; }
	cout << endl;
	
	// Perform full initialisation 
    /*if(!petap->Init(configfile.c_str()))
	{
        cout << "ERROR: MyPhysics Init failed!" << endl;
		return 0;
    }*/

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
            if(!petap->StartFile(file_in.c_str(), file_out.c_str())) cout << "ERROR: MyPhysics failed on file " << file_in << "!" << endl;
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

					files_found++;
                    // Run MyPhysics
                    if(!petap->StartFile(file_in.c_str(), file_out.c_str()))
                        cout << "ERROR: MyPhysics failed on file " << file_in << "!" << endl;

				}
			}
		}
	}
	if (files_found == 0)
	{
		cout << "ERROR: No GoAT files found!" << endl;
		return 0;
	}

	end = clock();
	cout << "Time required for execution: "
	<< (double)(end-start)/CLOCKS_PER_SEC
	<< " seconds." << "\n\n";

	return 0;
}

MyPhysics::MyPhysics()    :
    proton_eta(TString("eta")),
    hist_eta(TString("eta")),
    hist_eta_proton(TString("eta_proton")),
    proton_etap(TString("etap")),
    hist_etap(TString("etap"), kTRUE),
    hist_etap_proton(TString("etap_proton"), kTRUE)
{

}

MyPhysics::~MyPhysics()
{

}

Bool_t	MyPhysics::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    proton_eta.Clear();
    hist_eta.Clear();
    hist_eta_proton.Clear();
    proton_etap.Clear();
    hist_etap.Clear();
    hist_etap_proton.Clear();

    TraverseEntries(0, eta->GetNEntries());

    proton_eta.RandomSubtraction();
    hist_eta.RandomSubtraction();
    hist_eta_proton.RandomSubtraction();
    proton_etap.RandomSubtraction();
    hist_etap.RandomSubtraction();
    hist_etap_proton.RandomSubtraction();

    Write();
	return kTRUE;
}

void	MyPhysics::ProcessEvent()
{
    if((GetEventNumber() % 100000 == 0) && GetEventNumber()!=0) cout << "Event: "<< GetEventNumber() << " Total Etas found: " << hist_eta.GetNFound() << endl;

    if(eta->GetNParticles()>0)
    {
        if(protons->GetNParticles()>0)
        {
            if(proton_eta.ProcessEvent(eta->Particle(0), *protons, *tagger))
            {
                hist_eta_proton.ProcessEvent(*eta, *tagger);
                return;
            }
        }
        hist_eta.ProcessEvent(*eta, *tagger);
        return;
    }

    if(etap->GetNParticles()>0)
    {
        if(protons->GetNParticles()>0)
        {
            if(proton_etap.ProcessEvent(etap->Particle(0), *protons, *tagger))
            {
                hist_etap_proton.ProcessEvent(*etap, *tagger);
                return;
            }
        }
        hist_etap.ProcessEvent(*etap, *tagger);
        return;
    }
}


Bool_t 	MyPhysics::Write()
{
    file_out->cd();
    TDirectory* curDir  = gDirectory->GetDirectory("eta");
    if(!curDir)
    {
        file_out->cd();
        gDirectory->mkdir("eta");
        curDir  = file_out->GetDirectory("eta");
    }

    curDir->cd();
    curDir  = gDirectory->GetDirectory("NoProton");
    if(!curDir)
    {
        curDir  = file_out->GetDirectory("eta");
        curDir->cd();
        gDirectory->mkdir("NoProton");
        curDir  = curDir->GetDirectory("NoProton");
    }
    hist_eta.Write(*curDir);

    curDir->cd();
    curDir  = gDirectory->GetDirectory("WithProton");
    if(!curDir)
    {
        curDir  = file_out->GetDirectory("eta");
        curDir->cd();
        gDirectory->mkdir("WithProton");
        curDir  = curDir->GetDirectory("WithProton");
    }
    proton_eta.Write(*curDir);
    hist_eta_proton.Write(*curDir);



    file_out->cd();
    curDir  = gDirectory->GetDirectory("etap");
    if(!curDir)
    {
        file_out->cd();
        gDirectory->mkdir("etap");
        curDir  = file_out->GetDirectory("etap");
    }

    curDir->cd();
    curDir  = gDirectory->GetDirectory("NoProton");
    if(!curDir)
    {
        curDir  = file_out->GetDirectory("etap");
        curDir->cd();
        gDirectory->mkdir("NoProton");
        curDir  = curDir->GetDirectory("NoProton");
    }
    hist_etap.Write(*curDir);

    curDir->cd();
    curDir  = gDirectory->GetDirectory("WithProton");
    if(!curDir)
    {
        curDir  = file_out->GetDirectory("etap");
        curDir->cd();
        gDirectory->mkdir("WithProton");
        curDir  = curDir->GetDirectory("WithProton");
    }
    proton_etap.Write(*curDir);
    hist_etap_proton.Write(*curDir);



	return kTRUE;
}

#endif