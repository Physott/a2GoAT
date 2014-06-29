#ifndef __CINT__

#include "MyEtap.h"
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

    // Create instance of MyEtap class
    MyEtap* petap = new MyEtap;

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
    if(!petap->Init(configfile.c_str()))
	{
        cout << "ERROR: MyEtap Init failed!" << endl;
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
            if(!petap->StartFile(file_in.c_str(), file_out.c_str())) cout << "ERROR: MyEtap failed on file " << file_in << "!" << endl;
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
                    // Run MyEtap
                    if(!petap->StartFile(file_in.c_str(), file_out.c_str()))
                        cout << "ERROR: MyEtap failed on file " << file_in << "!" << endl;

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

MyEtap::MyEtap()    :
    time_raw(0),
    raw(TString("raw")),
    time_cutIM(0),
    cutIMevent(TString("cutIM")),
    time_cutMM(0),
    cutMMevent(TString("cutMM")),
    N_eta(0)
{
    cutIM[0][0] = 110;
    cutIM[0][1] = 155;
    cutIM[1][0] = 110;
    cutIM[1][1] = 155;
    cutIM[2][0] = 110;
    cutIM[2][1] = 155;

    cutMM[0] = 850;
    cutMM[1] = 1050;
}

MyEtap::~MyEtap()
{
    if(time_raw)       delete time_raw;
    if(time_cutIM)     delete time_cutIM;
    if(time_cutMM)     delete time_cutMM;
}

Bool_t	MyEtap::Init(const char* configfile)
{
    gROOT->cd();

    time_raw		= new TH1D("time_raw",		"time_raw",		1000,-500,500);
    time_cutIM		= new TH1D("time_cutIM",		"time_cutIM",		1000,-500,500);
    time_cutMM		= new TH1D("time_cutMM",		"time_cutMM",		1000,-500,500);

    raw.SetCuts(-10, 5, -515, -15, 15, 510);
    cutIMevent.SetCuts(-10, 5, -515, -15, 15, 510);
    cutMMevent.SetCuts(-10, 5, -515, -15, 15, 510);

    return kTRUE;
}

Bool_t	MyEtap::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseEntries(0, eta->GetNEntries());

    raw.RandomSubtraction();
    cutIMevent.RandomSubtraction();
    cutMMevent.RandomSubtraction();

    WriteHistograms();
	return kTRUE;
}

void	MyEtap::ProcessEvent()
{
    if(GetEventNumber() == 0) N_eta = 0;
    else if(GetEventNumber() % 100000 == 0) cout << "Event: "<< GetEventNumber() << " Total Etas found: " << N_eta << endl;

    Double_t    imSub[3];
    Double_t    misMass;
    bool        passIM;

    if(eta->GetNParticles()>0)
    {
        imSub[0]    = (eta->SubParticles(0, 0)+eta->SubParticles(0, 1)).M();
        imSub[1]    = (eta->SubParticles(0, 2)+eta->SubParticles(0, 3)).M();
        imSub[2]    = (eta->SubParticles(0, 4)+eta->SubParticles(0, 5)).M();
        if((imSub[0]>cutIM[0][0] && imSub[0]<cutIM[0][1]) && (imSub[1]>cutIM[1][0] && imSub[1]<cutIM[1][1]) && (imSub[2]>cutIM[2][0] && imSub[2]<cutIM[2][1]))
            passIM  = true;
        else
            passIM  = false;

        for(int i=0; i<tagger->GetNTagged(); i++)
        {
            time_raw->Fill(tagger->GetTagged_t(i));
            misMass = (tagger->GetVector(i)+TLorentzVector(0,0,0,MASS_PROTON) - eta->Meson(0)).M();
            raw.Fill(tagger->GetTagged_t(i), eta->Meson(0).M(), misMass);
            raw.FillSubMesons(tagger->GetTagged_t(i), imSub[0], imSub[1], imSub[2]);

            if(passIM)
            {
                time_cutIM->Fill(tagger->GetTagged_t(i));
                cutIMevent.Fill(tagger->GetTagged_t(i), eta->Meson(0).M(), (tagger->GetVector(i)+TLorentzVector(0,0,0,MASS_PROTON) - eta->Meson(0)).M());
                cutIMevent.FillSubMesons(tagger->GetTagged_t(i), imSub[0], imSub[1], imSub[2]);

                if(misMass>cutMM[0] && misMass<cutMM[1])
                {
                    time_cutMM->Fill(tagger->GetTagged_t(i));
                    cutMMevent.Fill(tagger->GetTagged_t(i), eta->Meson(0).M(), (tagger->GetVector(i)+TLorentzVector(0,0,0,MASS_PROTON) - eta->Meson(0)).M());
                    cutMMevent.FillSubMesons(tagger->GetTagged_t(i), imSub[0], imSub[1], imSub[2]);
                }
            }
        }
        N_eta++;
    }
}


Bool_t 	MyEtap::WriteHistograms()
{
	cout << "Writing histograms." << endl;

    if(!file_out) return kFALSE;

    file_out->cd();
    time_raw->Write();
    raw.Write(file_out);

    file_out->cd();
    time_cutIM->Write();
    cutIMevent.Write(file_out);

    file_out->cd();
    time_cutMM->Write();
    cutMMevent.Write(file_out);

	return kTRUE;
}

#endif
