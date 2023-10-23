#include "LoadRootFile.h"

LoadRootFile::LoadRootFile():Tool(){}


bool LoadRootFile::Initialise(std::string configfile, DataModel &data)
{
    if(configfile!="")  m_variables.Initialise(configfile);
    //m_variables.Print();

    m_data= &data;
    m_log= m_data->Log;

    if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
    if(!m_variables.Get("Path_In",m_data->TD.Path_In)) m_data->TD.Path_In="./";

    if(!m_variables.Get("Treename",Treename)) Treename="Nan";
    if(!m_variables.Get("Branchname",Branchname)) Branchname="Nan";

    file.open(m_data->TD.Path_In,ios_base::in); // Open the file
    if(!file.is_open()) 
    {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    return true;
}


bool LoadRootFile::Execute()
{
    if(m_data->TD.Stop==true)
    {
        std::cout << "Stopping loop now!" << std::endl;
        m_data->vars.Set("StopLoop",1);
        m_data->TD.Stop=false;
        m_data->TD.EndOfRun=false;
        return true;
    }

    m_data->TD.TTS_Map.clear();

    std::string current_file; 
    std::getline(file,current_file); 
    if(!m_variables.Get("Path_Out",m_data->TD.Path_Out)) m_data->TD.Path_Out=current_file;
    current_file += "Analysis.root";

    TFile* TTfile = TFile::Open(current_file.c_str());
    if(!TTfile || TTfile->IsZombie()) 
    {
        std::cerr << "Error opening file: " << current_file << std::endl;
        m_data->TD.Stop=true;
        return true;
    }

    TTree* current_tree;
    TTfile->GetObject(Treename.c_str(),current_tree);
    if(!current_tree) 
    {
        std::cerr << "Error accessing tree in file." << std::endl;
        TTfile->Close();
        return true;
    }      

    for(int i_ch=0; i_ch<MAX_NUM_CHANNEL; i_ch++)
    {
        std::string Branchname_current = Branchname + to_string(i_ch);
        std::cout<<"I want to read "<<Branchname_current<<std::endl;
        TBranch* branch = current_tree->GetBranch(Branchname_current.c_str());
        if (!branch) 
        {
            std::cerr << "Error accessing branch in tree." << std::endl;
            TTfile->Close();
            return false;
        }

        std::vector<float>* data = nullptr;
        current_tree->SetBranchAddress(Branchname_current.c_str(), &data);

        std::vector<float> dataVector;
        for (Long64_t i = 0; i < current_tree->GetEntries(); i++) 
        {
            std::cout<<"I want to read "<<i<<" from a total of "<<current_tree->GetEntries()<<std::endl;
            current_tree->GetEntry(i);
            if(data) 
            {
                std::cout<<data->size()<<std::endl;
                dataVector.push_back((*data)[0]);
            }else {
                std::cerr << "Error reading data from branch " << Branchname_current << std::endl;
                TTfile->Close();
                return false;
            }
        }
        std::cout<<"For "<<i_ch<<" "<<dataVector.size()<<" was written"<<std::endl;
        m_data->TD.TTS_Map.insert(std::pair<int,vector<float>>(i_ch,dataVector));
    }

    m_data->TD.EndOfRun=true;
    TTfile->Close();
    
    return true;
}


bool LoadRootFile::Finalise()
{
    file.close();
    return true;
}
