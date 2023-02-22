
#include <iostream>
#include <fstream>
#include <sstream>

#include "fasta.h"
#include "common.h"
#include "args.h"


void help(){
    
    std::cout << "#######################\nHelp for get4foldSites 13062022 \n#######################\n";
    std::cout << "Usage: get4foldSites -infile in.txt -outfile res.fas -iupac 0/1 -verbose 0/1" << std::endl;
    std::cout << "infile should contain full path of files to screen for 4 fold sites (1 file per line)." << std::endl;
    std::cout << "outfile will contain all 4-fold degenerate sites found." << std::endl;
    std::cout << "if -iupac 1, will consider positions with 'diploid' IUPAC codes." << std::endl;
    std::cout << "Note: tested on linux with g++ v5. Earlier version (4.8) did not work." << std::endl;
}



int main(int argc, const char * argv[]){
   

    // READ IN OPTIONS FOR THE RUN
    sargs myargs;
    try{
        myargs = args::getargs(argc, argv, std::vector<std::string> {"infile", "outfile"},
                               std::vector<std::string> {"iupac","verbose"}, std::vector<std::string>  {}, std::string {}, std::string {}); }
    catch (std::string e){
        std::cout << " Args failed: " << e << std::endl;
        help();
        exit(1);
    }
    
    std::string aInfile = myargs.args_string.at(0);
    std::string aOutfile =   myargs.args_string.at(1);
    bool aIupac = myargs.args_booleans.at(0);
    bool verbose = myargs.args_booleans.at(1);
    
    // GET NAMES OF FILES TO SCREEN
    std::string cline;
    std::vector<std::string> vFastaFiles;
    
    std::ifstream infile(aInfile);
    if(!infile.is_open()){
        std::cerr <<"<get4foldSites> ERROR: unable to open fasta list infile " << aInfile << std::endl;
        exit(1);
    }
    while(getline(infile, cline)){
        vFastaFiles.push_back(cline);
    }
    std::cout << "<get4foldSite> Read " << vFastaFiles.size() << " fasta file names from " << aInfile << std::endl;
    
    // GET ALL 4-FOLD SITES
    
    fasta fasta4Fold(1); // will hold resulting concatenated fasta with all 4-fold sites
    
    for (int ifile = 0; ifile < vFastaFiles.size(); ifile++ ) {
        fasta afasta(1,true);
        afasta.read_fasta_file(vFastaFiles.at(ifile));
        if( afasta.is_aligned() != 0 ){
            std::cerr << "<<get4foldSites> ERROR: file " << vFastaFiles.at(ifile) << " is not aligned, offending seq " << afasta.is_aligned() << std::endl;
            exit(1);
        }
        std::vector <int> vSitesToAdd;
        
        if(aIupac){
            fasta afasta_iupac(1, true);
            afasta_iupac.octo_from_hap(afasta);
            vSitesToAdd = afasta_iupac.get_4fold_sites0();
        }
        else{
            vSitesToAdd = afasta.get_4fold_sites0();
        }
        fasta afastaTemp4Fold(1);
        afastaTemp4Fold.set_num_inds(afasta.num_lines());
        afastaTemp4Fold.set_names(afasta.show_names());
        for(auto isite : vSitesToAdd){
            afastaTemp4Fold.append_site(afasta.show_site0(isite));
        }
        if(verbose){
            std::clog << "<get4foldSites> Adding " << vSitesToAdd.size()  <<" 4-fold degenerate sites from file " << vFastaFiles.at(ifile) << std::endl;
        }

        if(ifile == 0){
            fasta4Fold = afastaTemp4Fold;
        }
        else{
            fasta4Fold.free_concatenate(afastaTemp4Fold, false);
        }
    }
     fasta4Fold.write_to_file(aOutfile);
    
    std::clog << "<get4foldSites> Finished, read " << vFastaFiles.size() << " fasta files, printed "
    << fasta4Fold.num_bases() << " 4-fold degenerate sites to file " << aOutfile  <<  std::endl;
    
    return 0;
}

