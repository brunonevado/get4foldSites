
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <locale>

#include "fasta.h"

fasta::fasta ( int num_inds, bool coding ){
    names.reserve(num_inds);
    matrix.reserve(num_inds);
    

    if (coding){
        //only works with standard genetic_code.
        genetic_code.insert(std::make_pair("gct","A"));
        genetic_code.insert(std::make_pair("gcc","A"));
        genetic_code.insert(std::make_pair("gca","A"));
        genetic_code.insert(std::make_pair("gcg","A"));
        
        genetic_code.insert(std::make_pair("cgt","R"));
        genetic_code.insert(std::make_pair("cgc","R"));
        genetic_code.insert(std::make_pair("cga","R"));
        genetic_code.insert(std::make_pair("cgg","R"));
        genetic_code.insert(std::make_pair("aga","R"));
        genetic_code.insert(std::make_pair("agg","R"));
        
        genetic_code.insert(std::make_pair("aat","N"));
        genetic_code.insert(std::make_pair("aac","N"));
        
        genetic_code.insert(std::make_pair("gat","D"));
        genetic_code.insert(std::make_pair("gac","D"));
        
        genetic_code.insert(std::make_pair("tgt","C"));
        genetic_code.insert(std::make_pair("tgc","C"));
        
        genetic_code.insert(std::make_pair("caa","Q"));
        genetic_code.insert(std::make_pair("cag","Q"));
        
        genetic_code.insert(std::make_pair("gaa","E"));
        genetic_code.insert(std::make_pair("gag","E"));
        
        genetic_code.insert(std::make_pair("ggt","G"));
        genetic_code.insert(std::make_pair("ggc","G"));
        genetic_code.insert(std::make_pair("gga","G"));
        genetic_code.insert(std::make_pair("ggg","G"));
        
        genetic_code.insert(std::make_pair("cat","H"));
        genetic_code.insert(std::make_pair("cac","H"));
        
        genetic_code.insert(std::make_pair("att","I"));
        genetic_code.insert(std::make_pair("atc","I"));
        genetic_code.insert(std::make_pair("ata","I"));
        
        genetic_code.insert(std::make_pair("tta","L"));
        genetic_code.insert(std::make_pair("ttg","L"));
        genetic_code.insert(std::make_pair("ctt","L"));
        genetic_code.insert(std::make_pair("ctc","L"));
        genetic_code.insert(std::make_pair("cta","L"));
        genetic_code.insert(std::make_pair("ctg","L"));
        
        genetic_code.insert(std::make_pair("aaa","K"));
        genetic_code.insert(std::make_pair("aag","K"));
        
        genetic_code.insert(std::make_pair("atg","START"));
        
        genetic_code.insert(std::make_pair("ttt","F"));
        genetic_code.insert(std::make_pair("ttc","F"));
        
        genetic_code.insert(std::make_pair("cct","P"));
        genetic_code.insert(std::make_pair("ccc","P"));
        genetic_code.insert(std::make_pair("cca","P"));
        genetic_code.insert(std::make_pair("ccg","P"));
        
        genetic_code.insert(std::make_pair("tct","S"));
        genetic_code.insert(std::make_pair("tcc","S"));
        genetic_code.insert(std::make_pair("tca","S"));
        genetic_code.insert(std::make_pair("tcg","S"));
        genetic_code.insert(std::make_pair("agt","S"));
        genetic_code.insert(std::make_pair("agc","S"));
        
        genetic_code.insert(std::make_pair("act","T"));
        genetic_code.insert(std::make_pair("acc","T"));
        genetic_code.insert(std::make_pair("aca","T"));
        genetic_code.insert(std::make_pair("acg","T"));
        
        genetic_code.insert(std::make_pair("tgg","W"));
        
        genetic_code.insert(std::make_pair("tat","Y"));
        genetic_code.insert(std::make_pair("tac","Y"));
        
        genetic_code.insert(std::make_pair("gtt","V"));
        genetic_code.insert(std::make_pair("gtc","V"));
        genetic_code.insert(std::make_pair("gta","V"));
        genetic_code.insert(std::make_pair("gtg","V"));
        
        genetic_code.insert(std::make_pair("taa","STOP"));
        genetic_code.insert(std::make_pair("tga","STOP"));
        genetic_code.insert(std::make_pair("tag","STOP"));
    }
    
}


bool fasta::codon_is_4fold (const std::string & codon) const {
    // returns true if codon if 4-fold degenerate, otherwise false
    std::string codonA = codon.substr( 0, 2 ).append("a");
    std::string codonC = codon.substr( 0, 2 ).append("c");
    std::string codonG = codon.substr( 0, 2 ).append("g");
    std::string codonT = codon.substr( 0, 2 ).append("t");
    if(   genetic_code.at(codonA) != genetic_code.at(codon)
       || genetic_code.at(codonC) != genetic_code.at(codon)
       || genetic_code.at(codonG) != genetic_code.at(codon)
       || genetic_code.at(codonT) != genetic_code.at(codon)  ){
        return false;
    }
    return true;
}


void fasta::resize_matrix( unsigned int start, unsigned int end ){
    for( unsigned int line = 0; line < matrix.size(); line++ ){
        std::string sub = matrix.at(line).substr(start - 1, end - start + 1  );
        matrix.at(line) = sub;
    }
    
}


int fasta::is_aligned () const {
    int len0 = int( this->matrix.at(0).length() );
    for( int iline = 1; iline < this->num_lines(); iline++){
        if( this->matrix.at(iline).length() != len0 ){
            return (iline);
        }
    }
    return 0;
}

void fasta::info_to_stdout () const {
    std::cout << "Infile: " << this->input_file() << " (" << this->num_lines() << " sequences, " <<
    this->num_bases() << " bp long)\n";
}


void fasta::read_fasta_file( const std::string & infas){
    int cind = 0;
    matrix.clear();
    std::locale loc;
    std::string line;
    std::ifstream infile_fas (infas.c_str());
    if (infile_fas.is_open()){
        infile = infas;
        while ( ! infile_fas.eof() ) {
            getline(infile_fas, line);
            
            if( line[0] == '>' ){
                cind++;
                matrix.resize(cind);
                std::string name = line.substr(1);
                names.push_back(name);
            }
            else {
                matrix.at(cind-1).append( line );
            }
        }
    }
    else {
        std::cerr << "ERROR (read_fasta_file): Unable to open infile " << infas << "\n" ;
        exit(1);
    }
    infile_fas.close();
    // internal representation uses lower case for everything...
    for (unsigned int l = 0; l < matrix.size(); l++) {
        for (unsigned s = 0; s < matrix.at(l).length(); s++) {
            matrix.at(l).at(s) = tolower(matrix.at(l).at(s), loc);
        }
    }
}



void fasta::write_random( const std::string & out, const std::string & name, int len ) const {
    // Writes random sequence to file, probably not a very good way though
    std::ofstream outputFile;
    char bases[4] = { 'A', 'C', 'G', 'T' };
    outputFile.open(out.c_str());
    if( !outputFile.is_open() ){
        std::cerr << "ERROR (write_random): unable to open for output file " << out << "\n";
        exit(1);
    }
    
    outputFile << ">" << name << std::endl;
    
    for( int site = 0; site < len; site++){
        outputFile << bases[ rand() % 4 ];
    }
    outputFile << "\n";
    outputFile.close();
    
}


void fasta::write_to_file( const std::string & out, int append ) const {
    // Write alignment to file (fasta format)
    std::locale loc;
    std::ofstream outputFile;
    
    if( append == 0 ){
        outputFile.open(out.c_str());
    }
    else{
        outputFile.open(out.c_str(), std::ios::app );
    }
    if( !outputFile.is_open() ){
        std::cerr << "ERROR (write_to_file): unable to open for output file " << out << "\n";
        exit(1);
    }
    if( names.size() == matrix.size() ){
        for( unsigned int i = 0; i < matrix.size(); i++){
            outputFile << ">" << names.at(i) << std::endl;
            for (unsigned int site = 0; site < matrix.at(i).length(); site++) {
                outputFile << toupper(matrix.at(i).at(site),loc);
            }
            outputFile << std::endl;
        }
    }
    else{
        for( unsigned int i = 0; i < matrix.size(); i++){
            outputFile << ">i" << i+1 << std::endl;
            for (unsigned int site = 0; site < matrix.at(i).length(); site++) {
                outputFile << toupper(matrix.at(i).at(site),loc);
            }
        }
    }
    outputFile.close();
}

void fasta::octo_from_hap ( const fasta & original ){
    // assumes matrix of original is in frame and aligned
    // solves IUPAC codes into all possible combinations
    // resulting matrix will have 8x as many lines
    // only works for IUPAC "diploid" codes, R(AG) Y(CT) S(GC) W(AT) K(GT) M(AC)
    // e.g.:
    // AAA = AAA x 8
    // AAW = AAA x 4 + AAT x 4
    // AWW = AAA x 2 + AAT x 2 + ATA x 2 + ATT x 2
    // WWW = AAA + AAT + ATA + ATT + TAA + TAT + TTA + TTT
    // It does not add names: the resulting data is not intended to be written, but just to check
    // for 4 fold sites considering all possibilities with iupac codes
    matrix.resize(original.num_lines()*8);
    for (unsigned int iline = 0; iline < original.num_lines(); iline++) {
        for (unsigned int isite = 0; isite < original.num_bases(); isite+= 3 ) {
            std::vector <char> alleles1 = fromIUPAC( original.matrix.at(iline).at(isite) );
            std::vector <char> alleles2 = fromIUPAC( original.matrix.at(iline).at(isite+1) );
            std::vector <char> alleles3 = fromIUPAC( original.matrix.at(iline).at(isite+2) );
            
            this->matrix.at(iline).push_back(alleles1.at(0));
            this->matrix.at(iline).push_back(alleles2.at(0));
            this->matrix.at(iline).push_back(alleles3.at(0));
            
            this->matrix.at(iline + original.num_lines() ).push_back(alleles1.at(0));
            this->matrix.at(iline + original.num_lines() ).push_back(alleles2.at(0));
            this->matrix.at(iline + original.num_lines() ).push_back(alleles3.at(1));
            
            this->matrix.at(iline + 2*original.num_lines() ).push_back(alleles1.at(0));
            this->matrix.at(iline + 2*original.num_lines() ).push_back(alleles2.at(1));
            this->matrix.at(iline + 2*original.num_lines() ).push_back(alleles3.at(0));
            
            
            this->matrix.at(iline + 3*original.num_lines() ).push_back(alleles1.at(1));
            this->matrix.at(iline + 3*original.num_lines() ).push_back(alleles2.at(0));
            this->matrix.at(iline + 3*original.num_lines() ).push_back(alleles3.at(0));
            
            this->matrix.at(iline + 4*original.num_lines() ).push_back(alleles1.at(0));
            this->matrix.at(iline + 4*original.num_lines() ).push_back(alleles2.at(1));
            this->matrix.at(iline + 4*original.num_lines() ).push_back(alleles3.at(1));
            
            this->matrix.at(iline + 5*original.num_lines() ).push_back(alleles1.at(1));
            this->matrix.at(iline + 5*original.num_lines() ).push_back(alleles2.at(0));
            this->matrix.at(iline + 5*original.num_lines() ).push_back(alleles3.at(1));
            
            this->matrix.at(iline + 6*original.num_lines() ).push_back(alleles1.at(1));
            this->matrix.at(iline + 6*original.num_lines() ).push_back(alleles2.at(1));
            this->matrix.at(iline + 6*original.num_lines() ).push_back(alleles3.at(0));
            
            this->matrix.at(iline + 7*original.num_lines() ).push_back(alleles1.at(1));
            this->matrix.at(iline + 7*original.num_lines() ).push_back(alleles2.at(1));
            this->matrix.at(iline + 7*original.num_lines() ).push_back(alleles3.at(1));
        }
    }
}


void fasta::free_concatenate(const fasta &fasta2, bool verbose){
    // checks names instead of order
    // for sequences missing in fasta2 will add Ns (raises warning)
    // sequences present in fasta2 but not fasta1 are added with Ns irrespective of order
    // input fasta files should be aligned before running (does not check for this)
    
    unsigned int original_len = this->num_bases();
    for ( unsigned int iline = 0; iline < this->num_lines(); iline++ ) {
        bool found = false;
        for ( unsigned int iline2 = 0; iline2 < fasta2.num_lines(); iline2++ ) {
            if( this->names.at(iline) == fasta2.names.at(iline2)  ){
                found = true;
                this->matrix.at(iline).append(fasta2.matrix.at(iline2));
            }
        }
        if(!found){
            std::string seqN(fasta2.num_bases(), 'n');
            this->matrix.at(iline).append(seqN.c_str());
            if(verbose){
                std::clog << "<concatenateFasta2> WARNING: sequence " << this->names.at(iline) << " not found in file " << fasta2.infile << ", sequence included with Ns" << std::endl;
            }
        }
    }
    for ( unsigned int iline2 = 0; iline2 < fasta2.num_lines(); iline2++ ) {
        bool found = false;
        for( unsigned int iline = 0; iline < this->num_lines(); iline++ ) {
            if( this->names.at(iline) == fasta2.names.at(iline2)  ){
                found = true;
            }
        }
        if (!found){
            std::string missing( original_len, 'n');
            this->append_from_vector(missing.append( fasta2.matrix.at(iline2)), fasta2.names.at(iline2));
            if(verbose){
                std::clog << "<concatenateFasta2> WARNING: sequence " << fasta2.names.at(iline2) << " present in file " << fasta2.infile << " but not in " << this->infile << ". Sequence was included with Ns" << std::endl;
            }
        }
    }
}

std::vector <int> fasta::get_4fold_sites0( ) const {
    // returns vector of 0-index positions of all 4-fold sites
    // sequence should be aligned and in frame
    // skips any codons that have ambiguous/gaps
    std::vector <int> sites0;
    for (unsigned int isite = 0; isite < this->num_bases(); isite += 3) {
        std::string codon= "";
        bool addSite = false;  // skip if there are no fully resolved codons
        for (unsigned int iline = 0; iline < this->num_lines(); iline++) {
            codon =  matrix.at(iline).substr( isite, 3 );
            if(genetic_code.count(codon) == 0){
                continue;
            }
            if( codon_is_4fold(codon) ){
                addSite = true;
            }
            else{
                addSite = false;
                break;
            }
        }
        if(addSite){
            sites0.push_back(isite + 2);
        }
    }
    return sites0;
}


void fasta::new_fasta_from_inds ( const fasta & infasta, const std::vector <int> & index0 ){
    
    infile = infasta.infile;
    matrix.clear();
    names.clear();
    for (unsigned int line = 0; line < index0.size() ; line++ ) {
        matrix.push_back( infasta.matrix.at(index0.at(line)));
        names.push_back( infasta.names.at(index0.at(line)));
    }
}

std::string fasta::show_site0 (int isite) const {
    std::string site;
    for(int i = 0; i < this->num_lines(); i++ ){
        site.push_back(this->show_base(i, isite));
    }
    return site;
}

void fasta::append_site (const std::string & site){
    if(site.length() != this->num_lines() ){
        std::cerr << "ERROR (append_site): number of bases in site to append (" << site.length() << ") is not the same as number of sequences in data ("
        << this->num_lines() << ")." << std::endl;
    }
    for(int i = 0; i < site.size(); i++){
        this->matrix.at(i).push_back(site.at(i));
    }
}
