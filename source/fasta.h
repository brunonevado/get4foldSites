
#ifndef __bnFasta__
#define __bnFasta__

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>

#include "common.h"

class fasta {
    
private:
    std::unordered_map<std::string, std::string> genetic_code;
    std::vector <std::string> matrix;
    std::vector <std::string> names;
    std::string infile;
    bool codon_is_4fold (const std::string & codon) const;
    
public:
    fasta( int num_inds, bool coding = false );
    int num_lines () const {return int ( matrix.size() );}
    int num_bases () const {return int ( matrix[0].size() );}
    int num_names () const {return int ( names.size() );}
    int is_aligned() const;
    std::string input_file () const {return infile;}
    void info_to_stdout() const;
    char show_base( int line0, int site0 ) const { return matrix.at(line0).at(site0);} ;
    std::string show_seq ( int n ) const { return matrix.at(n);};
    std::vector <std::string> show_names() const { return names;};
    std::string show_site0 (int isite) const;
    
    void set_num_inds ( unsigned int value ) { matrix.resize(value); }
    void set_infile (const std::string & title) { infile = title; }
    void read_fasta_file ( const std::string & infile ) ;
    void append_from_vector ( const std::string & s, const std::string & name ) { matrix.push_back(s) ; names.push_back(name); }
    void new_fasta_from_inds ( const fasta & infasta, const std::vector <int> & index0 );
    void resize_matrix( unsigned int start, unsigned int end );
    void set_names ( std::vector < std::string > in ) {  names = in; };
    void append_site (const std::string & site);
    
    void octo_from_hap ( const fasta &original );
    std::vector <int> get_4fold_sites0( bool exclude4fold ) const;
    void free_concatenate ( const fasta & al2, bool verbose = false );
    
    void write_random ( const std::string & outfile, const std::string & name, int len ) const;
    void write_to_file ( const std::string & outfile, int append = 0 ) const;
    
    
    
};

#endif
