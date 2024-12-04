# concatenateFasta

get4foldSites:  returns 4-fold degenerate sites from fasta files.  
  
Author: B. Nevado  
  
Usage:  
Usage: get4foldSites -infile in.txt -outfile res.fas -iupac 0/1 -verbose 0/1  
    -infile: text file with list of input fasta files to process (1 file per line).  
    -outfile: output fasta file to write to.  
    -iupac: if set to 1 will consider positions with 'diploid' IUPAC codes. This may cause ambiguity for non-4-fold sites.    
    -exclude: if set to 1, will instead output sites that are NOT 4-fold degenerated.  
  
Output: 1 fasta file with all 4-fold sites (or all but 4-fold sites) found across all input files.  
    
Notes:  
    Tested on linux with g++ v5. Earlier version (4.8) did not work.  
    Set '-verbose 0' to turn off warnings.  
      
Installation (Linux):  
git clone https://github.com/brunonevado/get4foldSites.git  
cd get4foldSites  
make  
./get4foldSites  
  
Folders:  
/source: source code  
/examples: test/example files, run with get4foldSites -infile infiles.txt -outfile out.4fold.fas -iupac 1 -verbose 1  
