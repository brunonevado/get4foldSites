
#include "common.h"


std::vector <char> fromIUPAC(char in){
    std::vector <char> toreturn;
    
    if( in == 'a' || in == 'c' || in == 'g' || in == 't' || in == 'n' || in == '-' ){
        toreturn.push_back(in);
        toreturn.push_back(in);
        return toreturn;
    }
    else if ( in == 'm' ){
        toreturn.push_back('a');
        toreturn.push_back('c');
        return toreturn;
    }
    else if ( in == 'r' ){
        toreturn.push_back('a');
        toreturn.push_back('g');
        return toreturn;
    }
    else if ( in == 'w' ){
        toreturn.push_back('a');
        toreturn.push_back('t');
        return toreturn;
    }
    else if ( in == 's' ){
        toreturn.push_back('c');
        toreturn.push_back('g');
        return toreturn;
    }
    else if ( in == 'k' ){
        toreturn.push_back('g');
        toreturn.push_back('t');
        return toreturn;
    }
    else if ( in == 'y' ){
        toreturn.push_back('c');
        toreturn.push_back('t');
        return toreturn;
    }
    else{
        std::cerr << "ERROR (fromIUPAC): invalid code: " << in << "\n";
        exit(1);
    }
    
}

