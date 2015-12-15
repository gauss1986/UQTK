#include <string>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "Utilsave.h"

FILE* createfile(std::string name){
    //create files
    FILE* dump;
    if(!(dump=fopen(name.c_str(),"w"))){
        printf("Could not open file '%s'\n",name.c_str());
        exit(1);    
    }
    return dump;
}

void closefile(FILE* dump, std::string name){
    // close files
    if(fclose(dump)){
       printf("Could not close file '%s'\n",name.c_str());
       exit(1);    
    }
}
