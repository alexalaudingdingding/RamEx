// Aug.21, 2021
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Last update time: Nov 9, 2021
// Updated by Gongchao Jing

#ifndef _UTILITY_H
#define _UTILITY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include <stdlib.h>
#include <string.h>
#include <sys/dir.h>
#include <sys/stat.h>
#include "version.h"
#define BUFFER_SIZE 5000


using namespace std;

string Check_Env(){
    
    if (getenv("RamEx") == NULL){
                               
                               cerr << "Error: Please set the environment variable \"RamEx\" to the directory" << endl;
                               exit(0);
                               
                                   }
    
    string path =  getenv("RamEx");
    return path;
    
    //debug
    //return "/opt/tools/RamEX/";
    }

int Check_Path(const char * path, int type){
    
    if (strlen(path) < 1) return 0;
    
    DIR *pDir = opendir(path);
            
    if(pDir!=NULL){
                  closedir(pDir);                   
                  if (type == 0){
                     string command = "rm -rf ";
                     command += path;
                     system(command.c_str());
                     mkdir(path, 0755);                          
                     }
                  }                  
    else 
         mkdir(path, 0755);           
    return 0;
    
    }

bool Check_Path(const char * path){
    
    if (strlen(path) < 1) return false;
    
    DIR *pDir = opendir(path);
    
    if (pDir != NULL){
             closedir(pDir);
             return true;
             }
    
    return false;
    }

bool Check_File(const char * file){
    
    fstream infile(file, ifstream::in);
    
    if (!infile){
                 
                 cerr << "Error: Cannot open file: " << file << endl;
                 return false;
                 }
    
    infile.close();
    infile.clear();
    
    return true;
    
    }
void Run_With_Error(char * command, const char * error){
     
     string command_with_error;
     command_with_error = command;
     command_with_error +=" 2>>";
     command_with_error += error;
     system(command_with_error.c_str());
     }
#endif
