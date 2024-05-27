/*! \file FileIO.h
 * \brief File input/output operations.
 **/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <regex>
#include <cxxabi.h>
#include "Structs.h"
#ifndef __FILEIO_H
#define __FILEIO_H

#define NPWV    29
#define NFREQ   8301

template <typename T, typename U>
void readEtaATM(T **eta_array, U *pwv_atm, U *freq_atm);

#endif

template <typename T, typename U>
void readEtaATM(T **eta_array, U *pwv_atm, U *freq_atm) {
    
    pwv_atm->start = 0.1;
    pwv_atm->step = 0.1;
    pwv_atm->num = NPWV;

    freq_atm->start = 70.e9;
    freq_atm->step = 0.1e9;
    freq_atm->num = NFREQ;

    std::regex target("include\/FileIO\.h");
    std::string rel_loc = "resources/trans_data.dat";
    std::string abs_loc = std::regex_replace(__FILE__, target, rel_loc);

    *eta_array = new T[NPWV * NFREQ];
    
    std::string store;
    int status;
    //std::cout << abi::__cxa_demangle(typeid(store).name(), NULL, NULL, &status) << std::endl;

    std::ifstream myfile(abs_loc);
    std::string line;

    int line_nr = 0;
    int idx = 0;
    
    if (!myfile) std::cerr << "Could not open the file!" << std::endl;
    else {
        while(std::getline(myfile, line)) {
            std::istringstream iss(line);
            if(!line_nr) {
                line_nr++;
                continue;
            } 
            
            while(std::getline(iss, store, ' ')) {
                //std::cout << store << "test if space" << std::endl;
                if(!idx) {
                    idx++;
                    continue;
                } else if (store=="") {continue;}
                //std::cout << abi::__cxa_demangle(typeid(store).name(), NULL, NULL, &status) << std::endl;
                //iss.ignore();
                //std::cout << NFREQ * (idx-1) + (line_nr - 1) << " " << idx << " " << NFREQ * NPWV << std::endl; 
                (*eta_array)[NFREQ * (idx-1) + (line_nr - 1)] = std::stof(store);
                idx++;
                
            }
            line_nr++;
            idx = 0;
        }
        myfile.close();
        //std::cout << "here" << std::endl;
    }
}
