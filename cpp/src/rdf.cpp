/**
 * 
 * */

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>

using std::cout, std::endl;


int main (int argc, char** argv)
{
    // defining some default variables
    // parsing command line input
    char xyzfile;
    char outfile;
    float dr;

    for (int i=0; i < argc; i++){
        cout<<i<<" "<<argv[i]<<"\n";
    }
    int i = 1;

    while (i < argc) {
        if (strcmp(argv[i],"-xyz")) {
            xyzfile = *argv[i+1];
            cout<<argv[i]<<argv[i+1]<<"\n";
            cout<<xyzfile;
            i += 2;
        } 
        else if (strcmp(argv[i],"-o")) {
            outfile = *argv[i+1];
            cout<<argv[i]<<argv[i+1]<<"\n";
            i += 2;

        }
        else if (strcmp(argv[i],"-dr")) {
            dr = atof(argv[i+1]);
            cout<<argv[i]<<argv[i+1]<<"\n";
            i +=2;

        }
        else if (strcmp(argv[i],"-box")) {
            float box = atof(argv[i+1]);
            cout<<argv[i]<<argv[i+1]<<"\n";
            i +=2;
        }
        else {
            continue;
        }
    }

    cout<<xyzfile<<"\n";
    cout<<outfile<<"\n";
    cout<<dr<<"\n";



}