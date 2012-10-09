//AUTHOR: Helge Roider
//INSTITUTION: Max Planck Institute for Molecular Genetics - Berlin
//DATE: 16/07/2006
// $Id: TRAP.cpp 339 2012-09-25 09:16:14Z jbao $

//INPUT: 1. TRANSFAC MATRIX, 2. FASTA FILE; optional: 3. LAMBDA and 4.ln(R0)
//OUTPUT: returns, for a given sequence, <N> combined from forward and reverse strand and sequence length.

#include<stdio.h>
#include<cmath>
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
//#include <boost/multi_array.hpp>

using namespace std;

const int A = 0;
const int C = 1;
const int G = 2;
const int T = 3;

//typedef boost::multi_array<double, 2> array_type;
//typedef array_type::index array_index;
//typedef vector<vector<double> > My2DArray;

int main(int argc, char *argv[]){

  if(argc < 3){
    cout << "W A R N I N G !  Input parameter missing!\n";
    cout << "1. TRANSFAC MATRIX, 2. FASTA FILE; optional: 3. LAMBDA, 4.ln(R0)\n";
    exit(1);
  }

  //parameters
  double lambda;//lambda
  double Rmax;//R0 parameter
  if(argc < 5){ //if R0 and/or lambda are not defined use DEFAULT
    lambda = 0.7;
  }
  else{
    lambda = strtod(argv[3],NULL);
  }

  const int Pseudocount = 1;

  //GLOBAL YEAST BACKGROUND
  //const double gc_content = 0.35;
  //const double at_content = 0.65;
  const double gc_content = 0.5;
  const double at_content = 0.5;

//READ FASTA FILE
//----------------------------------------------------------------
  ifstream fasta(argv[2]);

  if(!fasta){
    cout << "FASTA file not opened\n";
  }
  
  string newID, seqID;
  int seqlength = 0;
  int sequenceheader = 0;
  int firstrun = 1; //indicates first header
  string bases; //complete sequence
  string newbases; //sequence in each row of fasta file

  while(!fasta.eof()){ //GO THROUGH FASTA FILE
    getline(fasta,newbases);

    if(newbases.substr(0,1) == ">"){ //NEW SEQUENCE HEADER
      sequenceheader = 1;
      newID = newbases.substr(1);
    }
    
    if(sequenceheader == 0){ //SEQUENCE LINE
      bases = bases + newbases;
      seqlength = bases.length();
      continue;
    }
    
    if(firstrun == 1){ //SKIP ANNOTATION FOR FIRST HEADER
      firstrun = 0;
      sequenceheader = 0;
      seqID = newID;
      continue;
    }

  }//loop over fasta file

  fasta.close();

        
    //RESET SEQUENCE VARIABLES
    //sequenceheader = 0;
    //seqID = newID;
    //bases = "";
    

//READ TRANSFAC FILE
//----------------------------------------------------------------
  ifstream transfac(argv[1]);

  if(!transfac){
    cout << "TRANSFAC file not opened\n";
  }

  double pwm[35][4], complement[35][4]; //matrices for forward and reverse strand
  //My2DArray pwm, complement; //matrices for forward and reverse strand
  string word[5]; //elements in row
  double value[4]; //count + pseudocount
  double max = 0; //consensus base count
  string row; //transfac file rows
  string delimiters = " \t"; //word seperators in each line
  int motiflength, reading, n;
  int start, end;
  string matrix_id;
  double info, sum, p;
  //array_type A(boost::extents[35][4]);
  //vector<My2DArray> all_pwm, all_complement;

  while(!transfac.eof()){
    getline(transfac,row);

    start = row.find_first_not_of(delimiters);
    int i = 0;

    while((start != string::npos)&&(i<5)){ //split row into tokens - word[]
      end = row.find_first_of(delimiters,start+1);
      if(end == string::npos){
	end = row.length();
      }
      word[i] = row.substr(start,end-start);
      i++;
      start = row.find_first_not_of(delimiters,end+1);
    }

    if((reading == 1)&&(word[0] == "XX")){ //end of matrix is reached
      reading = 0;
      motiflength++;

      for(int m = 0; m < motiflength; m++){ //create complement to matrix
        complement[motiflength-m-1][A] = pwm[m][T];
        complement[motiflength-m-1][C] = pwm[m][G];
        complement[motiflength-m-1][G] = pwm[m][C];
        complement[motiflength-m-1][T] = pwm[m][A];
        //vector<double> tmp;
        //tmp.push_back(pwm[m][T]);
        //tmp.push_back(pwm[m][G]);
        //tmp.push_back(pwm[m][C]);
        //tmp.push_back(pwm[m][A]);
        //complement.push_back(tmp);
      }

      //all_pwm.push_back(pwm);
      //all_complement.push_back(complement);
      
      //SET PARAMETERS Lambda and R0
      //----------------------------------------------------------------
        if(argc < 5){ //if R0 and/or lambda are not defined use DEFAULT
          Rmax = exp(0.584 * motiflength - 5.66);
        }
        else{
          Rmax = strtod(argv[4],NULL);
          Rmax = exp(Rmax);
        }
        cerr << "Lambda: " << lambda << "\t" << "R0: " << Rmax << "\t" << "MotifLength: " << motiflength << "\n";  

      // match current motif to seq
      double Pstored[seqlength - motiflength + 1];
      double P_combined = 0; //only palindrome correction
      double P_uncorrected = 0; //uncorrected expected count
      
      //LOOP OVER SEQUENCE
      int illegalBase;

      cout << matrix_id << "\t";

      for(int n = 0; n < seqlength - motiflength + 1; n++){ //LOOP OVER SEQUENCE
        double dE_forward = 0;
        double dE_compl = 0;

        //LOOP OVER MOTIF
        for(int m = 0; m < motiflength; m++){ //LOOP OVER MOTIF
	  
          int BASE;
          illegalBase = 0;
          switch(bases[n+m])
            {
            case 'A': 
              BASE = 0;break;
            case 'a':
              BASE = 0;break;
            case 'C':
              BASE = 1;break;
                case 'c':
                  BASE = 1;break;
            case 'G':
              BASE = 2;break;
                case 'g':
                  BASE = 2;break;
            case 'T':
              BASE = 3;break;
                case 't':
                  BASE = 3;break;	  
            case 'N':
              illegalBase = 1;break;
            case 'n':
              illegalBase = 1;break;
            default:
              cerr << "ILLEGAL CHARACTER " << bases[n+m] << " FOUND!\n Only A/a, T/t, G/g, C/c, N/n are allowed!\n"; illegalBase = 1;
            }

          dE_forward = dE_forward + pwm[m][BASE];
          dE_compl = dE_compl + complement[m][BASE];

        }//loop over motif
      

        //CALCULATE P(BOUND) FOR CURRENT SITE
        if(illegalBase == 0){
          double product = Rmax * exp(-1*dE_forward);
          double P_bound_F = product/(1 + product);

          product = Rmax * exp(-1*dE_compl);
          double P_bound_C = product/(1 + product);

          Pstored[n] = P_bound_F + (1 - P_bound_F) * P_bound_C; //correction for forward and reverse strand
          P_combined = P_combined + Pstored[n];
          //P_uncorrected = P_uncorrected + P_bound_F + P_bound_C;
          //cout<<P_bound_F<<"\t"<<P_bound_C<<"\t"<<Pstored[n]<<"\t"<<P_combined<<"\t"<<P_uncorrected<<"\n";
          cout << Pstored[n] << ",";
        }
        else{
          Pstored[n] = 0;
        }

      }//loop over sequence

      cout << "\t" << P_combined << "\t" << motiflength << "\t" <<
          info << "\n";
    
    }

    if(reading == 1){ //generate matrix
      n++;
      motiflength = n;

      //add pseudocounts
      value[0] = strtod( word[1].c_str(),NULL ) + Pseudocount;
      value[1] = strtod( word[2].c_str(),NULL ) + Pseudocount;
      value[2] = strtod( word[3].c_str(),NULL ) + Pseudocount;
      value[3] = strtod( word[4].c_str(),NULL ) + Pseudocount;
    
      // prob
      sum = 0;
      for (int i = 0; i < 4; ++i)
        sum += value[i];
      for (int i = 0; i < 4; ++i) {
        p = value[i] / sum; 
        info += -p * log2(p);
      }

      double maxAT, maxCG;

      if(value[0] > value[3]){
	maxAT = value[0];
      }
      else{
	maxAT = value[3];
      }

      if(value[1] > value[2]){
	maxCG = value[1];
      }
      else{
	maxCG = value[2];
      }

      if(maxAT > maxCG){
        pwm[n][A] = log(maxAT/value[0]) / lambda;
        pwm[n][C] = log((maxAT/at_content)*(gc_content/value[1])) / lambda;
        pwm[n][G] = log((maxAT/at_content)*(gc_content/value[2])) / lambda;
        pwm[n][T] = log(maxAT/value[3]) / lambda;
        //vector<double> tmp;
        //tmp.push_back(log(maxAT/value[0]) / lambda);
        //tmp.push_back(log((maxAT/at_content)*(gc_content/value[1])) / lambda);
        //tmp.push_back(log((maxAT/at_content)*(gc_content/value[2])) / lambda);
        //tmp.push_back(log(maxAT/value[3]) / lambda);
        //pwm.push_back(tmp);
      }
      else{
        pwm[n][A] = log((maxCG/gc_content)*(at_content/value[0])) / lambda;
        pwm[n][C] = log(maxCG/value[1]) / lambda;
        pwm[n][G] = log(maxCG/value[2]) / lambda;
        pwm[n][T] = log((maxCG/gc_content)*(at_content/value[3])) / lambda;
        //vector<double> tmp;
        //tmp.push_back(log((maxCG/gc_content)*(at_content/value[0])) / lambda);
        //tmp.push_back(log(maxCG/value[1]) / lambda);
        //tmp.push_back(log(maxCG/value[2]) / lambda);
        //tmp.push_back(log((maxCG/gc_content)*(at_content/value[3])) / lambda);
        //pwm.push_back(tmp);
      }
      if(maxAT == maxCG){
        pwm[n][A] = log(maxAT/value[0]) / lambda;
        pwm[n][C] = log(maxAT/value[1]) / lambda;
        pwm[n][G] = log(maxAT/value[2]) / lambda;
        pwm[n][T] = log(maxAT/value[3]) / lambda;
        //vector<double> tmp;
        //tmp.push_back(log(maxAT/value[0]) / lambda);
        //tmp.push_back(log(maxAT/value[1]) / lambda);
        //tmp.push_back(log(maxAT/value[2]) / lambda);
        //tmp.push_back(log(maxAT/value[3]) / lambda);
        //pwm.push_back(tmp);
      }
    }
    
    if(word[0] == "P0"){ //start of matrix
      n = -1;
      reading = 1;
      info = 0;
    }
    
    if(word[0] == "ID"){ //matrix ID
      matrix_id = word[1];
    }
  }
  
  transfac.close();





  /*
  //LAST SEQUENCE

  double Pstored[seqlength - motiflength + 1];
  double P_combined = 0; //only palindrome correction
  double P_uncorrected = 0; //uncorrected expected count

  //LOOP OVER SEQUENCE
  int illegalBase;

  for(int n = 0; n < seqlength - motiflength + 1; n++){ //LOOP OVER SEQUENCE
    double dE_forward = 0;
    double dE_compl = 0;

    //LOOP OVER MOTIF
    for(int m = 0; m < motiflength; m++){ //LOOP OVER MOTIF

      int BASE;
      illegalBase = 0;
      switch(bases[n+m])
	{
          case 'A':
            BASE = 0;break;
          case 'a':
            BASE = 0;break;
          case 'C':
            BASE = 1;break;
          case 'c':
            BASE = 1;break;
          case 'G':
            BASE = 2;break;
          case 'g':
            BASE = 2;break;
          case 'T':
            BASE = 3;break;
          case 't':
            BASE = 3;break;
          case 'N':
            illegalBase = 1;break;
          case 'n':
            illegalBase = 1;break;
          default:
            cerr << "ILLEGAL CHARACTER FOUND!\n Only A/a, T/t, G/g, C/c, N/n are allowed!\n"; illegalBase = 1;
	}

      dE_forward = dE_forward + pwm[m][BASE];
      dE_compl = dE_compl + complement[m][BASE];

    }//loop over motif

    //CALCULATE P(BOUND) FOR CURRENT SITE
    if(illegalBase == 0){
      double product = Rmax * exp(-1*dE_forward);
      double P_bound_F = product/(1 + product);

      product = Rmax * exp(-1*dE_compl);
      double P_bound_C = product/(1 + product);

      Pstored[n] = P_bound_F + (1 - P_bound_F) * P_bound_C; //correction for forward and reverse strand
      P_combined = P_combined + Pstored[n];
      //P_uncorrected = P_uncorrected + P_bound_F + P_bound_C;
      //cout<<P_bound_F<<"\t"<<P_bound_C<<"\t"<<Pstored[n]<<"\t"<<P_combined<<"\t"<<P_uncorrected<<"\n";
    }
    else{
      Pstored[n] = 0;
    }

  }//loop over sequence

  cout << matrix_id[0] << "\t" << P_combined << "\t" << seqlength << "\n";
  */

  return 0;
}
