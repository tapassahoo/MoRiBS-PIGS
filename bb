mc_estim.cc:  fstream fid;
mc_estim.cc:  fstream fid;
mc_estim.cc:  fstream fid;
mc_estim.cc:  fstream fid;
mc_estim.cc:  fstream fid;
mc_estim.cc:  fstream fid;
mc_estim.cc:   fstream fid(fdens.c_str(),mode);
mc_estim.cc:// stringstream stype;
mc_estim.cc:  fstream fid;
mc_estim.cc:  fstream fid;
mc_estim.cc:  fstream fid;
mc_estim.cc:  fstream fid;
mc_estim.cc:   fstream fid;
mc_input.cc:#include <fstream>
mc_input.cc:#include <sstream>
mc_input.cc:   ifstream inf(in_file,ios::in);
mc_input.cc:   fstream fid(file_name,mode);
mc_input.cc:   fstream fid(file_name,mode);
mc_input.cc://   fstream fid(file_name,mode | ios::binary);
mc_input.cc:   streamsize size;         
mc_input.cc:         fid.write((char *)&size,sizeof(streamsize));
mc_input.cc:         fid.read((char *)&size,sizeof(streamsize));
mc_input.cc:   fstream fid(file_name,mode | ios::binary);
mc_input.cc:   streamsize size;         
mc_input.cc:         fid.write((char *)&size,sizeof(streamsize));
mc_input.cc:         fid.read((char *)&size,sizeof(streamsize));
mc_input.cc:   fstream fid(file_name,mode | ios::binary);
mc_input.cc:// streamsize size;         
mc_input.cc:  fstream fid(file_name,mode);
mc_input.cc:  stringstream stype; 
mc_input.cc:// ifstream  infile(fileName,ios::in); 
mc_input.cc:void io_setout(fstream &out, int set_precision)
mc_main.cc:#include <sstream>
mc_main.cc:#include "rngstream.h"
mc_main.cc:fstream _feng;      // save accumulated energy
mc_main.cc:   if (!restart) // new run, generate new status, rnd() streams and config files     
mc_main.cc:   RandomIO(IORead,FRANDOM);  // load rnd streams 
mc_main.cc:                          stringstream bc;                // convert block # to string
mc_main.cc:                       stringstream bc;                // convert block # to string
mc_main.cc://  CHECKPOINT: save status, rnd streams and configs ------
mc_main.cc:  stringstream bc;                // convert block # to string
mc_main.cc:   fstream fid;
mc_piqmc.cc:#include "rngstream.h"
mc_poten.cc:#include <sstream>
mc_poten.cc:#include <fstream>
mc_poten.cc:    ifstream fid(fname.c_str(),ios::in);   
mc_poten.cc:    ifstream fid(fname.c_str(),ios::in);
mc_poten.cc:    ifstream fid(fname.c_str(),ios::in);   
mc_poten.cc:   stringstream time; time << NumbRotTimes;                  // number of time slices 
mc_poten.cc:   stringstream temp; temp << Temperature*Units.temperature; // temperature
mc_poten.cc:   ifstream fid(fden.c_str(),ios::in);
mc_poten.cc:   ifstream fid2(feng.c_str(),ios::in);
mc_poten.cc:   ifstream fid3(fesq.c_str(),ios::in);
mc_poten.cc:    stringstream time; time << NumbRotTimes;                  // number of time slices 
mc_poten.cc:    stringstream temp; temp << Temperature*Units.temperature; // temperature
mc_poten.cc:   ifstream fid(fname,ios::in);   
mc_poten.cc:   ifstream fid(fname,ios::in);   
mc_poten.cc:   ifstream fid(fname,ios::in);   
mc_poten.cc:   ifstream fid(fname,ios::in);
mc_randg.cc:int *STREAM[MAXGENS];  // rnd numbers streams
mc_randg.cc:  int max_streams=MAXGENS*max_mpi_procs;        // number of streams 
mc_randg.cc:     int streamnum = mpi_id*MAXGENS+ip;
mc_randg.cc:     STREAM[ip] = init_sprng(streamnum,max_streams,SEED,SPRNG_DEFAULT); // initialize stream 
mc_randg.cc://   print_sprng(stream);	
mc_randg.cc:void RandomIO(int tstatus, const char file_name[]) // check point for sprng streams
mc_randg.cc:             int size = pack_sprng(STREAM[ip],&bytes); // pack  a stream state 
mc_randg.cc:             fwrite(bytes,1,size,fp);                  // store a stream state
omprng.cc:#include "rngstream.h"
omprng.cc:#include <iostream>
rngstream.cc: * File:           RngStream.cpp for multiple streams of Random Numbers
rngstream.cc:#include <iostream>
rngstream.cc:#include "rngstream.h"
rngstream.cc:   /* Information on a stream. The arrays {Cg, Bg, Ig} contain the current
rngstream.cc:   state of the stream, the starting state of the current SubStream, and the
rngstream.cc:   starting state of the stream. This stream generates antithetic variates
rngstream.cc:void RngStream::ResetStartSubstream ()
rngstream.cc:void RngStream::ResetNextSubstream ()
rngstream.cc:    cout << "The current state of the Rngstream";
