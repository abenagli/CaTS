// g++ -Wall -o findLastProcess.exe findLastProcess.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cout << "insert fileName" << std::endl;
    return -1;
  }
  
  std::ifstream inFile(argv[1],std::ios::in);
  std::cout << ">>> reading file " << argv[1] << std::endl;
  
  std::string line;
  std::string line2;
  std::string str;
  int event = 0;
  int event_new = 0;
  
  while(1)
  {
    getline(inFile,line,'\n');
    if( !inFile.good() ) break;
    
    
    std::size_t found = line.find(std::string("G4Event"));
    if( found != std::string::npos )
    {
      ++event_new;
      continue;
    }
    
    if( event_new != event )
    {
      event = event_new;
      std::cout << "G4Event " << event << std::endl;
    }
    
    
    if( line.size() > 0 && line.at(0) == '*' )
    {
      std::string particleName,trackId,parentId;
      std::string step,vx,vy,vz,K,processName,fiStep;
      
      getline(inFile,line,'\n');
      std::stringstream ss(line);
      ss >> str >> str >> str >> str >> str >> particleName >> str >> str >> str >> trackId >> str >> str >> str >> parentId;
      
      getline(inFile,line,'\n'); // empty line
      getline(inFile,line,'\n'); // empty line
      getline(inFile,line,'\n'); // line with labels

      getline(inFile,line,'\n'); // line with first step
      std::stringstream ss2(line);
      ss2 >> step >> vx >> vy >> vz >> K;
      
      bool lastStep = false;
      while( !lastStep )
      {
        if( line.size() == 0 )
        {
          lastStep = true;
          continue;
        }
        
        std::size_t found = line.find(std::string("G4Event"));
        if( found != std::string::npos )
        {
          ++event_new;
          lastStep = true;
          continue;
        }        
        
        getline(inFile,line,'\n'); // line with other steps
        std::stringstream ss3(line);
        ss3 >> fiStep >> str >> str >> str >> str >> str >> str >> str >> str >> processName;
      }
      
      //if( particleName == "neutron," && trackId == "1,")
      {
        std::cout << ">>> particleName: " << std::fixed << std::setw(10) << particleName
                  << "   trackId: " << std::setw(4) << trackId
                  << "   parentId: " << std::setw(4) << parentId
                  << "   in. vertex: (" << std::setw(8) << vx << "," << std::setw(8) << vy << "," << std::setw(8) << vz << ")"
                  << "   in. K: " << std::setw(8) << K
                  << "   fi. process: " << processName
                  << "   fi. step: " << fiStep
                  << std::endl;
      }
    }
    
    //std::cout << line << std::endl;
    // if( line.at(0) == '*' )
    // {
    //   ss >> str >> str >> str >> str >> spill;
    //   WCRawData[spill] = new std::vector<int>;
    //   WCDecodedData[spill] = new std::vector<std::map<int,std::pair<std::vector<int>,std::vector<int> > > >;
    //   continue;
    // }
    
    // int token;
    // while( ss >> std::hex >> token )
    // {
    //   int val = int(token);
    // }
  }
}
