#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "datatypes/ub_BeamHeader.h"
#include "datatypes/ub_BeamData.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/filesystem.hpp>

using namespace gov::fnal::uboone::datatypes;
using namespace std;

int main(int argc, char** argv)
{
  if ( argc!=2 && argc!=3 ) {
    cout << argv[0]<<" file_name [device]"<<endl<<endl;
    cout << "\t [device] - optional; dumps min, max, avg and sum for that device"<<endl;
    cout << "\t          - for arrays dumps only the first element"<<endl;
    return 0;
  }

  stringstream ss(argv[1],ios_base::app | ios_base::out);
  ss<<".lock";
  int nlock=0;
  while ( boost::filesystem::exists(ss.str().c_str()) ) {
    sleep(1);
    nlock++;
    if (nlock==10) {
      cout<<"Error, will not read file until lock file ("<<ss.str()<<" is not removed. "
	  <<"Waited for 10s, no luck."<<endl;
      exit(0);
    }
  }

  bool printall=false;
  if (argc == 2) printall=true;
  std::ifstream ifs(argv[1],ios::in | ios::binary);
  if (!ifs.is_open()) {
    cerr << "Can't find "<<argv[1]<<endl;
    return 0;
  }

  ub_BeamHeader bhdr;
  ub_BeamData bdata;

  int nevent=0;
  int nevent_wdev=0;
  float sum=0;
  float min=999999;
  float max=-999999;
  float avg=0;

  boost::circular_buffer<uint32_t> rate(300);

  string devname="";

  uint32_t seconds; // GPS clock. Since Jan 1, 2012.
  uint16_t milli_seconds;

  uint32_t first_sec; // GPS clock. Since Jan 1, 2012.
  uint16_t first_msec;

  /*
  uint32_t last_sec; // GPS clock. Since Jan 1, 2012.
  uint16_t last_msec;

  //    printall=false;
  //ub_BeamHeader last_bhdr;
  */
  while (1) {
    try {
      boost::archive::binary_iarchive ia(ifs);
      ia>>bhdr;
      nevent++;
      seconds=bhdr.getSeconds();
      milli_seconds=bhdr.getMilliSeconds();
      rate.push_back(seconds);
      /*
  double diff=double(seconds-last_sec)+double(milli_seconds-last_msec)/1000.;
            if (diff<3 &&nevent>1) {
      	cout <<"diff = "<<diff<<endl;
      	cout <<last_bhdr;
      	cout <<bhdr;
	cout<<"-----"<<endl;
          }
	    //      last_bhdr.seconds=bhdr.seconds;
	    // last_bhdr.milli_seconds=bhdr.milli_seconds;
      last_bhdr=bhdr;
      last_sec=seconds;
      last_msec=milli_seconds;	   
      //if (0) cout <<diff<<endl; 
      */
      if (nevent==1) {
	first_sec=seconds;
	first_msec=milli_seconds;
      }
      if (printall) cout <<bhdr;
      for (int i=0;i<int(bhdr.getNumberOfDevices());i++) {
	ia>>bdata;
	if (printall) cout <<bdata;
	else if ((bdata.getDeviceName().find(argv[2])!=std::string::npos && devname=="")||
		 (bdata.getDeviceName()==devname) ) {
	  devname=bdata.getDeviceName();
	  sum+=bdata.getData()[0];
	  if ( min>bdata.getData()[0] ) min = bdata.getData()[0];
	  if ( max<bdata.getData()[0] ) max = bdata.getData()[0];
	  nevent_wdev++;
	}
      }
      if (printall) cout<<endl;
    } catch(boost::archive::archive_exception const& e) {

      if (printall) cout <<"Read "<<nevent<<" events"<<endl;
      else {
	avg=sum/float(nevent_wdev); 
	cout <<left
	     <<setw(12)<<"DEVICE"
	     <<setw(12)<<"NEVENTS"
	     <<setw(15)<<"MIN"
	     <<setw(15)<<"MAX"
	     <<setw(15)<<"AVG"
	     <<setw(15)<<"SUM"<<endl;
	cout.width(12);
	cout <<std::left
	     <<setw(12)<<devname
	     <<setw(12)<<nevent_wdev
	     <<std::setprecision(8)<<setw(15)<<min
	     <<std::setprecision(8)<<setw(15)<<max
	     <<std::setprecision(8)<<setw(15)<<avg
	     <<std::setprecision(8)<<setw(15)<<sum<<endl;
      }
      break;
    }
  }
  ifs.close();

  int n=0;
  uint64_t begtime=0;
  uint64_t endtime=0;
  for (auto itr=rate.begin();itr!=rate.end();itr++) {
    if (n==0) begtime=*itr;
    if (*itr!=0) {
      n++;
      endtime=*itr;
    }
  }

  cout <<"First event timestamp (sec msec) "<<first_sec<<" "<<first_msec<<endl;
  cout <<"Last event timestamp (sec msec)  "<<seconds<<" "<<milli_seconds<<endl;
  cout <<"Current rate (last "<<n<<" spills)    "<<float(n)/float(endtime-begtime)<<" Hz"<<endl;
}
