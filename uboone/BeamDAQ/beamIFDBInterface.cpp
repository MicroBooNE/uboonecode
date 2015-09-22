#include "beamIFDBInterface.h"
#include <iostream>

#include "httpResponse.h"

#include <messagefacility/MessageLogger/MessageLogger.h>

using namespace gov::fnal::uboone::beam;
using namespace std;
beamIFDBInterface::beamIFDBInterface()
{
  curl_global_init(CURL_GLOBAL_ALL);

  /* init the curl session */
  fCURLHandle = curl_easy_init();
  
  /* some servers don't like requests that are made without a user-agent 
     field, so we provide one 
  */
  curl_easy_setopt(fCURLHandle, CURLOPT_USERAGENT, "libcurl-agent/1.0");
  
  /* Enable redirection */
  curl_easy_setopt(fCURLHandle, CURLOPT_FOLLOWLOCATION, 1);
}

beamIFDBInterface::~beamIFDBInterface()
{
  /* cleanup curl stuff */
  curl_easy_cleanup(fCURLHandle);

  /* we're done with libcurl, so clean it up */
  curl_global_cleanup();
}

void beamIFDBInterface::GetData(const char *url,httpResponse* response)
{
  /* specify URL to get */
  curl_easy_setopt(fCURLHandle, CURLOPT_URL, url);
  
  /* send all data to this function  */
  curl_easy_setopt(fCURLHandle, CURLOPT_WRITEFUNCTION, &writeMemoryCallback);
  
  /* we pass our 'response' struct to the callback function */
  curl_easy_setopt(fCURLHandle, CURLOPT_WRITEDATA, (void *)response);

  /* set error buffer */
  static char error[1024];
  curl_easy_setopt(fCURLHandle, CURLOPT_ERRORBUFFER, error);

  /* get it! */
  int retcode=curl_easy_perform(fCURLHandle);
  if ( retcode != CURLE_OK ) {
    mf::LogError("")<<"Received error code "<<retcode<<" while retreiving data from server. Returned error buffer: "<<error;
    exit(EXIT_FAILURE);
  }

}


