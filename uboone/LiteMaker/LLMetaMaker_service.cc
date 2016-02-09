
#include "LLMetaMaker.h"

//-----------------------------------------------------------------------------------------
util::LLMetaMaker::LLMetaMaker(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
//----------------------------------------------------------------------------------------
{

  reconfigure(pset);

  reg.sPostOpenFile.watch     (this, &LLMetaMaker::postOpenFile     );

  reg.sPreBeginRun.watch      (this, &LLMetaMaker::preBeginRun      );
  reg.sPostBeginRun.watch     (this, &LLMetaMaker::postBeginRun     );

  reg.sPreProcessEvent.watch  (this, &LLMetaMaker::preProcessEvent  );
  reg.sPostProcessEvent.watch (this, &LLMetaMaker::postProcessEvent );

  reg.sPostBeginJob.watch     (this, &LLMetaMaker::postBeginJob );
  reg.sPostEndJob.watch       (this, &LLMetaMaker::postEndJob   );

}

//------------------------------------------------------------------
void util::LLMetaMaker::reconfigure(fhicl::ParameterSet const& pset)
//------------------------------------------------------------------
{
  _litefile_v = pset.get<std::vector<std::string> >("Files");
}

//------------------------------------------------------------------
void util::LLMetaMaker::postBeginJob()
//------------------------------------------------------------------
{
  std::cout << "This is after begin-job" << std::endl;
}

//------------------------------------------------------------------
void util::LLMetaMaker::postEndJob()
//------------------------------------------------------------------
{
  std::cout << "This is after end-job" << std::endl;
}

//------------------------------------------------------------
void util::LLMetaMaker::preProcessEvent(const art::Event& evt)
//------------------------------------------------------------
{
  std::cout << "This is before event process" << std::endl;
}

//------------------------------------------------------------
void util::LLMetaMaker::postProcessEvent(const art::Event& evt)
//------------------------------------------------------------
{
  std::cout << "This is after event process" << std::endl;
}

//------------------------------------------------------
void util::LLMetaMaker::preBeginRun(art::Run const& run)
//------------------------------------------------------
{
  std::cout << "This is prior to begin-run boundary" << std::endl;
}

//------------------------------------------------------
void util::LLMetaMaker::postBeginRun(art::Run const& run)
//------------------------------------------------------
{
  std::cout << "This is after to begin-run boundary" << std::endl;
}

//---------------------------------------------------------------
void util::LLMetaMaker::postOpenFile(const std::string& filename)
//---------------------------------------------------------------
{
  std::cout << "This is after open-file" << std::endl;
}


namespace util{

  DEFINE_ART_SERVICE(LLMetaMaker)

} // namespace util  

