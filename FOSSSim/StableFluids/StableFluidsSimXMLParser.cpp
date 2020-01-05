#include "StableFluidsSimXMLParser.h"

void StableFluidsSimXMLParser::loadExecutableSimulation( const std::string& file_name, scalar& max_time, int& velocity_pattern, std::string& description, std::string& scenetag, scalar &diff, scalar &visc, bool& use_advect, bool& use_proj )
{
  // Load the xml document
  std::vector<char> xmlchars;
  rapidxml::xml_document<> doc;
  loadXMLFile( file_name, xmlchars, doc );

  // Attempt to locate the root node
  rapidxml::xml_node<>* node = doc.first_node("scene");
  if( node == NULL ) 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse xml scene file. Failed to locate root <scene> node. Exiting." << std::endl;
    exit(1);
  }
  
  // Determine what simulation type this is (particle, rigid body, etc)
  std::string simtype;
  loadSimulationType( node, simtype, use_advect, use_proj );
  
  // Parse common state
  loadMaxTime( node, max_time );
  loadSceneDescriptionString( node, description );
  loadSceneTag( node, scenetag );
  
  // Parse the user-requested simulation type. The default is a particle simulation.
  if( simtype == "stable-fluids" )
  {
    loadVelocityPattern( node, velocity_pattern );
    loadDiffusion( node, diff );
    loadViscosity( node, visc );
  }
  else 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid simtype '" << simtype << "' specified. Valid options are 'particle-system' and 'rigid-body'. Exiting." << std::endl;
    exit(1);
  }
}

void StableFluidsSimXMLParser::loadXMLFile( const std::string& filename, std::vector<char>& xmlchars, rapidxml::xml_document<>& doc )
{
  // Attempt to read the text from the user-specified xml file
  std::string filecontents;
  if( !loadTextFileIntoString(filename,filecontents) )
  {
    std::cerr << "\033[31;1mERROR IN TWODSCENEXMLPARSER:\033[m XML scene file " << filename << ". Failed to read file." << std::endl;
    exit(1);
  }
  
  // Copy string into an array of characters for the xml parser
  for( int i = 0; i < (int) filecontents.size(); ++i ) xmlchars.push_back(filecontents[i]);
  xmlchars.push_back('\0');
  
  // Initialize the xml parser with the character vector
  doc.parse<0>(&xmlchars[0]);
}

bool StableFluidsSimXMLParser::loadTextFileIntoString( const std::string& filename, std::string& filecontents )
{
  // Attempt to open the text file for reading
  std::ifstream textfile(filename.c_str(),std::ifstream::in);
  if(!textfile) return false;
  
  // Read the entire file into a single string
  std::string line;
  while(getline(textfile,line)) filecontents.append(line);
  
  textfile.close();
  
  return true;
}

void StableFluidsSimXMLParser::loadSceneTag( rapidxml::xml_node<>* node, std::string& scenetag )
{
  assert( node != NULL );

  if( node->first_node("scenetag") )
  {
    if( node->first_node("scenetag")->first_attribute("tag") )
    {
      scenetag = node->first_node("scenetag")->first_attribute("tag")->value();
    }
    else
    {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value of tag attribute for scenetag. Value must be string. Exiting." << std::endl;
      exit(1);
    }
  }
}

void StableFluidsSimXMLParser::loadMaxTime( rapidxml::xml_node<>* node, scalar& max_t )
{
  assert( node != NULL );

  // Attempt to locate the duraiton node
  rapidxml::xml_node<>* nd = node->first_node("duration");
  if( nd == NULL ) 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No duration specified. Exiting." << std::endl;
    exit(1);
  }
  
  // Attempt to load the duration value
  rapidxml::xml_attribute<>* timend = nd->first_attribute("time"); 
  if( timend == NULL ) 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No duration 'time' attribute specified. Exiting." << std::endl;
    exit(1);
  }
  
  max_t = std::numeric_limits<scalar>::signaling_NaN();
  if( !stringutils::extractFromString(std::string(timend->value()),max_t) )
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'time' attribute for duration. Value must be numeric. Exiting." << std::endl;
    exit(1);
  }
}

void StableFluidsSimXMLParser::loadSceneDescriptionString( rapidxml::xml_node<>* node, std::string& description_string )
{
  assert( node != NULL );
  
  description_string = "No description specified.";
  
  // Attempt to locate the integrator node
  rapidxml::xml_node<>* nd = node->first_node("description");
  if( nd != NULL ) 
  {
    rapidxml::xml_attribute<>* typend = nd->first_attribute("text"); 
    if( typend == NULL ) 
    {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No text attribute specified for description. Exiting." << std::endl;
      exit(1);
    }
    description_string = typend->value();
  }
}

void StableFluidsSimXMLParser::loadVelocityPattern( rapidxml::xml_node<>* node, int& velocity_pattern )
{
  assert( node != NULL );

  // Attempt to locate the duraiton node
  rapidxml::xml_node<>* nd = node->first_node("velocityfieldpattern");
  if( nd == NULL ) 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No velocity field pattern specified. Exiting." << std::endl;
    exit(1);
  }
  
  // Attempt to load the duration value
  rapidxml::xml_attribute<>* velocitynd = nd->first_attribute("code"); 
  if( velocitynd == NULL ) 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No pattern 'code' attribute specified. Exiting." << std::endl;
    exit(1);
  }
  
  velocity_pattern = std::numeric_limits<int>::signaling_NaN();
  if( !stringutils::extractFromString(std::string(velocitynd->value()), velocity_pattern) )
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'code' attribute for velocity pattern. Value must be numeric. Exiting." << std::endl;
    exit(1);
  }

}

void StableFluidsSimXMLParser::loadSimulationType( rapidxml::xml_node<>* node, std::string& simtype, bool& use_advect, bool& use_proj )
{
  assert( node != NULL );
  rapidxml::xml_node<>* nd = node->first_node("simtype");
  
  if( node->first_node("simtype") ) {
    rapidxml::xml_node<>* nd = node->first_node("simtype");
    if( nd->first_attribute("type") ) simtype = nd->first_attribute("type")->value();
    
    rapidxml::xml_attribute<>* advnd = nd->first_attribute("advect");
    if( advnd != NULL )
    {
      if( !stringutils::extractFromString(std::string(advnd->value()),use_advect) )
      {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'advect' attribute for diffusion. Value must be boolean. Exiting." << std::endl;
        exit(1);
      }
    }
    
    rapidxml::xml_attribute<>* projnd = nd->first_attribute("project");
    if( projnd != NULL )
    {
      if( !stringutils::extractFromString(std::string(projnd->value()),use_proj) )
      {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'project' attribute for diffusion. Value must be boolean. Exiting." << std::endl;
        exit(1);
      }
    }
  }
}

void StableFluidsSimXMLParser::loadDiffusion( rapidxml::xml_node<>* node, scalar& diff )
{
  assert( node != NULL );

  // Attempt to locate the duraiton node
  rapidxml::xml_node<>* nd = node->first_node("diffusion");
  if( nd == NULL ) 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No diffusion specified. Exiting." << std::endl;
    exit(1);
  }
  
  // Attempt to load the diffusion value
  rapidxml::xml_attribute<>* diffnd = nd->first_attribute("diff"); 
  if( diffnd == NULL ) 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No diffusion 'diff' attribute specified. Exiting." << std::endl;
    exit(1);
  }
  
  diff = std::numeric_limits<scalar>::signaling_NaN();
  if( !stringutils::extractFromString(std::string(diffnd->value()),diff) )
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'diff' attribute for diffusion. Value must be numeric. Exiting." << std::endl;
    exit(1);
  }
}

void StableFluidsSimXMLParser::loadViscosity( rapidxml::xml_node<>* node, scalar& visc )
{
  assert( node != NULL );

  // Attempt to locate the duraiton node
  rapidxml::xml_node<>* nd = node->first_node("viscosity");
  if( nd == NULL ) 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No viscosity specified. Exiting." << std::endl;
    exit(1);
  }
  
  // Attempt to load the diffusion value
  rapidxml::xml_attribute<>* viscnd = nd->first_attribute("visc"); 
  if( viscnd == NULL ) 
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No viscosity 'visc' attribute specified. Exiting." << std::endl;
    exit(1);
  }
  
  visc = std::numeric_limits<scalar>::signaling_NaN();
  if( !stringutils::extractFromString(std::string(viscnd->value()),visc) )
  {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'visc' attribute for viscosity. Value must be numeric. Exiting." << std::endl;
    exit(1);
  }
}

