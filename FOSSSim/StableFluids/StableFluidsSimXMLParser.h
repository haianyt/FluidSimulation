#ifndef __STABLE_FLUIDS_SIM_XML_PARSER_H__
#define __STABLE_FLUIDS_SIM_XML_PARSER_H__

#include <Eigen/StdVector>

#include <iostream>
#include <fstream>
#include <limits>

#include "FOSSSim/MathDefines.h"
#include "FOSSSim/StringUtilities.h"

#include "rapidxml.hpp"

class StableFluidsSimXMLParser
{
public:
  
  void loadExecutableSimulation( const std::string& file_name, scalar& max_time, int& velocity_pattern, std::string& description, std::string& scenetag, scalar &diff, scalar &visc, bool& use_advect, bool& use_proj );

private:
  void loadXMLFile( const std::string& filename, std::vector<char>& xmlchars, rapidxml::xml_document<>& doc );
  bool loadTextFileIntoString( const std::string& filename, std::string& filecontents );
  void loadSceneTag( rapidxml::xml_node<>* node, std::string& scenetag );
  void loadMaxTime( rapidxml::xml_node<>* node, scalar& max_t );
  void loadSceneDescriptionString( rapidxml::xml_node<>* node, std::string& description_string );
  void loadVelocityPattern( rapidxml::xml_node<>* node, int& velocity_pattern );
  void loadSimulationType( rapidxml::xml_node<>* node, std::string& simtype, bool& use_advect, bool& use_proj );
  void loadViscosity( rapidxml::xml_node<>* node, scalar& visc );
  void loadDiffusion( rapidxml::xml_node<>* node, scalar& diff );


};

#endif