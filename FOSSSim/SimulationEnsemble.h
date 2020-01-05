#ifndef __SIMULATION_ENSEMBLE_H__
#define __SIMULATION_ENSEMBLE_H__

#include "MathDefines.h"
#include "RenderingUtilities.h"

class SimulationEnsemble
{
public:

  virtual void stepSystem( const scalar& dt ) = 0;

  /////////////////////////////////////////////////////////////////////////////
  // Mouse and keyboard input handlers

  virtual void keyboard( const unsigned char& key, const int& x, const int& y ) = 0;

  virtual void special( const int& key, const int& x, const int& y ) = 0;

  virtual void mouse( const int& button, const int& state, const int& x, const int& y ) = 0;

  virtual void motion( const int& x, const int& y ) = 0;

  /////////////////////////////////////////////////////////////////////////////
  // Functions for use by OpenGL and GLUT

  virtual void initializeOpenGL() = 0;
  
  virtual void reshape( const int& w, const int& h ) = 0;

  virtual void display() = 0;

  virtual unsigned int getGlutDisplayMode() const = 0;

  virtual const int& getWindowWidth() const = 0;

  virtual const int& getWindowHeight() const = 0;

  // TODO: I only keep this around for setting text color... hanlde this better
  // and scratch this method.
  virtual const renderingutils::Color& getBackgroundColor() const = 0;

  /////////////////////////////////////////////////////////////////////////////
  // Scene Benchmarking Functions

  virtual void updateSceneComparison() = 0;
  virtual void printErrorInformation( bool print_pass ) = 0;

  /////////////////////////////////////////////////////////////////////////////
  // Serialization Functions
  
  virtual void copyComparisonSceneToScene() = 0;
  virtual void serializeScene( std::ofstream& outputstream ) = 0;
  virtual void loadComparisonScene( std::ifstream& inputstream ) = 0;

  virtual void setVelocityPattern( int velocity_pattern ) = 0;
  virtual void setDiffusion( scalar diff ) = 0;
  virtual void setViscosity( scalar visc ) = 0;
  virtual void setUseAdvect( bool use_advect ) = 0;
  virtual void setUseProj( bool use_proj ) = 0;
};

#endif
