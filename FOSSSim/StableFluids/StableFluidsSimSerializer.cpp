#include "StableFluidsSimSerializer.h"

void StableFluidsSimSerializer::serializeScene( StableFluidsSim& scene, std::ofstream& outputstream ) const
{
  assert( outputstream.is_open() );
  
  int ndof_u = (scene.physicalRows()+2) * (scene.physicalCols()+1);
  int ndof_v = (scene.physicalRows()+1) * (scene.physicalCols()+2);

  scalar* udataDiffusion = scene.getHorizontalVelocitiesDiffusion().data();
  outputstream.write((char*)udataDiffusion,ndof_u*sizeof(scalar));
  scalar* vdataDiffusion = scene.getVerticalVelocitiesDiffusion().data();
  outputstream.write((char*)vdataDiffusion,ndof_v*sizeof(scalar));

  scalar* udataAdvect = scene.getHorizontalVelocitiesAdvect().data();
  outputstream.write((char*)udataAdvect,ndof_u*sizeof(scalar));
  scalar* vdataAdvect = scene.getVerticalVelocitiesAdvect().data();
  outputstream.write((char*)vdataAdvect,ndof_v*sizeof(scalar));

  scalar* udata = scene.getHorizontalVelocities().data();
  outputstream.write((char*)udata,ndof_u*sizeof(scalar));
  scalar* vdata = scene.getVerticalVelocities().data();
  outputstream.write((char*)vdata,ndof_v*sizeof(scalar));

  int ndof_d = (scene.physicalRows()+2) * (scene.physicalCols()+2);  
  scalar* datadensity = scene.getMarkerDensities().data();
  outputstream.write((char*)datadensity,ndof_d*sizeof(scalar));

  int ndof_f = (scene.physicalRows()+2) * (scene.physicalCols()+2);  
  bool* datahasfluid = scene.getHasFluid().data();
  outputstream.write((char*)datahasfluid,ndof_f*sizeof(bool));
    
  scalar diffusion = scene.getDiffusion();
  outputstream.write((char*)&diffusion, sizeof(scalar));
    
  scalar viscosity = scene.getViscosity();
  outputstream.write((char*)&viscosity, sizeof(scalar));
}

void StableFluidsSimSerializer::loadScene( StableFluidsSim& scene, std::ifstream& inputstream ) const
{
  assert( inputstream.is_open() );
  assert( !inputstream.eof() );

  int ndof_u = (scene.physicalRows()+2) * (scene.physicalCols()+1);
  int ndof_v = (scene.physicalRows()+1) * (scene.physicalCols()+2);

  scalar* udataDiffusion = scene.getHorizontalVelocitiesDiffusion().data();
  inputstream.read((char*)udataDiffusion,ndof_u*sizeof(scalar));
  scalar* vdataDiffusion = scene.getVerticalVelocitiesDiffusion().data();
  inputstream.read((char*)vdataDiffusion,ndof_v*sizeof(scalar));

  scalar* udataAdvect = scene.getHorizontalVelocitiesAdvect().data();
  inputstream.read((char*)udataAdvect,ndof_u*sizeof(scalar));
  scalar* vdataAdvect = scene.getVerticalVelocitiesAdvect().data();
  inputstream.read((char*)vdataAdvect,ndof_v*sizeof(scalar));

  scalar* udata = scene.getHorizontalVelocities().data();
  inputstream.read((char*)udata,ndof_u*sizeof(scalar));
  scalar* vdata = scene.getVerticalVelocities().data();
  inputstream.read((char*)vdata,ndof_v*sizeof(scalar));

  int ndof_d = (scene.physicalRows()+2) * (scene.physicalCols()+2);  
  scalar* datadensity = scene.getMarkerDensities().data();
  inputstream.read((char*)datadensity,ndof_d*sizeof(scalar));

  int ndof_f = (scene.physicalRows()+2) * (scene.physicalCols()+2); 
  bool* datahasfluid = scene.getHasFluid().data(); 
  inputstream.read((char*)datahasfluid,ndof_f*sizeof(bool));

  scalar diffusion = 0;
  inputstream.read((char*)&diffusion, sizeof(scalar));
  scene.setDiffusion(diffusion);
  
  scalar viscosity = 0;
  inputstream.read((char*)&viscosity, sizeof(scalar));
  scene.setViscosity(viscosity);
    
  if( inputstream.fail() )
  {
    std::cout << outputmod::startred << "Error in StableFluidsSimSerializer: " << outputmod::endred << "Failed to load timestep. Exiting." << std::endl;
    exit(1);
  }  
}
