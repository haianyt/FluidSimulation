#include "StableFluidsSimGrader.h"

StableFluidsSimGrader::StableFluidsSimGrader()
: m_accumulated_u_residual(0.0)
, m_max_u_residual(0.0)
, m_acceptable_accumulated_u_residual(1.0e-6)
, m_acceptable_max_u_residual(1.0e-10)
, m_accumulated_v_residual(0.0)
, m_max_v_residual(0.0)
, m_acceptable_accumulated_v_residual(1.0e-6)
, m_acceptable_max_v_residual(1.0e-10)
, m_accumulated_u_diffusion_residual(0.0)
, m_accumulated_v_diffusion_residual(0.0)
, m_accumulated_u_advect_residual(0.0)
, m_accumulated_v_advect_residual(0.0)
{}

void StableFluidsSimGrader::addToAccumulatedResidual( const StableFluidsSim& oracle_scene, const StableFluidsSim& testing_scene )
{
  assert( oracle_scene.physicalRows() == testing_scene.physicalRows() );
  assert( oracle_scene.physicalCols() == testing_scene.physicalCols() );
  assert( oracle_scene.getHorizontalVelocities().size() == testing_scene.getHorizontalVelocities().size() );
  assert( oracle_scene.getVerticalVelocities().size() == testing_scene.getVerticalVelocities().size() );
  
  const ArrayXs& oracle_u = oracle_scene.getHorizontalVelocities();
  const ArrayXs& testing_u = testing_scene.getHorizontalVelocities();
  const ArrayXs& oracle_v = oracle_scene.getVerticalVelocities();
  const ArrayXs& testing_v = testing_scene.getVerticalVelocities();
  

  for( int i = 0; i < oracle_scene.physicalRows() + 2; ++i )
  {
  	for( int j = 0; j < oracle_scene.physicalCols() + 1; ++j )
  	{
	  	scalar u_resid = fabs(oracle_u(i,j) - testing_u(i,j));
	    assert( u_resid >= 0.0 );
	    m_accumulated_u_residual += u_resid;
	    if( u_resid > m_max_u_residual ) m_max_u_residual = u_resid;
	  }
  }
  for( int i = 0; i < oracle_scene.physicalRows() + 1; ++i )
  {
  	for( int j = 0; j < oracle_scene.physicalCols() + 2; ++j )
  	{
	    scalar v_resid = fabs(oracle_v(i,j) - testing_v(i,j));;
	    assert( v_resid >= 0.0 );
	    m_accumulated_v_residual += v_resid;
	    if( v_resid > m_max_v_residual ) m_max_v_residual = v_resid;
  	}
  }

  // After diffusion
  const ArrayXs& oracle_uDiffusion = oracle_scene.getHorizontalVelocitiesDiffusion();
  const ArrayXs& testing_uDiffusion = testing_scene.getHorizontalVelocitiesDiffusion();
  const ArrayXs& oracle_vDiffusion = oracle_scene.getVerticalVelocitiesDiffusion();
  const ArrayXs& testing_vDiffusion = testing_scene.getVerticalVelocitiesDiffusion();

  for( int i = 0; i < oracle_scene.physicalRows() + 2; ++i )
  {
    for( int j = 0; j < oracle_scene.physicalCols() + 1; ++j )
    {
      scalar u_resid = fabs(oracle_uDiffusion(i,j) - testing_uDiffusion(i,j));
      assert( u_resid >= 0.0 );
      m_accumulated_u_diffusion_residual += u_resid;
    }
  }
  for( int i = 0; i < oracle_scene.physicalRows() + 1; ++i )
  {
    for( int j = 0; j < oracle_scene.physicalCols() + 2; ++j )
    {
      scalar v_resid = fabs(oracle_vDiffusion(i,j) - testing_vDiffusion(i,j));;
      assert( v_resid >= 0.0 );
      m_accumulated_v_diffusion_residual += v_resid;
    }
  }

  // After advection
  const ArrayXs& oracle_uAdvect = oracle_scene.getHorizontalVelocitiesAdvect();
  const ArrayXs& testing_uAdvect = testing_scene.getHorizontalVelocitiesAdvect();
  const ArrayXs& oracle_vAdvect = oracle_scene.getVerticalVelocitiesAdvect();
  const ArrayXs& testing_vAdvect = testing_scene.getVerticalVelocitiesAdvect();


  for( int i = 0; i < oracle_scene.physicalRows() + 2; ++i )
  {
    for( int j = 0; j < oracle_scene.physicalCols() + 1; ++j )
    {
      scalar u_resid = fabs(oracle_uAdvect(i,j) - testing_uAdvect(i,j));
      assert( u_resid >= 0.0 );
      m_accumulated_u_advect_residual += u_resid;
    }
  }

  for( int i = 0; i < oracle_scene.physicalRows() + 1; ++i )
  {
    for( int j = 0; j < oracle_scene.physicalCols() + 2; ++j )
    {
      scalar v_resid = fabs(oracle_vAdvect(i,j) - testing_vAdvect(i,j));;
      assert( v_resid >= 0.0 );
      m_accumulated_v_advect_residual += v_resid;
    }
  }
}

scalar StableFluidsSimGrader::getAccumulatedUResidual() const
{
  return m_accumulated_u_residual;
}

scalar StableFluidsSimGrader::getAccumulatedVResidual() const
{
  return m_accumulated_v_residual;
}

scalar StableFluidsSimGrader::getMaxUResidual() const
{
  return m_max_u_residual;
}

scalar StableFluidsSimGrader::getMaxVResidual() const
{
  return m_max_v_residual;
}

bool StableFluidsSimGrader::accumulatedUResidualPassed() const
{
  return m_accumulated_u_residual < m_acceptable_accumulated_u_residual;
}

bool StableFluidsSimGrader::accumulatedVResidualPassed() const
{
  return m_accumulated_v_residual < m_acceptable_accumulated_v_residual;
}

bool StableFluidsSimGrader::maxUResidualPassed() const
{
  return m_max_u_residual < m_acceptable_max_u_residual;
}

bool StableFluidsSimGrader::maxVResidualPassed() const
{
  return m_max_v_residual < m_acceptable_max_v_residual;
}

scalar StableFluidsSimGrader::getAccumulatedUAdvectResidual() const
{
  return m_accumulated_u_advect_residual;
}

scalar StableFluidsSimGrader::getAccumulatedVAdvectResidual() const
{
  return m_accumulated_v_advect_residual;
}

scalar StableFluidsSimGrader::getAccumulatedUDiffusionResidual() const
{
  return m_accumulated_u_diffusion_residual;
}

scalar StableFluidsSimGrader::getAccumulatedVDiffusionResidual() const
{
  return m_accumulated_v_diffusion_residual;
}

