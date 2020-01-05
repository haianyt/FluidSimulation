#ifndef __STABLE_FLUIDS_SIM_GRADER_H__
#define __STABLE_FLUIDS_SIM_GRADER_H__

#include "StableFluidsSim.h"

class StableFluidsSimGrader
{
public:
  
  StableFluidsSimGrader();
  
  void addToAccumulatedResidual( const StableFluidsSim& oracle_scene, const StableFluidsSim& testing_scene );

  scalar getAccumulatedUResidual() const;
  scalar getMaxUResidual() const;
  
  bool accumulatedUResidualPassed() const;
  bool maxUResidualPassed() const;

  scalar getAccumulatedVResidual() const;
  scalar getMaxVResidual() const;
  
  bool accumulatedVResidualPassed() const;
  bool maxVResidualPassed() const;  


  scalar getAccumulatedUDiffusionResidual() const;
  scalar getAccumulatedVDiffusionResidual() const;

  scalar getAccumulatedUAdvectResidual() const;
  scalar getAccumulatedVAdvectResidual() const;

  
private:
  scalar m_accumulated_u_residual;
  scalar m_max_u_residual;
  scalar m_acceptable_accumulated_u_residual;
  scalar m_acceptable_max_u_residual;

  scalar m_accumulated_v_residual;
  scalar m_max_v_residual;
  scalar m_acceptable_accumulated_v_residual;
  scalar m_acceptable_max_v_residual;

  // After diffusion
  scalar m_accumulated_u_diffusion_residual;
  scalar m_accumulated_v_diffusion_residual;
  // After advection
  scalar m_accumulated_u_advect_residual;
  scalar m_accumulated_v_advect_residual;
};

#endif
