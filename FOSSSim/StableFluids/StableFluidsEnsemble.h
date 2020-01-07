#ifndef __STABLE_FLUIDS_ENSEMBLE_H__
#define __STABLE_FLUIDS_ENSEMBLE_H__

#include "FOSSSim/SimulationEnsemble.h"
#include "FOSSSim/MathUtilities.h"
#include "FOSSSim/RenderingUtilities.h"

#include "StableFluidsSim.h"
#include "Colormap.h"
#include "StableFluidsSimGrader.h"
#include "StableFluidsSimSerializer.h"

// #include <CMD/CMD_Args.h>
// #include <UT/UT_Assert.h>
// #include <GU/GU_Detail.h>
// #include <GU/GU_PrimVolume.h>

class StableFluidsEnsemble : public SimulationEnsemble
{
public:
    StableFluidsEnsemble();

    virtual ~StableFluidsEnsemble();

    virtual void stepSystem(const scalar &dt);

    /////////////////////////////////////////////////////////////////////////////
    // Mouse and keyboard input handlers

    virtual void keyboard(const unsigned char &key, const int &x, const int &y);

    virtual void special(const int &key, const int &x, const int &y);

    virtual void mouse(const int &button, const int &state, const int &x, const int &y);

    virtual void motion(const int &x, const int &y);

    /////////////////////////////////////////////////////////////////////////////
    // Functions for use by OpenGL and GLUT

    virtual void initializeOpenGL();

    virtual void reshape(const int &w, const int &h);

    virtual void display();
    virtual void display_cell();
    virtual void display_particle();
    virtual void display_surface();

    virtual unsigned int getGlutDisplayMode() const;

    virtual const int &getWindowWidth() const;

    virtual const int &getWindowHeight() const;

    virtual const renderingutils::Color &getBackgroundColor() const;

    /////////////////////////////////////////////////////////////////////////////
    // Scene Benchmarking Functions

    virtual void updateSceneComparison();
    virtual void printErrorInformation(bool print_pass);

    /////////////////////////////////////////////////////////////////////////////
    // Serialization Functions

    virtual void copyComparisonSceneToScene();
    virtual void serializeScene(std::ofstream &outputstream);
    virtual void loadComparisonScene(std::ifstream &inputstream);

    virtual void setVelocityPattern(int velocity_pattern);
    virtual void setDiffusion(scalar diff);
    virtual void setViscosity(scalar visc);
    virtual void setUseAdvect(bool use_advect);
    virtual void setUseProj(bool use_proj);

    // virtual void save();

private:
    void addSphere(const int &row, const int &col, const int &R);
    ArrayXb generateSolidMap(int n);

    ArrayXb m_solids;
    StableFluidsSim *m_fluid_sim;

    renderingutils::Color m_bgcolor;

    Colormap m_colormap;

    int m_window_width;
    int m_window_height;

    bool m_left_drag;
    bool m_right_drag;

    int m_last_row;
    int m_last_col;

    bool m_display_debugging;
    bool m_display_surface;

    Eigen::Array<std::vector<int>, Eigen::Dynamic, Eigen::Dynamic> m_particle_cell_map;

    StableFluidsSim *m_comparison_fluid_sim;

    StableFluidsSimSerializer m_fluid_sim_serializer;

    StableFluidsSimGrader *m_fluid_sim_grader;

    // GU_Detail gdp;
    // GU_PrimVolume *volume;
    // int frame;
};

#endif
