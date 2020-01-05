#include "StableFluidsEnsemble.h"
//#define WITH_MARKERS
#ifdef WITH_MARKERS
#include "StableFluidsSimWithMarkers.h"
#endif

#define N (64)
#define PADDING (0)

StableFluidsEnsemble::StableFluidsEnsemble()
    : m_solids(generateSolidMap(N))
#ifdef WITH_MARKERS
      ,
      m_fluid_sim(new StableFluidsSimWithMarkers(N, N, m_solids))
#else
      ,
      m_fluid_sim(new StableFluidsSim(N, N))
#endif
      ,
      m_bgcolor(0.0, 0.0, 0.0), m_colormap(Colormap::MATLAB_JET, 256)
      //, m_window_width(4*m_fluid_sim->physicalCols())
      //, m_window_height(4*m_fluid_sim->physicalRows())
      ,
      m_window_width(500), m_window_height(500), m_left_drag(false), m_right_drag(false), m_last_row(0), m_last_col(0), m_display_debugging(false), m_display_surface(true), m_particle_cell_map(N + 1, N + 1), frame(0)
{
#ifdef WITH_MARKERS
    m_comparison_fluid_sim = new StableFluidsSimWithMarkers(N, N, m_solids);
#else
    m_comparison_fluid_sim = new StableFluidsSim(N, N);
#endif
    assert(m_comparison_fluid_sim != NULL);
    m_fluid_sim_grader = new StableFluidsSimGrader;

    volume = (GU_PrimVolume *)GU_PrimVolume::build(&gdp);
}

StableFluidsEnsemble::~StableFluidsEnsemble()
{
}

void StableFluidsEnsemble::stepSystem(const scalar &dt)
{
    m_fluid_sim->stepSystem(dt);
    frame++;
}

void StableFluidsEnsemble::keyboard(const unsigned char &key, const int &x, const int &y)
{
    if (key == '=' || key == '+')
    {
        m_colormap.incrementColormap();
    }
    if (key == '-' || key == '_')
    {
        m_colormap.decrementColormap();
    }
    if (key == 'r' || key == 'R')
    {
        m_fluid_sim->clear();
    }
    if (key >= '0' && key <= '9')
    {
        m_fluid_sim->setPrescribedVelocity(key - '0');
    }
    if (key == 'v' || key == 'V')
    {
        m_fluid_sim->setVerbose(!m_fluid_sim->verbose());
    }
    //    if( key == 'd' || key == 'D' )
    //    {
    //        if (m_display_surface)
    //        {
    //            m_display_surface = false;
    //            m_display_debugging = false;
    //        } else
    //        {
    //            if (!m_display_debugging)
    //            {
    //                m_display_debugging = true;
    //            } else
    //            {
    //                m_display_surface = true;
    //            }
    //        }
    //    }
    if (key == 'd' || key == 'D')
    {
        m_display_debugging = !m_display_debugging;
    }
    if (key == 'f' || key == 'F')
    {
        m_display_surface = !m_display_surface;
    }
}

void StableFluidsEnsemble::special(const int &key, const int &x, const int &y)
{
}

void StableFluidsEnsemble::mouse(const int &button, const int &state, const int &x, const int &y)
{
    const int &row = (y * m_fluid_sim->physicalRows()) / m_window_height;
    const int &col = (x * m_fluid_sim->physicalCols()) / m_window_width;

    // Record that the left mouse button was pressed (begin left drag)
    if (!m_right_drag && button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        m_left_drag = true;
        m_last_row = row;
        m_last_col = col;
    }
    // Record that the left mouse button was released (end left drag)
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    {
        m_left_drag = false;
        m_last_row = 0;
        m_last_col = 0;
    }

    // Record that the right mouse button was pressed (begin right drag)
    if (!m_right_drag && button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
    {
        m_right_drag = true;
        m_last_row = row;
        m_last_col = col;
    }
    // Record that the right mouse button was released (end right drag)
    if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP)
    {
        addSphere(row, col, 4);
        m_right_drag = false;
        m_last_row = 0;
        m_last_col = 0;
    }
}

void StableFluidsEnsemble::motion(const int &x, const int &y)
{
    const int &row = (y * m_fluid_sim->physicalRows()) / m_window_height;
    const int &col = (x * m_fluid_sim->physicalCols()) / m_window_width;

    // If a left button drag is in progress
    if (m_left_drag)
    {
        // Points returned by glut could be far apart, so connect with
        //  an unrboken line.
        std::vector<int> xpts;
        std::vector<int> ypts;
        renderingutils::bressenhamLine(m_last_row, m_last_col, row, col, xpts, ypts);
        assert(xpts.size() == ypts.size());

        m_last_row = row;
        m_last_col = col;

        // If there is only one point, we can't define a direction for the force
        if (xpts.size() <= 1)
            return;

        // Compute a vector connecting the two points returned by glut
        Vector2s dir;
        dir << ((scalar)xpts.back()) - ((scalar)xpts.front()), ((scalar)ypts.back()) - ((scalar)ypts.front());
        scalar norm = dir.norm();
        if (norm == 0)
            return;
        dir /= norm;
        // Make the vector magnitude 2
        dir *= 2.0;

        // Add scaled endpoint vectors to each point on the line
        ArrayXs &u = m_fluid_sim->getVerticalVelocities();
        ArrayXs &v = m_fluid_sim->getHorizontalVelocities();
        for (std::vector<int>::size_type i = 0; i < xpts.size(); ++i)
        {
            if (xpts[i] < 0)
                continue;
            if (xpts[i] >= m_fluid_sim->physicalRows())
                continue;
            if (ypts[i] < 0)
                continue;
            if (ypts[i] >= m_fluid_sim->physicalCols())
                continue;

            u(xpts[i] + 1, ypts[i] + 1) += dir.x();
            v(xpts[i] + 1, ypts[i] + 1) += dir.y();
        }
    }

    // If a right button drag is in progress
    if (m_right_drag)
    {
        // Points returned by glut could be far apart, so connect with
        //  an unrboken line.
        std::vector<int> xpts;
        std::vector<int> ypts;
        renderingutils::bressenhamLine(m_last_row, m_last_col, row, col, xpts, ypts);
        assert(xpts.size() == ypts.size());

        // Place a sphere at each point on the line
        for (std::vector<int>::size_type i = 0; i < xpts.size(); ++i)
            addSphere(xpts[i], ypts[i], 4);

        m_last_row = row;
        m_last_col = col;
    }
}

void StableFluidsEnsemble::initializeOpenGL()
{
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glClearColor((GLdouble)m_bgcolor.r(), (GLdouble)m_bgcolor.g(), (GLdouble)m_bgcolor.b(), (GLdouble)1.0);
}

void StableFluidsEnsemble::reshape(const int &w, const int &h)
{
    m_window_width = w;
    m_window_height = h;

    glViewport(0, 0, m_window_width, m_window_height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float padding = PADDING;
    gluOrtho2D(0 - padding, m_fluid_sim->physicalRows() + padding, 0 - padding, m_fluid_sim->physicalCols() + padding);

    // Ensure that the window is square
    if (m_window_width != m_window_height)
    {
        glutReshapeWindow(std::max(m_window_width, m_window_height), std::max(m_window_width, m_window_height));
    }
    // Ensure that the window is at least as big as the number of fluid cells
    if (m_window_width < m_fluid_sim->physicalRows())
    {
        glutReshapeWindow(m_fluid_sim->physicalRows(), m_fluid_sim->physicalRows());
    }
}

void StableFluidsEnsemble::display()
{
    glClear(GL_COLOR_BUFFER_BIT);

#ifdef WITH_MARKERS
    if (m_display_surface)
    {
        display_surface();
    }
    else
    {
        if (m_display_debugging)
            display_cell();
        display_particle();
    }
#else
    save();
    display_cell();
#endif
}

void StableFluidsEnsemble::save()
{
    UT_VoxelArrayWriteHandleF handle = volume->getVoxelWriteHandle();
    const ArrayXs &marker_density = m_fluid_sim->getMarkerDensities();
    handle->size(N, N, N);

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
            {
                handle->setValue(i, j, k, marker_density(i+1,j+1)*5);
            }

    std::ofstream myfile;
    myfile.open(std::to_string(frame) + ".bgeo");
    gdp.save(myfile, 1, NULL);
    myfile.close();
}

void StableFluidsEnsemble::display_cell()
{
    const ArrayXs &marker_density = m_fluid_sim->getMarkerDensities();
#ifdef WITH_MARKERS
    const ArrayXb &fluids = (dynamic_cast<StableFluidsSimWithMarkers *>(m_fluid_sim))->has_fluid();
#else
    ArrayXb fluids(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
        for (int j = 0; j < N + 2; j++)
            fluids(i, j) = true;
#endif

    for (int row = 0; row < m_fluid_sim->physicalRows(); ++row)
        for (int col = 0; col < m_fluid_sim->physicalCols(); ++col)
        {
            // Compute the color of the current grid cell
            scalar r, g, b;
            m_colormap.getColorByDensity(marker_density(row + 1, col + 1), r, g, b);
            if (m_solids(row + 1, col + 1))
                glColor3d(1, 1, 1);
            else if (fluids(row + 1, col + 1))
                glColor3d((GLdouble)r, (GLdouble)g, (GLdouble)b);
            else
                glColor3d(0, 0, 0);

            // Compute the left and right endpoints of the cell
            const scalar &xl = col;
            const scalar &xr = col + 1;
            assert(xl <= xr);

            // Compute the bottom and top endpoints of the cell
            const scalar &yb = m_fluid_sim->physicalRows() - row - 1;
            const scalar &yu = m_fluid_sim->physicalRows() - row;
            assert(yb <= yu);

            glBegin(GL_QUADS);
            glVertex2d((GLdouble)xl, (GLdouble)yb);
            glVertex2d((GLdouble)xr, (GLdouble)yb);
            glVertex2d((GLdouble)xr, (GLdouble)yu);
            glVertex2d((GLdouble)xl, (GLdouble)yu);
            glEnd();
        }

    if (m_display_debugging)
    {
        // fluid boundary
        glColor4f(1.0, 1.0, 1.0, 1.0);
        glBegin(GL_LINE_STRIP);
        glVertex2d(0, 0);
        glVertex2d(N, 0);
        glVertex2d(N, N);
        glVertex2d(0, N);
        glVertex2d(0, 0);
        glEnd();

        // grid
        glColor4f(0.5, 0.5, 0.5, 1.0);
        glBegin(GL_LINES);
        for (int i = 1; i < N; i++)
        {
            glVertex2d(0, i);
            glVertex2d(N, i);
            glVertex2d(i, 0);
            glVertex2d(i, N);
        }
        glEnd();

        // velocity
        scalar s = 1.0;
        glBegin(GL_LINES);
        for (int i = 1; i <= N; i++)
            for (int j = 0; j <= N; j++)
            {
                glColor4f(0.0, 0.6, 0.0, 1.0);
                glVertex2d(j, (N + 1 - i) - 0.5);
                glVertex2d(j + m_fluid_sim->getHorizontalVelocities()(i, j) * s, (N + 1 - i) - 0.5);
            }
        for (int i = 0; i <= N; i++)
            for (int j = 1; j <= N; j++)
            {
                glColor4f(0.6, 0.0, 0.0, 1.0);
                glVertex2d(j - 0.5, (N - i));
                glVertex2d(j - 0.5, (N - i) - m_fluid_sim->getVerticalVelocities()(i, j) * s);
            }
        glEnd();
    }
}

void StableFluidsEnsemble::display_particle()
{
    const ArrayXs &marker_density = m_fluid_sim->getMarkerDensities();
#ifdef WITH_MARKERS
    const std::vector<Vector2s> &particles = (dynamic_cast<StableFluidsSimWithMarkers *>(m_fluid_sim))->particles();
#else
    const std::vector<Vector2s> &particles = std::vector<Vector2s>();
#endif

    for (size_t i = 0; i < particles.size(); i++)
    {
        int row = (int)(particles[i].x() * N + 0.5);
        int col = (int)(particles[i].y() * N + 0.5);
        assert(row >= 0);
        assert(row <= N + 1);
        assert(col >= 0);
        assert(col <= N + 1);
        //std::cout << row << ", " << col << ", " << particles[i].x() << std::endl;

        // Compute the color of the current grid cell
        scalar r, g, b;
        m_colormap.getColorByDensity(marker_density(row, col), r, g, b);

        // Compute the left and right endpoints of the cell
        const scalar radius = (m_display_debugging ? 0.1 : 1.0);
        const scalar xc = particles[i].y() * N - 0.5;
        const scalar xl = xc - radius;
        const scalar xr = xc + radius;

        // Compute the bottom and top endpoints of the cell
        const scalar yc = m_fluid_sim->physicalRows() - (particles[i].x() * N - 0.5);
        const scalar yb = yc - radius;
        const scalar yt = yc + radius;

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, m_display_debugging ? GL_ONE : GL_ONE_MINUS_SRC_ALPHA);

        const int nv = 12;
        glBegin(GL_TRIANGLE_FAN);
        glColor4d(r, g, b, 1);
        glVertex2d(xc, yc);
        glColor4d(r, g, b, m_display_debugging ? 1 : 0);
        glVertex2d(xc + radius, yc);
        for (int i = 0; i < nv; i++)
        {
            glVertex2d(xc + radius * cos((i + 1) * 2 * 3.141593 / nv), yc + radius * sin((i + 1) * 2 * 3.141593 / nv));
        }
        glEnd();
        glDisable(GL_BLEND);
    }

    if (m_display_debugging)
    {
        // fluid boundary
        glColor4f(1.0, 1.0, 1.0, 1.0);
        glBegin(GL_LINE_STRIP);
        glVertex2d(0, 0);
        glVertex2d(N, 0);
        glVertex2d(N, N);
        glVertex2d(0, N);
        glVertex2d(0, 0);
        glEnd();

        // grid
        glColor4f(0.5, 0.5, 0.5, 1.0);
        glBegin(GL_LINES);
        for (int i = 1; i < N; i++)
        {
            glVertex2d(0, i);
            glVertex2d(N, i);
            glVertex2d(i, 0);
            glVertex2d(i, N);
        }
        glEnd();

        // velocity
        scalar s = 1.0;
        glBegin(GL_LINES);
        for (int i = 1; i <= N; i++)
            for (int j = 0; j <= N; j++)
            {
                glColor4f(0.0, 0.6, 0.0, 1.0);
                glVertex2d(j, (N + 1 - i) - 0.5);
                glVertex2d(j + m_fluid_sim->getHorizontalVelocities()(i, j) * s, (N + 1 - i) - 0.5);
            }
        for (int i = 0; i <= N; i++)
            for (int j = 1; j <= N; j++)
            {
                glColor4f(0.6, 0.0, 0.0, 1.0);
                glVertex2d(j - 0.5, (N - i));
                glVertex2d(j - 0.5, (N - i) - m_fluid_sim->getVerticalVelocities()(i, j) * s);
            }
        glEnd();
    }
}

void StableFluidsEnsemble::display_surface()
{
    const ArrayXs &marker_density = m_fluid_sim->getMarkerDensities();
#ifdef WITH_MARKERS
    const std::vector<Vector2s> &particles = (dynamic_cast<StableFluidsSimWithMarkers *>(m_fluid_sim))->particles();
    const ArrayXb &fluids = (dynamic_cast<StableFluidsSimWithMarkers *>(m_fluid_sim))->has_fluid();
#else
    const std::vector<Vector2s> &particles = std::vector<Vector2s>();
    ArrayXb fluids(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
        for (int j = 0; j < N + 2; j++)
            fluids(i, j) = true;
#endif

    scalar r = 0.5;
    scalar h = r * 3;
    scalar hcull = h * 2;
    scalar tao = (1 - r * r / h / h) * (1 - r * r / h / h) * (1 - r * r / h / h);

    // extract the free surface from particles
    // reference: R. Bridson 2008, "Fluid Simulation for Computer Graphics" p84.
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            m_particle_cell_map(i, j).clear();

    for (size_t p = 0; p < particles.size(); p++)
    {
        int i = (int)(particles[p].x() * N - 0.5);
        int j = (int)(particles[p].y() * N - 0.5);
        assert(i >= 0 && i < N);
        assert(j >= 0 && j < N);
        m_particle_cell_map(i, j).push_back(p);
        m_particle_cell_map(i + 1, j).push_back(p);
        m_particle_cell_map(i, j + 1).push_back(p);
        m_particle_cell_map(i + 1, j + 1).push_back(p);
    }

    ArrayXs implicit_surface(N + 1, N + 1);
    for (int i = 0; i <= N; i++)
    {
        for (int j = 0; j <= N; j++)
        {
            Vector2s node_pos = Vector2s(i + 0.5, j + 0.5);

            scalar is = 0;
            for (size_t p = 0; p < m_particle_cell_map(i, j).size(); p++)
            {
                Vector2s par_pos = particles[m_particle_cell_map(i, j)[p]] * N;
                scalar dsq = (node_pos - par_pos).dot(node_pos - par_pos) / h / h;
                is += (dsq > 1 ? 0 : (1 - dsq) * (1 - dsq) * (1 - dsq));
            }
            implicit_surface(i, j) = is - tao;
        }
    }

    glBegin(GL_TRIANGLES);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            scalar r, g, b;
            m_colormap.getColorByDensity(marker_density(i + 1, j + 1), r, g, b);
            if (m_solids(i + 1, j + 1))
                glColor3d(1, 1, 1);
            else
                glColor3d((GLdouble)r, (GLdouble)g, (GLdouble)b);

            // the positive pattern of the four corners of a cell
            int positive_pattern = 0;
            scalar lt = implicit_surface(i, j);
            scalar rt = implicit_surface(i, j + 1);
            scalar lb = implicit_surface(i + 1, j);
            scalar rb = implicit_surface(i + 1, j + 1);
            positive_pattern = (lt > 0 ? 0x8 : 0x0) + (rt > 0 ? 0x4 : 0x0) + (lb > 0 ? 0x2 : 0x0) + (rb > 0 ? 0x1 : 0x0);
            if (m_solids(i + 1, j + 1))
                positive_pattern = 0xF;

            scalar xl = j;
            scalar xr = j + 1;
            scalar yt = N - i;
            scalar yb = N - (i + 1);

            scalar xt = xl - lt / (rt - lt);
            scalar xb = xl - lb / (rb - lb);
            scalar yl = yt + lt / (lb - lt);
            scalar yr = yt + rt / (rb - rt);

            switch (positive_pattern)
            {
            case 0x0: // 0000
                // nothing to do
                break;
            case 0x1: // 0001
                // bottom-right triangle
                glVertex2d(xr, yr);
                glVertex2d(xb, yb);
                glVertex2d(xr, yb);
                break;
            case 0x2: // 0010
                // bottom-left triangle
                glVertex2d(xb, yb);
                glVertex2d(xl, yl);
                glVertex2d(xl, yb);
                break;
            case 0x3: // 0011
                // bottom half
                glVertex2d(xr, yr);
                glVertex2d(xl, yl);
                glVertex2d(xr, yb);
                glVertex2d(xl, yl);
                glVertex2d(xr, yb);
                glVertex2d(xl, yb);
                break;
            case 0x4: // 0100
                // top-right triangle
                glVertex2d(xt, yt);
                glVertex2d(xr, yr);
                glVertex2d(xr, yt);
                break;
            case 0x5: // 0101
                // right half
                glVertex2d(xt, yt);
                glVertex2d(xb, yb);
                glVertex2d(xr, yb);
                glVertex2d(xr, yb);
                glVertex2d(xr, yt);
                glVertex2d(xt, yt);
                break;
            case 0x6: // 0110
                // top-right triangle and bottom-left triangle
                glVertex2d(xt, yt);
                glVertex2d(xr, yr);
                glVertex2d(xr, yt);
                glVertex2d(xb, yb);
                glVertex2d(xl, yl);
                glVertex2d(xl, yb);
                break;
            case 0x7: // 0111
                // complement of top-left triangle
                glVertex2d(xr, yb);
                glVertex2d(xl, yl);
                glVertex2d(xl, yb);

                glVertex2d(xr, yb);
                glVertex2d(xt, yt);
                glVertex2d(xl, yl);

                glVertex2d(xr, yb);
                glVertex2d(xr, yt);
                glVertex2d(xt, yt);
                break;
            case 0x8: // 1000
                // top-left triangle
                glVertex2d(xl, yl);
                glVertex2d(xt, yt);
                glVertex2d(xl, yt);
                break;
            case 0x9: // 1001
                // top-left triangle, and bottom-right triangle
                glVertex2d(xl, yl);
                glVertex2d(xt, yt);
                glVertex2d(xl, yt);

                glVertex2d(xr, yr);
                glVertex2d(xb, yb);
                glVertex2d(xr, yb);
                break;
            case 0xA: // 1010
                // left half
                glVertex2d(xb, yb);
                glVertex2d(xt, yt);
                glVertex2d(xl, yt);

                glVertex2d(xl, yt);
                glVertex2d(xl, yb);
                glVertex2d(xb, yb);
                break;
            case 0xB: // 1011
                // complement of top-right triangle
                glVertex2d(xl, yb);
                glVertex2d(xr, yb);
                glVertex2d(xr, yr);

                glVertex2d(xl, yb);
                glVertex2d(xr, yr);
                glVertex2d(xt, yt);

                glVertex2d(xl, yb);
                glVertex2d(xt, yt);
                glVertex2d(xl, yt);
                break;
            case 0xC: // 1100
                // top half
                glVertex2d(xl, yl);
                glVertex2d(xr, yr);
                glVertex2d(xr, yt);

                glVertex2d(xr, yt);
                glVertex2d(xl, yt);
                glVertex2d(xl, yl);
                break;
            case 0xD: // 1101
                // complement to bottom-left triangle
                glVertex2d(xr, yt);
                glVertex2d(xl, yt);
                glVertex2d(xl, yl);

                glVertex2d(xr, yt);
                glVertex2d(xl, yl);
                glVertex2d(xb, yb);

                glVertex2d(xr, yt);
                glVertex2d(xb, yb);
                glVertex2d(xr, yb);
                break;
            case 0xE: // 1110
                // complement to bottom-right triangle
                glVertex2d(xl, yt);
                glVertex2d(xl, yb);
                glVertex2d(xb, yb);

                glVertex2d(xl, yt);
                glVertex2d(xb, yb);
                glVertex2d(xr, yr);

                glVertex2d(xl, yt);
                glVertex2d(xr, yr);
                glVertex2d(xr, yt);
                break;
            case 0xF: // 1111
                // whole cell
                glVertex2d(xl, yb);
                glVertex2d(xr, yb);
                glVertex2d(xr, yt);

                glVertex2d(xr, yt);
                glVertex2d(xl, yt);
                glVertex2d(xl, yb);
                break;
            }
        }
    }
    glEnd();

    if (m_display_debugging)
    {
        // implicit surface values
        glBegin(GL_QUADS);
        for (int i = 0; i <= N; i++)
        {
            for (int j = 0; j <= N; j++)
            {
                if (implicit_surface(i, j) > 0)
                    glColor4d(0, 1, 0, 1);
                else
                    glColor4d(1, 0, 0, 1);

                scalar r = 0.2;
                glVertex2d(j - r, N - i - r);
                glVertex2d(j + r, N - i - r);
                glVertex2d(j + r, N - i + r);
                glVertex2d(j - r, N - i + r);
            }
        }
        glEnd();
    }

    // fluid boundary
    glColor4f(1.0, 1.0, 1.0, 1.0);
    glBegin(GL_LINE_STRIP);
    glVertex2d(0, 0);
    glVertex2d(N, 0);
    glVertex2d(N, N);
    glVertex2d(0, N);
    glVertex2d(0, 0);
    glEnd();

    if (m_display_debugging)
    {
        // grid
        glColor4f(0.5, 0.5, 0.5, 1.0);
        glBegin(GL_LINES);
        for (int i = 1; i < N; i++)
        {
            glVertex2d(0, i);
            glVertex2d(N, i);
            glVertex2d(i, 0);
            glVertex2d(i, N);
        }
        glEnd();
    }
}

unsigned int StableFluidsEnsemble::getGlutDisplayMode() const
{
    return GLUT_DOUBLE | GLUT_RGBA;
}

const int &StableFluidsEnsemble::getWindowWidth() const
{
    return m_window_width;
}

const int &StableFluidsEnsemble::getWindowHeight() const
{
    return m_window_height;
}

const renderingutils::Color &StableFluidsEnsemble::getBackgroundColor() const
{
    return m_bgcolor;
}

void StableFluidsEnsemble::addSphere(const int &row, const int &col, const int &R)
{
    const scalar &scalarR = (scalar)R;

    ArrayXs &marker_densities = m_fluid_sim->getMarkerDensities();

    for (int i = -R; i <= R; ++i)
        for (int j = -R; j <= R; ++j)
        {
            // If the sphere extends off the grid, nothing to do
            if (row + i < 0)
                continue;
            if (row + i >= m_fluid_sim->physicalRows())
                continue;
            if (col + j < 0)
                continue;
            if (col + j >= m_fluid_sim->physicalCols())
                continue;

            assert(row + i >= 0);
            assert(row + i < m_fluid_sim->physicalRows());
            assert(col + j >= 0);
            assert(col + j < m_fluid_sim->physicalCols());

            // Compute the squared distance from the center
            const scalar &x = (scalar)i;
            const scalar &y = (scalar)j;
            scalar zz = scalarR * scalarR - x * x - y * y;
            if (zz < 0.0)
                continue;
            // Compute and normalize the distance
            scalar z = sqrt(zz) / R;
            assert(z >= 0.0);
            assert(z <= 1.0);

            // If the update will increase the value, do it
            marker_densities(row + i + 1, col + j + 1) = std::max(z, marker_densities(row + i + 1, col + j + 1));
        }
}

ArrayXb StableFluidsEnsemble::generateSolidMap(int n)
{
    ArrayXb solids = ArrayXb::Zero(n + 2, n + 2);
    scalar cx = 0.5;
    scalar cy = 0.6;
    scalar r = 0.2;
    for (int i = 0; i < n + 2; i++)
    {
        for (int j = 0; j < n + 2; j++)
        {
            if (i == 0 || i == n + 1 || j == 0 || j == n + 1)
                solids(i, j) = true;

#ifdef WITH_MARKERS
            scalar x = (scalar)j / n;
            scalar y = (scalar)i / n;
            if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= r * r)
                solids(i, j) = true;
#endif
        }
    }

    return solids;
}

// Scene Benchmarking Functions

void StableFluidsEnsemble::updateSceneComparison()
{
    assert(m_fluid_sim != NULL);
    assert(m_comparison_fluid_sim != NULL);
    assert(m_fluid_sim_grader != NULL);

    // position/velocity residuals
    m_fluid_sim_grader->addToAccumulatedResidual(*m_fluid_sim, *m_comparison_fluid_sim);
}

void StableFluidsEnsemble::printErrorInformation(bool print_pass)
{
    assert(m_fluid_sim_grader != NULL);

    scalar accumulatedUResidual = m_fluid_sim_grader->getAccumulatedUResidual();
    scalar accumulatedVResidual = m_fluid_sim_grader->getAccumulatedVResidual();
    scalar maxUResidual = m_fluid_sim_grader->getMaxUResidual();
    scalar maxVResidual = m_fluid_sim_grader->getMaxVResidual();

    scalar accumulatedUDiffusionResidual = m_fluid_sim_grader->getAccumulatedUDiffusionResidual();
    scalar accumulatedVDiffusionResidual = m_fluid_sim_grader->getAccumulatedVDiffusionResidual();
    scalar accumulatedUAdvectResidual = m_fluid_sim_grader->getAccumulatedUAdvectResidual();
    scalar accumulatedVAdvectResidual = m_fluid_sim_grader->getAccumulatedVAdvectResidual();

    if (!print_pass)
    {
        std::cout << "Accumulated U Residual after diffusion: " << accumulatedUDiffusionResidual << std::endl;
        std::cout << "Accumulated V Residual after diffusion: " << accumulatedVDiffusionResidual << std::endl;
        std::cout << "Accumulated U Residual after advection: " << accumulatedUAdvectResidual << std::endl;
        std::cout << "Accumulated V Residual after advection: " << accumulatedVAdvectResidual << std::endl;
        std::cout << "Accumulated U Residual: " << accumulatedUResidual << std::endl;
        std::cout << "Accumulated V Residual: " << accumulatedVResidual << std::endl;
        std::cout << "Maximum U Residual: " << maxUResidual << std::endl;
        std::cout << "Maximum V Residual: " << maxVResidual << std::endl;

        std::cout << outputmod::startpink << "FOSSSim message: " << outputmod::endpink;
        std::cout << "Simulation did not run to completion. No pass/fail status will be reported." << std::endl;
    }
    else
    {
        bool accumulatedUResidualAcceptable = m_fluid_sim_grader->accumulatedUResidualPassed();
        bool accumulatedVResidualAcceptable = m_fluid_sim_grader->accumulatedVResidualPassed();
        bool maxUResidualAcceptable = m_fluid_sim_grader->maxUResidualPassed();
        bool maxVResidualAcceptable = m_fluid_sim_grader->maxVResidualPassed();

        bool total_success = accumulatedUResidualAcceptable && accumulatedVResidualAcceptable && maxUResidualAcceptable && maxVResidualAcceptable;

        std::cout << "Accumulated U Residual after diffusion: " << accumulatedUDiffusionResidual << std::endl;
        std::cout << "Accumulated V Residual after diffusion: " << accumulatedVDiffusionResidual << std::endl;
        std::cout << "Accumulated U Residual after advection: " << accumulatedUAdvectResidual << std::endl;
        std::cout << "Accumulated V Residual after advection: " << accumulatedVAdvectResidual << std::endl;
        std::cout << "Accumulated U Residual: " << accumulatedUResidual << "   " << (accumulatedUResidualAcceptable ? "Passed." : "Failed.") << std::endl;
        std::cout << "Accumulated V Residual: " << accumulatedVResidual << "   " << (accumulatedVResidualAcceptable ? "Passed." : "Failed.") << std::endl;
        std::cout << "Maximum U Residual: " << maxUResidual << "   " << (maxUResidualAcceptable ? "Passed." : "Failed.") << std::endl;
        std::cout << "Maximum V Residual: " << maxVResidual << "   " << (maxVResidualAcceptable ? "Passed." : "Failed.") << std::endl;

        std::cout << "Overall success: " << (total_success ? "Passed." : "Failed.") << std::endl;

        std::cout << outputmod::startpink << "FOSSSim message: " << outputmod::endpink << "Saved residual to residual.txt." << std::endl;

        std::ofstream residfile;
        residfile.open("residual.txt");
        residfile << accumulatedUResidualAcceptable << "   " << accumulatedUResidual << std::endl;
        residfile << accumulatedVResidualAcceptable << "   " << accumulatedVResidual << std::endl;
        residfile << maxUResidualAcceptable << "   " << maxUResidual << std::endl;
        residfile << maxVResidualAcceptable << "   " << maxVResidual << std::endl;
        residfile.close();
    }
}

/////////////////////////////////////////////////////////////////////////////
// Serialization Functions

void StableFluidsEnsemble::copyComparisonSceneToScene()
{
    assert(m_fluid_sim != NULL);
    assert(m_comparison_fluid_sim != NULL);
#ifdef WITH_MARKERS
    ((StableFluidsSimWithMarkers *)m_fluid_sim)->copyState(*((StableFluidsSimWithMarkers *)m_comparison_fluid_sim));
#else
    m_fluid_sim->copyState(*m_comparison_fluid_sim);
#endif
}

void StableFluidsEnsemble::serializeScene(std::ofstream &outputstream)
{
    assert(m_fluid_sim != NULL);
    m_fluid_sim_serializer.serializeScene(*m_fluid_sim, outputstream);
}

void StableFluidsEnsemble::loadComparisonScene(std::ifstream &inputstream)
{
    // #ifdef WITH_MARKERS
    // m_comparison_fluid_sim = new StableFluidsSimWithMarkers(N,N,m_solids);
    // #else
    // m_comparison_fluid_sim = new StableFluidsSim(N,N));
    // #endif
    assert(m_comparison_fluid_sim != NULL);
    // m_fluid_sim_grader = new StableFluidsSimGrader;
    m_fluid_sim_serializer.loadScene(*m_comparison_fluid_sim, inputstream);
}

void StableFluidsEnsemble::setVelocityPattern(int velocity_pattern)
{
    assert(m_fluid_sim != NULL);
    m_fluid_sim->setPrescribedVelocity(velocity_pattern);
}

void StableFluidsEnsemble::setDiffusion(scalar diff)
{
    assert(m_fluid_sim != NULL);
    m_fluid_sim->setDiffusion(diff);
}

void StableFluidsEnsemble::setViscosity(scalar visc)
{
    assert(m_fluid_sim != NULL);
    m_fluid_sim->setViscosity(visc);
}

void StableFluidsEnsemble::setUseAdvect(bool use_advect)
{
    assert(m_fluid_sim != NULL);
    m_fluid_sim->setUseAdvect(use_advect);
}

void StableFluidsEnsemble::setUseProj(bool use_proj)
{
    assert(m_fluid_sim != NULL);
    m_fluid_sim->setUseProj(use_proj);
}
