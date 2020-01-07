#include "StableFluidsSim.h"
#include <Eigen/LU>

//#define VERBOSE (0)

StableFluidsSim::StableFluidsSim( const int& rows, const int& cols, const scalar& diff, const scalar& visc, bool use_advect, bool use_proj)
: m_diff(diff)
, m_visc(visc)
, m_use_advect(use_advect)
, m_use_project(use_proj)
, m_N(rows)
, m_d(m_N + 2, m_N + 2, m_N + 2)
, m_u(m_N + 2, m_N + 1, m_N + 2)
, m_v(m_N + 1, m_N + 2, m_N + 2)
, m_w(m_N + 2, m_N + 2, m_N + 1)
, m_uAfterDiffusion(m_N + 2, m_N + 1, m_N + 2)
, m_vAfterDiffusion(m_N + 1, m_N + 2, m_N + 2)
, m_wAfterDiffusion(m_N + 2, m_N + 2, m_N + 1)
, m_uAfterAdvect(m_N + 2, m_N + 1, m_N + 2)
, m_vAfterAdvect(m_N + 1, m_N + 2, m_N + 2)
, m_wAfterAdvect(m_N + 2, m_N + 2, m_N + 1)
, VERBOSE(false)
, m_all_ones(m_N + 2, m_N + 2, m_N + 2)
, frame(0)
{
  assert(rows==cols);

  clear();
  volume = (GU_PrimVolume *)GU_PrimVolume::build(&gdp);

}

StableFluidsSim::~StableFluidsSim()
{

}

////////////////////////////////////////////////////////////////////////////////////
void StableFluidsSim::SWAP(ArrayXs *& x1, ArrayXs *& x2)
{
  ArrayXs * t = x1;
  x1 = x2;
  x2 = t;
}

void StableFluidsSim::add_source(int N, ArrayXs * x, ArrayXs * x0, scalar dt)
{
  // nothing to do. adding marker/velocity is handled in StableFluidsEnsemble.
  *x = *x0;
}

void StableFluidsSim::diffuseD(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
    // assert((*x0 == *x0).all());

    scalar a = diff * dt * N * N;
    *x = *x0;

    for (int k = 0; k < 30; k++) {
      for(int z = 1; z<=N; z++)
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) { // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
                // STUDENTS: You will certainly need code here, do diffuse for ([1, N], [1, N])
                // Gauss-Seidel relaxation
                (*x)(i, j, z) = ((*x0)(i, j, z) + a * (
                  (*x)(i-1, j, z) + (*x)(i+1, j, z) + 
                  (*x)(i, j-1, z) + (*x)(i, j+1, z) +
                  (*x)(i, j, z-1) + (*x)(i, j, z+1)
                  )) / (1 + 6 * a);
            }
        }
    }
}

void StableFluidsSim::diffuseU(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
    // assert((*x0 == *x0).all());

    scalar a = diff * dt * N * N;
    *x = *x0;

    for (int k = 0; k < 30; k++) {
      for (int z = 1; z <= N; z++)
        for (int i = 1; i <= N; i++) {
            for (int j = 0; j <= N; j++) { // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
                // STUDENTS: You will certainly need code here, do diffuse for ([1, N], [0, N]), note the case when (j == 0) or (j == N) need special treatment
                if (j == 0)
                    (*x)(i, j, z) = ((*x0)(i, j, z) + a * ((*x)(i-1, j, z) + (*x)(i+1, j, z) + 
                                                    (*x)(i,j,z-1) + (*x)(i,j,z+1) +
                                                    (*x)(i, j+1, z))) / (1 + 5 * a);
                else if (j == N)
                    (*x)(i, j, z) = ((*x0)(i, j, z) + a * ((*x)(i-1, j, z) + (*x)(i+1, j, z) + 
                                                    (*x)(i,j,z-1) + (*x)(i,j,z+1) +
                                                    (*x)(i, j-1, z))) / (1 + 5 * a);
                else
                    (*x)(i, j, z) = ((*x0)(i, j, z) + a * ((*x)(i-1, j, z) + (*x)(i+1, j, z) + 
                                                    (*x)(i, j, z-1) + (*x)(i, j, z+1) +
                                                    (*x)(i, j-1, z) + (*x)(i, j+1, z))) / (1 + 6 * a);
            }
        }
    }
}

void StableFluidsSim::diffuseV(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
    // assert((*x0 == *x0).all());

    scalar a = diff * dt * N * N;
    *x = *x0;

    for (int k = 0; k < 30; k++) {
      for (int z = 1; z <= N; z++)
        for (int i = 0; i <= N; i++) {
            for (int j = 1; j <= N; j++) { // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
                // STUDENTS: You will certainly need code here, do diffuse for ([1, N], [0, N]), note the case when (j == 0) or (j == N) need special treatment
                if (i == 0)
                    (*x)(i, j, z) = ((*x0)(i, j, z) + a * ((*x)(i+1, j, z) + 
                                                    (*x)(i, j, z-1) + (*x)(i, j, z+1) +
                                                    (*x)(i, j-1, z) + (*x)(i, j+1, z))) / (1 + 5 * a);
                else if (i == N)
                    (*x)(i, j, z) = ((*x0)(i, j, z) + a * ((*x)(i-1, j, z) +
                                                    (*x)(i, j, z-1) + (*x)(i, j, z+1) +
                                                    (*x)(i, j-1, z) + (*x)(i, j+1, z))) / (1 + 5 * a);
                else
                    (*x)(i, j, z) = ((*x0)(i, j, z) + a * ((*x)(i-1, j, z) + (*x)(i+1, j, z) + 
                                                    (*x)(i, j, z-1) + (*x)(i, j, z+1) +
                                                    (*x)(i, j-1, z) + (*x)(i, j+1, z))) / (1 + 6 * a);
            }
        }
    }
}

void StableFluidsSim::diffuseW(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
    // assert((*x0 == *x0).all());

    scalar a = diff * dt * N * N;
    *x = *x0;

    for (int k = 0; k < 30; k++) {
      for (int z = 0; z <= N; z++)
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) { // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
                // STUDENTS: You will certainly need code here, do diffuse for ([1, N], [0, N]), note the case when (j == 0) or (j == N) need special treatment
                if (z == 0)
                    (*x)(i, j, z) = ((*x0)(i, j, z) + a * ((*x)(i-1, j, z) + (*x)(i+1, j, z) + 
                                                    (*x)(i, j, z+1) +
                                                    (*x)(i, j-1, z) + (*x)(i, j+1, z))) / (1 + 5 * a);
                else if (z == N)
                    (*x)(i, j, z) = ((*x0)(i, j, z) + a * ((*x)(i-1, j, z) + (*x)(i+1, j, z) + 
                                                    (*x)(i, j, z-1) +
                                                    (*x)(i, j-1, z) + (*x)(i, j+1, z))) / (1 + 5 * a);
                else
                    (*x)(i, j, z) = ((*x0)(i, j, z) + a * ((*x)(i-1, j, z) + (*x)(i+1, j, z) + 
                                                    (*x)(i, j, z-1) + (*x)(i, j, z+1) +
                                                    (*x)(i, j-1, z) + (*x)(i, j+1, z))) / (1 + 6 * a);
            }
        }
    }
}

void StableFluidsSim::advectD(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, ArrayXs * w, scalar dt)
{
    // assert((*x0 == *x0).all());
    // assert((*u == *u).all());
    // assert((*v == *v).all());

    // STUDENTS: You will certainly need code here, advect for ([1, N], [1, N])
    for (int z = 1; z <= N; z++)
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            scalar ii = i - dt * N * interpolateV(v, i, j, z);
            scalar jj = j - dt * N * interpolateU(u, i, j, z);
            scalar zz = z - dt * N * interpolateW(w, i, j, z);
            (*x)(i, j, z) = interpolateD(x0, ii, jj, zz);
        }
    }
}

// Real time fluid dynamics for games page 8
scalar StableFluidsSim::interpolateD(ArrayXs * d, scalar i, scalar j, scalar z)
{
    // STUDENTS: You will certainly need code here, note the indices should be CLAMP-ed to [0, m_N], since we have to use (i + 1) and (j + 1)
    int i0 = CLAMP((int)i, 0, m_N);
    int j0 = CLAMP((int)j, 0, m_N);
    int z0 = CLAMP((int)z, 0, m_N);

    int i1 = i0 + 1, j1 = j0 + 1, z1 = z0 + 1;

    scalar s1 = CLAMP(i - i0, 0, 1);
    scalar s0 = 1 - s1;
    scalar t1 = CLAMP(j - j0, 0, 1);
    scalar t0 = 1 - t1;
    scalar r1 = CLAMP(z - z0, 0, 1);
    scalar r0 = 1 - r1;

    scalar a0 = s0 * (t0 * (*d)(i0, j0, z0) + t1 * (*d)(i0, j1, z0)) + s1 * (t0 * (*d)(i1, j0, z0) + t1 * (*d)(i1, j1, z0));
    scalar a1 = s0 * (t0 * (*d)(i0, j0, z1) + t1 * (*d)(i0, j1, z1)) + s1 * (t0 * (*d)(i1, j0, z1) + t1 * (*d)(i1, j1, z1));

    return a0 * r0 + a1 * r1;
}

scalar StableFluidsSim::interpolateU(ArrayXs * u, scalar i, scalar j, scalar z)
{
    // STUDENTS: You will certainly need code here, note the i index should be CLAMP-ed to [0, m_N], while j index should be CLAMP-ed to [0, m_N-1], since we have to use (i + 1) and (j + 1)
    int i0 = CLAMP((int)i, 0, m_N);
    int j0 = CLAMP((int)(j - 0.5), 0, m_N - 1);
    int z0 = CLAMP((int)z, 0, m_N);

    int i1 = i0 + 1, j1 = j0 + 1, z1 = z0 + 1;

    scalar s1 = CLAMP(i - i0, 0, 1);
    scalar s0 = 1 - s1;
    scalar t1 = CLAMP(j - 0.5 - j0, 0, 1);
    scalar t0 = 1 - t1;
    scalar r1 = CLAMP(z - z0, 0, 1);
    scalar r0 = 1 - r1;

    scalar a0 = s0 * (t0 * (*u)(i0, j0, z0) + t1 * (*u)(i0, j1, z0)) + s1 * (t0 * (*u)(i1, j0, z0) + t1 * (*u)(i1, j1, z0));
    scalar a1 = s0 * (t0 * (*u)(i0, j0, z1) + t1 * (*u)(i0, j1, z1)) + s1 * (t0 * (*u)(i1, j0, z1) + t1 * (*u)(i1, j1, z1));

    return a0 * r0 + a1 * r1;
}

scalar StableFluidsSim::interpolateV(ArrayXs * v, scalar i, scalar j, scalar z)
{
    // STUDENTS: You will certainly need code here
    int i0 = CLAMP((int)(i - 0.5), 0, m_N - 1);
    int j0 = CLAMP((int)j, 0, m_N);
    int z0 = CLAMP((int)z, 0, m_N);
    int i1 = i0 + 1, j1 = j0 + 1, z1 = z0 + 1;

    scalar s1 = CLAMP(i - 0.5 - i0, 0, 1);
    scalar s0 = 1 - s1;
    scalar t1 = CLAMP(j- j0, 0, 1);
    scalar t0 = 1 - t1;
    scalar r1 = CLAMP(z - z0, 0, 1);
    scalar r0 = 1 - r1;

    scalar a0 = s0 * (t0 * (*v)(i0, j0, z0) + t1 * (*v)(i0, j1, z0)) + s1 * (t0 * (*v)(i1, j0, z0) + t1 * (*v)(i1, j1, z0));
    scalar a1 = s0 * (t0 * (*v)(i0, j0, z1) + t1 * (*v)(i0, j1, z1)) + s1 * (t0 * (*v)(i1, j0, z1) + t1 * (*v)(i1, j1, z1));

    return a0 * r0 + a1 * r1;
}

scalar StableFluidsSim::interpolateW(ArrayXs * w, scalar i, scalar j, scalar z)
{
    // STUDENTS: You will certainly need code here
    int i0 = CLAMP((int)i, 0, m_N);
    int j0 = CLAMP((int)j, 0, m_N);
    int z0 = CLAMP((int)(z-0.5), 0, m_N-1);
    int i1 = i0 + 1, j1 = j0 + 1, z1 = z0 + 1;

    scalar s1 = CLAMP(i -i0, 0, 1);
    scalar s0 = 1 - s1;
    scalar t1 = CLAMP(j- j0, 0, 1);
    scalar t0 = 1 - t1;
    scalar r1 = CLAMP(z - 0.5 - z0, 0, 1);
    scalar r0 = 1 - r1;

    scalar a0 = s0 * (t0 * (*w)(i0, j0, z0) + t1 * (*w)(i0, j1, z0)) + s1 * (t0 * (*w)(i1, j0, z0) + t1 * (*w)(i1, j1, z0));
    scalar a1 = s0 * (t0 * (*w)(i0, j0, z1) + t1 * (*w)(i0, j1, z1)) + s1 * (t0 * (*w)(i1, j0, z1) + t1 * (*w)(i1, j1, z1));

    return a0 * r0 + a1 * r1;
}

void StableFluidsSim::advectU(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, ArrayXs * w, scalar dt)
{
    // assert((*x0 == *x0).all());
    // assert((*u == *u).all());
    // assert((*v == *v).all());
    for (int z = 1; z <= N; z++)
    for (int i = 1; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            // STUDENTS: You will certainly need code here,
            // add the origin of U grid to the coordinate before sampling, for example, sample at (i + 0, j + 0.5) when you need backtracing the old velocity at (i, j)
            // now you have the backward-traced velocity, minus it from the current position (i + 0, j + 0.5), then sample the velocity again.
            scalar ii = i - dt * N * interpolateV(v, i, j + 0.5, z);
            scalar jj = j + 0.5 - dt * N * interpolateU(u, i, j + 0.5, z);
            scalar zz = z - dt * N * interpolateW(w, i, j + 0.5, z);
            (*x)(i, j, z) = interpolateU(x0, ii, jj, zz);
        }
    }
}

void StableFluidsSim::advectV(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, ArrayXs * w, scalar dt)
{
    // assert((*x0 == *x0).all());
    // assert((*u == *u).all());
    // assert((*v == *v).all());
    for (int z = 1; z <= N; z++)
    for (int i = 0; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            // STUDENTS: You will certainly need code here
            scalar ii = i + 0.5 - dt * N * interpolateV(v, i + 0.5, j, z);
            scalar jj = j - dt * N * interpolateU(u, i + 0.5, j, z);
            scalar zz = z - dt * N * interpolateW(w, i + 0.5, j, z);
            (*x)(i, j, z) = interpolateV(x0, ii, jj, zz);
        }
    }
}

void StableFluidsSim::advectW(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, ArrayXs * w, scalar dt)
{
    // assert((*x0 == *x0).all());
    // assert((*u == *u).all());
    // assert((*v == *v).all());
    for (int z = 0; z <= N; z++)
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            // STUDENTS: You will certainly need code here
            scalar ii = i - dt * N * interpolateV(v, i , j, z+0.5);
            scalar jj = j - dt * N * interpolateU(u, i , j, z+0.5);
            scalar zz = z + 0.5 - dt * N * interpolateW(w, i , j, z+0.5);
            (*x)(i, j, z) = interpolateW(x0, ii, jj, zz);
        }
    }
}

void StableFluidsSim::project(int N, ArrayXs * u, ArrayXs * v, ArrayXs * w, ArrayXs * u0, ArrayXs * v0, ArrayXs * w0)
{
    if (VERBOSE) std::cout << "u0: " << std::endl << *u0 << std::endl << std::endl;
    if (VERBOSE) std::cout << "v0: " << std::endl << *v0 << std::endl << std::endl;

    ArrayXs div(N + 2, N + 2, N+2);
    ArrayXs p(N + 2, N + 2, N+2);
    div.setZero();
    p.setZero();
    scalar h = 1.0 / N;

    // STUDENTS: You will certainly need code here

    // set solid boundary conditions, 0 the most top and bottom row / left and right column of u0, v0
    // [(0, 0), (N + 1, 0)], [(0, N ), (N + 1, N)], [(0, 0), (0, N)], and [(N + 1, 0), (N + 1, N)] for u
    // [(0, 0), (0, N + 1)], [(N, 0),(N, N + 1)], [(0, 0), (N, 0)], and [(0, N + 1), (N, N + 1)] for v)
    
      std::cout << "set boundary" << std::endl;
    for(int j = 0; j <= N+1; j++)
    for(int i = 0; i <= N+1; i++) {
        (*u0)(i, 0, j) = 0;
        (*u0)(i, N, j) = 0;
        (*v0)(0, i, j) = 0;
        (*v0)(N, i, j) = 0;
        (*w0)(i, j, 0) = 0;
        (*w0)(i, j, N) = 0;
    }

    for(int j = 0; j <= N+1; j++)
    for(int i = 0; i <= N; i++) {
        (*u0)(0, i, j) = 0;
        (*u0)(N + 1, i, j) = 0;
        (*u0)(j, i, 0) = 0;
        (*u0)(j, i, N + 1) = 0;

        (*v0)(i, 0, j) = 0;
        (*v0)(i, N + 1, j) = 0;
        (*v0)(i, j, 0) = 0;
        (*v0)(i, j, N + 1) = 0;

        (*w0)(j, 0, i) = 0;
        (*w0)(j, N + 1, i) = 0;
        (*w0)(0, j, i) = 0;
        (*w0)(N + 1, j, i) = 0;

    }

    // (*u0)(N + 1, 0) = 0;
    // (*u0)(N + 1, N) = 0;
    // (*v0)(0, N + 1) = 0;
    // (*v0)(N, N + 1) = 0;

  std::cout << "cal div" << std::endl;

    for (int z = 1; z <= N; z++)
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            // compute divergence of the velocity field, note the divergence field is available from ([1, N], [1, N])
            // equation from page 3
            div(i, j, z) = ((*v0)(i, j, z) - (*v0)(i - 1, j, z) + 
                        (*u0)(i, j, z) - (*u0)(i, j - 1, z) +
                        (*w0)(i, j, z) - (*w0)(i, j, z - 1)) / h;
        }
    }

  std::cout << "gauss seidel" << std::endl;
    for (int k = 0; k < 30; k++) {
      for (int z = 1; z <= N; z++)
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) { // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
                // solve for pressure inside the region ([1, N], [1, N])
                // equation from page 3, account for corner cases
                scalar d = div(i, j,z) * h * h;
                scalar num = 0;
                if(i > 1) {
                    num += 1;
                    d -= p(i - 1, j, z);
                }
                if(i < N) {
                    num += 1;
                    d -= p(i + 1, j, z);
                }
                if(j > 1) {
                    num += 1;
                    d -= p(i, j - 1, z);
                }
                if(j < N) {
                    num += 1;
                    d -= p(i, j + 1, z);
                }
                if(z > 1) {
                    num += 1;
                    d -= p(i, j, z-1);
                }
                if(z < N) {
                    num += 1;
                    d -= p(i, j, z+1);
                }
                p(i, j,z) = -d / num;
            }
        }
    }

    (*u) = (*u0);
    (*v) = (*v0);

  std::cout << "apply pressure" << std::endl;
    for(int z = 1; z <= N; z++)
    for(int i = 1; i <= N; i++)
    for(int j = 1; j <= N; j++){
      if(j!=N)
        (*u)(i, j, z) = (*u0)(i, j, z) - (p(i, j+1, z) - p(i, j, z)) / h;
      if(i!=N)
        (*v)(i, j, z) = (*v0)(i, j, z) - (p(i+1, j, z) - p(i, j, z)) / h;
      if(z!=N)
        (*w)(i, j, z) = (*w0)(i, j, z) - (p(i, j, z+1) - p(i, j, z)) / h;
    }


    // for (int i = 1; i <= N; i++) {
    //     for (int j = 1; j < N; j++) {
    //         // apply pressure to correct velocities ([1, N], [1, N)) for u, ([1, N), [1, N]) for v
    //         // equation from page 6
    //         (*u)(i, j) = (*u0)(i, j) - (p(i, j+1) - p(i, j)) / h;
    //         // switch i, j because boundary is ([1, N), [1, N]) for v
    //         (*v)(j, i) = (*v0)(j, i) - (p(j+1, i) - p(j, i)) / h;
    //     }
    // }
}

void StableFluidsSim::dens_step(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, ArrayXs * w, scalar diff, scalar dt)
{
  // IMPORTANT: DO NOT MODIFY THIS CODE!

  ArrayXs * outu = u;
  ArrayXs * outv = v;
  ArrayXs * outw = w;

  add_source(N, x, x0, dt);

  SWAP(x0, x);
  diffuseD(N, x, x0, diff, dt);

  SWAP(x0, x);
  advectD(N, x, x0, u, v, w, dt);

  // if (outu != u)
    *outu = *u;
  // if (outv != v)
    *outv = *v;
  // if (outw != w)
    *outw = *w;
}

void StableFluidsSim::vel_step(int N, ArrayXs * u, ArrayXs * v, ArrayXs * w, ArrayXs * u0, ArrayXs * v0, ArrayXs * w0, scalar visc, scalar dt)
{
  // IMPORTANT: DO NOT MODIFY THIS CODE!

  ArrayXs * outu = u;
  ArrayXs * outv = v;
  ArrayXs * outw = w;

  add_source(N, u, u0, dt);
  add_source(N, v, v0, dt);
  add_source(N, w, w0, dt);

  std::cout << "diffuse vel" << std::endl;

  if(visc > 0.0) {
    SWAP(u0, u);
    SWAP(v0, v);
    SWAP(w0, w);
    diffuseU(N, u, u0, visc, dt);
    diffuseV(N, v, v0, visc, dt);
    diffuseW(N, w, w0, visc, dt);
  }

  std::cout << "project" << std::endl;
  if(m_use_project) {
    SWAP(u0, u);
    SWAP(v0, v);
    SWAP(w0, w);
    project(N, u, v, w, u0, v0, w0);
  }

  m_uAfterDiffusion.setZero();
  m_vAfterDiffusion.setZero();
  m_wAfterDiffusion.setZero();
  m_uAfterDiffusion = *u;
  m_vAfterDiffusion = *v;
  m_wAfterDiffusion = *w;

  std::cout << "advect" << std::endl;
  if(m_use_advect) {
    SWAP(u0, u);
    SWAP(v0, v);
    SWAP(w0, w);
      std::cout << "advect u" << std::endl;
    advectU(N, u, u0, u0, v0, w0, dt);
      std::cout << "advect v" << std::endl;
    advectV(N, v, v0, u0, v0, w0, dt);
      std::cout << "advect w" << std::endl;
    advectW(N, w, w0, u0, v0, w0, dt);
  }
std::cout << "after advect" << std::endl;

  m_uAfterAdvect.setZero();
  m_vAfterAdvect.setZero();
  m_wAfterAdvect.setZero();
  m_uAfterAdvect = *u;
  m_vAfterAdvect = *v;
  m_wAfterAdvect = *w;

  std::cout << "project" << std::endl;
  if(m_use_project) {
    SWAP(u0, u);
    SWAP(v0, v);
    SWAP(w0, w);
    project(N, u, v, w, u0, v0, w0);
  }

  // if (outu != u)
    *outu = *u;
  // if (outv != v)
    *outv = *v;
    *outw = *w;

}

void StableFluidsSim::stepSystem( const scalar& dt)
{
  // IMPORTANT: DO NOT MODIFY THIS CODE!
  if (VERBOSE) std::cout << "step" << std::endl;
  ArrayXs new_d(m_N + 2, m_N + 2, m_N + 2);
  ArrayXs new_u(m_N + 2, m_N + 1, m_N + 2);
  ArrayXs new_v(m_N + 1, m_N + 2, m_N + 2);
  ArrayXs new_w(m_N + 2, m_N + 2, m_N + 1);

  new_d.setZero();
  new_u.setZero();
  new_v.setZero();
  new_w.setZero();


  //sourcing
  std::cout << "sourcing" << std::endl;

  Vector3s c((m_N+2)/2, 10, (m_N+2)/2);
  scalar r = 6;
  for (int z = 1; z <= m_N; z++)
  for (int i = 1; i <= m_N; i++)
  for (int j = 1; j <= m_N; j++){
    Vector3s p(i,j,z);
    if((p-c).norm() < r){
      m_d(i,j,z) += 20 * dt;
      m_u(i,j,z) = 5;
    }
    m_d(i,j,z) = CLAMP(m_d(i,j,z),0,1);
  }

std::cout << "step vel" << std::endl;
  vel_step(m_N, &new_u, &new_v, &new_w, &m_u, &m_v, &m_w, m_visc, dt);
  std::cout << "step dens" << std::endl;
  dens_step(m_N, &new_d, &m_d, &new_u, &new_v, &new_w, m_diff, dt);

  m_d = new_d;
  m_u = new_u;
  m_v = new_v;
  m_w = new_w;

  frame++;
  save(frame);
}

const ArrayXs& StableFluidsSim::getMarkerDensities() const
{
  return m_d;
}

ArrayXs& StableFluidsSim::getMarkerDensities()
{
  return m_d;
}

ArrayXs& StableFluidsSim::getHorizontalVelocities()
{
  return m_u;
}

const ArrayXs& StableFluidsSim::getHorizontalVelocities() const
{
  return m_u;
}

ArrayXs& StableFluidsSim::getVerticalVelocities()
{
  return m_v;
}

const ArrayXs& StableFluidsSim::getVerticalVelocities() const
{
  return m_v;
}

int StableFluidsSim::physicalRows() const
{
  return m_N;
}

int StableFluidsSim::physicalCols() const
{
  return m_N;
}

void StableFluidsSim::clear()
{
  m_d.setZero();
  m_u.setZero();
  m_v.setZero();
  m_w.setZero();
  m_all_ones.setOnes();
}

void StableFluidsSim::setPrescribedVelocity(int p)
{
  // IMPORTANT: DO NOT MODIFY THIS CODE!

  switch (p)
  {
    case 0:
      break;
    case 1:
      for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
        for (int j = m_N * 0.2; j <= m_N * 0.8; j++)
          m_u(i, j) = m_v(i, j) = 0.8;
      break;
    case 2:
      for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
        m_u(i, i) = m_v(i, i) = 0.8;
      break;
    case 3:
      for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
        m_v(i, m_N * 0.2 + 1) = -1.6,
        m_v(i, m_N * 0.8 + 1) = -1.6;
      break;
    case 4:
      for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
        m_u(m_N, i) = 2.4;
      break;
    case 5:
      for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
        for (int j = m_N * 0.83; j <= m_N; j++)
          m_u(j, i) = 2.4;
      break;
    case 8:
      for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
        m_v(i, m_N / 2) = 0.8;
      break;
    case 9:
      for (int i = m_N * 0.2; i <= m_N * 0.8; i++)
        m_v(i, m_N / 2) = -0.8;
      break;
  }
}


void StableFluidsSim::copyState( const StableFluidsSim& otherscene )
{
  // IMPORTANT: DO NOT MODIFY THIS CODE!
  m_diff = otherscene.m_diff;
  m_visc = otherscene.m_visc;
  m_N = otherscene.m_N;
  m_d = otherscene.m_d;
  m_u = otherscene.m_u;
  m_v = otherscene.m_v;
}

ArrayXs& StableFluidsSim::getHorizontalVelocitiesAdvect()
{
  return m_uAfterAdvect;
}

const ArrayXs& StableFluidsSim::getHorizontalVelocitiesAdvect() const
{
  return m_uAfterAdvect;
}

ArrayXs& StableFluidsSim::getVerticalVelocitiesAdvect()
{
  return m_vAfterAdvect;
}

const ArrayXs& StableFluidsSim::getVerticalVelocitiesAdvect() const
{
  return m_vAfterAdvect;
}

ArrayXs& StableFluidsSim::getHorizontalVelocitiesDiffusion()
{
  return m_uAfterDiffusion;
}

const ArrayXs& StableFluidsSim::getHorizontalVelocitiesDiffusion() const
{
  return m_uAfterDiffusion;
}

ArrayXs& StableFluidsSim::getVerticalVelocitiesDiffusion()
{
  return m_vAfterDiffusion;
}

const ArrayXs& StableFluidsSim::getVerticalVelocitiesDiffusion() const
{
  return m_vAfterDiffusion;
}

ArrayXb& StableFluidsSim::getHasFluid()
{
  return m_all_ones;
}

const ArrayXb& StableFluidsSim::getHasFluid() const
{
  return m_all_ones;
}

void StableFluidsSim::setDiffusion(scalar diff)
{
  m_diff = diff;
}

void StableFluidsSim::setViscosity(scalar visc)
{
  m_visc = visc;
}

void StableFluidsSim::setUseAdvect( bool use_advect )
{
  m_use_advect = use_advect;
}

void StableFluidsSim::setUseProj( bool use_proj )
{
  m_use_project = use_proj;
}

scalar StableFluidsSim::getDiffusion()
{
  return m_diff;
}

scalar StableFluidsSim::getViscosity()
{
  return m_visc;
}



void StableFluidsSim::save(int framenum)
{
    UT_VoxelArrayWriteHandleF handle = volume->getVoxelWriteHandle();
    // ArrayXs &marker_density = m_d;
    handle->size(m_N, m_N, m_N);

    for (int z = 0; z < m_N; ++z)
    for (int i = 0; i < m_N; ++i)
        for (int j = 0; j < m_N; ++j) 
            {
              // if(m_d(i+1,j+1,z+1)>0){
              //   std::cout << m_d(i+1,j+1,z+1) << std::endl;
              // }
                handle->setValue(i, j, z, m_d(i+1,j+1,z+1));
            }

    std::ofstream myfile;
    myfile.open(std::to_string(framenum) + ".bgeo");
    gdp.save(myfile, 1, NULL);
    myfile.close();
}
