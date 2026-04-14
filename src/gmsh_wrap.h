/* gmsh_wrap.h
 * Forward declarations for extern "C" gmsh wrappers (gmsh_wrap.cpp).
 * Force-included into poisson2d_mex.c compilation via -include flag.
 */
#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

void gmsh_init_(void);
void gmsh_add_polygon_(int64_t* n, double* pts_x, double* pts_y,
                       double* mesh_size, int64_t* loop_tag);
void gmsh_generate_(int64_t* algorithm, int64_t* n_nodes, int64_t* n_tris);
void gmsh_fetch_(int64_t* n_nodes, int64_t* n_tris,
                 double* nodes, double* tris);
void gmsh_finalize_(void);

#ifdef __cplusplus
}
#endif
