// gmsh_wrap.cpp
// extern "C" wrappers for mwrap-generated poisson2d_mex.c.
// Static cache holds generated mesh between gmsh_generate_ and gmsh_fetch_.

#include "../extern/gmsh/api/gmsh.h"
#include <vector>
#include <cstdint>

static std::vector<int>     s_loops;
static std::vector<double>  s_nodes;    // interleaved x,y,z per node
static std::vector<int64_t> s_tris;     // interleaved v0,v1,v2 per tri (1-based)
static int64_t s_n_nodes = 0;
static int64_t s_n_tris  = 0;

extern "C" {

// gmsh_init_ ---------------------------------------------------------------
void gmsh_init_(void) {
    if (!gmsh::isInitialized()) gmsh::initialize();
    gmsh::model::add("mex_model");
    gmsh::option::setNumber("General.Terminal", 0);
    s_loops.clear();
    s_nodes.clear();
    s_tris.clear();
    s_n_nodes = 0;
    s_n_tris  = 0;
}

// gmsh_add_polygon_ --------------------------------------------------------
// pts_x, pts_y: length-n coordinate arrays (split for mwrap compatibility).
// mesh_size: characteristic length.
// loop_tag: Gmsh curve loop tag returned to caller.
void gmsh_add_polygon_(int64_t* n, double* pts_x, double* pts_y,
                       double* mesh_size, int64_t* loop_tag) {
    int64_t np = *n;
    std::vector<int> p_tags;
    p_tags.reserve(np);
    for (int64_t i = 0; i < np; ++i)
        p_tags.push_back(
            gmsh::model::geo::addPoint(pts_x[i], pts_y[i], 0.0, *mesh_size));
    std::vector<int> l_tags;
    l_tags.reserve(np);
    for (int64_t i = 0; i < np; ++i)
        l_tags.push_back(
            gmsh::model::geo::addLine(p_tags[i], p_tags[(i+1) % np]));
    int lt = gmsh::model::geo::addCurveLoop(l_tags);
    s_loops.push_back(lt);
    *loop_tag = (int64_t)lt;
}

// gmsh_generate_ -----------------------------------------------------------
// Runs 2D mesh, caches result, returns counts so caller can pre-allocate
// before calling gmsh_fetch_.
void gmsh_generate_(int64_t* algorithm, int64_t* n_nodes, int64_t* n_tris) {
    if (s_loops.empty()) { *n_nodes = 0; *n_tris = 0; return; }

    gmsh::model::geo::addPlaneSurface(s_loops);
    gmsh::model::geo::synchronize();
    gmsh::option::setNumber("Mesh.Algorithm", (double)*algorithm);
    gmsh::model::mesh::generate(2);

    // Extract nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, paramCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, paramCoord);
    s_n_nodes = (int64_t)nodeTags.size();
    s_nodes.resize(3 * s_n_nodes);
    for (int64_t i = 0; i < s_n_nodes; ++i) {
        s_nodes[3*i+0] = coord[3*i+0];
        s_nodes[3*i+1] = coord[3*i+1];
        s_nodes[3*i+2] = coord[3*i+2];
    }

    // Build tag -> 1-based index map
    std::vector<int64_t> idx_map(nodeTags.back() + 1, 0);
    for (int64_t i = 0; i < s_n_nodes; ++i)
        idx_map[nodeTags[i]] = i + 1;

    // Extract triangles (element type 2 = 3-node triangle)
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, nodeTagsPerElem;
    gmsh::model::mesh::getElements(elemTypes, elemTags, nodeTagsPerElem, 2);
    s_n_tris = 0;
    s_tris.clear();
    for (std::size_t t = 0; t < elemTypes.size(); ++t) {
        if (elemTypes[t] == 2) {
            const auto& tn = nodeTagsPerElem[t];
            s_n_tris = (int64_t)(tn.size() / 3);
            s_tris.resize(3 * s_n_tris);
            for (int64_t i = 0; i < s_n_tris; ++i) {
                s_tris[3*i+0] = idx_map[tn[3*i+0]];
                s_tris[3*i+1] = idx_map[tn[3*i+1]];
                s_tris[3*i+2] = idx_map[tn[3*i+2]];
            }
            break;
        }
    }

    *n_nodes = s_n_nodes;
    *n_tris  = s_n_tris;
}

// gmsh_fetch_ --------------------------------------------------------------
// Copies cached mesh into caller-allocated MATLAB column-major arrays:
//   nodes: n_nodes x 3  (col 0=x, col 1=y, col 2=z)
//   tris:  n_tris  x 3  (col 0=v0, col 1=v1, col 2=v2, 1-based)
void gmsh_fetch_(int64_t* n_nodes, int64_t* n_tris,
                 double* nodes, double* tris) {
    int64_t nn = *n_nodes;
    int64_t nt = *n_tris;
    for (int64_t i = 0; i < nn; ++i) {
        nodes[i]        = s_nodes[3*i+0];
        nodes[i + nn]   = s_nodes[3*i+1];
        nodes[i + 2*nn] = s_nodes[3*i+2];
    }
    for (int64_t i = 0; i < nt; ++i) {
        tris[i]        = (double)s_tris[3*i+0];
        tris[i + nt]   = (double)s_tris[3*i+1];
        tris[i + 2*nt] = (double)s_tris[3*i+2];
    }
}

// gmsh_finalize_ -----------------------------------------------------------
void gmsh_finalize_(void) {
    gmsh::model::remove();
    gmsh::finalize();
    s_loops.clear();
    s_nodes.clear();
    s_tris.clear();
    s_n_nodes = 0;
    s_n_tris  = 0;
}

} // extern "C"
