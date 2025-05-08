//
// Testing/debugging Dirichlet boundary code
//
// A single domain of a single hex in this one, surrounded by multiple other hexes
//

#include <iostream>
#include <vector>
#include <list>
#include <array>
#include <stdexcept>

#include <sm/vec>
#include <sm/hexgrid>

#include <morph/ReadCurves.h>
#include <morph/tools.h>
#include <morph/ColourMap.h>
#include <morph/ShapeAnalysis.h>
#include <morph/Visual.h>
#include <morph/PolygonVisual.h>

int main()
{
    int rtn = 0;
    try {
        sm::hexgrid hg(0.2, 1, 0);

        hg.setBoundaryOnOuterEdge();

        std::cout << hg.extent() << std::endl;

        std::cout << "Number of hexes in grid:" << hg.num() << std::endl;
        std::cout << "Last vector index:" << hg.lastVectorIndex() << std::endl;

        // Make up a variable.
        std::vector<float> f (hg.num(), 0.1f);

        // Set values in the variable so that it's an identity variable.
        auto hi = hg.hexen.begin();
        auto hi2 = hi;
        while (hi->has_nse()) {
            while (hi2->has_ne()) {
                f[hi2->vi] = 0.2f;
                hi2 = hi2->ne;
            }
            f[hi2->vi] = 0.2f;
            hi2 = hi->nse;
            hi = hi->nse;
        }
        f[hi2->vi] = 0.2f;

        hi = hg.hexen.begin()->nw;
        hi2 = hi;
        while (hi->has_nse()) {
            while (hi2->has_nw()) {
                f[hi2->vi] = 0.4f;
                hi2 = hi2->nw;
            }
            f[hi2->vi] = 0.4f;
            hi2 = hi->nse;
            hi = hi->nse;
        }
        f[hi2->vi] = 0.4f;

        hi = hg.hexen.begin();
        f[hi->vi] = 0.3f;

        f[hi->ne->vi] = 0.55f;

        f[hi->nw->vi] = 0.35f;

        // The code to actually test:
        std::list<morph::DirichVtx<float>> vertices;
        std::list<morph::DirichDom<float>> domains =
        morph::ShapeAnalysis<float>::dirichlet_vertices (&hg, f, vertices);

        // There should be precise number of vertices
        unsigned int reqd = 31;
        if (vertices.size() != reqd) { rtn -= 1; }

        // Expecting one domain
        if (domains.size() != 3) { rtn -= 1; }

        morph::Visual v(1600, 1000, "Dirichlet code");
        v.lightingEffects();
        sm::vec<float, 3> offset = { 0.0f, 0.0f, 0.0f };
        sm::vec<float, 3> offset2 = offset;
        offset2 += {0,0,0.002f};
        std::array<float,3> cl_b = morph::ColourMap<float>::jetcolour (0.78);
        float sz = hg.hexen.front().d;
        for (auto h : hg.hexen) {
            std::array<float,3> cl_a = morph::ColourMap<float>::jetcolour (f[h.vi]);
            std::array<float,3> p = h.position();
            sm::vec<float,3> pv = { p[0], p[1], p[2] };
            sm::vec<float,3> vtx = pv;
            vtx += sm::vec<float, 3>({1,0,0});
            auto pvp = std::make_unique<morph::PolygonVisual<>> (offset, pv, vtx, sz/1.8f, 0.002f, cl_a, 6);
            v.bindmodel (pvp);
            pvp->finalize();
            v.addVisualModel (pvp);
            if (h.boundaryHex()) {
                auto pvp2 = std::make_unique<morph::PolygonVisual<>> (offset2, pv, vtx, sz/12.0f, 0.002f, cl_b, 6);
                v.bindmodel (pvp2);
                pvp2->finalize();
                v.addVisualModel (pvp2);
            }
        }

        std::array<float,3> cl_c = morph::ColourMap<float>::jetcolour (0.98);
        for (auto verti : vertices) {
            sm::vec<float,3> posn = verti.v.plus_one_dim (0.002);
            sm::vec<float,3> vtx = posn + sm::vec<float, 3>({1,0,0});
            auto pvp = std::make_unique<morph::PolygonVisual<>> (offset, posn, vtx, sz/8.0f, 0.002f, cl_c, 60);
            v.bindmodel (pvp);
            pvp->finalize();
            v.addVisualModel (pvp);
        }

        offset += { 0, 0, 0.004 };
        std::array<float,3> cl_d = morph::ColourMap<float>::jetcolour (0.7);
        std::array<float,3> cl_e = morph::ColourMap<float>::jetcolour (0.01);
        for (auto dom_outer : domains) {
            for (auto dom_inner : dom_outer.vertices) {
                // Draw the paths
                for (auto path : dom_inner.pathto_next) {
                    sm::vec<float,3> posn = path.plus_one_dim (0.0);
                    sm::vec<float,3> vtx = posn + sm::vec<float, 3>({1,0,0});
                    auto pvp = std::make_unique<morph::PolygonVisual<>> (offset, posn, vtx, sz/16.0f, 0.002f, cl_d, 6);
                    v.bindmodel (pvp);
                    pvp->finalize();
                    v.addVisualModel (pvp);
                }
                for (auto path : dom_inner.pathto_neighbour) {
                    sm::vec<float,3> posn = path.plus_one_dim (0.0);
                    sm::vec<float,3> vtx = posn + sm::vec<float, 3>({1,0,0});
                    auto pvp = std::make_unique<morph::PolygonVisual<>> (offset, posn, vtx, sz/16.0f, 0.002f, cl_e, 6);
                    v.bindmodel (pvp);
                    pvp->finalize();
                    v.addVisualModel (pvp);
                }
            }
        }

        // Draw small hex at boundary centroid.
        sm::vec<float,3> centroid = hg.boundaryCentroid.plus_one_dim();
        sm::vec<float,3> centroidv = centroid + sm::vec<float,3> ({ 0.0f, 1.0f, 0.0f });
        auto pvp = std::make_unique<morph::PolygonVisual<>> (sm::vec<float>({0,0,0}), centroid, centroidv, sz/16.0f, 0.01f, sm::vec<float>({0,0,1}), 10);
        v.bindmodel (pvp);
        pvp->finalize();
        v.addVisualModel (pvp);
        // red hex at zero
        auto pvp2 = std::make_unique<morph::PolygonVisual<>> (sm::vec<float>( {0,0,0.01f}), sm::vec<float>({0,0,0}), sm::vec<float>({0,1,0}), sz/20.0f, 0.01f, sm::vec<float>({1,0,0}), 8);
        v.bindmodel (pvp2);
        pvp2->finalize();
        v.addVisualModel (pvp2);

        while (v.readyToFinish() == false) {
            v.waitevents (0.018);
            v.render();
        }

    } catch (const std::exception& e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        std::cerr << "Current working directory: " << morph::tools::getPwd() << std::endl;
        rtn = -1;
    }
    return rtn;
}
