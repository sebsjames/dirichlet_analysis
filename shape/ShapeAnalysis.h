/*!
 * \file
 *
 * \author Seb James
 * \date 2019
 */
#pragma once

#include <vector>
#include <list>
#include <set>
#include <map>
#include <limits>
#include <stdexcept>

#include <sm/vec>
#include <sm/hex>
#include <sm/hexgrid>

#include <morph/DirichDom.h>
#include <morph/DirichVtx.h>

namespace morph {

    /*!
     * Rotational direction Rotn::Clock or Rotn::Anticlock
     */
    enum class Rotn {
        Unknown,
        Clock,
        Anticlock
    };

    /*!
     * A helper class, containing pattern analysis code to analyse patterns within hexgrids.
     */
    template <typename Flt>
    class ShapeAnalysis
    {
    public:
        static constexpr bool dbg = false;

        /*!
         * Obtain the contours (as a vector of list<sm::hex>) in the scalar fields f, where threshold is
         * crossed.
         */
        static std::vector<std::list<sm::hex> > get_contours (sm::hexgrid* hg,
                                                          std::vector<std::vector<Flt> >& f,
                                                          Flt threshold) {

            unsigned int nhex = hg->num();
            unsigned int N = f.size();

            std::vector<std::list<sm::hex> > rtn;
            // Initialise
            for (unsigned int li = 0; li < N; ++li) {
                std::list<sm::hex> lh;
                rtn.push_back (lh);
            }

            Flt maxf = -1e7;
            Flt minf = +1e7;
            for (auto h : hg->hexen) {
                if (h.onBoundary() == false) {
                    for (unsigned int i = 0; i<N; ++i) {
                        if (f[i][h.vi] > maxf) { maxf = f[i][h.vi]; }
                        if (f[i][h.vi] < minf) { minf = f[i][h.vi]; }
                    }
                }
            }
            Flt scalef = 1.0 / (maxf-minf);

            // Re-normalize
            std::vector<std::vector<Flt> > norm_f;
            norm_f.resize (N);
            for (unsigned int i=0; i<N; ++i) {
                norm_f[i].resize (nhex, 0.0);
            }

            for (unsigned int i = 0; i<N; ++i) {
                for (unsigned int h=0; h<nhex; h++) {
                    norm_f[i][h] = (f[i][h] - minf) * scalef;
                }
            }

            // Collate
            for (unsigned int i = 0; i<N; ++i) {
                for (auto h : hg->hexen) {
                    if (h.onBoundary() == false) {
                        if (norm_f[i][h.vi] >= threshold) {
                            if ( (h.has_ne() && norm_f[i][h.ne->vi] < threshold)
                                 || (h.has_nne() && norm_f[i][h.nne->vi] < threshold)
                                 || (h.has_nnw() && norm_f[i][h.nnw->vi] < threshold)
                                 || (h.has_nw() && norm_f[i][h.nw->vi] < threshold)
                                 || (h.has_nsw() && norm_f[i][h.nsw->vi] < threshold)
                                 || (h.has_nse() && norm_f[i][h.nse->vi] < threshold) ) {
                                rtn[i].push_back (h);
                            }
                        }
                    } else { // h.onBoundary() is true
                        if (norm_f[i][h.vi] >= threshold) {
                            rtn[i].push_back (h);
                        }
                    }
                }
            }

            return rtn;
        }

        /*!
         * Like get_contours, but returns a full hexgrid's worth of Flts instead of
         * lists of hexes.
         */
        static std::vector<Flt> get_contour_map (sm::hexgrid* hg,
                                                 std::vector<std::vector<Flt> >& f,
                                                 Flt threshold) {

            unsigned int nhex = hg->num();
            unsigned int N = f.size();

            std::vector<Flt> rtn (nhex, 0.0);

            Flt maxf = -1e7;
            Flt minf = +1e7;
            for (auto h : hg->hexen) {
                if (h.onBoundary() == false) {
                    for (unsigned int i = 0; i<N; ++i) {
                        if (f[i][h.vi] > maxf) { maxf = f[i][h.vi]; }
                        if (f[i][h.vi] < minf) { minf = f[i][h.vi]; }
                    }
                }
            }
            Flt scalef = 1.0 / (maxf-minf);

            // Re-normalize
            std::vector<std::vector<Flt> > norm_f;
            norm_f.resize (N);
            for (unsigned int i=0; i<N; ++i) {
                norm_f[i].resize (nhex, 0.0);
            }

            for (unsigned int i = 0; i<N; ++i) {
                for (unsigned int h=0; h<nhex; h++) {
                    norm_f[i][h] = (f[i][h] - minf) * scalef;
                }
            }

            // Collate
            for (unsigned int i = 0; i<N; ++i) {
                for (auto h : hg->hexen) {
                    if (h.onBoundary() == false) {
                        if (norm_f[i][h.vi] >= threshold) {
                            if ( (h.has_ne() && norm_f[i][h.ne->vi] < threshold)
                                 || (h.has_nne() && norm_f[i][h.nne->vi] < threshold)
                                 || (h.has_nnw() && norm_f[i][h.nnw->vi] < threshold)
                                 || (h.has_nw() && norm_f[i][h.nw->vi] < threshold)
                                 || (h.has_nsw() && norm_f[i][h.nsw->vi] < threshold)
                                 || (h.has_nse() && norm_f[i][h.nse->vi] < threshold) ) {
                                rtn[h.vi] = (Flt)i/(Flt)N;
                            }
                        }
                    } else { // h.onBoundary() is true
                        if (norm_f[i][h.vi] >= threshold) {
                            rtn[h.vi] = (Flt)i/(Flt)N;
                        }
                    }
                }
            }

            return rtn;
        }

        //! Like get_contour_map, but no pre-normalizing and sets contours to the flag value
        //! (used by SPW in SOM model analysis steps)
        static std::vector<Flt> get_contour_map_flag_nonorm (sm::hexgrid* hg,std::vector<Flt> & f, Flt threshold, Flt flagVal) {
            unsigned int nhex = hg->num();
            std::vector<Flt> rtn (nhex, 0.0);
            for (auto h : hg->hexen) {
                if (h.onBoundary() == false) {
                    if (f[h.vi] >= threshold) {
                        if ((h.has_ne() && f[h.ne->vi] < threshold)
                            || (h.has_nne() && f[h.nne->vi] < threshold)
                            || (h.has_nnw() && f[h.nnw->vi] < threshold)
                            || (h.has_nw() && f[h.nw->vi] < threshold)
                            || (h.has_nsw() && f[h.nsw->vi] < threshold)
                            || (h.has_nse() && f[h.nse->vi] < threshold) ) {
                            rtn[h.vi] = flagVal;
                        }
                    }
                }
            }
            return rtn;
        }

        //! Like get_contour_map, but for N vector<Flt>s in @f, return N+1 positive values in the
        //! return vector. Good for plotting contours with ColourMapType::RainbowZeroBlack or
        //! ColourMapType::RainbowZeroWhite
        static std::vector<Flt> get_contour_map_nozero (sm::hexgrid* hg,
                                                        std::vector<std::vector<Flt> >& f,
                                                        Flt threshold) {

            unsigned int nhex = hg->num();
            unsigned int N = f.size();

            std::vector<Flt> rtn (nhex, 0.0);

            Flt maxf = -1e7;
            Flt minf = +1e7;
            for (auto h : hg->hexen) {
                if (h.onBoundary() == false) {
                    for (unsigned int i = 0; i<N; ++i) {
                        if (f[i][h.vi] > maxf) { maxf = f[i][h.vi]; }
                        if (f[i][h.vi] < minf) { minf = f[i][h.vi]; }
                    }
                }
            }
            Flt scalef = 1.0 / (maxf-minf);

            // Re-normalize
            std::vector<std::vector<Flt> > norm_f;
            norm_f.resize (N);
            for (unsigned int i=0; i<N; ++i) {
                norm_f[i].resize (nhex, 0.0);
            }

            for (unsigned int i = 0; i<N; ++i) {
                for (unsigned int h=0; h<nhex; h++) {
                    norm_f[i][h] = (f[i][h] - minf) * scalef;
                }
            }

            // Collate
            for (unsigned int i = 0; i<N; ++i) {
                for (auto h : hg->hexen) {
                    if (h.onBoundary() == false) {
                        if (norm_f[i][h.vi] >= threshold) {
                            if ( (h.has_ne() && norm_f[i][h.ne->vi] < threshold)
                                 || (h.has_nne() && norm_f[i][h.nne->vi] < threshold)
                                 || (h.has_nnw() && norm_f[i][h.nnw->vi] < threshold)
                                 || (h.has_nw() && norm_f[i][h.nw->vi] < threshold)
                                 || (h.has_nsw() && norm_f[i][h.nsw->vi] < threshold)
                                 || (h.has_nse() && norm_f[i][h.nse->vi] < threshold) ) {
                                rtn[h.vi] = (Flt)(i+1)/(Flt)(N+1); // only this line...
                            }
                        }
                    } else { // h.onBoundary() is true
                        if (norm_f[i][h.vi] >= threshold) {
                            rtn[h.vi] = (Flt)(i+1)/(Flt)(N+1); // ...and this line differ from get_contour_map()
                        }
                    }
                }
            }

            return rtn;
        }

        /*!
         * Take a set of variables, @f, for the given hexgrid @hg. Return a vector of Flts (again,
         * based on the hexgrid @hg) which marks each hex with the outer index of the @f which has
         * highest value in that hex, scaled and converted to a float.
         */
        static std::vector<Flt>
        dirichlet_regions (sm::hexgrid* hg, std::vector<std::vector<Flt> >& f) {
            unsigned int N = f.size();

            // Single variable to return
            std::vector<Flt> rtn (f[0].size(), 0.0);

            // Mark regions first.
            for (auto h : hg->hexen) {

                Flt maxf = -1e7;
                for (unsigned int i = 0; i<N; ++i) {
                    if (f[i][h.vi] > maxf) {
                        maxf = f[i][h.vi];
                        Flt fi = 0.0f;
                        fi = (Flt)i;
                        rtn[h.vi] = (fi/N);
                    }
                }
            }

            return rtn;
        }

        /*!
         * @regions is a vector of a size specified in the hexgrid @hg containing N
         * unique IDs. For each unique ID, compute the centroid of all the hexes
         * having that ID. Return a map keyed by ID, containing the coordinates of
         * each centroid.
         */
        static std::map<Flt, sm::vec<Flt, 2> >
        region_centroids (sm::hexgrid* hg, const std::vector<Flt>& regions) {
            std::map<Flt, sm::vec<Flt, 2> > centroids;
            std::map<Flt, Flt> counts;
            for (unsigned int h = 0; h<regions.size(); ++h) {
                centroids[regions[h]][0] += hg->d_x[h];
                centroids[regions[h]][1] += hg->d_y[h];
                counts[regions[h]] += 1.0;
            }
            auto ci = centroids.begin();
            while (ci != centroids.end()) {
                ci->second /= counts[ci->first];
                ++ci;
            }
            return centroids;
        }

        /*!
         * A method to test the hex give by @h, which must live on the hexgrid pointed to by @hg, to
         * see if it is a Dirichlet vertex. If so, a vertex should be created in @vertices.
         */
        static void
        vertex_test (sm::hexgrid* hg, std::vector<Flt>& f,
                     std::list<sm::hex>::iterator h, std::list<DirichVtx<Flt> >& vertices) {

            // For each hex, examine its neighbours, counting number of different neighbours.
            std::set<Flt> n_ids;
            n_ids.insert (f[h->vi]);
            for (unsigned int ni = 0; ni < 6; ++ni) {
                if (h->has_neighbour(ni)) {
                    n_ids.insert (f[h->get_neighbour(ni)->vi]);
                }
            }
            unsigned int n_ids_sz = n_ids.size();

            // n_ids is WRONG. Need to test for n_ids that are adjoining each vertex. This number
            // will be 1, 2 or 3.
            if constexpr (dbg) { std::cout << h->outputRG() << ": boundaryHex:"
                                           << h->boundaryHex() << ", n_ids.size():" << n_ids_sz << std::endl; }

            if (n_ids_sz >= 2) {

                // Then there's the possibility of a vertex on this hex.
                if constexpr (dbg) { std::cout << h->outputRG() << ": Possibility of boundary or internal vertex...\n"; }

                if (h->boundaryHex() == true) { // 1. Test for boundary vertices. n_ids size >= 2.

                    if constexpr (dbg) { std::cout << h->outputRG() << ": Possibility of boundary vertex...\n"; }

                    // Here, I need to set a vertex where two hexes join and we're on the
                    // boundary. This provides information to set the angles to discover the best
                    // center for each domain (see Honda 1983).
                    for (int ni = 0; ni < 6; ++ni) { // ni==0 is neighbour E. 1 is neighbour NE, etc.

                        // If there's a neighbour in direction ni and it has a different ID:
                        if (h->has_neighbour(ni) && f[h->get_neighbour(ni)->vi] != f[h->vi]) {

                            // Examine which direction DOESN'T have a neighbour and that will
                            // determine which hex vertex is the domain vertex.

                            // The first non-identical ID
                            int nii = (ni+1)%6;
                            if (!h->has_neighbour(nii)) {
                                // Then vertex direction is "vertex direction ni"
                                if constexpr (dbg) { std::cout << "vtx bh 1\n"; }
                                DirichVtx<Flt> a(
                                    h->get_vertex_coord(ni),
                                    hg->getd(),
                                    f[h->vi],
                                    sm::vec<Flt, 2>({Flt{-1}, f[h->get_neighbour(ni)->vi]}),
                                    h);
                                a.onBoundary = true;
                                vertices.push_back (a);

                            } else {
                                nii = ni>0?ni-1:5;
                                if (!h->has_neighbour(nii)) {
                                    // Then vertex dirn is "vertex direction (ni-1) or 5", i.e. nii.
                                    if constexpr (dbg) { std::cout << "vtx bh 2\n"; }
                                    DirichVtx<Flt> a (
                                        h->get_vertex_coord(nii),
                                        hg->getd(),
                                        f[h->vi],
                                        sm::vec<Flt, 2>({f[h->get_neighbour(ni)->vi], Flt{-1}}),
                                        h);
                                    a.onBoundary = true;
                                    vertices.push_back (a);
                                }
                            }
                        }
                    }
                } // end boundaryHex specific extra test.

                // 2. Test for internal vertices with >2 (i.e. 3) different types in self &
                // neighbouring hexes, so now work out which of the Hex's vertices is the vertex of
                // the domain. Note that an internal vertex CAN appear on a boundaryHex as long as
                // n_ids is >= 3.
                if (n_ids_sz >= 3) {
                    if constexpr (dbg) { std::cout << h->outputRG() << ": Possibility of internal vertex...\n"; }
                    for (int ni = 0; ni < 6; ++ni) { // ni==0 is neighbour east. 1 is neighbour NE, etc.

                        // If there's a neighbour in direction ni and that neighbour has different ID:
                        if (h->has_neighbour(ni) && f[h->get_neighbour(ni)->vi] != f[h->vi]) {

                            // The first non-identical ID
                            Flt f1 = f[h->get_neighbour(ni)->vi];
                            int nii = (ni+1)%6;

                            // Test all six vertices, no breaking.
                            if (h->has_neighbour(nii)
                                && f[h->get_neighbour(nii)->vi] != f[h->vi]
                                && f[h->get_neighbour(nii)->vi] != f1 // f1 already tested != f[h->vi]
                                ) {
                                // Then vertex is "vertex ni"
                                if constexpr (dbg) { std::cout << "vtx ih 1\n"; }
                                vertices.push_back (
                                    DirichVtx<Flt>(
                                        h->get_vertex_coord(ni),
                                        hg->getd(),
                                        f[h->vi],
                                        sm::vec<Flt, 2>({f[h->get_neighbour(nii)->vi], f[h->get_neighbour(ni)->vi]}),
                                        h)
                                    );
                            }
                        }
                    }
                }
            }
        }

        /*!
         * Walk an edge between two domains. Common code used by walk_to_neighbour and walk_to_next.
         *
         * @f The map of identities for the hexgrid @hg
         *
         * @v The starting Dirichlet vertex at which the edge starts
         *
         * @path A return variable into which the coordinates of the path will be written
         *
         * @edgedoms The identities of the two domains that the edge should trace a path between.
         *
         * @next_neighb_dom A variable in which to return the next neighbour domain ID that is
         * found.
         */
        static sm::vec<Flt, 2>
        walk_common (std::vector<Flt>& f,
                     DirichVtx<Flt>& v,
                     std::list<sm::vec<Flt, 2>>& path,
                     sm::vec<Flt, 2>& edgedoms,
                     Flt& next_neighb_dom) {

            //WALK ("** Called to walk from " << v.v << " or " << v.hi->outputRG() << ". edgedoms: "<< edgedoms);

            // Really, we only have coordinates to return.
            sm::vec<Flt, 2> next_one = { std::numeric_limits<Flt>::max(), std::numeric_limits<Flt>::max() };

            // Used later
            int i = 0;
            int j = 0;

            // Walk the edge, with hexit pointing to the hexes on the edgedoms[0]
            // side. _Initially_, point hexit at the hex that's on the inside of the domain for
            // which v is a Dirichlet vertex - v.hi. At least, this is what you do when walking OUT
            // to a neighbour vertex that's part of another domain.
            std::list<sm::hex>::iterator hexit = v.hi;
            // point hexit_neighb to the hexes on the edgedoms[1] side
            std::list<sm::hex>::iterator hexit_neighb = v.hi;
            // The first hex, inside the domain.
            std::list<sm::hex>::iterator hexit_first = v.hi;
            // Temporary hex pointers
            std::list<sm::hex>::iterator hexit_next = v.hi;
            std::list<sm::hex>::iterator hexit_last = v.hi;

            // Set true when we find the partner vertex.
            bool partner_found = false;

            sm::vec<Flt, 2> v_init = v.v;

            /*
             * A Find hexit_first
             */

            // Test if the initial hex itself is on the side of the edge (as it will be when walking
            // from one domain vertex to the next domain vertex). This code is side-stepped when
            // walking OUT to a neighbour vertex that's part of another domain. This looks a bit
            // hacky, but I do have to find out which hex to start from when walking along the edge
            // and so some sort of code like this has to go somewhere.
            if (f[hexit_first->vi] == edgedoms[0]) {
                //WLK2 ("Hex AT " << hexit_first->outputRG() << " has f=" << f[hexit_first->vi]);
                // Then I need to find out which of the other hexes should swap to hexit_first.
                for (i = 0; i<6; ++i) {
                    if (hexit_first->has_neighbour(i)) {
                        // For each neighbour to hexit_first, check its centre is one long-radius
                        // from v.v. Only two will fulfil this criterion.
                        Flt x_ = hexit_first->get_neighbour(i)->x - v_init[0];
                        Flt y_ = hexit_first->get_neighbour(i)->y - v_init[1];
                        Flt distance = sqrt (x_*x_ + y_*y_);
                        if constexpr (dbg) { std::cout << "vertex to hex-centre distance: " << distance << std::endl; }
                        if constexpr (dbg) { std::cout << "hexgrid long radius distance: " << hexit_first->getLR() << std::endl; }
                        bool correct_distance = distance-hexit_first->getLR() < 0.001 ? true : false;
                        if (correct_distance) {
                            if (f[hexit_first->get_neighbour(i)->vi] != edgedoms[1]
                                && f[hexit_first->get_neighbour(i)->vi] != edgedoms[0]) {
                                //WLK2 ("Found the true hexit_first");
                                hexit = hexit_first;
                                hexit_first = hexit_first->get_neighbour(i);
                                break;
                            }
                        }
                    }
                }
            }

            // Now the main walking algorithm

            // Include loopcount to catch any infinite loops to debug their reason for existence.
            unsigned int loopcount = 0;

            while (!partner_found) {

                if (loopcount++ > 1000) {
                    //WALK ("UH OH: Breaking out of while loop after " << loopcount << " iterations");
                    next_one = v_init;
                    partner_found = true;
                    break;
                }

                //WLK2 ("===== while loop. partner_found==false ======");

                /*
                 * B Find the iterator hexit.
                 */

                // Find the initial direction of the edge and the hex containing edgedoms[0]:
                int hexit_first_dirn = std::numeric_limits<int>::max();
                for (i = 0; i<6; ++i) {
                    if constexpr (dbg) { std::cout << "i=" << i << ", Comparing coordinates: "
                                                   << hexit_first->get_vertex_coord(i)
                                                   << " and v_init = " << v_init << std::endl; }

                    if (hexit_first->compare_vertex_coord(i, v_init) == true) {

                        // Then the neighbours are either side of vertex direction i.
                        //WLK2 ("initial *vertex* in direction " << i << "/"
                        //      << Hex::vertex_name(i) << " wrt current hexit_first: "
                        //      << hexit_first->outputRG());

                        if (hexit_first->has_neighbour ((i+1)%6)) {
                            //WLK2 ("Hex adjoining " << hexit_first->outputRG()
                            //      << " to " << Hex::neighbour_pos((i+1)%6)
                            //      << " has f=" << f[hexit_first->get_neighbour ((i+1)%6)->vi] << ", "
                            //      << hexit_first->get_neighbour ((i+1)%6)->outputRG());

                            // The next hex to be pointed to by hexit is the one with
                            // f==edgedoms[0]
                            hexit = hexit_first->get_neighbour ((i+1)%6);

                            if (f[hexit->vi] == edgedoms[0]) {
                                hexit_first_dirn = (i+1)%6;
                                //WLK2 ("Good, hex in direction "
                                //      << Hex::neighbour_pos(hexit_first_dirn) << ", which is "
                                //      << hexit->outputRG()
                                //      << ", has ID = edgedoms[0] = " << edgedoms[0]);
                                break;
                            }
                            //else {
                            //    WLK2 ("Hex in direction " << Hex::neighbour_pos((i+1)%6)
                            //          << " has ID!=edgedoms[0] = " << edgedoms[0]);
                            //}
                        } // else No neighbour in direction " << ((i+1)%6) << " of a vertex hex.

                        if (hexit_first->has_neighbour ((i>0)?(i-1):5)) {
                            //WLK2 ("Hex adjoining " << hexit_first->outputRG()
                            //      << " to " << Hex::neighbour_pos((i>0)?(i-1):5)
                            //      << " has f=" << f[hexit_first->get_neighbour ((i>0)?(i-1):5)->vi]
                            //      << ", " << hexit_first->get_neighbour ((i>0)?(i-1):5)->outputRG());

                            // The next hex to be pointed to by hexit is the one with
                            // f==edgedoms[0]
                            hexit = hexit_first->get_neighbour (i>0?i-1:5);

                            if (f[hexit->vi] == edgedoms[0]) {
                                hexit_first_dirn = i>0?i-1:5;
                                //WLK2 ("Good, hex in direction "
                                //      << Hex::neighbour_pos(hexit_first_dirn) << ", which is "
                                //      << hexit->outputRG()
                                //      << ", has ID = edgedoms[0] = " << edgedoms[0]
                                //      << " (edgedoms[1]=" << edgedoms[1] << ")");
                                break;
                            } else {
                                //WLK2 ("Hex in direction " << Hex::neighbour_pos((i>0)?(i-1):5)
                                //      << " has ID!=edgedoms[0] = " << edgedoms[0]);
                            }
                        } else {
                            // If we get here, then neither hex to each side of the initial hexes
                            // were on the edge. That means that the edge has length 2 vertices
                            // only.
                            hexit_first_dirn = (i>0)?(i-1):5;
                            //WLK2 ("Okay, hex in direction "
                            //      << Hex::neighbour_pos(hexit_first_dirn) << ", which is "
                            //      << hexit->outputRG()
                            //      << ", has neither ID. Only 2 vertices in the edge.");
                            break;
                        }
                    }
                }
                //WLK2 ("After determining initial direction of edge, i=" << i
                //      << " or vertex dirn: " << Hex::vertex_name(i));
                //WLK2 ("...and direction to edgedoms[0] Hex is " << hexit_first_dirn
                //      << " or hex dirn: " << Hex::neighbour_pos(hexit_first_dirn));

                /*
                 * C Find the iterator hexit_neighb.
                 */

                // Now point hexit_neighb at the hex_first containing edgedoms[1]_first Look at
                // hex neighbours in directions i+1 and i-1 from hexit_first.
                bool found_second = false;
                // dirn from hexit_first to hexit_neighb which has edgedoms[1] identity
                int hexit_second_dirn = std::numeric_limits<int>::max();

                // Use *hex identity* to find neighbour, along with hex location with respect to
                // the initial vertex, v.
                j = (hexit_first_dirn+1)%6;
                if (hexit_first->has_neighbour (j)) {
                    // If we have a neighbour, then check if it's on the other side of the edge;
                    // i.e. that the initial vertex v.v lies between the neighbour
                    // hexit->get_neighbour(j) and hexit.
                    hexit_neighb = hexit_first->get_neighbour(j);
                    //WLK2 ("hexit_neighb in j+1 dirn, which should be over the edge has f="
                    //      << f[hexit_neighb->vi] << ", " << hexit_neighb->outputRG()
                    //      << ". Comparing with edgedoms[1]=" << edgedoms[1]);
                    // Ok, there IS a neighbour in this direction.

                    // Distance to vertex
                    float d_to_v = hexit_neighb->distanceFrom (v_init);
                    // Is distance to vertex nearly equal to one long hex radius?
                    bool oneLR = abs(d_to_v - hexit_neighb->getLR()) < (hexit_neighb->d/100);

                    if (f[hexit_neighb->vi] == edgedoms[1] && oneLR) {
                        //WLK2 ("1 hexit_neighb found in 'j+1' dirn based on identity.");
                        found_second = true;
                        hexit_second_dirn = j;
                    }
                } // else no neighbour in 'j+1' dirn
                if (!found_second) {
                    j = (hexit_first_dirn>0)?(hexit_first_dirn-1):5;
                    if (hexit_first->has_neighbour (j)) {
                        hexit_neighb = hexit_first->get_neighbour(j);
                        //WLK2 ("hexit_neighb in j-1 dirn, which should be over the edge has f="
                        //      << f[hexit_neighb->vi] << ", " << hexit_neighb->outputRG());
                        // Is distance to vertex nearly equal to one long hex radius?
                        float d_to_v = hexit_neighb->distanceFrom (v_init);
                        bool oneLR = abs(d_to_v - hexit_neighb->getLR()) < (hexit_neighb->d/100);
                        if (f[hexit_neighb->vi] == edgedoms[1] && oneLR) {
                            //WLK2 ("2 hexit_neighb found in 'j-1' dirn based on identity.");
                            found_second = true;
                            hexit_second_dirn = j;
                        }
                    } // else no neighbour in 'j-1' dirn
                }

                if (!found_second) {
                    // Failed to find the second hex associated with the initial vertex!
                    return v_init;
                }

                //WALK ("hexit_first: "<< hexit_first->outputRG()
                //      << ", hexit: " << hexit->outputRG()
                //      << ", hexit_neighb: " << hexit_neighb->outputRG());

                /*
                 * D Determine rotations around the iterator hexit to make.
                 */

                // Can now say whether the edgedoms are in clockwise or anti-clockwise order.
                Rotn rot = Rotn::Unknown;
                int hex_hex_neighb_dirn = std::numeric_limits<int>::max();
                if (hexit_second_dirn == ((hexit_first_dirn>0)?(hexit_first_dirn-1):5)) {
                    // Rotation of edgedoms[0] to edgedoms[1] is clockwise around hexit_first.

                    // Direction from hexit to hexit_neighb is the hexit anti-direction + 1
                    hex_hex_neighb_dirn = (((hexit_first_dirn+3)%6)+1)%6;

                    // Rotation from hexit_first to hexit_neighb is therefore ANTI-clockwise.
                    //WLK2 ("Rotate anticlockwise around hexit");
                    rot = Rotn::Anticlock;

                } else if (hexit_second_dirn  == (hexit_first_dirn+1)%6) {
                    // Rotation of edgedoms[0] to edgedoms[1] is anti-clockwise around
                    // hexit_first

                    // Direction from hexit to hexit_neighb is the hexit anti-direction - 1
                    hex_hex_neighb_dirn = (hexit_first_dirn+3)%6;
                    hex_hex_neighb_dirn = hex_hex_neighb_dirn>0?(hex_hex_neighb_dirn-1):5;

                    // Rotation from hexit_first to hexit_neighb is therefore CLOCKWISE.
                    //WLK2 ("Rotate clockwise around hexit");
                    rot = Rotn::Clock;

                } // else rot == Rotn::Unknown

                // Now hexit and hexit_neighb hexes straddle the edge that I want to walk along, and
                // I know which way to rotate around hexit to find all the edge vertices that
                // surround hexit.

                // It should be that case that hexit_neighb ==
                // hexit->get_neighbour(hex_hex_neighb_dirn);
                if (hexit_neighb != hexit->get_neighbour(hex_hex_neighb_dirn)) {
                    // Return?
                    if constexpr (dbg) { std::cout << "hexit_neighb is in the WRONG direction, return v_init now\n"; }
                    next_one = v_init;
                    return next_one;
                }

                // Here we have, in hexit, a hex with value edgedoms[0]. Find the neighbour hex
                // with value edgedoms[1] and add two vertices to v.edge accordingly.

                // Rotate all the way around each "inner edge hex", starting from
                // hex_hex_neighb_dirn.  Rotation may be clockwise or anticlockwise.
                int last_j;
                if (rot == Rotn::Anticlock) {
                    last_j = hex_hex_neighb_dirn>0?(hex_hex_neighb_dirn-1):5;
                } else {
                    last_j = (hex_hex_neighb_dirn+1)%6;
                }

                // This for loop rotates around each inner edge hexit:
                hexit_last = hexit_neighb;
                for (j  = hex_hex_neighb_dirn;
                     j != last_j;
                     j  = (rot == Rotn::Anticlock)?((j+1)%6):(j>0?j-1:5) ) {

                    //WLK2 ("Inner hex loop, j=" << j << ", last_j=" << last_j);

                    if (!hexit->has_neighbour (j)) {
                        //WLK2 ("No neighbour in direction " << Hex::neighbour_pos(j));
                        // so edge ends
                        v_init = hexit->get_vertex_coord (j>0?j-1:5);
                        //WALK ("Edge ends; push_back very last vertex coordinate "
                        //      << (j>0?j-1:5) << " for the path: ("
                        //      << v_init[0] << "," << v_init[1] << ")");
                        path.push_back (v_init);
                        next_neighb_dom = -1.0;
                        partner_found = true;
                        next_one = v_init;
                        break;
                    } else {
                        // If we have a neighbour, then check if it's on the other side of the edge.
                        hexit_next = hexit->get_neighbour (j);

                        if constexpr (dbg) { std::cout << "hexit_next, " << hexit_next->outputRG()
                                                       << ", which should be over the edge has f=" << f[hexit_next->vi] << std::endl; }
                        if (f[hexit_next->vi] == edgedoms[1]) {
                            // hexit_next has identity edgedoms[1], so add vertex j
                            v_init = hexit->get_vertex_coord (j>0?j-1:5);
                            //WALK ("push_back vertex coordinate " << (j>0?j-1:5) << " for the path: ("
                            //      << v_init[0] << "," << v_init[1] << ")");
                            path.push_back (v_init);
                            // Update hexit_last
                            hexit_last = hexit_next;

                        } else {
                            //WLK2 ("This neighbour does not have identity = edgedoms[1] = "
                            //             << edgedoms[1]);
                            if (f[hexit_next->vi] == edgedoms[0]) {
                                //WLK2 ("This neighbour DOES have identity = edgedoms[0] = "
                                //      << edgedoms[0]);

                                v_init = hexit->get_vertex_coord (j>0?j-1:5);

                                // This is the time to cycle around the hexes
                                //WALK ("Step to next hex");
                                //WLK2 ("Set hexit_first to " << hexit->outputRG());
                                hexit_first = hexit;
                                //WLK2 ("Set hexit to " << hexit_next->outputRG());
                                hexit = hexit_next;
                                //WLK2 ("Set hexit_neighb to hexit_last: " << hexit_last->outputRG());
                                // hexit_last is the last neighbour with identity edgedoms[1]
                                hexit_neighb = hexit_last;
                                //WLK2 ("break out of loop as hexits updated.");
                                break;
                            } else {
                                //WLK2 ("Neither of the edgedom identities, must be end of the edge");
                                v_init = hexit->get_vertex_coord (j>0?j-1:5);
                                //WALK ("Edge ends; push_back final vertex coordinate " << (j>0?j-1:5)
                                //      << " for the path: ("
                                //      << v_init[0] << "," << v_init[1] << ")");
                                path.push_back (v_init);
                                next_one = v_init;
                                next_neighb_dom = f[hexit_next->vi];
                                //WLK2 ("Set next_neighb_dom to " << next_neighb_dom
                                //      << " with edgedoms "
                                //      << edgedoms[0] << "," << edgedoms[1]);
                                partner_found = true;
                                break;
                            }
                        }
                    }
                }
            //WLK2 ("Finished loop around inner hex");
            } // end while !partner_found

            //WALK ("*** returning. loopcount is " << loopcount);
            return next_one;
        }

        /*!
         * Walk out to the next vertex from vertx @v on hexgrid @hg for which identities are in @f.
         *
         * @hg The hexgrid on which the action takes place
         *
         * @f The identity variable.
         *
         * @v The Dirichlet Vertex that we're walking from, and into which we're going to write the
         * path to the next neighbour
         *
         * @next_neighb_dom The identity of the next domain neighbour when we find the vertex
         * neighbour.
         */
        static sm::vec<Flt, 2>
        walk_to_next (std::vector<Flt>& f, DirichVtx<Flt>& v, Flt& next_neighb_dom)
        {
            // Starting from hex v.hi, find neighbours whos f values are v.f/v.neighb[0]. Record
            // (in v.path_to_next) a series of coordinates that make up the path between that vertex
            // and the next vertex in the domain.
            sm::vec<Flt, 2> edgedoms = { v.f, v.neighb[0] };
            return walk_common (f, v, v.pathto_next, edgedoms, next_neighb_dom);
        }

        /*!
         * Walk out to a neighbour from vertex @v.
         *
         * @f The identity variable.
         *
         * @v The Dirichlet Vertex that we're walking from, and into which we're going to write the
         * path to the next neighbour
         *
         * @next_neighb_dom The identity of the next domain neighbour when we find the vertex
         * neighbour.
         */
        static sm::vec<Flt, 2>
        walk_to_neighbour (std::vector<Flt>& f, DirichVtx<Flt>& v, Flt& next_neighb_dom) {

            // Don't set neighbours for the edge vertices (though edge vertices *can be set* as
            // neighbours for other vertices).
            if (v.neighb[0] == -1.0f || v.neighb[1] == -1.0f) {
                return sm::vec<Flt, 2>({Flt{0}, Flt{0}});
            }
            sm::vec<Flt, 2> edgedoms = v.neighb;
            return walk_common (f, v, v.pathto_neighbour, edgedoms, next_neighb_dom);
        }

        /*!
         * Given an iterator into the list of DirichVtxs @vertices, find the next vertex in the
         * domain, along with the vertex neighbours, and recurse until @domain has been populated
         * with all the vertices that define it.
         *
         * Return true for success, false for failure, and leave dv pointing to the next vertex in
         * vertices so that @domain can be stored, reset and the next Dirichlet domain can be found.
         */
        static bool
        process_domain (sm::hexgrid* hg, std::vector<Flt>& f,
                        typename std::list<DirichVtx<Flt>>::iterator dv,
                        std::list<DirichVtx<Flt>>& vertices,
                        DirichDom<Flt>& domain,
                        DirichVtx<Flt> first_vtx) {

            // Domain ID is set in dv as dv->f;
            DirichVtx<Flt> v = *dv;

            // On the first call, first_vtx should have been set to vertices.end()
            if (first_vtx.unset()) {
                // Mark the first vertex in our domain
                if constexpr (dbg) { std::cout << "Mark first vertex at " << v << std::endl; }
                first_vtx = v;
            } // else "Don't update first_vtx as it was already set."

            // Find the neighbour of this vertex, if possible. Can't
            // do this if it's a boundary vertex, but nothing happens
            // in that case.
            Flt next_neighb_dom = std::numeric_limits<Flt>::max();

            sm::vec<Flt, 2> neighb_vtx = walk_to_neighbour (f, v, next_neighb_dom);
            if constexpr (dbg) { std::cout << "walk_to_neighbour returned with vertex " << neighb_vtx << std::endl; }
            v.vn = neighb_vtx;

            // Walk to the next vertex
            next_neighb_dom = std::numeric_limits<Flt>::max();
            sm::vec<Flt, 2> next_vtx = walk_to_next (f, v, next_neighb_dom);
            if constexpr (dbg) { std::cout << "starting from " << v.v << ", walk_to_next returned with vertex " << next_vtx << std::endl; }

            if constexpr (dbg) { std::cout << "Closing dv at " << dv->v << " with f=" << dv->f << " [" << dv->neighb << "]\n"; }
            dv->closed = true;
            domain.vertices.push_back (v);

            typename std::list<DirichVtx<Flt>>::iterator dv2 = vertices.begin();
            if (first_vtx.compare (next_vtx) == false) {
                // Find a dv which matches next_vtx.
                bool matched_next_vertex = false;
                if constexpr (dbg) { std::cout << "Search vertices for " << next_vtx << std::endl; }
                while (dv2 != vertices.end()) {
                    // Instead of: dv2 = next_vtx;, do:
                    if (dv2->closed == false && dv2->compare (next_vtx) == true) {
                        // vertex has correct coordinate. Check it has correct neighbours.
                        if constexpr (dbg) {
                            std::cout << "Coordinate match. Is dv2->f == v.f? " << dv2->f <<" == "<< v.f << "?\n";
                            std::cout << "Is dv2->neighb[0] == next_neighb_dom? " << dv2->neighb[0]
                                      << " == " << next_neighb_dom << "?\n";
                            std::cout << "Is dv2->neighb[1] == v.neighb[0]? " << dv2->neighb[1]
                                      << " == " << v.neighb[0] << "?\n";
                            std::cout << "(v.neighb[1] = " << v.neighb[1] << ")\n";
                        }
                        if (dv2->f == v.f
                            && dv2->neighb[1] == v.neighb[0]
                            && dv2->neighb[0] == next_neighb_dom) {
                            // Match for current dv2
                            if constexpr (dbg) { std::cout << "Match\n"; }
                            matched_next_vertex = true;
                            dv = dv2;
                            break;
                        } // else no match
                    }
                    ++dv2;
                }
                if (!matched_next_vertex) {
                    if constexpr (dbg) { std::cout << "Failed to find a match for the next_vtx which walk_to_next found. Return false.\n"; }
                    return false;
                }

            } else {
                if constexpr (dbg) { std::cout << "walk_to_next() arrived back at the first vertex. Return true\n"; }
                return true;
            }

            // Now move on to the next vertex in the domain, re-calling process_domain
            // recursively, or exiting if we got to the start of the domain perimeter. We
            // shouldn't get anywhere close to the recursion limit in this system.
            if (dv->onBoundary == false) {
                if constexpr (dbg) { std::cout << "Recursively call process_domain...\n"; }
                return process_domain (hg, f, dv, vertices, domain, first_vtx);
            }

            if constexpr (dbg) { std::cout << "Arrived at a boundary vertex. Return false (because we didn't find a full domain\n"; }
            return false;
        }

        /*!
         * Determine the locations of the vertices on a Hex grid which are surrounded by three
         * different values of @f. @f is indexed by the hexgrid @hg. Return a list containing lists
         * of the vertices, each of which define a domain.
         */
        static std::list<DirichDom<Flt>>
        dirichlet_vertices (sm::hexgrid* hg, std::vector<Flt>& f, std::list<DirichVtx<Flt>>& vertices) {

            // 1. Go though and find a list of all vertices, in no particular order.  This will
            // lead to duplications because >1 domain for a given ID, f, is possible early in
            // simulations. From this list, I can find vertex sets, whilst deleting from the list
            // until it is empty, and know that I will have discovered all the domain vertex sets.
            // list<DirichVtx<Flt> > vertices;
            std::list<sm::hex>::iterator h = hg->hexen.begin();
            while (h != hg->hexen.end()) {
                vertex_test (hg, f, h, vertices);
                // Move on to the next hex in hexen
                ++h;
            }

            // 2. Delete from the list<DirichVtx> and construct a list<list<DirichVtx>> of all the
            // domains. The list<DirichVtx> for a single domain should be ordered so that the
            // perimeter of the domain is traversed. I have to do Dirichlet domain boundary walks
            // to achieve this (to disambiguate between vertices from separate, but same-ID
            // domains).
            std::list<DirichDom<Flt>> dirich_domains;
            typename std::list<DirichVtx<Flt>>::iterator dv = vertices.begin();
            //unsigned int domcount = 0;
            while (dv != vertices.end() /* && domcount++ < 3 */) {
                DirichDom<Flt> one_domain;
                DirichVtx<Flt> first_vtx;
                if (dv->hi->boundaryHex() == true) {
                    if constexpr (dbg) { std::cout << "Don't process hexes on the boundary\n"; }
                    dv->closed = true;
                    dv++;
                } else {
                    bool success = process_domain (hg, f, dv, vertices, one_domain, first_vtx);
                    dv++;
                    if (success) {
                        // Set the identity, f of the domain
                        one_domain.f = one_domain.vertices.front().f;
                        if constexpr (dbg) { std::cout << "Found outline of a domain (ID " << one_domain.f << ")\n"; }

                        // Calculate the area of the domain.
                        one_domain.compute_area (hg, f);
                        one_domain.compute_edge_deviation();
                        // Add the domain
                        dirich_domains.push_back (one_domain);
                    } // process_domain failed to find the outline of a domain
                }
            }

            return dirich_domains;
        }

        //! Count number of instances of the value @val in the vector of values @vec.
        static unsigned int
        count_up (const std::vector<Flt>& vec, const Flt val) {
            unsigned int count = 0;
            for (auto v : vec) { count += (v == val) ? 1 : 0; }
            return count;
        }

        /*!
         * Take a list of Dirichlet domains and compute a metric for the Dirichlet-ness
         * of the vertices after Honda1983. Return the overall Honda measure, and place
         * all delta_j values in \a delta_j.
         */
        static Flt
        dirichlet_analyse (std::list<DirichDom<Flt>>& doms, std::vector<sm::vec<Flt, 2>>& d_centres)
        {
            std::map<Flt, Flt> delta_js;
            return ShapeAnalysis::dirichlet_analyse (doms, d_centres, delta_js);
        }

        static Flt
        dirichlet_analyse (std::list<DirichDom<Flt>>& doms, std::vector<sm::vec<Flt, 2>>& d_centres,
                           std::map<Flt, Flt>& delta_js)
        {
            Flt sum_delta_j = 0.0;
            Flt sum_areas = 0.0;
            d_centres.clear();
            delta_js.clear();
            auto di = doms.begin();
            while (di != doms.end()) {
                sm::vec<Flt, 2> P;
                // metric returns is delta_j, the mean_sos_per_vertex
                Flt delta_j = di->dirichlet_analyse_single_domain (P);
                //cout << "Domain delta_j = " << delta_j << endl;
                sum_delta_j += delta_j;
                delta_js[di->f] = delta_j;
                d_centres.push_back (P);
                // Sum up area too.
                sum_areas += di->area;
                ++di;
            }

            // The Ns cancel out of the equation given in Honda1983 as "For practical calculation."
            // on p196.
            return sum_delta_j/sum_areas;
        }

    }; // ShapeAnalysis

} // namespace morph
