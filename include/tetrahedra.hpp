/*
* tetrahedra.hpp
*
* ---------------------------------------------------------------------
* Copyright (C) 2022 Matthew (matthewoots at gmail.com)
*
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License
*  as published by the Free Software Foundation; either version 2
*  of the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
* ---------------------------------------------------------------------
*/

#ifndef TETRAHEDRA_H
#define TETRAHEDRA_H

#include <iostream>
#include <string>
#include <chrono>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

using namespace Eigen;
using namespace std;
using namespace std::chrono;

namespace geometry
{
    class tetrahedra
    {
        public:

        /** @brief Description of the faces using 3 edges **/
        struct face
        {
            int v_i[3]; // vertices
        };

        /** @brief Description of a tetrahedron **/
        struct tetrahedron
        {
            Eigen::Vector3d c; // centroid
            Eigen::Vector3d v[4]; // vertices
            Eigen::Vector3d n[4]; // normals
            face f[4]; // faces
        };

        inline tetrahedron construct_tetrahedron(Eigen::Vector3d *vert)
        {
            tetrahedron t;
            
            for (int i = 0; i < 4; i++)
                t.v[i] = vert[i];
            
            // time_point<std::chrono::system_clock> start = system_clock::now();
            t.c = get_centroid(vert);
            // auto t_c = duration<double>(system_clock::now() - start).count() * 1000;
            // printf("%sget_centroid %lf ms %s\n", KBLU, t_c, KNRM);

            // start = system_clock::now();
            construct_faces(t);
            // auto t_f = duration<double>(system_clock::now() - start).count() * 1000;
            // printf("%sconstruct_faces %lf ms %s\n", KBLU, t_f, KNRM);

            return t;
        }


        private:

        /** @brief 
         * @param 
        **/
        inline Eigen::Vector3d get_centroid(Eigen::Vector3d *v)
        {
            return (v[0] + v[1] + v[2]+ v[3]) / 4;
        }

        /** @brief 
         * @param 
         * @param
         * @param
        **/
        inline void construct_faces(tetrahedron& t)
        {
            // (face 0) = 0, 1, 2 (face 1) = 1, 2, 3 (face 2) = 2, 3, 0 (face 3) = 3, 0, 1
            Eigen::MatrixXi m(4,3);
            m << 0, 1, 2,
                1, 2, 3,
                2, 3, 0,
                3, 0, 1;
            
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 4; j++)
                    t.f[j].v_i[i] = m(j,i);
            
            // Resolve the windings problem and add the normals to the object
            for (int i = 0; i < 4; i++)
            {
                int x[3] = {t.f[i].v_i[0], t.f[i].v_i[1], t.f[i].v_i[2]};
                // https://stackoverflow.com/questions/40454789/computing-face-normals-and-winding
                Eigen::Vector3d v0 = t.v[x[0]] - t.v[x[2]];
                Eigen::Vector3d v1 = t.v[x[1]] - t.v[x[2]];
                Eigen::Vector3d c_p = Eigen::Vector3d{
                    v0[1] * v1[2] - v0[2] * v1[1],
                    v0[2] * v1[0] - v0[0] * v1[2],
                    v0[0] * v1[1] - v0[1] * v1[0]
                };
                t.n[i] = c_p.normalized();

                // Check whether normal is facing inwards or outwards
                // https://stackoverflow.com/a/40465409
                Eigen::Vector3d v_c = t.v[x[0]] - t.c;
                double d = t.n[i].dot(v_c);

                // the normal is inward so negate (each component of) the normal to get an outward normal.
                if (d < 0)
                {
                    int tmp = t.f[i].v_i[0];
                    t.f[i].v_i[0] = t.f[i].v_i[1];
                    t.f[i].v_i[1] = tmp;

                    t.n[i] = -t.n[i];

                    // v0 = t.v[t.f[i].v_i[0]] - t.v[x[2]];
                    // v1 = t.v[t.f[i].v_i[1]] - t.v[x[2]];
                    // c_p = Eigen::Vector3d{
                    //     v0[1] * v1[2] - v0[2] * v1[1],
                    //     v0[2] * v1[0] - v0[0] * v1[2],
                    //     v0[0] * v1[1] - v0[1] * v1[0]
                    // };
                    // t.n[i] = c_p.normalized();
                }
            }
        }

    };
}

#endif