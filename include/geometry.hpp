/*
* geometry.hpp
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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <string>
#include <chrono>
#include <cmath>
#include <vector>
#include <float.h>
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
    class helper
    {
        private:
            Eigen::Vector3d find_furthest_point(
                Eigen::Vector3d dir, Eigen::Vector3d *vert)
            {
                Eigen::Vector3d p;
                double m_d = -FLT_MAX;
        
                for (int i = 0; i < 4; i++) 
                {
                    double dist = vert[i].dot(dir);
                    if (dist > m_d) 
                    {
                        m_d = dist;
                        p = vert[i];
                    }
                }
        
                return p;
            }

            Eigen::Vector3d find_support(
                tetrahedra::tetrahedron t_0, tetrahedra::tetrahedron t_1, Eigen::Vector3d dir)
            {        
                return find_furthest_point(dir, t_0.v) - find_furthest_point(-dir, t_1.v);
            }

            bool point_inside_tetrahedron(
                tetrahedra::tetrahedron t, Eigen::Vector3d p)
            {
                for (auto &normal : t.n)
                {
                    Eigen::Vector3d p_c = p - t.c;
                    double d = normal.dot(p_c);

                    if (d < 0)
                        return false;
                }
                return true;
            }

            // The algorithm is based on 
            // "David Eberly, 'Distance Between Point and Triangle in 3D',
            // Geometric Tools, LLC, (1999)"
            // http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
            // https://math.stackexchange.com/questions/140195/estimate-distance-between-point-and-triangle
            void closed_point_to_triangle_gt(
                Eigen::Vector3d& d_p, double& m_d, Eigen::Vector3d *v, Eigen::Vector3d query)
            {
                // [dist,PP0] = pointTriangleDistance(TRI,P)
                // https://www.mathworks.com/matlabcentral/fileexchange/22857-distance-between-a-point-and-a-triangle-in-3d

                Eigen::Vector3d B = v[2];
                Eigen::Vector3d E0 = v[0] - B;
                Eigen::Vector3d E1 = v[1] - B;

                Eigen::Vector3d D = B - query;
                double a = E0.dot(E0);
                double b = E0.dot(E1);
                double c = E1.dot(E1);
                double d = E0.dot(D);
                double e = E1.dot(D);
                double f = D.dot(D);

                double det = a*c - b*b;
                double s = b*e - c*d;
                double t = b*d - a*e;
                
                double sqrDistance = 0;

                if ((s+t) <= det)
                {
                    if (s < 0)
                    {
                        // Start of region 4
                        if (t < 0)
                        {
                            if (d < 0)
                            {
                                t = 0;
                                if (-d >= a)
                                {
                                    s = 1;
                                    sqrDistance = a + 2*d + f;
                                }
                                else
                                {
                                    s = -d/a;
                                    sqrDistance = d*s + f;
                                }
                            }
                            else
                            {
                                s = 0;
                                if (e >= 0)
                                {
                                    t = 0;
                                    sqrDistance = f;
                                }
                                else
                                {
                                    if (-e >= c)
                                    {
                                        t = 1;
                                        sqrDistance = c + 2*e + f;
                                    }
                                    else
                                    {
                                        t = -e/c;
                                        sqrDistance = e*t + f;
                                    }
                                }
                            }
                        }
                        // End of region 4
                        // Start of region 3
                        else
                        {
                            s = 0;
                            if (e >= 0)
                            {
                                t = 0;
                                sqrDistance = f;
                            }
                            else
                            {
                                if (-e >= c)
                                {
                                    t = 1;
                                    sqrDistance = c + 2*e +f;
                                }
                                else
                                {
                                    t = -e/c;
                                    sqrDistance = e*t + f;
                                }
                            }
                        }
                        // End of region 3
                    }
                    else
                    {
                        // Start of region 5
                        if (t < 0)
                        {
                            t = 0;
                            if (d >= 0)
                            {
                                s = 0;
                                sqrDistance = f;
                            }
                            else
                            {
                                if (-d >= a)
                                {
                                    s = 1;
                                    sqrDistance = a + 2*d + f;
                                }
                                else
                                {
                                    s = -d/a;
                                    sqrDistance = d*s + f;
                                }
                            }
                        }
                        // End of region 5
                        // Start of region 0
                        else
                        {
                            double invDet = 1/det;
                            s = s * invDet;
                            t = t * invDet;
                            sqrDistance = 
                                s * (a*s + b*t + 2*d) + t * (b*s + c*t + 2*e) + f;
                        }
                        // End of region 0
                    }
                }
                else
                {
                    // Start of region 2
                    if (s < 0)
                    {
                        double tmp0 = b + d;
                        double tmp1 = c + e;
                        if (tmp1 > tmp0) // minimum on edge s+t=1
                        {
                            double numer = tmp1 - tmp0;
                            double denom = a - 2*b + c;
                            if (numer >= denom)
                            {
                                s = 1;
                                t = 0;
                                sqrDistance = a + 2*d + f;
                            }
                            else
                            {
                                s = numer/denom;
                                t = 1-s;
                                sqrDistance = 
                                    s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                            }
                        }
                        else // minimum on edge s=0
                        {
                            s = 0;
                            if (tmp1 <= 0)
                            {
                                t = 1;
                                sqrDistance = c + 2*e + f;
                            }
                            else
                            {
                                if (e >= 0)
                                {
                                    t = 0;
                                    sqrDistance = f;
                                }
                                else
                                {
                                    t = -e/c;
                                    sqrDistance = e*t + f;
                                }
                            }
                        }
                    }
                    // End of region 2
                    else
                    {
                        // Start of region 6
                        if (t < 0)
                        {
                            double tmp0 = b + e;
                            double tmp1 = a + d;
                            if (tmp1 > tmp0)
                            {
                                double numer = tmp1 - tmp0;
                                double denom = a-2*b+c;
                                if (numer >= denom)
                                {
                                    t = 1;
                                    s = 0;
                                    sqrDistance = c + 2*e + f;
                                }
                                else
                                {
                                    t = numer/denom;
                                    s = 1 - t;
                                    sqrDistance = 
                                        s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                                }
                            }
                            else  
                            {
                                t = 0;
                                if (tmp1 <= 0)
                                {
                                    s = 1;
                                    sqrDistance = a + 2*d + f;
                                }
                                else
                                {
                                    if (d >= 0)
                                    {
                                        s = 0;
                                        sqrDistance = f;
                                    }
                                    else
                                    {
                                        s = -d/a;
                                        sqrDistance = d*s + f;
                                    }
                                }
                            }
                        }
                        // End of region 6
                        // Start of region 1
                        else
                        {
                            double numer = c + e - b - d;
                            if (numer <= 0)
                            {
                                s = 0;
                                t = 1;
                                sqrDistance = c + 2*e + f;
                            }
                            else
                            {
                                double denom = a - 2*b + c;
                                if (numer >= denom)
                                {
                                    s = 1;
                                    t = 0;
                                    sqrDistance = a + 2*d + f;
                                }
                                else
                                {
                                    s = numer/denom;
                                    t = 1-s;
                                    sqrDistance = 
                                        s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                                }
                            }
                        }
                        // End of region 1
                    }
                }

                // account for numerical round-off error
                if (sqrDistance < 0)
                    sqrDistance = 0;

                m_d = sqrt(sqrDistance);
                d_p = B + s * E0 + t * E1;
            }
                


        public:

            vector<Eigen::Vector3d> construct_minkowski_difference(
                tetrahedra::tetrahedron t_0, tetrahedra::tetrahedron t_1)
            {
                vector<Eigen::Vector3d> list;

                std::vector<tetrahedra::tetrahedron> t_list;
                t_list.push_back(t_0);
                t_list.push_back(t_1);

                for (int n = 0; n < 2; n++)
                {
                    for (auto &vertices : t_list[n].v)
                    {
                        Eigen::Vector3d dir = vertices - t_list[n].c;
                        dir = dir / dir.norm();
                        Eigen::Vector3d supp = find_support(t_list[n%2], t_list[(n+1)%2], dir);
                        std::cout << "direction " << dir.transpose() << 
                            " support " << supp.transpose() << " (" <<
                            n%2 << "-" << (n+1)%2 << ")" << std::endl;
                        bool not_added = true;
                        if (!list.empty())
                        {
                            for (auto &vert : list)
                            {
                                if ((supp - vert).norm() < 1E-8)
                                {
                                    not_added = false;
                                    break;
                                }
                            }
                        }

                        if (not_added)
                            list.push_back(supp);
                    }
                }

                return list;
            }

            bool point_poly_collision(
                tetrahedra::tetrahedron t, Eigen::Vector3d p, 
                Eigen::Vector3d& c_d, double& seperation)
            {
                // Check if point is inside the tetrahedron
                bool point_inside = point_inside_tetrahedron(t, p);
                seperation = FLT_MAX;

                for (auto &face : t.f)
                {
                    Eigen::Vector3d v[3]; // vertices
                    v[0] = t.v[face.v_i[0]];
                    v[1] = t.v[face.v_i[1]];
                    v[2] = t.v[face.v_i[2]];
                    Eigen::Vector3d closest;
                    double distance; 
                    closed_point_to_triangle_gt(closest, distance, v, p);

                    if (abs(distance) < seperation)
                    {
                        seperation = abs(distance);
                        c_d = closest;
                    }
                }

                return point_inside;
            }

    };
}

#endif