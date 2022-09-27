/*
* single_construction.cpp
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

#include "tetrahedra.hpp"
#include "geometry.hpp"
#include "matplotlibcpp.h"
#include <iostream>
#include <vector>
#include <random>
#include <thread>

#define KNRM  "\033[0m"
#define KRED  "\033[31m"
#define KGRN  "\033[32m"
#define KYEL  "\033[33m"
#define KBLU  "\033[34m"
#define KMAG  "\033[35m"
#define KCYN  "\033[36m"
#define KWHT  "\033[37m"

using namespace geometry;
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds

namespace plt = matplotlibcpp;

/** @brief linspace function with given min and max **/
vector<double> linspace(double min, double max, double n)
{
    vector<double> linspaced;
    double delta = (max - min) / (n - 1.0);
    linspaced.push_back(min);
    
    for (int i = 1; i < (int)n; i++)
    {
        linspaced.push_back(linspaced[i-1] + delta);
    }

    return linspaced;
}

int main()
{
    geometry::tetrahedra t;
    geometry::helper h;
    time_point<std::chrono::system_clock> start;

    std::vector<geometry::tetrahedra::tetrahedron> solids;

    int number_of_solids = 2;
    std::random_device dev;
    std::mt19937 generator(dev());

    std::uniform_real_distribution<double> dis(-2.0, 2.0);
    std::uniform_real_distribution<double> fac(0.5, 1.0);

    Eigen::Vector3d random_point = Eigen::Vector3d(
        dis(generator), dis(generator), dis(generator)/2);
    printf("random_point %lf %lf %lf\n", random_point.x(), random_point.y(), random_point.z());

    vector<double> s;
    vector<Eigen::Vector3d> c_p;

    for (int i = 0; i < number_of_solids; i++)
    {
        Eigen::Vector3d vert[4];
        Eigen::Vector3d offset = Eigen::Vector3d(
            dis(generator), dis(generator), dis(generator)/2);
        vert[0] = fac(generator) * Eigen::Vector3d(0, 0, 1) + offset;
        vert[1] = fac(generator) * Eigen::Vector3d(1, 0, -1) + offset;
        vert[2] = fac(generator) * Eigen::Vector3d(-1, -1, -1) + offset;
        vert[3] = fac(generator) * Eigen::Vector3d(-1, 1, -1) + offset;
        
        for (auto &vertices : vert)
            printf("obj(%d) vertices %lf %lf %lf\n", i, vertices.x(), vertices.y(), vertices.z());

        start = system_clock::now();
        geometry::tetrahedra::tetrahedron t_obj = t.construct_tetrahedron(vert);
        auto t_t = duration<double>(system_clock::now() - start).count() * 1000;
        printf("%sconstruct_tetrahedron %lf ms %s\n", KRED, t_t, KNRM);

        solids.push_back(t_obj);

        Eigen::Vector3d c_d;
        double seperation;

        start = system_clock::now();
        bool col = h.point_poly_collision(t_obj, random_point, c_d, seperation);
        auto p_c = duration<double>(system_clock::now() - start).count() * 1000;
        printf("%spoint_poly_collision %lf ms %s\n", KGRN, p_c, KNRM);
        
        printf("col(%s) distance %lf closest_point %lf %lf %lf\n", 
            col ? "y" : "n", seperation, c_d.x(), c_d.y(), c_d.z());
        
        c_p.push_back(c_d);
        s.push_back(seperation);
    }

    // vector<Eigen::Vector3d> v_m = h.construct_minkowski_difference(solids[0], solids[1]);
    // printf("minkowski vertices (%ld)\n", v_m.size());

    /** @brief Visualization **/
    // 6 pairs 1-2, 2-3, 3-4, 4-1, 1-3, 2-4
    Eigen::MatrixXi m(6,2);
    m << 0, 1,
        1, 2,
        2, 3,
        3, 0,
        0, 2,
        1, 3;
    
    int division = 50;
    
    std::vector<double> x, y, z;
    std::vector<double> r_x, r_y, r_z;
    std::vector<double> x_c_f, y_c_f, z_c_f;
    for (int n = 0; n < number_of_solids; n++)
    {
        geometry::tetrahedra::tetrahedron t_obj = solids[n];

        /** @brief Push back edge lines xyz **/
        for (int i = 0; i < 6; i++)
        {
            Eigen::Vector3d v0 = t_obj.v[m(i,0)];
            Eigen::Vector3d v1 = t_obj.v[m(i,1)];
            std::vector<double> t_x = linspace(v0.x(), v1.x(), division);
            std::vector<double> t_y = linspace(v0.y(), v1.y(), division);
            std::vector<double> t_z = linspace(v0.z(), v1.z(), division);
            x.insert(std::end(x), std::begin(t_x), std::end(t_x));
            y.insert(std::end(y), std::begin(t_y), std::end(t_y));
            z.insert(std::end(z), std::begin(t_z), std::end(t_z));

            printf("(%d-%d) x (%.3lf to %.3lf), y (%.3lf to %.3lf), z (%.3lf to %.3lf)\n", 
                m(i,0), m(i,1), v0.x(), v1.x(), v0.y(), v1.y(), v0.z(), v1.z());
        }

        /** @brief Push back centroid xyz **/
        x.push_back(t_obj.c.x());
        y.push_back(t_obj.c.y());
        z.push_back(t_obj.c.z());

        std::vector<double> t_r_x = linspace(random_point.x(), c_p[n].x(), division);
        std::vector<double> t_r_y = linspace(random_point.y(), c_p[n].y(), division);
        std::vector<double> t_r_z = linspace(random_point.z(), c_p[n].z(), division);
        r_x.insert(std::end(r_x), std::begin(t_r_x), std::end(t_r_x));
        r_y.insert(std::end(r_y), std::begin(t_r_y), std::end(t_r_y));
        r_z.insert(std::end(r_z), std::begin(t_r_z), std::end(t_r_z));

        /** @brief Push back random xyz **/
        r_x.push_back(random_point.x());
        r_y.push_back(random_point.y());
        r_z.push_back(random_point.z());
        
        /** @brief Push back center of faces **/
        for (int i = 0; i < 4; i++)
        {
            Eigen::Vector3d v0 = t_obj.v[t_obj.f[i].v_i[0]];
            Eigen::Vector3d v1 = t_obj.v[t_obj.f[i].v_i[1]];
            Eigen::Vector3d v2 = t_obj.v[t_obj.f[i].v_i[2]];
            
            Eigen::Vector3d c_f = (v0 + v1 + v2) / 3;
            Eigen::Vector3d c_ff = c_f + t_obj.n[i] * 0.5;
            std::vector<double> t_x = linspace(c_f.x(), c_ff.x(), 10);
            std::vector<double> t_y = linspace(c_f.y(), c_ff.y(), 10);
            std::vector<double> t_z = linspace(c_f.z(), c_ff.z(), 10);
            x_c_f.insert(std::end(x_c_f), std::begin(t_x), std::end(t_x));
            y_c_f.insert(std::end(y_c_f), std::begin(t_y), std::end(t_y));
            z_c_f.insert(std::end(z_c_f), std::begin(t_z), std::end(t_z));

            printf("(%d) (%.3lf %.3lf %.3lf)\n", 
                i, t_obj.n[i].x(), t_obj.n[i].y(), t_obj.n[i].z());
        }

    }

    // std::vector<double> x_m, y_m, z_m;
    // std::vector<double> x_m_c, y_m_c, z_m_c;
    // for (int n = 0; n < (int)v_m.size(); n++)
    // {
    //     x_m.push_back(v_m[n].x());
    //     y_m.push_back(v_m[n].y());
    //     z_m.push_back(v_m[n].z());
    // }
    // /** @brief Push back origin xyz **/
    // x_m_c.push_back(0.0);
    // y_m_c.push_back(0.0);
    // z_m_c.push_back(0.0);

    // Plot the tetrahedrons in 3d space
    plt::scatter(x, y, z, 0.4, {}, 1);
    plt::scatter(x_c_f, y_c_f, z_c_f, 1.5, {}, 1);
    plt::scatter(r_x, r_y, r_z, 2.0, {}, 1);
    plt::xlabel("x/m");
    plt::ylabel("y/m");
    plt::set_zlabel("z/m");
    // Plot the Minkowski difference
    // plt::scatter(x_m, y_m, z_m, 1.0, {}, 2);
    // plt::scatter(x_m_c, y_m_c, z_m_c, 1.5, {}, 2);
    // plt::xlabel("x/m");
    // plt::ylabel("y/m");
    // plt::set_zlabel("z/m");

    plt::show();

    return 0;
}