/**********************************************************************
 *
 * EMCC_cone-sphere.geo
 *
 * Copyright (C) 2012 Idesbald Van den Bosch
 *
 * This file is part of Puma-EM.
 * 
 * Puma-EM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Puma-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Puma-EM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Suggestions/bugs : <vandenbosch.idesbald@gmail.com>
 *
 **********************************************************************/

lc = 0.0154548127642;
inch = 0.0254;
R = 2.947*inch;
cos = Cos(7.0/180.0*Pi);
Point(1) = {0,0,0,lc};
Point(2) = {23.821 * inch,0,0,lc/2.0};
Point(3) = {0.0, R*cos, 0, lc};
Point(4) = {-0.359*inch,0.0,0,0.1};
Point(5) = {-0.359*inch,2.947*inch,0,0.1};
Point(6) = {-(2.947+0.359)*inch,0,0,lc};
Line(1) = {2,3};
Circle(2) = {3,4,6};
Extrude Line {1, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {3, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {6, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {9, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {12, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {15, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {18, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {21, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {2, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {27, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {30, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {33, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {36, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {39, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {42, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {45, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
