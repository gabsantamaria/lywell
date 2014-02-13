/**********************************************************************
 *
 * EMCC_cone-gap-sphere.geo
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

lc = 0.0276842236587;
inch = 0.0254;
R = 2.947*inch;
cos = Cos(7.0/180.0*Pi);
Point(1) = {0,0,0,lc};
Point(2) = {23.821 * inch,0,0,lc/2.0};
Point(3) = {0.0, R*cos, 0, lc};
Point(4) = {-0.359*inch,0.0,0,lc};
Point(5) = {-0.359*inch,2.947*inch,0,lc};
Point(6) = {-(2.947+0.359)*inch,0,0,lc};
Line(1) = {2,3};
Point(7) = {0,2.697*inch,0,lc};
Point(8) = {-0.359*inch,2.697*inch,0,lc};
Line(2) = {3,7};
Line(3) = {7,8};
Line(4) = {8,5};
Circle(5) = {5,4,6};
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{1,2,3,4,5};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{6,9,13,17,21};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{24,27,31,35,39};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{42,45,49,53,57};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{60,63,67,71,75};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{78,81,85,89,93};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{96,99,103,107,111};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{114,117,121,125,129};
}
