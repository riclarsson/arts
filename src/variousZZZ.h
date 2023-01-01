/* Copyright (C) 2021 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

#ifndef variousZZZ_h
#define variousZZZ_h

#include "gridded_fields.h"


/** Checks if a vector is a valid line-of-sight (direction vector)

    The function gives an error message if this is not the case.

    @param[in] name  Name of the variable, if there is an error
    @param[in] los   The vector

    @author Patrick Eriksson 
    @date   2022-12-31
 */
void chk_rte_los(const String& name,
                 ConstVectorView los);

/** Checks if a vector is a valid position

    The function gives an error message if this is not the case.

    @param[in] name  Name of the variable, if there is an error
    @param[in] pos   The vector

    @author Patrick Eriksson 
    @date   2022-12-31
 */
void chk_rte_pos(const String& name,
                 ConstVectorView pos);

/** Performs all needed checks of refellipsoid

    @param[in]   refellipsoid   As the WSV with the same name.

    @author Patrick Eriksson 
    @date   2021-07-30
 */
void chk_refellipsoidZZZ(ConstVectorView refellipsoid);

/** Checks if a matrix is a valid sensor_pos

    @param[in] name        Name of the variable, if there is an error
    @param[in] sensor_pos  The matrix

    @author Patrick Eriksson 
    @date   2021-07-30
 */
void chk_sensor_pos(const String& name,
                    ConstMatrixView sensor_pos);


/** Checks if a matrix is a valid sensor_los

    @param[in] name        Name of the variable, if there is an error
    @param[in] sensor_pos  The matrix

    @author Patrick Eriksson 
    @date   2021-07-30
 */
void chk_sensor_los(const String& name,
                    ConstMatrixView sensor_los);

/** Checks if two matrices are valid combination of sensor_posLlos

    @param[in] name1       Name of the pos variable, if there is an error
    @param[in] sensor_pos  A first matrix
    @param[in] name2       Name of the los variable, if there is an error
    @param[in] sensor_los  A second matrix

    @author Patrick Eriksson 
    @date   2021-07-30
 */
void chk_sensor_poslos(const String& name1,
                       ConstMatrixView sensor_pos,
                       const String& name2,
                       ConstMatrixView sensor_los);

/** Calculates the geometrical length to the surface

   A negative length is returned if the geomtrical path has no intersection
   with the surface.

   The function also check that the observation position is actually above the
   surface. 

   @param[in] rte_pos           As the WSV with the same name.
   @param[in] rte_los           As the WSV with the same name.
   @param[in] ecef              rte_pos in ECEF.
   @param[in] decef             rte_los in ECEF.
   @param[in] atmosphere_dim    As the WSV with the same name.
   @param[in] refellipsoid      As the WSV with same name.
   @param[in] surface_elevation As the WSV with same name.
   @param[in] surface_search_accuracy See WSM IntersectionGeometricalWithSurface.
              surface_search_safe     See WSM IntersectionGeometricalWithSurface.       

   @return   Length to the surface.

   @author Patrick Eriksson
   @date   2021-08-06
 */
Numeric find_crossing_with_surface_z(const Vector rte_pos,
                                     const Vector rte_los,
                                     const Vector ecef,
                                     const Vector decef,
                                     const Index& atmosphere_dim,
                                     const Vector& refellipsoid,
                                     const GriddedField2& surface_elevation,
                                     const Numeric& surface_search_accuracy,
                                     const Index& surface_search_safe);

/** Returns surface elevation at lat and lon of a position

    See *surface_elevation* for allowed longitude grids and applied
    extrapolation.

   @param[in] pos               A position vector.
   @param[in] surface_elevation As the WSV with the same name.

   @return   Elevation.

   @author Patrick Eriksson
   @date   2022-12-31
 */
Numeric surface_z_at_pos(const Vector& pos,
                         const GriddedField2& surface_elevation);

#endif  // variousZZZ_h