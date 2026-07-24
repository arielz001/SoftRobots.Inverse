/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
*                       Plugin SoftRobots.Inverse                             *
*                                                                             *
* Authors: Christian Duriez, Eulalie Coevoet, Yinoussa Adagolodjo             *
* (c) 2023 INRIA                                                              *
******************************************************************************/
#define SOFTROBOTS_INVERSE_CameraProjectionEffector_CPP
#include <SoftRobots.Inverse/component/config.h>
#include <sofa/core/ObjectFactory.h>
#include <SoftRobots.Inverse/component/constraint/CameraProjectionEffector.inl>

namespace softrobotsinverse::constraint
{

using namespace sofa::defaulttype;
using sofa::core::ConstraintParams;

////////////////////////////////////////////    FACTORY    //////////////////////////////////////////////
using namespace sofa::helper;

// Registering the component in the SOFA ObjectFactory
int CameraProjectionEffectorClass = sofa::core::RegisterObject(
        "CameraProjectionEffector constrains a 3D rigid model to match target 2D image ellipse parameters "
        "[cx, cy, major, minor, angle] via camera projection.")
    .add< CameraProjectionEffector<Rigid3Types> >(true) // Rigid3Types is the primary default template
    .add< CameraProjectionEffector<Vec3Types> >()
;

////////////////////////////////////////////////////////////////////////////////////////////////////////

// Explicit template instantiation
template class SOFA_SOFTROBOTS_INVERSE_API CameraProjectionEffector<sofa::defaulttype::Rigid3Types>;
template class SOFA_SOFTROBOTS_INVERSE_API CameraProjectionEffector<sofa::defaulttype::Vec3Types>;

} // namespace softrobotsinverse::constraint