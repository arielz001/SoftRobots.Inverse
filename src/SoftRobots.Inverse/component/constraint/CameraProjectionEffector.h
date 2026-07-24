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
#pragma once

#include <SoftRobots.Inverse/component/behavior/Effector.h>
#include <SoftRobots.Inverse/component/config.h>

#include <sofa/type/Vec.h>
#include <Eigen/Dense>

namespace softrobotsinverse::constraint
{
using softrobotsinverse::behavior::Effector;

/**
 * @class CameraProjectionEffector
 * @brief Effector component for SoftRobots.Inverse that projects 3D rigid poses 
 *        into a 2D image ellipse space [cx, cy, major, minor, angle].
 */
template< class DataTypes >
class CameraProjectionEffector : public Effector<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(CameraProjectionEffector, DataTypes), SOFA_TEMPLATE(Effector, DataTypes));

    typedef typename sofa::core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename DataTypes::VecCoord            VecCoord;
    typedef typename DataTypes::Coord               Coord;
    typedef typename DataTypes::Deriv               Deriv;
    typedef typename DataTypes::Real                Real;

public:
    CameraProjectionEffector(MechanicalState* object = nullptr);
    ~CameraProjectionEffector() override;

    /////////////// Inherited from Effector ///////////////////
    void init() override;

    void getConstraintViolation(const sofa::core::ConstraintParams* cParams,
                                sofa::linearalgebra::BaseVector *resV,
                                const sofa::linearalgebra::BaseVector *Jdx) override;
    ///////////////////////////////////////////////////////////

    // Goal data in 2D image space: [cx, cy, major, minor, angle]
    sofa::Data<VecCoord>            d_effectorGoal;

    // Camera intrinsic & 3D geometric parameters
    sofa::Data<sofa::type::Vec2>    d_focalLength;      ///< Focal length [fx, fy] in pixels
    sofa::Data<sofa::type::Vec2>    d_principalPoint;   ///< Principal point [u0, v0] in pixels
    sofa::Data<double>              d_ellipseRadius;    ///< Real 3D radius of the target ellipse in meters

    void setTargetDefaultValue();
    void resizeData();

    ////////////////////////// Inherited attributes ////////////////////////////
    using softrobots::behavior::SoftRobotsConstraint<DataTypes>::m_state;
    using softrobots::behavior::SoftRobotsConstraint<DataTypes>::d_componentState;
    using Effector<DataTypes>::d_indices;
    using Effector<DataTypes>::d_constraintIndex;
    using Effector<DataTypes>::getTarget;
    ///////////////////////////////////////////////////////////////////////////

private:
    /**
     * @brief Projects a 3D rigid pose into a 2D image ellipse [cx, cy, major, minor, angle].
     */
    Eigen::Vector5d computeProjectedEllipse(const Coord& pose3D);
};

#if !defined(SOFTROBOTS_INVERSE_CameraProjectionEFFECTOR_CPP)
extern template class SOFA_SOFTROBOTS_INVERSE_API CameraProjectionEffector<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace softrobotsinverse::constraint