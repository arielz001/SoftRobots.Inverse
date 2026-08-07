
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
* This plugin is distributed under the GNU AGPL v3 (Affero General            *
* Public License) license.                                                    *
*                                                                             *
* Authors: Christian Duriez, Eulalie Coevoet, Yinoussa Adagolodjo             *
*                                                                             *
* (c) 2023 INRIA                                                              *
*                                                                             *
* Contact information: https://project.inria.fr/softrobot/contact/            *
******************************************************************************/

#pragma once

#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/logging/Messaging.h>

#include <SoftRobots.Inverse/component/constraint/CameraProjectionPointsEffector.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace softrobotsinverse::constraint
{

using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::core::objectmodel::ComponentState;


template<class DataTypes>
CameraProjectionPointsEffector<DataTypes>::CameraProjectionPointsEffector(MechanicalState* object)
    : Effector<DataTypes>(object)
    , softrobots::constraint::CameraProjectionPointsModel<DataTypes>(object)
    , d_effectorGoal(initData(&d_effectorGoal, "effectorGoal",
                    "Desired 3D positions or target for the effector."))
    , d_pointParameters(initData(&d_pointParameters, "pointParameters",
                    "Target 2D point coordinates detected in image [u_target, v_target]."))
    , d_cameraPosition(initData(&d_cameraPosition, sofa::type::Vec3d(0.0, 0.0, -13.2), "cameraPosition",
                    "Position of the camera in world coordinates [x, y, z]."))
{
}


template<class DataTypes>
CameraProjectionPointsEffector<DataTypes>::~CameraProjectionPointsEffector()
{
}


template<class DataTypes>
void CameraProjectionPointsEffector<DataTypes>::init()
{
    softrobots::constraint::CameraProjectionPointsModel<DataTypes>::init();

    if(!d_effectorGoal.isSet())
    {
        msg_warning(this) << "TargetPosition not defined. Default value assigned (" << Coord() << ").";
        setTargetDefaultValue();
    }

    if(d_indices.getValue().size() != d_effectorGoal.getValue().size())
        resizeData();
}

template<class DataTypes>
void CameraProjectionPointsEffector<DataTypes>::setTargetDefaultValue()
{
    WriteAccessor<sofa::Data<VecCoord>> defaultTarget = d_effectorGoal;
    defaultTarget.resize(1);
    defaultTarget[0] = Coord();
}

template<class DataTypes>
void CameraProjectionPointsEffector<DataTypes>::resizeData()
{
    if(d_indices.getValue().size() < d_effectorGoal.getValue().size())
    {
        msg_warning(this) << "Indices size is lower than target size, some targets will not be considered.";
    }
    else
    {
        msg_warning(this) << "Indices size is larger than target size. Launch resize process.";
        WriteAccessor<sofa::Data<sofa::type::vector<unsigned int>>> indices = d_indices;
        indices.resize(d_effectorGoal.getValue().size());
    }
}


// PROJECTION 3D FROM 2D
Eigen::Vector3d ProjectionPosition(double u_point, double v_point,
                                   double cx_camera, double cy_camera, 
                                   double x_focal_length, double y_focal_length,
                                   double estimated_Z)
{
    double X = (u_point - cx_camera) * estimated_Z / x_focal_length;
    double Y = (v_point - cy_camera) * estimated_Z / y_focal_length;

    return Eigen::Vector3d(X, Y, estimated_Z);
}


// PROJECTION 2D FROM 3D 

Eigen::Vector2d calculateProjectedPoint(
    double x, double y, double z,
    const sofa::type::Vec2d& focalLength,
    const sofa::type::Vec2d& principalPoint,
    const sofa::type::Vec3d& cameraPos)
{
    // Coordenadas relativas a la cámara
    double x_rel = x - cameraPos[0];
    double y_rel = y - cameraPos[1];
    double z_rel = z - cameraPos[2];

    if (std::abs(z_rel) < 1e-6)
        z_rel = (z_rel >= 0) ? 1e-6 : -1e-6;

    // Proyección Pinhole 2D (u, v)
    double u = focalLength[0] * (x_rel / z_rel) + principalPoint[0];
    double v = focalLength[1] * (y_rel / z_rel) + principalPoint[1]; 

    return Eigen::Vector2d(u, v);
}



template<class DataTypes>
void CameraProjectionPointsEffector<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* cParams,
                                                                 sofa::linearalgebra::BaseVector *resV,
                                                                 const sofa::linearalgebra::BaseVector *Jdx)
{
    if (d_componentState.getValue() != ComponentState::Valid)
        return;

    SOFA_UNUSED(cParams);

    ReadAccessor<sofa::Data<VecCoord>> x = m_state->readPositions();
    
    const auto& pointParams         = sofa::helper::getReadAccessor(d_pointParameters); 
    const sofa::type::Vec3d cameraPosition = d_cameraPosition.getValue();
    const auto& weight              = sofa::helper::getReadAccessor(d_weight);
    const auto& indices             = sofa::helper::getReadAccessor(d_indices);
    const auto& constraintIndex     = sofa::helper::getReadAccessor(d_constraintIndex);

    // intrinisics
    const sofa::type::Vec2d fLength = d_focalLength.getValue();
    const sofa::type::Vec2d pPoint  = d_principalPoint.getValue();

    // target 2D point in image coordinates
    Eigen::Vector2d P_target;
    P_target << pointParams[0], pointParams[1];

    unsigned int index = 0;

    for (unsigned int i = 0; i < indices.size(); i++)
    {
        const auto& coord = x[indices[i]];//current position

        double x_pos = coord[0];
        double y_pos = coord[1];
        double z_pos = coord[2];

        // 2d projection
        Eigen::Vector2d P_current = calculateProjectedPoint(
            x_pos, y_pos, z_pos, 
            fLength, 
            pPoint,
            cameraPosition
        );

        // difference
        Eigen::Vector2d P_diff = P_current - P_target;

        msg_info(this) << "P_diff (u, v): " << P_diff.transpose();

        // jacobian 
        for (sofa::Size j = 0; j < 2; j++)
        {
            Real dfree = Jdx->element(constraintIndex + index) + P_diff[j] * weight[j];
            resV->set(constraintIndex + index, dfree);
            index++;
        }
    }
}

} // namespace softrobotsinverse::constraint