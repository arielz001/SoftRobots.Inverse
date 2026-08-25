
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
#include <Eigen/Geometry>

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
    , d_pointCenter(initData(&d_pointCenter, "pointCenter",
                    "Target point center detected in image [cx, cy, a, b, angle]."))
    , d_cameraPosition(initData(&d_cameraPosition, sofa::type::Vec3d(0.0, 0.0, -13.2), "cameraPosition",
                    "Position of the camera in world coordinates [x, y, z]."))
    , d_cameraOrientation(initData(&d_cameraOrientation, sofa::type::Vec3d(0.0, 0.0, 0.0), "cameraOrientation",
                    "Orientation of the camera in Euler angles [pitch, yaw, roll] in degrees."))
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



// Eigen::Matrix<double, 2, 1> calculateProjectedPoint(
//     double x, double y, double z,
//     const sofa::type::Vec2d& focalLength,
//     const sofa::type::Vec2d& principalPoint,
//     const sofa::type::Vec3d& cameraPos)
// {
//     // relative coordinates of the camera 
//     double x_rel = x - cameraPos[0];
//     double y_rel = y - cameraPos[1];
//     double z_rel = z - cameraPos[2];
//     if (std::abs(z_rel) < 1e-3)
//         {
//             z_rel = (z_rel >= 0) ? 1e-3 : -1e-3;
//         }

//     // std::cout << "z_rel: " << z_rel << std::endl;
//     // 2d projection
//     double u = focalLength[0] * (x_rel / z_rel) + principalPoint[0];
//     double v = focalLength[1] * (y_rel / z_rel) + principalPoint[1]; 


//     Eigen::Matrix<double, 2, 1> point;
//     point << u, v;

//     return point;
// }


Eigen::Matrix<double, 2, 1> calculateProjectedPoint(
    double x, double y, double z,
    const sofa::type::Vec2d& focalLength,
    const sofa::type::Vec2d& principalPoint,
    const sofa::type::Vec3d& cameraPos,
    const sofa::type::Vec3d& cameraOrientation)
{
    double pitch = cameraOrientation[0] * M_PI / 180.0;
    double yaw   = cameraOrientation[1] * M_PI / 180.0;
    double roll  = cameraOrientation[2] * M_PI / 180.0; 

    Eigen::Matrix3d Rx, Ry, Rz;
    Rx = Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitX());
    Ry = Eigen::AngleAxisd(yaw,   Eigen::Vector3d::UnitY());
    Rz = Eigen::AngleAxisd(roll,  Eigen::Vector3d::UnitZ());

    Eigen::Matrix3d R = Rz * Ry * Rx;

    Eigen::Vector3d P_world(x, y, z);
    Eigen::Vector3d C(cameraPos[0], cameraPos[1], cameraPos[2]);
    Eigen::Vector3d P_rel = P_world - C;

    Eigen::Vector3d P_cam = R * P_rel;

    double x_cam = P_cam.x();
    double y_cam = P_cam.y();
    double z_cam = P_cam.z();

    if (std::abs(z_cam) < 1e-5) z_cam = 1e-5;

    double u = focalLength[0] * (x_cam / z_cam) + principalPoint[0];
    double v = focalLength[1] * (y_cam / z_cam) + principalPoint[1];

    Eigen::Matrix<double, 2, 1> point;
    point << u, v;

    return point;
}




template<class DataTypes>
void CameraProjectionPointsEffector<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* cParams,
                                                                 sofa::linearalgebra::BaseVector *resV,
                                                                 const sofa::linearalgebra::BaseVector *Jdx)
{
    if (d_componentState.getValue() != ComponentState::Valid)
        return;
 
    SOFA_UNUSED(cParams);

    // with this we can acces to the parameters
    ReadAccessor<sofa::Data<VecCoord>> x = m_state->readPositions();
    
    const auto& pointCenter         = sofa::helper::getReadAccessor(d_pointCenter); 
    // const double realRadius         = d_radiusEllipse.getValue();
    const sofa::type::Vec3d cameraPosition = d_cameraPosition.getValue(); // acess to camera position
    const sofa::type::Vec3d cameraOrientation = d_cameraOrientation.getValue(); // acess to camera orientation
    const auto& weight              = sofa::helper::getReadAccessor(d_weight);
    const auto& indices             = sofa::helper::getReadAccessor(d_indices);
    const auto& constraintIndex     = sofa::helper::getReadAccessor(d_constraintIndex);

    // calib parameters
    const sofa::type::Vec2d fLength = d_focalLength.getValue();
    const sofa::type::Vec2d pPoint  = d_principalPoint.getValue();

    // this is the objective to follow 
    // in this case we detect with other code and our camera the ellipse parameters 
    Eigen::Matrix<double, 2, 1> Point_target;
    Point_target << pointCenter[0], pointCenter[1];

    unsigned int index = 0;


    for (unsigned int i = 0; i < indices.size(); i++)
    {
        const auto& coord = x[indices[i]];  // actual position in sofa ! 

        double x_pos = coord[0];
        double y_pos = coord[1];
        double z_pos = coord[2];

        // quaternions 
        // Eigen::Quaterniond q(coord[6], coord[3], coord[4], coord[5]);
        // Eigen::Matrix3d R = q.toRotationMatrix();

        // 2d projection of the 3d position and orientation of the effector
        Eigen::Matrix<double, 2, 1> Point_current = calculateProjectedPoint(
            x_pos, y_pos, z_pos, 
            fLength, 
            pPoint,
            cameraPosition,
            cameraOrientation
        );

        Eigen::Matrix<double, 2, 1> Point_diff = Point_current - Point_target;


        msg_info(this) << "Point_diff: " << Point_diff.transpose();

        for (sofa::Size j = 0; j < 2; j++)
        {
            // Real dfree = Jdx->element(constraintIndex + index) + Point_diff[j] * weight[j];
            Real dfree = Jdx->element(index) + Point_diff[j] * weight[j];
            resV->set(constraintIndex + index, dfree);
            index++;
        }

        // resV->set(constraintIndex + index, 0.0);
        // index++;
    }
}

} 

