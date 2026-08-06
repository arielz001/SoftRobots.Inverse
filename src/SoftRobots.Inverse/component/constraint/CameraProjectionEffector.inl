
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

#include <SoftRobots.Inverse/component/constraint/CameraProjectionEffector.h>
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
CameraProjectionEffector<DataTypes>::CameraProjectionEffector(MechanicalState* object)
    : Effector<DataTypes>(object)
    , softrobots::constraint::CameraProjectionModel<DataTypes>(object)
    , d_effectorGoal(initData(&d_effectorGoal, "effectorGoal",
                    "Desired 3D positions or target for the effector."))
    , d_ellipseParameters(initData(&d_ellipseParameters, "ellipseParameters",
                    "Target ellipse parameters detected in image [cx, cy, a, b, angle]."))
    , d_cameraPosition(initData(&d_cameraPosition, sofa::type::Vec3d(0.0, 0.0, -13.2), "cameraPosition",
                    "Position of the camera in world coordinates [x, y, z]."))
{
}


template<class DataTypes>
CameraProjectionEffector<DataTypes>::~CameraProjectionEffector()
{
}


template<class DataTypes>
void CameraProjectionEffector<DataTypes>::init()
{
    softrobots::constraint::CameraProjectionModel<DataTypes>::init();

    if(!d_effectorGoal.isSet())
    {
        msg_warning(this) << "TargetPosition not defined. Default value assigned (" << Coord() << ").";
        setTargetDefaultValue();
    }

    if(d_indices.getValue().size() != d_effectorGoal.getValue().size())
        resizeData();
}

template<class DataTypes>
void CameraProjectionEffector<DataTypes>::setTargetDefaultValue()
{
    WriteAccessor<sofa::Data<VecCoord>> defaultTarget = d_effectorGoal;
    defaultTarget.resize(1);
    defaultTarget[0] = Coord();
}

template<class DataTypes>
void CameraProjectionEffector<DataTypes>::resizeData()
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


Eigen::Vector3d ProjectionPosition(double cx_ellipse, double cy_ellipse, 
                                   double a_ellipse, double b_ellipse, double angle_ellipse,
                                   double cx_camera, double cy_camera, 
                                   double x_focal_length, double y_focal_length,
                                   double real_radius)
{
    if (a_ellipse <= 1e-6)
    {
        return Eigen::Vector3d::Zero();
    }

    double a = a_ellipse;

    // axis z distnace
    double Z = (x_focal_length * real_radius) / a;

    // axis (x, y) distance
    double X = (cx_ellipse - cx_camera) * Z / x_focal_length;
    double Y = (cy_ellipse - cy_camera) * Z / y_focal_length;

    SOFA_UNUSED(b_ellipse);
    SOFA_UNUSED(angle_ellipse);

    return Eigen::Vector3d(X, Y, Z);
}


Eigen::Matrix<double, 5, 1> calculateProjectedEllipse(
    double x, double y, double z,
    const Eigen::Matrix3d& R,
    double radius,
    const sofa::type::Vec2d& focalLength,
    const sofa::type::Vec2d& principalPoint,
    const sofa::type::Vec3d& cameraPos)
{
    // coords relative to camera postiicon given in sofa
    double x_rel = x - cameraPos[0];
    double y_rel = y - cameraPos[1];
    double z_rel = z - cameraPos[2];

    // projection to 2d
    double u = focalLength[0] * (x_rel / z_rel) + principalPoint[0];
    double v = focalLength[1] * (y_rel / z_rel) + principalPoint[1]; 

    // normal vector of the effector in camera coordinates
    Eigen::Vector3d normal = R.col(2); 
    double cos_tilt = std::abs(normal(2));
    if (cos_tilt < 1e-3) cos_tilt = 1e-3;

    // semi-axes of the projected ellipse
    double semi_a = focalLength[0] * (2*radius / z_rel);  
    double semi_b = semi_a * cos_tilt;           

    // angle of the semi-major axis (+90° to align with the normal vector)
    double alpha_rad = std::atan2(-normal(1), normal(0)) + (M_PI / 2.0);
    double alpha_deg = alpha_rad * (180.0 / M_PI);

    // normalization of the angle to the range [0, 180)
    while (alpha_deg < 0.0) alpha_deg += 180.0;
    while (alpha_deg >= 180.0) alpha_deg -= 180.0;

    Eigen::Matrix<double, 5, 1> ellipse;
    ellipse << u, v, semi_a, semi_b, alpha_deg;

    std::cout << "Simulated Projection: " << ellipse.transpose() << std::endl;
    // std::cout << "relative position: " << x_rel << ", " << y_rel << ", " << z_rel << std::endl;
    // std::cout << "camera position: " << cameraPos[0] << ", " << cameraPos[1] << ", " << cameraPos[2] << std::endl;
    // std::cout << "xyz position: " << x << ", " << y << ", " << z  << std::endl << '\n\n';


    return ellipse;
}

template<class DataTypes>
void CameraProjectionEffector<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* cParams,
                                                                 sofa::linearalgebra::BaseVector *resV,
                                                                 const sofa::linearalgebra::BaseVector *Jdx)
{
    if (d_componentState.getValue() != ComponentState::Valid)
        return;

    SOFA_UNUSED(cParams);

    // with this we can acces to the parameters
    ReadAccessor<sofa::Data<VecCoord>> x = m_state->readPositions();
    
    const auto& ellipseParams       = sofa::helper::getReadAccessor(d_ellipseParameters); 
    const double realRadius         = d_radiusEllipse.getValue();
    const sofa::type::Vec3d cameraPosition = d_cameraPosition.getValue(); // acess to camera position
    const auto& weight              = sofa::helper::getReadAccessor(d_weight);
    const auto& indices             = sofa::helper::getReadAccessor(d_indices);
    const auto& constraintIndex     = sofa::helper::getReadAccessor(d_constraintIndex);

    // calib parameters
    const sofa::type::Vec2d fLength = d_focalLength.getValue();
    const sofa::type::Vec2d pPoint  = d_principalPoint.getValue();

    // this is the objective to follow 
    // in this case we detect with other code and our camera the ellipse parameters 
    Eigen::Matrix<double, 5, 1> E_target;
    E_target << ellipseParams[0], ellipseParams[1], ellipseParams[2], ellipseParams[3], ellipseParams[4];

    unsigned int index = 0;


    for (unsigned int i = 0; i < indices.size(); i++)
    {
        const auto& coord = x[indices[i]];  // actual position in sofa ! 

        double x_pos = coord[0];
        double y_pos = coord[1];
        double z_pos = coord[2];

        // quaternions 
        Eigen::Quaterniond q(coord[6], coord[3], coord[4], coord[5]);
        Eigen::Matrix3d R = q.toRotationMatrix();

        // 2d projection of the 3d position and orientation of the effector
        Eigen::Matrix<double, 5, 1> E_current = calculateProjectedEllipse(
            x_pos, y_pos, z_pos, 
            R, 
            realRadius, 
            fLength, 
            pPoint,
            cameraPosition
        );

        Eigen::Matrix<double, 5, 1> E_diff = E_current - E_target;

        // normalization
        while (E_diff[4] > 90.0)  E_diff[4] -= 180.0;
        while (E_diff[4] < -90.0) E_diff[4] += 180.0;

        // this is the projection error 
        msg_info(this) << "E_diff: " << E_diff.transpose();

        for (sofa::Size j = 0; j < 5; j++)
        {
            Real dfree = Jdx->element(constraintIndex + index) + E_diff[j] * weight[j];
            resV->set(constraintIndex + index, dfree);
            index++;
        }

        resV->set(constraintIndex + index, 0.0);
        index++;
    }
}

} 

