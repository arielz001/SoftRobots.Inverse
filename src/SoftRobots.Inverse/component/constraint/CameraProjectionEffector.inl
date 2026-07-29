
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

// ------------------------------------------------------------------------------
// Constructor & Destructor
// ------------------------------------------------------------------------------

// template<class DataTypes>
// CameraProjectionEffector<DataTypes>::CameraProjectionEffector(MechanicalState* object)
//     : Effector<DataTypes>(object)
//     , softrobots::constraint::CameraProjectionModel<DataTypes>(object) // Llama a la clase base
//     , d_effectorGoal(initData(&d_effectorGoal, "effectorGoal",
//                     "Desired 3D positions or target for the effector."))
//     , d_ellipseParameters(initData(&d_ellipseParameters, "ellipseParameters",
//                     "Target ellipse parameters detected in image [cx, cy, a, b, angle]."))
// {
// }

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

    // Cálculo de la profundidad Z
    double Z = (x_focal_length * real_radius) / a;

    // Desproyección del centro
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
    double camera_z = cameraPos[2]; // Posición Z de la cámara

    // 1. Calcular la distancia REAL desde la lente de la cámara hasta el disco
    double z_rel = z - camera_z; // 115.0 - (-13.2) = 128.2 mm
    
    double safe_z = std::abs(z_rel);
    if (safe_z < 10.0) safe_z = 10.0; 

    // 2. Proyección usando z_rel
    double u = focalLength[0] * (x / safe_z) + principalPoint[0];
    double v = focalLength[1] * (y / safe_z) + principalPoint[1];

    // 3. Normal e inclinación
    Eigen::Vector3d normal = R.col(2); 
    double cos_tilt = std::abs(normal(2));
    if (cos_tilt < 1e-3) cos_tilt = 1e-3;

    // 4. Semiejes calculados con la profundidad relativa correcta
    double semi_a = focalLength[0] * (radius / safe_z);  
    double semi_b = semi_a * cos_tilt;           

    // 5. Ángulo
    double alpha = std::atan2(normal(1), normal(0));

    Eigen::Matrix<double, 5, 1> ellipse;
    ellipse << u, v, semi_a, semi_b, alpha;
    return ellipse;
}

// ------------------------------------------------------------------------------
// Método principal: Cálculo de la Violación de Restricciones (Solver Inverso)
// ------------------------------------------------------------------------------
template<class DataTypes>
void CameraProjectionEffector<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* cParams,
                                                                 sofa::linearalgebra::BaseVector *resV,
                                                                 const sofa::linearalgebra::BaseVector *Jdx)
{
    if (d_componentState.getValue() != ComponentState::Valid)
        return;

    SOFA_UNUSED(cParams);

    // 1. Acceso a las posiciones actuales de SOFA y parámetros
    ReadAccessor<sofa::Data<VecCoord>> x = m_state->readPositions();
    
    const auto& ellipseParams   = sofa::helper::getReadAccessor(d_ellipseParameters); 
    const double realRadius     = d_radiusEllipse.getValue();
    const auto& cameraPosition = sofa::helper::getReadAccessor(d_cameraPosition);
    const auto& weight          = sofa::helper::getReadAccessor(d_weight);
    const auto& indices         = sofa::helper::getReadAccessor(d_indices);
    const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);

    // Extraer calibración de la cámara directamente desde las variables Data
    sofa::type::Vec2d fLength = d_focalLength.getValue();
    sofa::type::Vec2d pPoint  = d_principalPoint.getValue();

    // Elipse objetivo [cx, cy, a, b, angle]
    Eigen::Matrix<double, 5, 1> E_target;
    E_target << ellipseParams[0], ellipseParams[1], ellipseParams[2], ellipseParams[3], ellipseParams[4];

    unsigned int index = 0;

    // 2. Iterar sobre cada nodo/efector configurado
    for (unsigned int i = 0; i < indices.size(); i++)
    {
        const auto& coord = x[indices[i]]; // Posición 3D y orientación actual en SOFA

        double x_pos = coord[0];
        double y_pos = coord[1];
        double z_pos = coord[2];

        // Extraer orientación (Cuaternión w, x, y, z)
        Eigen::Quaterniond q(coord[6], coord[3], coord[4], coord[5]);
        Eigen::Matrix3d R = q.toRotationMatrix();

        // 3. Proyectar la elipse 2D según la postura actual del robot en 3D
        Eigen::Matrix<double, 5, 1> E_current = calculateProjectedEllipse(
            x_pos, y_pos, z_pos, 
            R, 
            realRadius, 
            fLength, 
            pPoint,
            cameraPosition
        );

        // 4. Calcular la diferencia (error) entre lo proyectado y lo deseado
        Eigen::Matrix<double, 5, 1> E_diff = E_current - E_target;
        // std::cout << "E_diff: " << E_diff << std::endl;
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

} // namespace softrobotsinverse::constraint

