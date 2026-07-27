
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

template<class DataTypes>
CameraProjectionEffector<DataTypes>::CameraProjectionEffector(MechanicalState* object)
    : Effector<DataTypes>(object)
    , softrobots::constraint::CameraProjectionModel<DataTypes>(object) // Llama a la clase base
    , d_effectorGoal(initData(&d_effectorGoal, "effectorGoal",
                    "Desired 3D positions or target for the effector."))
    , d_ellipseParameters(initData(&d_ellipseParameters, "ellipseParameters",
                    "Target ellipse parameters detected in image [cx, cy, a, b, angle]."))
{
}

template<class DataTypes>
CameraProjectionEffector<DataTypes>::~CameraProjectionEffector()
{
}

// ------------------------------------------------------------------------------
// Métodos de inicialización
// ------------------------------------------------------------------------------

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

// ------------------------------------------------------------------------------
// Función Auxiliar: Desproyección (2D Ellipse -> 3D Position)
// ------------------------------------------------------------------------------

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
    const sofa::type::Vec2d& principalPoint)
{
    // 1. Centro de la elipse en píxeles (Proyección perspectiva Pinhole)
    double u = focalLength[0] * (x / z) + principalPoint[0];
    double v = focalLength[1] * (y / z) + principalPoint[1];

    // 2. Extraer la normal del disco 3D (tercera columna de la matriz de rotación)
    Eigen::Vector3d normal = R.col(2); 

    // inclination
    double cos_tilt = std::abs(normal(2));
    if (cos_tilt < 1e-3) cos_tilt = 1e-3; // Evitar división por cero si está de canto

    // Semiejes
    double semi_a = focalLength[0] * (radius / z); 
    double semi_b = semi_a * cos_tilt;           

    // 5. Angle
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
            pPoint
        );

        // 4. Calcular la diferencia (error) entre lo proyectado y lo deseado
        Eigen::Matrix<double, 5, 1> E_diff = E_current - E_target;

        // 5. Asignar el residuo al vector resV (5 dimensiones por efector)
        for (sofa::Size j = 0; j < 5; j++)
        {
            Real dfree = Jdx->element(constraintIndex + index) + E_diff[j] * weight[j];
            resV->set(constraintIndex + index, dfree);
            index++;
        }
    }
}
} // namespace softrobotsinverse::constraint


// template<class DataTypes>
// void CameraProjectionEffector<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* cParams,
//                                                          sofa::linearalgebra::BaseVector *resV,
//                                                          const sofa::linearalgebra::BaseVector *Jdx)
// {
//     if(d_componentState.getValue() != ComponentState::Valid)
//         return;

//     SOFA_UNUSED(cParams);
//     // const auto& PosSensor = sofa::helper::getReadAccessor(d_PosSensor);
//     // const auto& mum = sofa::helper::getReadAccessor(d_mum);
//     const auto& EllipseParameters = sofa::helper::getReadAccessor(d_EllipseParameters);
//     const auto& CameraParameters = sofa::helper::getReadAccessor(d_CameraParameters);

    
//     ReadAccessor<sofa::Data<VecCoord> > x = m_state->readPositions();
//     ReadAccessor<sofa::Data<VecCoord> > effectorGoal = d_effectorGoal;



// // -----------------------------------------------------------------------------------------------------------------------------

//     // std::cout << "Pos_Sensor: " << PosSensor << std::endl;
//     double PosSensor_x = PosSensor[0][0] ;
//     double PosSensor_y = PosSensor[0][1] ;
//     double PosSensor_z = PosSensor[0][2] ;
//     double PosIman_x = PosSensor[0][3] ;
//     double PosIman_y = PosSensor[0][4] ;
//     double PosIman_z = PosSensor[0][5] ;
//     // std::cout << "PosSensor_x: " << PosSensor_x << std::endl;
//     // std::cout << "PosSensor_y: " << PosSensor_y << std::endl;
//     // std::cout << "PosSensor_z: " << PosSensor_z << std::endl;
//     // std::cout << "PosIman_x: " << PosIman_x << std::endl;
//     // std::cout << "PosIman_y: " << PosIman_y << std::endl;
//     // std::cout << "PosIman_z: " << PosIman_z << std::endl;

//     // Acceder directamente al Data de 'x'
//     auto& data = *x;  // Desreferenciamos el ReadAccessor para obtener el Data
//     // Acceder al vector que contiene los valores de las coordenadas (ya que es un vector de Vec<2, double>)
//     auto& vec = data;  // 'data' es un vector de Vec<2, double>
//     Eigen::Vector3d B_calculada;  // Variable global o de ámbito extendido
//     Eigen::Vector3d B_calculada_X;  // Variable global o de ámbito extendido
//     Eigen::Vector3d B_calculada_Y;  // Variable global o de ámbito extendido

//     for (const auto& coord : vec) {
//         // Acceder a todas las componentes de cada Vec<2, double> y mostrar las tres componentes si es posible
//         // msg_warning() << "Coord: (" << coord[0] << ", " << coord[1] << ", " << coord[2] << ")";  // Asumiendo que 'Vec<2, double>' tiene tres componentes
//         Eigen::Quaterniond MiR(coord[6], coord[3], coord[4], coord[5]);  // (w, x, y, z)
//         Eigen::Matrix3d rotation_matrix = MiR.toRotationMatrix();
//         // std::cout << "Matriz de rotación:\n" << rotation_matrix << std::endl;
//         // Extraer la última columna
//         const double mu_x = rotation_matrix(0, 2);  // Elemento (0, 2)
//         const double mu_y = rotation_matrix(1, 2);  // Elemento (1, 2)
//         const double mu_z = rotation_matrix(2, 2);  // Elemento (2, 2)

//         // std::cout << "mu_x: " << mu_x << std::endl;
//         // std::cout << "mu_y: " << mu_y << std::endl;
//         // std::cout << "mu_z: " << mu_z << std::endl;

//         const double ajuste_x = PosSensor_x ;
//         const double ajuste_y = PosSensor_y;
//         const double ajuste_z = PosSensor_z;
// //        const double ajuste_z = 2.7 - 15; // Este ajuste era para corroborar calculo en c++ con el de python
//         // std::cout << "coord : " << coord << std::endl;
//         B_calculada = Calculo_B_z(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2]-ajuste_z,mum[0][2],mu_x,mu_y,mu_z);
//         B_calculada_X = Calculo_B_x(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2]-ajuste_z,mum[0][0],mu_x,mu_y,mu_z);
//         B_calculada_Y = Calculo_B_y(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2]-ajuste_z,mum[0][1],mu_x,mu_y,mu_z);

//         // std::cout << "mum Antes de campo CameraProjectionEffector.inl: " << mum[0][0] << std::endl;
//         std::cout << "Campo magnético B_calculado c++ X: " << B_calculada_X[0] << std::endl;
//         std::cout << "Campo magnético B_calculado c++ Y: " << B_calculada_Y[1] << std::endl;
//         std::cout << "Campo magnético B_calculado c++ Z: " << B_calculada[2] << std::endl;
//     }


//     double B_calculada_x = B_calculada_X[0]; // Accede al primer elemento 
//     double B_calculada_y = B_calculada_Y[1]; // Accede al segundo elemento
//     double B_calculada_z = B_calculada[2]; // Accede al tercer elemento 
//     // std::cout << "B_calculada_x: " << B_calculada_x << std::endl;
//     // std::cout << "B_calculada_y: " << B_calculada_y << std::endl;
//     // std::cout << "B_calculada_z: " << B_calculada_z << std::endl;

// // -----------------------------------------------------------------------------------------------------------------------------


//     const auto& useDirections = sofa::helper::getReadAccessor(d_useDirections);
//     const auto& directions = sofa::helper::getReadAccessor(d_directions);
//     const auto& Jacobian = sofa::helper::getReadAccessor(d_Jacobian);
//     const auto& weight = sofa::helper::getReadAccessor(d_weight);
//     const auto& indices = sofa::helper::getReadAccessor(d_indices);
//     sofa::Index sizeIndices = indices.size();
//     const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);

//     int index = 0;
//     for (unsigned int i=0; i<sizeIndices; i++)
//     {
//         Coord pos = x[indices[i]]; //con pos calculé B_calculado
//         Coord goalPos = getTarget(effectorGoal[i],pos);

//         // double PosSensor_y =
//         // double PosSensor_z =

//         // std::cout << "Posicion Iman Sofa: " << pos << std::endl; //con pos calculé B_calculado
//         // std::cout << "B_Calculada: "< << B_calculada.transpose() << std::endl;

//         double B_Goal_x = effectorGoal[i][0]; // Accede al primer elemento
//         double B_Goal_y = effectorGoal[i][1]; // Accede al segundo elemento
//         double B_Goal_z = effectorGoal[i][2]; // Accede al tercer elemento
//         // std::cout << "GoalPosx: " << B_Goal_x << std::endl;
//         // std::cout << "GoalPosy: " << B_Goal_y << std::endl;
//         // std::cout << "GoalPosz: " << B_Goal_z << std::endl;

//         double B_diff_x = B_Goal_x - B_calculada_x; 
//         double B_diff_y = B_Goal_y - B_calculada_y;    
//         double B_diff_z = B_Goal_z - B_calculada_z; 
//         // std::cout << "B_diff_x: " << B_diff_x << std::endl;
//         // std::cout << "B_diff_y: " << B_diff_y << std::endl;
//         // std::cout << "B_diff_z: " << B_diff_z << std::endl;


//         Eigen::Vector3d vec(-B_diff_x, -B_diff_y, -B_diff_z);
//         std::cout << "vector diferencia: " << vec << std::endl;
//         pos[0] = B_calculada_x;  // Asignar un nuevo valor
//         pos[1] = B_calculada_y;  // Asignar un nuevo valor
//         pos[2] = B_calculada_z;  // Extraer el valor escalar
//         pos[3] = 0;
//         pos[4] = 0;
//         pos[5] = 0;  
//         pos[6] = 1;  
// //        std::cout << "Pos2: " << pos << std::endl; //con pos calculé B_calculado
//         // std::cout << "goalPos: " << goalPos << std::endl;

//         Deriv d = DataTypes::coordDifference(pos,goalPos);
// //        std::cout << "d: " << d << std::endl; //con pos calculé B_calculado
//         // calcular diferencia entre B calculado sobre x y B medido en el sensor

//         for(sofa::Size j=0; j<3; j++)
//             {
//                 // Real dfree = Jdx->element(index) + d*directions[j]*weight[j];
//                 // std::cout << "Jacobian: " << Jacobian << std::endl; //con pos calculé B_calculado
//                 // std::cout << "Jdx->element " << Jdx->element(index) << std::endl;
//                 Real dfree = Jdx->element(index) + vec[j];
//                 // std::cout << "dfree " << dfree << std::endl;
//                 resV->set(constraintIndex+index, dfree);
//                 index++;
//             }
//     }
// }




// } // namespace
