
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

#include <SoftRobots.Inverse/component/constraint/MagnetEffector.h>

// ------------------------------------------------------------------------------

#include <Eigen/Dense>
#include <Eigen/Geometry>  // Para cuaterniones

using namespace std;
using namespace Eigen;

// ------------------------------------------------------------------------------


namespace softrobotsinverse::constraint
{

using sofa::helper::ReadAccessor ;
using sofa::helper::WriteAccessor ;
using sofa::core::objectmodel::ComponentState;

template<class DataTypes>
MagnetEffector<DataTypes>::MagnetEffector(MechanicalState* object)
    : Effector<DataTypes>(object)
    , softrobots::constraint::MagnetModel<DataTypes>(object)
    , d_effectorGoal(initData(&d_effectorGoal,"effectorGoal",
                    "Desired positions. \n"
                    "If the size does not match with the size of indices, \n"
                    "one will resize considerering the smallest one."))
{
}

template<class DataTypes>
MagnetEffector<DataTypes>::~MagnetEffector()
{
}

template<class DataTypes>
void MagnetEffector<DataTypes>::init()
{
    softrobots::constraint::MagnetModel<DataTypes>::init();

    if(!d_effectorGoal.isSet())
    {
        msg_warning(this) <<"TargetPosition not defined. Default value assigned  ("<<Coord()<<").";
        setTargetDefaultValue();
    }

    if(d_indices.getValue().size() != d_effectorGoal.getValue().size())
        resizeData();
}

template<class DataTypes>
void MagnetEffector<DataTypes>::setTargetDefaultValue()
{
    WriteAccessor<sofa::Data<VecCoord> > defaultTarget = d_effectorGoal;
    defaultTarget.resize(1);
    defaultTarget[0] = Coord();
}

template<class DataTypes>
void MagnetEffector<DataTypes>::resizeData()
{
    if(d_indices.getValue().size() < d_effectorGoal.getValue().size())
    {
        msg_warning(this)<<"Indices size is lower than target size, some targets will not be considered.";
    }
    else
    {
        msg_warning(this) <<"Indices size is larger than target size. Launch resize process.";
        WriteAccessor<sofa::Data<sofa::type::vector<unsigned int> > > indices = d_indices;
        indices.resize(d_effectorGoal.getValue().size());
    }
}





// --------------------------------------------------------------------------------------------------------------------------------------------







Eigen::Vector3d Calculo_B_z(double x,double y,double z,double mum, double mu_hat_x,double mu_hat_y,double mu_hat_z)
{    
    // Definición de variables
    Vector3d Distancia_r(x, y, z);
    //std::cerr << "Valor de Distancia_r: " << Distancia_r << std::endl;
    double length_r = Distancia_r.norm();  // Magnitud del vector
    //std::cerr << "Valor de length_r: " << length_r << std::endl;
    Vector3d r_hat = Distancia_r / length_r; // Vector unitario
    // cout << "r_hat: " << r_hat.transpose() << endl;
    Vector3d Mu_hat(mu_hat_x, mu_hat_y,mu_hat_z);  // Dirección del momento magnético (vector unitario)
    // cout << "Mu_hat : " << Mu_hat.transpose() << endl;
    // Multiplicación por mu (magnitud)
    Eigen::Vector3d mu = Mu_hat * mum;
    // Producto tensorial r_hat * r_hat^T
    Eigen::Matrix3d AAA = 3 * (r_hat * r_hat.transpose()) - Eigen::Matrix3d::Identity();
    // Multiplicamos AAA por Mu_hat y luego por mu
    Eigen::Vector3d numerador = AAA * mu;
    double denominador = 4 * M_PI * std::pow(std::abs(length_r), 3);
    // Campo magnético B (vector)
    Eigen::Vector3d Campo_Magnetico_resultado = numerador / denominador;
    Campo_Magnetico_resultado *= 1e15;  // Multiplicamos por 10^12 (en unidades apropiadas)
    return Campo_Magnetico_resultado;
}
Eigen::Vector3d Calculo_B_x(double x,double y,double z,double mum, double mu_hat_x,double mu_hat_y,double mu_hat_z)
{    
    // Definición de variables
    Vector3d Distancia_r(x, y, z);
    //std::cerr << "Valor de Distancia_r: " << Distancia_r << std::endl;
    double length_r = Distancia_r.norm();  // Magnitud del vector
    //std::cerr << "Valor de length_r: " << length_r << std::endl;
    Vector3d r_hat = Distancia_r / length_r; // Vector unitario
    // cout << "r_hat: " << r_hat.transpose() << endl;
    Vector3d Mu_hat(mu_hat_x, mu_hat_y,mu_hat_z);  // Dirección del momento magnético (vector unitario)
    // cout << "Mu_hat : " << Mu_hat.transpose() << endl;
    // Multiplicación por mu (magnitud)
    Eigen::Vector3d mu = Mu_hat * mum;
    // Producto tensorial r_hat * r_hat^T
    Eigen::Matrix3d AAA = 3 * (r_hat * r_hat.transpose()) - Eigen::Matrix3d::Identity();
    // Multiplicamos AAA por Mu_hat y luego por mu
    Eigen::Vector3d numerador = AAA * mu;
    double denominador = 4 * M_PI * std::pow(std::abs(length_r), 3);
    // Campo magnético B (vector)
    Eigen::Vector3d Campo_Magnetico_resultado = numerador / denominador;
    Campo_Magnetico_resultado *= 1e15;  // Multiplicamos por 10^12 (en unidades apropiadas)
    return Campo_Magnetico_resultado;
}
Eigen::Vector3d Calculo_B_y(double x,double y,double z,double mum, double mu_hat_x,double mu_hat_y,double mu_hat_z)
{    
    // Definición de variables
    Vector3d Distancia_r(x, y, z);
    //std::cerr << "Valor de Distancia_r: " << Distancia_r << std::endl;
    double length_r = Distancia_r.norm();  // Magnitud del vector
    //std::cerr << "Valor de length_r: " << length_r << std::endl;
    Vector3d r_hat = Distancia_r / length_r; // Vector unitario
    // cout << "r_hat: " << r_hat.transpose() << endl;
    Vector3d Mu_hat(mu_hat_x, mu_hat_y,mu_hat_z);  // Dirección del momento magnético (vector unitario)
    // cout << "Mu_hat : " << Mu_hat.transpose() << endl;
    // Multiplicación por mu (magnitud)
    Eigen::Vector3d mu = Mu_hat * mum;
    // Producto tensorial r_hat * r_hat^T
    Eigen::Matrix3d AAA = 3 * (r_hat * r_hat.transpose()) - Eigen::Matrix3d::Identity();
    // Multiplicamos AAA por Mu_hat y luego por mu
    Eigen::Vector3d numerador = AAA * mu;
    double denominador = 4 * M_PI * std::pow(std::abs(length_r), 3);
    // Campo magnético B (vector)
    Eigen::Vector3d Campo_Magnetico_resultado = numerador / denominador;
    Campo_Magnetico_resultado *= 1e15;  // Multiplicamos por 10^12 (en unidades apropiadas)
    return Campo_Magnetico_resultado;
}


// --------------------------------------------------------------------------------------------------------------------------------------------




template<class DataTypes>
void MagnetEffector<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* cParams,
                                                         sofa::linearalgebra::BaseVector *resV,
                                                         const sofa::linearalgebra::BaseVector *Jdx)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    SOFA_UNUSED(cParams);
    const auto& PosSensor = sofa::helper::getReadAccessor(d_PosSensor);
    const auto& mum = sofa::helper::getReadAccessor(d_mum);
    ReadAccessor<sofa::Data<VecCoord> > x = m_state->readPositions();
    ReadAccessor<sofa::Data<VecCoord> > effectorGoal = d_effectorGoal;



// -----------------------------------------------------------------------------------------------------------------------------

    // std::cout << "Pos_Sensor: " << PosSensor << std::endl;
    double PosSensor_x = PosSensor[0][0] ;
    double PosSensor_y = PosSensor[0][1] ;
    double PosSensor_z = PosSensor[0][2] ;
    double PosIman_x = PosSensor[0][3] ;
    double PosIman_y = PosSensor[0][4] ;
    double PosIman_z = PosSensor[0][5] ;
    // std::cout << "PosSensor_x: " << PosSensor_x << std::endl;
    // std::cout << "PosSensor_y: " << PosSensor_y << std::endl;
    // std::cout << "PosSensor_z: " << PosSensor_z << std::endl;
    // std::cout << "PosIman_x: " << PosIman_x << std::endl;
    // std::cout << "PosIman_y: " << PosIman_y << std::endl;
    // std::cout << "PosIman_z: " << PosIman_z << std::endl;

    // Acceder directamente al Data de 'x'
    auto& data = *x;  // Desreferenciamos el ReadAccessor para obtener el Data
    // Acceder al vector que contiene los valores de las coordenadas (ya que es un vector de Vec<2, double>)
    auto& vec = data;  // 'data' es un vector de Vec<2, double>
    Eigen::Vector3d B_calculada;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_calculada_X;  // Variable global o de ámbito extendido
    Eigen::Vector3d B_calculada_Y;  // Variable global o de ámbito extendido

    for (const auto& coord : vec) {
        // Acceder a todas las componentes de cada Vec<2, double> y mostrar las tres componentes si es posible
        // msg_warning() << "Coord: (" << coord[0] << ", " << coord[1] << ", " << coord[2] << ")";  // Asumiendo que 'Vec<2, double>' tiene tres componentes
        Eigen::Quaterniond MiR(coord[6], coord[3], coord[4], coord[5]);  // (w, x, y, z)
        Eigen::Matrix3d rotation_matrix = MiR.toRotationMatrix();
        // std::cout << "Matriz de rotación:\n" << rotation_matrix << std::endl;
        // Extraer la última columna
        const double mu_x = rotation_matrix(0, 2);  // Elemento (0, 2)
        const double mu_y = rotation_matrix(1, 2);  // Elemento (1, 2)
        const double mu_z = rotation_matrix(2, 2);  // Elemento (2, 2)

        // std::cout << "mu_x: " << mu_x << std::endl;
        // std::cout << "mu_y: " << mu_y << std::endl;
        // std::cout << "mu_z: " << mu_z << std::endl;

        const double ajuste_x = PosSensor_x ;
        const double ajuste_y = PosSensor_y;
        const double ajuste_z = PosSensor_z;
//        const double ajuste_z = 2.7 - 15; // Este ajuste era para corroborar calculo en c++ con el de python
        // std::cout << "coord : " << coord << std::endl;
        B_calculada = Calculo_B_z(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2]-ajuste_z,mum[0][2],mu_x,mu_y,mu_z);
        B_calculada_X = Calculo_B_x(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2]-ajuste_z,mum[0][0],mu_x,mu_y,mu_z);
        B_calculada_Y = Calculo_B_y(coord[0]- ajuste_x,coord[1]- ajuste_y,coord[2]-ajuste_z,mum[0][1],mu_x,mu_y,mu_z);

        // std::cout << "mum Antes de campo magneteffector.inl: " << mum[0][0] << std::endl;
        std::cout << "Campo magnético B_calculado c++ X: " << B_calculada_X[0] << std::endl;
        std::cout << "Campo magnético B_calculado c++ Y: " << B_calculada_Y[1] << std::endl;
        std::cout << "Campo magnético B_calculado c++ Z: " << B_calculada[2] << std::endl;
    }


    double B_calculada_x = B_calculada_X[0]; // Accede al primer elemento 
    double B_calculada_y = B_calculada_Y[1]; // Accede al segundo elemento
    double B_calculada_z = B_calculada[2]; // Accede al tercer elemento 
    // std::cout << "B_calculada_x: " << B_calculada_x << std::endl;
    // std::cout << "B_calculada_y: " << B_calculada_y << std::endl;
    // std::cout << "B_calculada_z: " << B_calculada_z << std::endl;

// -----------------------------------------------------------------------------------------------------------------------------







    const auto& useDirections = sofa::helper::getReadAccessor(d_useDirections);
    const auto& directions = sofa::helper::getReadAccessor(d_directions);
    const auto& Jacobian = sofa::helper::getReadAccessor(d_Jacobian);
    const auto& weight = sofa::helper::getReadAccessor(d_weight);
    const auto& indices = sofa::helper::getReadAccessor(d_indices);
    sofa::Index sizeIndices = indices.size();
    const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);

    int index = 0;
    for (unsigned int i=0; i<sizeIndices; i++)
    {
        Coord pos = x[indices[i]]; //con pos calculé B_calculado
        Coord goalPos = getTarget(effectorGoal[i],pos);

        // double PosSensor_y =
        // double PosSensor_z =

        // std::cout << "Posicion Iman Sofa: " << pos << std::endl; //con pos calculé B_calculado
        // std::cout << "B_Calculada: "< << B_calculada.transpose() << std::endl;

        double B_Goal_x = effectorGoal[i][0]; // Accede al primer elemento
        double B_Goal_y = effectorGoal[i][1]; // Accede al segundo elemento
        double B_Goal_z = effectorGoal[i][2]; // Accede al tercer elemento
        // std::cout << "GoalPosx: " << B_Goal_x << std::endl;
        // std::cout << "GoalPosy: " << B_Goal_y << std::endl;
        // std::cout << "GoalPosz: " << B_Goal_z << std::endl;

        double B_diff_x = B_Goal_x - B_calculada_x; 
        double B_diff_y = B_Goal_y - B_calculada_y;    
        double B_diff_z = B_Goal_z - B_calculada_z; 
        // std::cout << "B_diff_x: " << B_diff_x << std::endl;
        // std::cout << "B_diff_y: " << B_diff_y << std::endl;
        // std::cout << "B_diff_z: " << B_diff_z << std::endl;


        Eigen::Vector3d vec(-B_diff_x, -B_diff_y, -B_diff_z);
        std::cout << "vector diferencia: " << vec << std::endl;
        pos[0] = B_calculada_x;  // Asignar un nuevo valor
        pos[1] = B_calculada_y;  // Asignar un nuevo valor
        pos[2] = B_calculada_z;  // Extraer el valor escalar
        pos[3] = 0;
        pos[4] = 0;
        pos[5] = 0;  
        pos[6] = 1;  
//        std::cout << "Pos2: " << pos << std::endl; //con pos calculé B_calculado
        // std::cout << "goalPos: " << goalPos << std::endl;

        Deriv d = DataTypes::coordDifference(pos,goalPos);
//        std::cout << "d: " << d << std::endl; //con pos calculé B_calculado
        // calcular diferencia entre B calculado sobre x y B medido en el sensor

        for(sofa::Size j=0; j<3; j++)
            {
                // Real dfree = Jdx->element(index) + d*directions[j]*weight[j];
                // std::cout << "Jacobian: " << Jacobian << std::endl; //con pos calculé B_calculado
                // std::cout << "Jdx->element " << Jdx->element(index) << std::endl;
                Real dfree = Jdx->element(index) + vec[j];
                // std::cout << "dfree " << dfree << std::endl;
                resV->set(constraintIndex+index, dfree);
                index++;
            }
    }
}




} // namespace
