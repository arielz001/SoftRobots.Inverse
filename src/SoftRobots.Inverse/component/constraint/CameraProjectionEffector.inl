
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

// ------------------------------------------------------------------------------

#include <Eigen/Dense>
#include <Eigen/Geometry>  // To use Quaternions

using namespace std;
using namespace Eigen;

// ------------------------------------------------------------------------------


namespace softrobotsinverse::constraint
{

// Effector for SoftRobots.Inverse that projects 3D poses (Rigid3d) 
//    into 2D image ellipse space [cx, cy, major, minor, angle].


using sofa::helper::ReadAccessor ;
using sofa::helper::WriteAccessor ;
using sofa::core::objectmodel::ComponentState;


//===============================================================
// [-] ORIGINAL FILE
// template<class DataTypes>
// CameraProjectionEffector<DataTypes>::CameraProjectionEffector(MechanicalState* object)
//     : Effector<DataTypes>(object)
//     , softrobots::constraint::CameraProjectionModel<DataTypes>(object)
//     , d_effectorGoal(initData(&d_effectorGoal,"effectorGoal",
//                     "Desired positions. \n"
//                     "If the size does not match with the size of indices, \n"
//                     "one will resize considerering the smallest one."))
//===============================================================


//===============================================================
// [+] MODIFICATION
// this is for get the real camera parameters, you can change it your program.
// Initializes the effector with camera intrinsic parameters and default ellipse target data.
template<class DataTypes>
CameraProjectionEffector<DataTypes>::CameraProjectionEffector(MechanicalState* object)
    : Effector<DataTypes>(object)
    , d_effectorGoal(initData(&d_effectorGoal, "effectorGoal",
                    "Desired ellipse parameters in image space: [cx, cy, major, minor, angle]."))
    , d_focalLength(initData(&d_focalLength, sofa::type::Vec2(600.0, 600.0), "focalLength", 
                    "Focal length [fx, fy] in pixels."))
    , d_principalPoint(initData(&d_principalPoint, sofa::type::Vec2(320.0, 240.0), "principalPoint", 
                    "Principal point [u0, v0] in pixels."))
    , d_diskRadius(initData(&d_ellipseRadius, 0.05, "ellipseRadius", 
                    "Real 3D radius of the target ellipse in meters."))
{}
//===============================================================

// Destructor
template<class DataTypes>
CameraProjectionEffector<DataTypes>::~CameraProjectionEffector()
{}


//===============================================================
//[-] ORIGINAL FILE
// template<class DataTypes>
// void CameraProjectionEffector<DataTypes>::init()
// {
//     softrobots::constraint::CameraProjectionModel<DataTypes>::init();

//     if(!d_effectorGoal.isSet())
//     {
//         msg_warning(this) <<"TargetPosition not defined. Default value assigned  ("<<Coord()<<").";
//         setTargetDefaultValue();
//     }

//     if(d_indices.getValue().size() != d_effectorGoal.getValue().size())
//         resizeData();
// }
//===============================================================



//===============================================================
// [+] MODIFICATION
// Calls base class initialization and validates data array sizes.
template<class DataTypes>
void CameraProjectionEffector<DataTypes>::init()
{
    // initialize the softrobots.inverse class
    Effector<DataTypes>::init();

    if(!d_effectorGoal.isSet())
    {
        msg_warning(this) << "Target position not defined. Default value assigned.";
        setTargetDefaultValue();
    }

    if(d_indices.getValue().size() != d_effectorGoal.getValue().size())
        resizeData();
}

//===============================================================
// [+] ORIGINAL FILE
// Sets default values for the target ellipse array when not explicitly defined.
template<class DataTypes>
void CameraProjectionEffector<DataTypes>::setTargetDefaultValue()
{
    WriteAccessor<sofa::Data<VecCoord> > defaultTarget = d_effectorGoal;
    defaultTarget.resize(1);
    defaultTarget[0] = Coord();
}

//===============================================================
// [+] ORIGINAL FILE
// Resizes the target ellipse array to match the size of the indices array.
template<class DataTypes>
void CameraProjectionEffector<DataTypes>::resizeData()
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
//===============================================================


// ________________________________________________________________________________
// ________________________________________________________________________________


/**
 * @brief Projects a 3D rigid pose into a 2D image ellipse [cx, cy, major, minor, angle].
 * @param pose3D The current 3D position and orientation (quaternion) from SOFA.
 * @return Eigen::Vector5d Vector containing [cx, cy, major, minor, angle] to be projected in the image space.
 */
template<class DataTypes>
Eigen::Vector5d CameraProjectionEffector<DataTypes>::computeProjectedEllipse(const Coord& pose3D)
{
    // 1. Extract 3D position and orientation quaternion (w, x, y, z)
    Eigen::Vector3d center3D(pose3D[0], pose3D[1], pose3D[2]);
    Eigen::Quaterniond orientation(pose3D[6], pose3D[3], pose3D[4], pose3D[5]);
    Eigen::Matrix3d R = orientation.toRotationMatrix();

    // 2. Retrieve camera intrinsic parameters and 3D target dimensions configured in SOFA
    const auto fx = d_focalLength.getValue()[0];
    const auto fy = d_focalLength.getValue()[1];
    const auto u0 = d_principalPoint.getValue()[0];
    const auto v0 = d_principalPoint.getValue()[1];
    const double radius = d_ellipseRadius.getValue(); // 3D radius of the target ellipse

    // 3. Define local 3D ellipse boundary vectors
    Eigen::Vector3d localEdgeU(radius, 0.0, 0.0);
    Eigen::Vector3d localEdgeV(0.0, radius, 0.0);

    // 4. Transform local ellipse points to camera 3D space
    Eigen::Vector3d centerCam = center3D;
    Eigen::Vector3d edgeUCam  = R * localEdgeU + center3D;
    Eigen::Vector3d edgeVCam  = R * localEdgeV + center3D;

    // Prevent division by zero or negative depth
    double Zc = std::max(centerCam.z(), 1e-4);
    double Zu = std::max(edgeUCam.z(), 1e-4);
    double Zv = std::max(edgeVCam.z(), 1e-4);

    // 5. Project 3D points to 2D image coordinates (Pinhole model)
    Eigen::Vector2d p_center(fx * (centerCam.x() / Zc) + u0, fy * (centerCam.y() / Zc) + v0);
    Eigen::Vector2d p_u(fx * (edgeUCam.x() / Zu) + u0,       fy * (edgeUCam.y() / Zu) + v0);
    Eigen::Vector2d p_v(fx * (edgeVCam.x() / Zv) + u0,       fy * (edgeVCam.y() / Zv) + v0);

    // 6. Compute 2D projected ellipse parameters
    Eigen::Vector2d u = p_u - p_center;
    Eigen::Vector2d v = p_v - p_center;

    double cx = p_center.x();
    double cy = p_center.y();
    double major = 2.0 * u.norm();
    double minor = 2.0 * v.norm();
    double angle = std::atan2(u.y(), u.x()) * (180.0 / M_PI); // Angle in degrees

    return Eigen::Vector5d(cx, cy, major, minor, angle);
}



// ________________________________________________________________________________
// ________________________________________________________________________________


template<class DataTypes>
void CameraProjectionEffector<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* cParams,
                                                                 sofa::linearalgebra::BaseVector *resV,
                                                                 const sofa::linearalgebra::BaseVector *Jdx)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return;

    SOFA_UNUSED(cParams);

    ReadAccessor<sofa::Data<VecCoord>> x = m_state->readPositions();
    ReadAccessor<sofa::Data<VecCoord>> effectorGoal = d_effectorGoal;

    const auto& indices = sofa::helper::getReadAccessor(d_indices);
    sofa::Index sizeIndices = indices.size();
    const auto& constraintIndex = sofa::helper::getReadAccessor(d_constraintIndex);

    int index = 0;
    for (unsigned int i = 0; i < sizeIndices; i++)
    {
        Coord currentPose3D = x[indices[i]];

        // projection current 3D pose to simulated 2D image ellipse [cx, cy, major, minor, angle]
        Eigen::Vector5d calculatedEllipse = computeProjectedEllipse(currentPose3D);

        // Goal defined by effectorGoal
        Eigen::Vector5d goalEllipse(effectorGoal[i][0],
                                   effectorGoal[i][1],
                                   effectorGoal[i][2],
                                   effectorGoal[i][3],
                                   effectorGoal[i][4]);

        // Error = Goal - Current
        Eigen::Vector5d ellipseError = goalEllipse - calculatedEllipse;

        // Assign constraint violations (5 dimensions per node) to the SOFA residual vector
        for(sofa::Size j = 0; j < 5; j++)
        {
            Real dfree = Jdx->element(index) - ellipseError[j];
            resV->set(constraintIndex + index, dfree);
            index++;
        }
    }
}


} // namespace