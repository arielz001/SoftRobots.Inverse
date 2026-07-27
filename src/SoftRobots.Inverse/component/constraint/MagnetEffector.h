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

#include <SoftRobots/component/constraint/model/MagnetModel.h>
#include <SoftRobots.Inverse/component/behavior/Effector.h>

#include <SoftRobots.Inverse/component/config.h>

namespace softrobotsinverse::constraint
{
using softrobotsinverse::behavior::Effector ;

/**
 * The "MagnetEffector" component is used to constrain one or several points of a model
 * to reach desired positions, by acting on chosen actuator(s).
 * Description can be found at:
 * https://softrobotscomponents.readthedocs.io
*/
template< class DataTypes >
class MagnetEffector : public Effector<DataTypes>, public softrobots::constraint::MagnetModel<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(MagnetEffector,DataTypes), SOFA_TEMPLATE(Effector,DataTypes));

    typedef typename sofa::core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename DataTypes::VecCoord            VecCoord;
    typedef typename DataTypes::Coord               Coord;
    typedef typename DataTypes::Deriv               Deriv;
    typedef typename DataTypes::Real                Real;

public:
    MagnetEffector(MechanicalState* object = nullptr);
    ~MagnetEffector() override;

    /////////////// Inherited from Effector ////////////

    void init() override;

    void getConstraintViolation(const sofa::core::ConstraintParams* cParams ,
                                sofa::linearalgebra::BaseVector *resV,
                                const sofa::linearalgebra::BaseVector *Jdx) override;
    ///////////////////////////////////////////////////////////////

    sofa::Data<VecCoord>                                d_effectorGoal;

    void setTargetDefaultValue();
    void resizeData();

    ////////////////////////// Inherited attributes ////////////////////////////
    using softrobots::constraint::MagnetModel<DataTypes>::d_indices ;
    using softrobots::constraint::MagnetModel<DataTypes>::d_directions ;
    using softrobots::constraint::MagnetModel<DataTypes>::d_Jacobian ;
    using softrobots::constraint::MagnetModel<DataTypes>::d_PosSensor ;
    using softrobots::constraint::MagnetModel<DataTypes>::d_mum ;
    using softrobots::constraint::MagnetModel<DataTypes>::d_useDirections ;
    using softrobots::constraint::MagnetModel<DataTypes>::d_constraintIndex ;
    using softrobots::constraint::MagnetModel<DataTypes>::d_weight ;
    using softrobots::behavior::SoftRobotsConstraint<DataTypes>::m_state ;
    using softrobots::behavior::SoftRobotsConstraint<DataTypes>::d_componentState ;
    using Effector<DataTypes>::getTarget ;
    ///////////////////////////////////////////////////////////////////////////
    
    

};

#if !defined(SOFTROBOTS_INVERSE_MAGNETEFFECTOR_CPP)
extern template class SOFA_SOFTROBOTS_INVERSE_API MagnetEffector<sofa::defaulttype::Vec1Types>;
extern template class SOFA_SOFTROBOTS_INVERSE_API MagnetEffector<sofa::defaulttype::Vec2Types>;
extern template class SOFA_SOFTROBOTS_INVERSE_API MagnetEffector<sofa::defaulttype::Vec3Types>;
extern template class SOFA_SOFTROBOTS_INVERSE_API MagnetEffector<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace


