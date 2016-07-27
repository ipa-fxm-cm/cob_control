/*!
 *****************************************************************
 * \file
 *
 * \note
 *   Copyright (c) 2014 \n
 *   Fraunhofer Institute for Manufacturing Engineering
 *   and Automation (IPA) \n\n
 *
 *****************************************************************
 *
 * \note
 *   Project name: care-o-bot
 * \note
 *   ROS stack name: cob_control
 * \note
 *   ROS package name: cob_twist_controller
 *
 * \author
 *   Author: Felix Messmer, email: Felix.Messmer@ipa.fraunhofer.de
 *
 * \date Date of creation: April, 2014
 *
 * \brief
 *   This package provides the implementation of an inverse kinematics solver.
 *
 ****************************************************************/

#include <ros/ros.h>
#include <eigen_conversions/eigen_kdl.h>
#include <kdl/chainfksolvervel_recursive.hpp>

#include "cob_twist_controller/inverse_differential_kinematics_solver.h"
#include <cob_twist_controller/inverse_jacobian_calculations/inverse_jacobian_calculation.h>

/**
 * Solve the inverse kinematics problem at the first order differential level.
 */
int InverseDifferentialKinematicsSolver::CartToJnt(const JointStates& joint_states,
                                                   const geometry_msgs::Pose pose,
                                                   const KDL::Twist& twist,
                                                   const KDL::Twist& v_in,
                                                   KDL::JntArray& qdot_chain_out,
                                                   KDL::JntArray& qdot_base_out)
{
    // ROS_INFO_STREAM("joint_states.current_q_: " << joint_states.current_q_.rows());
    int8_t retStat = -1;

    /// Let the ChainJntToJacSolver calculate the jacobian "jac_chain" for the current joint positions "q_in"
    KDL::Jacobian jac_chain(chain_.getNrOfJoints());
    jnt2jac_.JntToJac(joint_states.current_q_, jac_chain);
    // ROS_INFO_STREAM("jac_chain.rows: " << jac_chain.rows() << ", jac_chain.columns: " << jac_chain.columns());

    JointStates joint_states_full = this->kinematic_extension_->adjustJointStates(joint_states, pose, twist);
    // ROS_INFO_STREAM("joint_states_full.current_q_: " << joint_states_full.current_q_.rows());

    /// append columns to Jacobian in order to reflect additional DoFs of kinematical extension
    KDL::Jacobian jac_full = this->kinematic_extension_->adjustJacobian(jac_chain);
    // ROS_INFO_STREAM("jac_full.rows: " << jac_full.rows() << ", jac_full.columns: " << jac_full.columns());

    KDL::Jacobian jac_base;
    Vector6d_t v_in_vec;
    Eigen::MatrixXd qdot_out_chain_vec;
    PInvBySVD pinv_calc;
    Eigen::MatrixXd qdots_base_out;
    Eigen::VectorXd singularValues;
    Eigen::MatrixXd pinv_base;
    Eigen::MatrixXd svd_U;

    tf::twistKDLToEigen(v_in, v_in_vec);

    if(params_.kinematic_extension == BASE_ACTIVE)
    {
        // Seperate Base and Chain jacobian
        Eigen::MatrixXd jac_base_matrix = jac_full.data.block<6,6>(0,6);
        jac_base.data = jac_base_matrix;


        // Calculate Singular directions
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(jac_chain.data, Eigen::ComputeThinU | Eigen::ComputeThinV);
        singularValues = svd.singularValues();
        pinv_base = pinv_calc.calculate(params_, jac_base.data);
        svd_U = svd.matrixU();
    }

    retStat = constraint_solver_factory_.calculateJointVelocities(jac_chain.data,
                                                                  v_in_vec,
                                                                  joint_states,
                                                                  qdot_out_chain_vec);


    KDL::JntArray qdot_out_base(jac_base.columns());
    KDL::Twist t_bl;

    if(params_.kinematic_extension == BASE_ACTIVE)
    {
        // Decide whether base should be moved or not


        if(singularValues(singularValues.size()-1) < 0.1)
        {
            qdots_base_out = pinv_base * svd_U.block<6,1>(0,svd_U.cols()-1);

            Eigen::MatrixXd last_U = svd_U.block<6,1>(0,svd_U.cols()-1);

            ROS_WARN_STREAM("u: " << last_U(0,0) << ", " << last_U(0,1) << ", " << last_U(0,2));


            /// convert output
            for (int i = 0; i < jac_base.columns(); i++)
            {
                qdot_out_base(i) = qdots_base_out(i,0);
            }
        }
        else
        {
            for (int i = 0; i < jac_base.columns(); i++)
            {
                qdot_out_base(i) = 0.0;
            }
        }

        tf::StampedTransform cb_transform_bl;
        KDL::Frame  cb_frame_bl;

        /// get required transformations
        try
        {
            ros::Time now = ros::Time(0);
            tf_listener_.waitForTransform(params_.chain_base_link, "base_link", now, ros::Duration(0.5));
            tf_listener_.lookupTransform(params_.chain_base_link, "base_link", now, cb_transform_bl);
        }
        catch (tf::TransformException& ex)
        {
            ROS_ERROR("%s", ex.what());
        }
        cb_frame_bl.p = KDL::Vector(cb_transform_bl.getOrigin().x(), cb_transform_bl.getOrigin().y(), cb_transform_bl.getOrigin().z());
        cb_frame_bl.M = KDL::Rotation::Quaternion(cb_transform_bl.getRotation().x(), cb_transform_bl.getRotation().y(), cb_transform_bl.getRotation().z(), cb_transform_bl.getRotation().w());

        KDL::Twist t_cb;

        t_cb.vel = KDL::Vector(qdot_out_base(0), qdot_out_base(1), qdot_out_base(2));
        t_cb.rot = KDL::Vector(qdot_out_base(3), qdot_out_base(4), qdot_out_base(5));

        t_bl = cb_frame_bl * t_cb;

    }

    /// convert output
    KDL::JntArray qdot_out_chain(jac_chain.columns());
    for (int i = 0; i < jac_chain.columns(); i++)
    {
        qdot_out_chain(i) = qdot_out_chain_vec(i);
    }

    /// limiters shut be applied here in order to be able to consider the additional DoFs within "AllLimit", too
    qdot_out_chain = this->limiters_->enforceLimits(qdot_out_chain, joint_states.current_q_);

    /// process result for kinematical extension
//    this->kinematic_extension_->processResultExtension(qdot_out_chain);

    geometry_msgs::Twist base_vel_msg;

    double gain = 0.01;

    base_vel_msg.linear.x = gain*t_bl.vel[0];
    base_vel_msg.linear.y = gain*t_bl.vel[1];
    base_vel_msg.linear.z = gain*t_bl.vel[2];
    base_vel_msg.angular.x = gain*t_bl.rot[0];
    base_vel_msg.angular.y = gain*t_bl.rot[1];
    base_vel_msg.angular.z = gain*t_bl.rot[2];

    ROS_WARN_STREAM("twist: " << base_vel_msg);
    base_vel_pub_.publish(base_vel_msg);

    /// then qdot_out should be resized to contain only the chain_qdot_out's again
    for (int i = 0; i < jac_chain.columns(); i++)
    {
        qdot_chain_out(i) = qdot_out_chain(i);
    }

//    for (int i = 0; i < jac_base.columns(); i++)
//    {
//        qdot_base_out(i) = qdot_out_base(i);
//    }

    return retStat;
}

void InverseDifferentialKinematicsSolver::resetAll(TwistControllerParams params)
{
    this->params_ = params;

    this->kinematic_extension_.reset(KinematicExtensionBuilder::createKinematicExtension(this->params_));
    this->limiter_params_ = this->kinematic_extension_->adjustLimiterParams(this->params_.limiter_params);

    this->limiters_.reset(new LimiterContainer(this->limiter_params_));
    this->limiters_->init();

    this->task_stack_controller_.clearAllTasks();
    if (0 != this->constraint_solver_factory_.resetAll(this->params_, this->limiter_params_))  // params member as reference!!! else process will die!
    {
        ROS_ERROR("Failed to reset IDK constraint solver after dynamic_reconfigure.");
    }
}
