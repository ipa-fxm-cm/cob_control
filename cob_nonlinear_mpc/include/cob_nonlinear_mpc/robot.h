/*!
 *****************************************************************
 * \file
 *
 * \note
 *   Copyright (c) 2017 \n
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
 *   ROS package name: cob_nonlinear_mpc
 *
 * \author
 *   Author: Bruno Brito  email: Bruno.Brito@ipa.fraunhofer.de
 *   Christoph Mark, email: Christoph.Mark@ipa.fraunhofer.de
 *
 * \date Date of creation: April, 2017
 *
 * \brief
 *
 *
 ****************************************************************/

#ifndef ROBOT_H
#define ROBOT_H

#include <ros/ros.h>

#include <sensor_msgs/JointState.h>
#include <geometry_msgs/Twist.h>
#include <std_msgs/Float64MultiArray.h>
#include <nav_msgs/Odometry.h>

#include <urdf/model.h>

#include <kdl_parser/kdl_parser.hpp>
#include <kdl/jntarray.hpp>
#include <kdl/jntarrayvel.hpp>
#include <kdl/frames.hpp>

struct Robot
{
    int dof;
    urdf::Model urdf;
    KDL::Chain kinematic_chain;
    std::vector<KDL::Joint> joints;
    std::vector<KDL::Frame> joint_frames;
    std::vector<KDL::Segment> forward_kinematics;
    bool base_active_;
    std::unordered_map<std::string, std::vector<std::string> > self_collision_map_;
    std::string root_frame;
};

#endif  // NONLINEAR_MPC_H
