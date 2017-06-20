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
#include <string>
#include <vector>
#include <limits>
#include <ros/ros.h>

#include <cob_nonlinear_mpc/cob_nonlinear_mpc.h>

#include <kdl_conversions/kdl_msg.h>
#include <tf/transform_datatypes.h>

#include <cob_srvs/SetString.h>
#include <string>
#include <Eigen/Dense>
#include <casadi/core/function/sx_function.hpp>
#include <stdlib.h>

bool CobNonlinearMPC::initialize()
{
    ros::NodeHandle nh_nmpc("nmpc");
    ros::NodeHandle nh_nmpc_constraints("nmpc/constraints");

    ForwardKinematics fk;
    BoundingVolume bv;

    // JointNames
    if (!nh_.getParam("joint_names", joint_names))
    {
        ROS_ERROR("Parameter 'joint_names' not set");
        return false;
    }

    if(!nh_.getParam("root_frame", root_frame_))
    {
        ROS_ERROR("Parameter 'root_frame' not set");
        return false;
    }

    robot_.root_frame = root_frame_;

    // nh_nmpc
    int num_shooting_nodes;
    if (!nh_nmpc.getParam("shooting_nodes", num_shooting_nodes))
    {
        ROS_ERROR("Parameter 'num_shooting_nodes_' not set");
        return false;
    }
    double time_horizon;
    if (!nh_nmpc.getParam("time_horizon", time_horizon))
    {
        ROS_ERROR("Parameter 'time_horizon' not set");
        return false;
    }
    int state_dim;
    if (!nh_nmpc.getParam("state_dim", state_dim))
    {
        ROS_ERROR("Parameter 'state_dim' not set");
        return false;
    }
    int control_dim;
    if (!nh_nmpc.getParam("control_dim", control_dim))
    {
        ROS_ERROR("Parameter 'control_dim' not set");
        return false;
    }

    mpc_ctr_.reset(new MPC(num_shooting_nodes,time_horizon ,state_dim,control_dim));
    if (!nh_nmpc.getParam("base/base_active", robot_.base_active_))
    {
        ROS_ERROR("Parameter 'base/base_active' not set");
        return false;
    }

    // nh_nmpc_constraints
    vector<double> state_path_constraints_min;
    if (!nh_nmpc_constraints.getParam("state/path_constraints/min", state_path_constraints_min))
    {
        ROS_ERROR("Parameter 'state/path_constraints/min' not set");
        return false;
    }
    vector<double> state_path_constraints_max;
    if (!nh_nmpc_constraints.getParam("state/path_constraints/max", state_path_constraints_max))
    {
        ROS_ERROR("Parameter 'state/path_constraints/max' not set");
        return false;
    }
    mpc_ctr_->set_path_constraints(state_path_constraints_min,state_path_constraints_max);

    vector<double> state_terminal_constraints_min;
    if (!nh_nmpc_constraints.getParam("state/terminal_constraints/min", state_terminal_constraints_min))
    {
        ROS_ERROR("Parameter 'state/terminal_constraints/min' not set");
        return false;
    }
    vector<double> state_terminal_constraints_max;
    if (!nh_nmpc_constraints.getParam("state/terminal_constraints/max", state_terminal_constraints_max))
    {
        ROS_ERROR("Parameter 'state/terminal_constraints/max' not set");
        return false;
    }
    mpc_ctr_->set_state_constraints(state_terminal_constraints_min,state_terminal_constraints_max);
    vector<double> input_constraints_min;
    if (!nh_nmpc_constraints.getParam("input/input_constraints/min", input_constraints_min))
    {
        ROS_ERROR("Parameter 'input/input_constraints/min' not set");
        return false;
    }
    vector<double> input_constraints_max;
    if (!nh_nmpc_constraints.getParam("input/input_constraints/max", input_constraints_max))
    {
        ROS_ERROR("Parameter 'input/input_constraints/max' not set");
        return false;
    }
    mpc_ctr_->set_input_constraints(input_constraints_min,input_constraints_max);

    if (!nh_nmpc_constraints.getParam("min_distance", min_dist))
    {
            ROS_ERROR("Parameter 'input/input_constraints/min' not set");
            return false;
    }

    // Chain
    if (!nh_.getParam("chain_base_link", chain_base_link_))
    {
        ROS_ERROR("Parameter 'chain_base_link' not set");
        return false;
    }

    if (!nh_.getParam("chain_tip_link", chain_tip_link_))
    {
        ROS_ERROR("Parameter 'chain_tip_link' not set");
        return false;
    }

    if (!nh_.getParam("self_collision_matrix", scm_))
    {
        ROS_ERROR_STREAM("Parameter '" << nh_.getNamespace() << "/self_collision_matrix' not set");
        //return false;
    }

    for (XmlRpc::XmlRpcValue::iterator it = scm_.begin(); it != scm_.end(); ++it)
    {
        std::vector<std::string> empty_vec;
        robot_.self_collision_map_[it->first] = empty_vec;
        ROS_ASSERT(it->second.getType() == XmlRpc::XmlRpcValue::TypeArray);
        for (int j=0; j < it->second.size(); ++j)
        {
            ROS_ASSERT(it->second[j].getType() == XmlRpc::XmlRpcValue::TypeString);
            robot_.self_collision_map_[it->first].push_back(it->second[j]);
        }
    }

    XmlRpc::XmlRpcValue bvb;

    vector<double> bvb_positions, bvb_radius;

    if(robot_.base_active_)
    {
        if (!nh_.getParam("bounding_volume_base", bvb))
        {
            ROS_ERROR("Parameter 'bounding_volume_base' not set");
            return false;
        }

        for (XmlRpc::XmlRpcValue::iterator it = bvb.begin(); it != bvb.end(); ++it)
        {
            if(it->first == "position")
            {
                for (int j=0; j < it->second.size(); ++j)
                {
                    ROS_ASSERT(it->second[j].getType() == XmlRpc::XmlRpcValue::TypeDouble);
                    bvb_positions.push_back(it->second[j]);
                }
            }
            else if(it->first == "bv_radius")
            {
                for (int j=0; j < it->second.size(); ++j)
                {
                    ROS_ASSERT(it->second[j].getType() == XmlRpc::XmlRpcValue::TypeDouble);
                    bvb_radius.push_back(it->second[j]);
                }
            }
            else
            {
                ROS_ERROR("Wrong bounding volume format");
            }
        }
    }

    bv.setBVBpositions(bvb_positions);
    bv.setBVBradius(bvb_radius);

    // Parse Robot model
    this->process_KDL_tree();

    SX u = SX::sym("u", control_dim);  // control
    SX x = SX::sym("x", state_dim); // states

    // Calculate symbolic forward kinematics
    fk.setX(x);
    fk.setU(u);
    fk.symbolic_fk(robot_);

    // Bounding Volumes
    bv.setRobot(robot_);
    bv.setForwardKinematic(fk);
    bv.generate_bounding_volumes();

    // MPC stuff
    mpc_ctr_->setBoundingVolumes(bv);
    mpc_ctr_->setForwardKinematics(fk);
    mpc_ctr_->init();

    joint_state_ = KDL::JntArray(robot_.kinematic_chain.getNrOfJoints());
    jointstate_sub_ = nh_.subscribe("joint_states", 1, &CobNonlinearMPC::jointstateCallback, this);

    if(robot_.base_active_)
    {
        odometry_state_ = KDL::JntArray(3);
        odometry_sub_ = nh_.subscribe("base/odometry", 1, &CobNonlinearMPC::odometryCallback, this);
        base_vel_pub_ = nh_.advertise<geometry_msgs::Twist>("base/command", 1);
    }

    frame_tracker_sub_ = nh_.subscribe("command_pose", 1, &CobNonlinearMPC::FrameTrackerCallback, this);
    pub_ = nh_.advertise<std_msgs::Float64MultiArray>("joint_group_velocity_controller/command", 1);

    ROS_WARN_STREAM(nh_.getNamespace() << "/NMPC...initialized!");
    return true;
}

void CobNonlinearMPC::FrameTrackerCallback(const geometry_msgs::Pose::ConstPtr& msg)
{
    KDL::JntArray state = getJointState();

    ros::Time time = ros::Time::now();

    Eigen::MatrixXd qdot = mpc_ctr_->mpc_step(*msg, state);

    ros::Time time_new = ros::Time::now();

    ROS_INFO_STREAM("mpc time: " << (time_new - time).toSec());

    geometry_msgs::Twist base_vel_msg;
    std_msgs::Float64MultiArray vel_msg;

    if(robot_.base_active_)
    {
        base_vel_msg.linear.x = qdot(0);
        base_vel_msg.linear.y = qdot(1);
        base_vel_msg.linear.z = 0;
        base_vel_msg.angular.x = 0;
        base_vel_msg.angular.y = 0;
        base_vel_msg.angular.z = qdot(2);

        base_vel_pub_.publish(base_vel_msg);

        for (unsigned int i = 3; i < joint_state_.rows()+3; i++)
        {
            vel_msg.data.push_back(qdot(i));
        }
        pub_.publish(vel_msg);

    }
    else{
        for (unsigned int i = 0; i < joint_state_.rows(); i++)
        {
            vel_msg.data.push_back(qdot(i));
        }
        pub_.publish(vel_msg);
    }
}


void CobNonlinearMPC::jointstateCallback(const sensor_msgs::JointState::ConstPtr& msg)
{
    KDL::JntArray q_temp = joint_state_;

    int count = 0;

    for (uint16_t j = 0; j < joint_state_.rows(); j++)
    {
        for (uint16_t i = 0; i < msg->name.size(); i++)
        {
            if (strcmp(msg->name[i].c_str(), joint_names[j].c_str()) == 0)
            {
                q_temp(j) = msg->position[i];
                count++;
                break;
            }
        }
    }
    joint_state_ = q_temp;
}

void CobNonlinearMPC::odometryCallback(const nav_msgs::Odometry::ConstPtr& msg)
{
    KDL::JntArray temp = odometry_state_;
    KDL::Frame odom_frame_bl;
    tf::StampedTransform odom_transform_bl;

    temp(0) = msg->pose.pose.position.x;
    temp(1) = msg->pose.pose.position.y;
    //CONVERT FROM QUATERNION TO JOINT ANGLE ROTATION
    double ysqr = msg->pose.pose.orientation.y * msg->pose.pose.orientation.y;
    double t3 = +2.0 * (msg->pose.pose.orientation.w * msg->pose.pose.orientation.z + msg->pose.pose.orientation.x * msg->pose.pose.orientation.y);
    double t4 = +1.0 - 2.0 * (ysqr + msg->pose.pose.orientation.z * msg->pose.pose.orientation.z);
    temp(2) = std::atan2(t3, t4);

    odometry_state_ = temp;
}


KDL::JntArray CobNonlinearMPC::getJointState()
{
    KDL:: JntArray state(joint_state_.rows() + odometry_state_.rows());

    for(int i = 0; i < odometry_state_.rows(); i++)
    {
        state(i) = odometry_state_(i);
    }

    for(int i = 0 ; i < joint_state_.rows(); i++)
    {
        state(i+odometry_state_.rows()) = this->joint_state_(i);
    }

    for(unsigned int i=0; i < state.rows();i++)
    {
        mpc_ctr_->x0_min.push_back(state(i));
        mpc_ctr_->x0_max.push_back(state(i));
        mpc_ctr_->x_init.push_back(state(i));
    }

    return state;
}
