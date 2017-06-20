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

#ifndef COB_NONLINEAR_MPC_COB_NONLINEAR_MPC_H
#define COB_NONLINEAR_MPC_COB_NONLINEAR_MPC_H

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

#include <tf/transform_listener.h>
#include <tf/tf.h>
#include <visualization_msgs/MarkerArray.h>

#include <ctime>
#include <casadi/casadi.hpp>

#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include <cob_nonlinear_mpc/forward_kinematics.h>
#include <cob_nonlinear_mpc/bounding_volumes.h>
#include <cob_nonlinear_mpc/nonlinear_mpc.h>


using namespace casadi;
using namespace std;

class CobNonlinearMPC
{
private:
    ros::NodeHandle nh_;
    std::vector<std::string> transformation_names_;
    std::vector<std::string> transformation_names_base_;
    std::vector<std::string> joint_names;
    std::string root_frame_;

    std::string chain_base_link_;
    std::string chain_tip_link_;

    Robot robot_;

    boost::shared_ptr<MPC> mpc_ctr_;

    tf::TransformListener tf_listener_;

    ros::Subscriber jointstate_sub_;
    ros::Subscriber odometry_sub_;
    ros::Subscriber frame_tracker_sub_;

    ros::Publisher base_vel_pub_;
    ros::Publisher pub_;


    KDL::JntArray joint_state_;
    KDL::JntArray odometry_state_;
    KDL::Tree robot_tree_;

    double min_dist;



    XmlRpc::XmlRpcValue scm_;
public:
    CobNonlinearMPC()
    {
    }
    ~CobNonlinearMPC(){

    }


    bool initialize();

    void jointstateCallback(const sensor_msgs::JointState::ConstPtr& msg);
    void odometryCallback(const nav_msgs::Odometry::ConstPtr& msg);

    void FrameTrackerCallback(const geometry_msgs::Pose::ConstPtr& msg);
    KDL::JntArray getJointState();

    bool process_KDL_tree();
};


#endif  // COB_NONLINEAR_MPC_COB_NONLINEAR_MPC_H
