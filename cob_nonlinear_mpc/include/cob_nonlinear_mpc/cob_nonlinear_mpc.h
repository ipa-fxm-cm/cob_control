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
 *   ROS stack name: cob_driver
 * \note
 *   ROS package name: cob_nonlinear_mpc
 *
 * \author
 *   Author: Christoph Mark, email: Christoph.Mark@ipa.fraunhofer.de
 *
 * \date Date of creation: April, 2017
 *
 * \brief
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

using namespace casadi;
using namespace std;


struct DH
{
    double alpha;
    double a;
    std::string d;
    std::string theta;

    std::string type;
};

struct T_BVH
{
    SX T;
    SX BVH_p;

    bool constraint = false;
};


class CobNonlinearMPC
{
private:
    ros::NodeHandle nh_;
    std::vector<std::string> transformation_names_;
    std::vector<std::string> transformation_names_base_;
    std::vector<std::string> joint_names;

    std::string chain_base_link_;
    std::string chain_tip_link_;

    int dof;

    KDL::Chain chain_;

    DH dh_param;
    std::vector<DH> dh_params;
    std::vector<DH> dh_params_base_;

    std::vector<SX> transformation_vector_dual_quat_;

    std::vector<SX> bvh_points_;

    std::vector<T_BVH> transform_vec_bvh_;
    std::vector<SX> transform_base_vec_;
    int state_dim_;
    int control_dim_;
    int discretization_steps_;
    double time_horizon_;
    double h_;

    int p_order_;
    int num_finite_elements_;
    // Coefficients of the collocation equation
    vector<vector<double> > C_;
    // Coefficients of the continuity equation
    vector<double> D_;

    // Coefficients of the quadrature function
    vector<double> B_;
    visualization_msgs::MarkerArray marker_array_;

    // Declare variables
    SX u_;
    SX x_;
    SX fk_;
    SX fk_link4_;
    std::vector<SX> fk_vector_;
    SX fk_base_;

    SX fk_dual_quat_;
    SX p_;
    tf::TransformListener tf_listener_;
    SX qdot_;

    // Constraints
    vector<double> state_path_constraints_min_, state_path_constraints_max_;
    vector<double> state_terminal_constraints_min_, state_terminal_constraints_max_;
    vector<double> input_constraints_min_, input_constraints_max_;

    ros::Subscriber jointstate_sub_;
    ros::Subscriber odometry_sub_;
    ros::Subscriber pose_sub_;

    ros::Publisher base_vel_pub_;
    ros::Publisher pub_;
    ros::Publisher marker_pub_;

    KDL::JntArray joint_state_;
    KDL::JntArray odometry_state_;
    vector<vector<double>> u_open_loop_;
    vector<vector<double>> x_open_loop_;

    bool base_active_;


public:
    CobNonlinearMPC()
    {
//            u_init_ = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        marker_pub_ = nh_.advertise<visualization_msgs::MarkerArray>("nmpc/bvh", 1);
    }
    ~CobNonlinearMPC(){}


    bool initialize();

    void jointstateCallback(const sensor_msgs::JointState::ConstPtr& msg);
    void odometryCallback(const nav_msgs::Odometry::ConstPtr& msg);

    void poseCallback(const geometry_msgs::Pose::ConstPtr& msg);
    KDL::JntArray getJointState();

    Eigen::MatrixXd mpc_step(const KDL::JntArray& state,
                             SX objective);

    SX get_objective(const geometry_msgs::Pose pose,
                       const KDL::JntArray& state);

    void visualizeBVH(const geometry_msgs::Point point, double radius, int id);

    SX quaternion_product(SX q1, SX q2);
    SX dual_quaternion_product(SX q1, SX q2);
};


#endif  // COB_NONLINEAR_MPC_COB_NONLINEAR_MPC_H
