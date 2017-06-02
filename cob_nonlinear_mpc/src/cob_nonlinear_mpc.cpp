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
 *   Author: Christoph Mark, email: Christoph.Mark@ipa.fraunhofer.de
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
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <cob_srvs/SetString.h>
#include <string>
#include <Eigen/Dense>
#include <casadi/core/function/sx_function.hpp>
#include <stdlib.h>

bool CobNonlinearMPC::initialize()
{
    ros::NodeHandle nh_nmpc("nmpc");
    ros::NodeHandle nh_nmpc_dh("nmpc/dh");
    ros::NodeHandle nh_nmpc_base_dh("nmpc/base/dh");
    ros::NodeHandle nh_nmpc_constraints("nmpc/constraints");


    // nh_nmpc
    if (!nh_nmpc.getParam("transformations", transformation_names_))
    {
        ROS_ERROR("Parameter 'transformation_names' not set");
        return false;
    }

    if (!nh_nmpc.getParam("shooting_nodes", num_shooting_nodes_))
    {
        ROS_ERROR("Parameter 'num_shooting_nodes_' not set");
        return false;
    }

    if (!nh_nmpc.getParam("time_horizon", time_horizon_))
    {
        ROS_ERROR("Parameter 'time_horizon' not set");
        return false;
    }

    if (!nh_nmpc.getParam("state_dim", state_dim_))
    {
        ROS_ERROR("Parameter 'state_dim' not set");
        return false;
    }
    if (!nh_nmpc.getParam("control_dim", control_dim_))
    {
        ROS_ERROR("Parameter 'control_dim' not set");
        return false;
    }
    if (!nh_nmpc.getParam("base/base_active", base_active_))
    {
        ROS_ERROR("Parameter 'base/base_active' not set");
        return false;
    }
    if (!nh_nmpc.getParam("base/transformations", transformation_names_base_))
    {
        ROS_ERROR("Parameter 'base/transformations' not set");
        return false;
    }

    // nh_nmpc_constraints
    if (!nh_nmpc_constraints.getParam("state/path_constraints/min", state_path_constraints_min_))
    {
        ROS_ERROR("Parameter 'state/path_constraints/min' not set");
        return false;
    }
    if (!nh_nmpc_constraints.getParam("state/path_constraints/max", state_path_constraints_max_))
    {
        ROS_ERROR("Parameter 'state/path_constraints/max' not set");
        return false;
    }

    if (!nh_nmpc_constraints.getParam("state/terminal_constraints/min", state_terminal_constraints_min_))
    {
        ROS_ERROR("Parameter 'state/terminal_constraints/min' not set");
        return false;
    }
    if (!nh_nmpc_constraints.getParam("state/terminal_constraints/max", state_terminal_constraints_max_))
    {
        ROS_ERROR("Parameter 'state/terminal_constraints/max' not set");
        return false;
    }

    if (!nh_nmpc_constraints.getParam("input/input_constraints/min", input_constraints_min_))
    {
        ROS_ERROR("Parameter 'input/input_constraints/min' not set");
        return false;
    }
    if (!nh_nmpc_constraints.getParam("input/input_constraints/max", input_constraints_max_))
    {
        ROS_ERROR("Parameter 'input/input_constraints/max' not set");
        return false;
    }

    if(base_active_)
    {
        for(unsigned int i = 0; i < transformation_names_base_.size(); i++)
        {
            DH param;

            if (!nh_nmpc_base_dh.getParam(transformation_names_base_.at(i)+"/type", param.type))
            {
                ROS_ERROR("Parameter 'type' not set");
                return false;
            }

            if (!nh_nmpc_base_dh.getParam(transformation_names_base_.at(i)+"/theta", param.theta))
            {
                ROS_ERROR("Parameter 'theta' not set");
                return false;
            }
            if (!nh_nmpc_base_dh.getParam(transformation_names_base_.at(i)+"/d", param.d))
            {
                ROS_ERROR("Parameter 'd' not set");
                return false;
            }

            if (!nh_nmpc_base_dh.getParam(transformation_names_base_.at(i)+"/a", param.a))
            {
                ROS_ERROR("Parameter 'a' not set");
                return false;
            }

            if (!nh_nmpc_base_dh.getParam(transformation_names_base_.at(i)+"/alpha", param.alpha))
            {
                ROS_ERROR("Parameter 'alpha' not set");
                return false;
            }

            dh_params_base_.push_back(param);
        }
    }

    // nh_nmpc_dh
    for(unsigned int i = 0; i < transformation_names_.size(); i++)
    {
        DH param;

        if (!nh_nmpc_dh.getParam(transformation_names_.at(i)+"/type", param.type))
        {
            ROS_ERROR("Parameter 'type' not set");
            return false;
        }

        if (!nh_nmpc_dh.getParam(transformation_names_.at(i)+"/theta", param.theta))
        {
            ROS_ERROR("Parameter 'theta' not set");
            return false;
        }
        if (!nh_nmpc_dh.getParam(transformation_names_.at(i)+"/d", param.d))
        {
            ROS_ERROR("Parameter 'd' not set");
            return false;
        }

        if (!nh_nmpc_dh.getParam(transformation_names_.at(i)+"/a", param.a))
        {
            ROS_ERROR("Parameter 'a' not set");
            return false;
        }

        if (!nh_nmpc_dh.getParam(transformation_names_.at(i)+"/alpha", param.alpha))
        {
            ROS_ERROR("Parameter 'alpha' not set");
            return false;
        }

        dh_params.push_back(param);
    }

    if (!nh_.getParam("self_collision_matrix", scm_))
    {
        ROS_ERROR("Parameter 'self_collision_matrix' not set");
        return false;
    }

    for (XmlRpc::XmlRpcValue::iterator it = scm_.begin(); it != scm_.end(); ++it)
    {
        std::vector<std::string> empty_vec;
        self_collision_map_[it->first] = empty_vec;
        ROS_ASSERT(it->second.getType() == XmlRpc::XmlRpcValue::TypeArray);
        for (int j=0; j < it->second.size(); ++j)
        {
            ROS_ASSERT(it->second[j].getType() == XmlRpc::XmlRpcValue::TypeString);
            self_collision_map_[it->first].push_back(it->second[j]);
        }
    }

    if(base_active_)
    {
        if (!nh_.getParam("bounding_volume_base", bvb_))
        {
            ROS_ERROR("Parameter 'bounding_volume_base' not set");
            return false;
        }

        for (XmlRpc::XmlRpcValue::iterator it = bvb_.begin(); it != bvb_.end(); ++it)
        {
            if(it->first == "position")
            {
                for (int j=0; j < it->second.size(); ++j)
                {
                    ROS_ASSERT(it->second[j].getType() == XmlRpc::XmlRpcValue::TypeDouble);
                    bvb_positions_.push_back(it->second[j]);
                }
            }
            else if(it->first == "bv_radius")
            {
                for (int j=0; j < it->second.size(); ++j)
                {
                    ROS_ASSERT(it->second[j].getType() == XmlRpc::XmlRpcValue::TypeDouble);
                    bvb_radius_.push_back(it->second[j]);
                }
            }
            else
            {
                ROS_ERROR("Wrong bounding volume format");
            }
        }
    }

    ROS_WARN_STREAM(bvb_positions_);
    ROS_WARN_STREAM(bvb_radius_);

    // Casadi symbolics
    u_ = SX::sym("u", state_dim_);  // control
    x_ = SX::sym("x", control_dim_); // states

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

    /// parse robot_description and generate KDL chains
    KDL::Tree my_tree;
    if (!kdl_parser::treeFromParam("/robot_description", my_tree))
    {
        ROS_ERROR("Failed to construct kdl tree");
        return false;
    }

    my_tree.getChain(chain_base_link_, chain_tip_link_, chain_);
    if (chain_.getNrOfJoints() == 0)
    {
        ROS_ERROR("Failed to initialize kinematic chain");
        return false;
    }
    ROS_INFO_STREAM("Number of joints:" << chain_.getNrOfJoints());
    ROS_INFO_STREAM("Number of segments:" << chain_.getNrOfSegments());
    dof = joint_names.size();

    /// parse robot_description and set velocity limits
    urdf::Model model;
    if (!model.initParam("/robot_description"))
    {
        ROS_ERROR("Failed to parse urdf file for JointLimits");
        return false;
    }
    ROS_WARN("Robot Description loaded...");
    std::vector<double> joint_params_;
    urdf::Vector3 position;

    std::vector<KDL::Joint> joints;
    KDL::Frame F;
    double roll,pitch,yaw;
    for (uint16_t i = 0; i < chain_.getNrOfSegments(); i++)
    {
        joints.push_back(chain_.getSegment(i).getJoint());
        ROS_INFO_STREAM("Chain segment "<< chain_.getSegment(i).getName());
        F=chain_.getSegment(i).getFrameToTip();
        F.M.GetRPY(roll,pitch,yaw);
        ROS_INFO_STREAM("Chain frame "<< " X: " << F.p.x()<< " Y: " << F.p.y()<< " Z: "<<F.p.z());
        ROS_INFO_STREAM("Chain frame "<< " ROLL: " << roll<< " PITCH: " << pitch<< " YAW: "<<yaw);

    }

    // JointNames
    std::vector<KDL::Vector> joint_origins;
    for (uint16_t i = 0; i < joints.size(); i++)
    {
        joint_origins.push_back(joints[i].JointOrigin());
        ROS_INFO_STREAM("Joint name "<< joints[i].getName()<< " type: " <<joints[i].getType() << " origin: " << joint_origins[i].x());
        ROS_INFO_STREAM("Joint origin "<< " X: " << joint_origins[i].x()<< " Y: " << joint_origins[i].y()<< " Z: " << joint_origins[i].z());

    }

    KDL::Vector pos;
    KDL::Rotation rot;
    std::vector<KDL_frame> joint_frames;
    std::vector<KDL_frame> F_previous;


    for (uint16_t i = 0; i < chain_.getNrOfSegments(); i++)
    {
        KDL_frame f;
        if(joints[i].getType()==8)
        {
            ROS_INFO("Fixed joint");
            ROS_INFO_STREAM("Chain segment "<< chain_.getSegment(i).getName());
            if(i==0)
            {
                f.frame = chain_.getSegment(i).getFrameToTip();
                f.link = chain_.getSegment(i).getName();
                F_previous.push_back(f);
            }
            else
            {
                f.frame = F_previous.at(i-1).frame*chain_.getSegment(i).getFrameToTip();
                f.link = chain_.getSegment(i).getName();
                F_previous.push_back(f);
            }
//            ROS_INFO_STREAM("Joint position "<< " X: " << F_previous.at(i).p.x()<< " Y: " << F_previous.at(i).p.y()<< " Z: " << F_previous.at(i).p.z());
            rot=F_previous.at(i).frame.M;
//            ROS_WARN("Rotation matrix %f %f %f \n %f %f %f \n %f %f %f \n",rot(0,0),rot(0,1),rot(0,2),rot(1,0),rot(1,1),rot(1,2),rot(2,0),rot(2,1),rot(2,2));
//            ROS_INFO_STREAM("Joint position of transformation"<< " X: " << F_previous.at(i).p.x()<< " Y: " << F_previous.at(i).p.y()<< " Z: " << F_previous.at(i).p.z());
        }

        if(joints[i].getType()==0)
        {
            KDL_frame f;
            ROS_INFO("Rotational joint");
            ROS_INFO_STREAM("Joint name "<< chain_.getSegment(i).getJoint().getName());
            f.frame = F_previous.at(i-1).frame*chain_.getSegment(i).getFrameToTip();
            f.link = chain_.getSegment(i).getName();

            F_previous.push_back(f);
            pos=F_previous.at(i).frame.p;
            if(joint_frames.size()==0)
            {
//                ROS_INFO("FIRST JOINT");
                joint_frames.push_back(F_previous.at(i));
                rot=F_previous.at(i).frame.M;
                pos=F_previous.at(i).frame.p;
            }
            else
            {
                f.frame = chain_.getSegment(i).getFrameToTip();
                f.link = chain_.getSegment(i).getName();

                joint_frames.push_back(f);
                rot=chain_.getSegment(i).getFrameToTip().M;
                pos=chain_.getSegment(i).getFrameToTip().p;
            }

//
//            ROS_INFO_STREAM("Joint position "<< " X: " << pos.x()<< " Y: " << pos.y()<< " Z: " << pos.z());
//            ROS_WARN("Rotation matrix %f %f %f \n %f %f %f \n %f %f %f \n",rot(0,0),rot(0,1),rot(0,2),rot(1,0),rot(1,1),rot(1,2),rot(2,0),rot(2,1),rot(2,2));
//            ROS_INFO_STREAM("Joint position of transformation"<< " X: " << pos.x()<< " Y: " << pos.y()<< " Z: " << pos.z());            //F_previous.p= pos;
        }
    }

    SX T = SX::sym("T",4,4);

    //Base config
    T(0,0) = cos(x_(2)); T(0,1) = -sin(x_(2));  T(0,2) = 0.0; T(0,3) = x_(0);
    T(1,0) = sin(x_(2)); T(1,1) = cos(x_(2));   T(1,2) = 0.0; T(1,3) = x_(1);
    T(2,0) = 0.0;        T(2,1) = 0.0;          T(2,2) = 1.0; T(2,3) = 0;
    T(3,0) = 0.0;        T(3,1) = 0.0;          T(3,2) = 0.0; T(3,3) = 1.0;

    fk_base_ = T;

    for(int i=0;i<joint_frames.size();i++)
    {
        rot=joint_frames.at(i).frame.M;
        pos=joint_frames.at(i).frame.p;
        ROS_WARN("Rotation matrix %f %f %f \n %f %f %f \n %f %f %f \n",joint_frames.at(i).frame(0,0),joint_frames.at(i).frame(0,1),joint_frames.at(i).frame(0,2),joint_frames.at(i).frame(1,0),joint_frames.at(i).frame(1,1),joint_frames.at(i).frame(1,2),joint_frames.at(i).frame(2,0),joint_frames.at(i).frame(2,1),joint_frames.at(i).frame(2,2));
        ROS_INFO_STREAM("Joint position of transformation"<< " X: " << pos.x()<< " Y: " << pos.y()<< " Z: " << pos.z());
        T(0,0) = rot(0,0)*cos(x_(i+3))+rot(0,1)*sin(x_(i+3));
        T(0,1) = -rot(0,0)*sin(x_(i+3))+rot(0,1)*cos(x_(i+3));
        T(0,2) = rot(0,2); T(0,3) = pos.x();
        T(1,0) = rot(1,0)*cos(x_(i+3))+rot(1,1)*sin(x_(i+3));
        T(1,1) = -rot(1,0)*sin(x_(i+3))+rot(1,1)*cos(x_(i+3));
        T(1,2) = rot(1,2); T(1,3) = pos.y();
        T(2,0) = rot(2,0)*cos(x_(i+3))+rot(2,1)*sin(x_(i+3));
        T(2,1) = -rot(2,0)*sin(x_(i+3))+rot(2,1)*cos(x_(i+3));
        T(2,2) = rot(2,2); T(2,3) = pos.z();
        T(3,0) = 0.0; T(3,1) = 0.0; T(3,2) = 0.0; T(3,3) = 1.0;

        T_BVH p;
        p.link = joint_frames.at(i).link;
        p.T = T;

        if(pos.z() > 0.1)
        {
            T(2,3) = T(2,3)/2;
            p.constraint = true;
        }
        else
        {
            p.constraint = false;
        }
        p.BVH_p = T;
        transform_vec_bvh_.push_back(p);
    }

    // Get Endeffector FK
    for(int i=0; i< transform_vec_bvh_.size(); i++)
    {
        if(base_active_)
        {
            if(i==0)
            {   ROS_WARN("BASE IS ACTIVE");
                fk_ = mtimes(fk_base_,transform_vec_bvh_.at(i).T);
            }
            else
            {
                fk_ = mtimes(fk_,transform_vec_bvh_.at(i).T);
            }
        }
        else
        {
            if(i==0)
            {
                fk_ = transform_vec_bvh_.at(i).T;
            }
            else
            {
                fk_ = mtimes(fk_,transform_vec_bvh_.at(i).T);
            }
        }
        fk_vector_.push_back(fk_); // stacks up multiplied transformation until link n
    }

    // Get bounding volume forward kinematics
    for(int i=0; i<transform_vec_bvh_.size(); i++)
    {
        T_BVH bvh = transform_vec_bvh_.at(i);
        std::vector<SX> bvh_arm;
        if(i-1<0)
        {
            SX transform = mtimes(fk_vector_.at(i),bvh.T);
            SX tmp = SX::vertcat({transform(0,3), transform(1,3), transform(2,3)});
            bvh_arm.push_back(tmp);
            bvh_matrix[bvh.link].push_back(bvh_arm);

            if(bvh.constraint)
            {
                bvh_arm.clear();
                tmp.clear();
                transform = mtimes(fk_vector_.at(i),bvh.BVH_p);
                tmp = SX::vertcat({transform(0,3), transform(1,3), transform(2,3)});
                bvh_arm.push_back(tmp);
                bvh_matrix[bvh.link].push_back(bvh_arm);
            }
        }
        else
        {
            bvh_arm.clear();
            SX transform = mtimes(fk_vector_.at(i-1),bvh.T);
            SX tmp = SX::vertcat({transform(0,3), transform(1,3), transform(2,3)});
            bvh_arm.push_back(tmp);
            bvh_matrix[bvh.link].push_back(bvh_arm);
            bvh_arm.clear();

            if(bvh.constraint)
            {
                tmp.clear();
                transform = mtimes(fk_vector_.at(i-1),bvh.BVH_p);
                tmp = SX::vertcat({transform(0,3), transform(1,3), transform(2,3)});
                bvh_arm.push_back(tmp);
                bvh_matrix[bvh.link].push_back(bvh_arm);
            }
        }
    }
    if(base_active_)
    {
        for(int i=0; i<bvb_positions_.size(); i++)
        {
            std::vector<SX> base_bvh;
            SX tmp = SX::vertcat({fk_base_(0,3), fk_base_(1,3), fk_base_(2,3)+bvb_positions_.at(i)});
            base_bvh.push_back(tmp);
            bvh_matrix["body"].push_back(base_bvh);
        }

    }


    vector<double> tmp;
    for(int k=0; k < num_shooting_nodes_; ++k)
    {
        tmp.clear();
        for(int i=0; i < control_dim_; ++i)
        {
            tmp.push_back(0);
        }
        u_open_loop_.push_back(tmp);
    }

    for(int k=1; k <= num_shooting_nodes_; ++k)
    {
        tmp.clear();

        for(int i=0; i < state_dim_; ++i)
        {
            tmp.push_back(0);
        }

        x_open_loop_.push_back(tmp);
    }

    joint_state_ = KDL::JntArray(7);
    odometry_state_ = KDL::JntArray(3);
    jointstate_sub_ = nh_.subscribe("joint_states", 1, &CobNonlinearMPC::jointstateCallback, this);
    odometry_sub_ = nh_.subscribe("base/odometry", 1, &CobNonlinearMPC::odometryCallback, this);
    pose_sub_ = nh_.subscribe(nh_.getNamespace()+"/command_pose", 1, &CobNonlinearMPC::poseCallback, this);

    base_vel_pub_ = nh_.advertise<geometry_msgs::Twist>("base/command", 1);
    pub_ = nh_.advertise<std_msgs::Float64MultiArray>(nh_.getNamespace()+"/joint_group_velocity_controller/command", 1);

    ROS_WARN_STREAM(nh_.getNamespace() << "/NMPC...initialized!");
    return true;
}

void CobNonlinearMPC::poseCallback(const geometry_msgs::Pose::ConstPtr& msg)
{
    KDL::JntArray state = getJointState();

    Eigen::MatrixXd qdot = mpc_step(*msg, state);

    geometry_msgs::Twist base_vel_msg;
    std_msgs::Float64MultiArray vel_msg;

    base_vel_msg.linear.x = qdot(0);
    base_vel_msg.linear.y = qdot(1);
    base_vel_msg.linear.z = 0;
    base_vel_msg.angular.x = 0;
    base_vel_msg.angular.y = 0;
    base_vel_msg.angular.z = qdot(2);

    base_vel_pub_.publish(base_vel_msg);


    for (unsigned int i = 3; i < 10; i++)
    {
        vel_msg.data.push_back((double)qdot(i));
    }
    pub_.publish(vel_msg);
}


void CobNonlinearMPC::jointstateCallback(const sensor_msgs::JointState::ConstPtr& msg)
{

    KDL::JntArray q_temp = joint_state_;

    joint_names = {"arm_left_1_joint", "arm_left_2_joint", "arm_left_3_joint", "arm_left_4_joint", "arm_left_5_joint", "arm_left_6_joint", "arm_left_7_joint"};
//    joint_names = {"arm_1_joint","arm_2_joint","arm_3_joint","arm_4_joint","arm_5_joint","arm_6_joint","arm_7_joint"};

    int count = 0;

    for (uint16_t j = 0; j < 7; j++)
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
    temp(2) = msg->pose.pose.orientation.z;

    odometry_state_ = temp;
}


KDL::JntArray CobNonlinearMPC::getJointState()
{
    KDL:: JntArray tmp(joint_state_.rows() + odometry_state_.rows());
//    KDL:: JntArray tmp(joint_state_.rows());

//    tmp = this->odometry_state_;

    for(int i = 0; i < odometry_state_.rows(); i++)
    {
        tmp(i) = odometry_state_(i);
    }

    for(int i = 0 ; i < joint_state_.rows(); i++)
    {
        tmp(i+odometry_state_.rows()) = this->joint_state_(i);
    }

    return tmp;
}


Eigen::MatrixXd CobNonlinearMPC::mpc_step(const geometry_msgs::Pose pose,
                                          const KDL::JntArray& state)
{
    // Distance to obstacle
    double min_dist = 0.4;

    // Bounds and initial guess for the control
    vector<double> u_min =  input_constraints_min_;
    vector<double> u_max  = input_constraints_max_;

    // Bounds and initial guess for the state
    vector<double> x0_min;
    vector<double> x0_max;
    vector<double> x_init;
    for(unsigned int i=0; i < state.rows();i++)
    {
        x0_min.push_back(state(i));
        x0_max.push_back(state(i));
        x_init.push_back(state(i));
    }

    x_open_loop_.at(0) = x_init;

    vector<double> x_min  = state_path_constraints_min_;
    vector<double> x_max  = state_path_constraints_max_;
    vector<double> xf_min = state_terminal_constraints_min_;
    vector<double> xf_max = state_terminal_constraints_max_;

    // ODE right hand side and quadrature
    SX qdot = SX::vertcat({u_});

    // Current Quaternion and Position Vector.
    double kappa = 0.001; // Small regulation term for numerical stability for the NLP

    SX q_c = SX::vertcat({
        0.5 * sqrt(fk_(0,0) + fk_(1,1) + fk_(2,2) + 1.0 + kappa),
        0.5 * (sign((fk_(2,1) - fk_(1,2)))) * sqrt(fk_(0,0) - fk_(1,1) - fk_(2,2) + 1.0 + kappa),
        0.5 * (sign((fk_(0,2) - fk_(2,0)))) * sqrt(fk_(1,1) - fk_(2,2) - fk_(0,0) + 1.0 + kappa),
        0.5 * (sign((fk_(1,0) - fk_(0,1)))) * sqrt(fk_(2,2) - fk_(0,0) - fk_(1,1) + 1.0 + kappa)
    });

    SX p_c = SX::vertcat({fk_(0,3), fk_(1,3), fk_(2,3)});

    // Desired Goal-pose
    SX x_d = SX::vertcat({pose.position.x, pose.position.y, pose.position.z});
    SX q_d = SX::vertcat({pose.orientation.w, pose.orientation.x, pose.orientation.y, pose.orientation.z});

    // Prevent collision with Base_link
    SX barrier;
    SX dist;

    std::unordered_map<std::string, std::vector<std::string> >::iterator it_scm;

    int counter = 0;
    double bv_radius;

    for( it_scm = self_collision_map_.begin(); it_scm != self_collision_map_.end(); it_scm++)
    {
        std::vector<string> scm_collision_links = it_scm->second;
        for(int i=0; i<scm_collision_links.size(); i++)
        {
            ROS_WARN_STREAM(it_scm->first);
            vector<vector<SX>> p1_mat = bvh_matrix[it_scm->first];
            vector<vector<SX>> p2_mat = bvh_matrix[scm_collision_links.at(i)];

            for(int k=0; k<p1_mat.size(); k++)
            {
                if(it_scm->first == "body")
                {
                    bv_radius = bvb_radius_.at(k);
                }
                else
                {
                    bv_radius = 0.1;
                }

                vector<SX> p1_vec = p1_mat.at(k);
                for(int m=0; m<p2_mat.size(); m++)
                {
                    vector<SX> p2_vec = p2_mat.at(m);

                    SX p1 = SX::vertcat({p1_vec.at(0)});
                    SX p2 = SX::vertcat({p2_vec.at(0)});
                    dist = dot(p1 - p2, p1 - p2);

                    if(counter == 0)
                    {
                        barrier = exp((bv_radius - sqrt(dist))/0.01);
                        counter = 1;
                    }
                    else
                    {
                        barrier += exp((bv_radius - sqrt(dist))/0.01);
                    }
                }
            }
        }
    }

    // Get orientation error
    SX q_c_inverse = SX::vertcat({q_c(0), -q_c(1), -q_c(2), -q_c(3)});
    SX e_quat= quaternion_product(q_c_inverse,q_d);
    SX error_attitute = SX::vertcat({ e_quat(1), e_quat(2), e_quat(3)});

    // L2 norm of the control signal
    SX R = 1*SX::vertcat({100, 100, 100, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1});
    SX energy = dot(sqrt(R)*u_,sqrt(R)*u_);

    // L2 norm of the states
    SX S = 0.1*SX::vertcat({0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1});
    SX motion = dot(sqrt(S)*x_,sqrt(S)*x_);

    // Objective
    SX L = 10*dot(p_c-x_d,p_c-x_d) + energy + 10 * dot(error_attitute,error_attitute) + barrier;

    // Create Euler integrator function
    Function F = create_integrator(state_dim_, control_dim_, time_horizon_, num_shooting_nodes_, qdot, x_, u_, L);

    // Total number of NLP variables
    int NV = state_dim_*(num_shooting_nodes_+1) + control_dim_*num_shooting_nodes_;

    // Declare variable vector for the NLP
    MX V = MX::sym("V",NV);

    // NLP variable bounds and initial guess
    vector<double> v_min,v_max,v_init;

    // Offset in V
    int offset=0;

    // State at each shooting node and control for each shooting interval
    vector<MX> X, U;

    for(unsigned int k=0; k<num_shooting_nodes_; ++k)
    {
        // Local state
        X.push_back( V.nz(Slice(offset,offset+state_dim_)));

        vector<double> x_ol;
        x_ol = x_open_loop_.at(k);

        if(k==0)
        {
            v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
            v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
        }
        else
        {
            v_min.insert(v_min.end(), x_min.begin(), x_min.end());
            v_max.insert(v_max.end(), x_max.begin(), x_max.end());
        }
        v_init.insert(v_init.end(), x_ol.begin(), x_ol.end());
        offset += state_dim_;

        // Local control via shift initialization
        U.push_back( V.nz(Slice(offset,offset+control_dim_)));
        v_min.insert(v_min.end(), u_min.begin(), u_min.end());
        v_max.insert(v_max.end(), u_max.begin(), u_max.end());

        vector<double> u_ol;
        u_ol = u_open_loop_.at(k);

        v_init.insert(v_init.end(), u_ol.begin(), u_ol.end());
        offset += control_dim_;
    }

    vector<double> x_ol;
    x_ol = x_open_loop_.at(num_shooting_nodes_-1);
    // State at end
    X.push_back(V.nz(Slice(offset,offset+state_dim_)));
    v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
    v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
    v_init.insert(v_init.end(), x_ol.begin(), x_ol.end());
    offset += state_dim_;

    // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
    casadi_assert(offset==NV);

    // Objective function
    MX J = 0;

    //Constraint function and bounds
    vector<MX> g;

    // Loop over shooting nodes
    for(unsigned int k=0; k<num_shooting_nodes_; ++k)
    {
        // Create an evaluation node
        MXDict I_out = F( MXDict{ {"x0", X[k]}, {"p", U[k]} });

        // Save continuity constraints
        g.push_back( I_out.at("xf") - X[k+1] );

        // Add objective function contribution
        J += I_out.at("qf");
    }

    // NLP
    MXDict nlp = {{"x", V}, {"f", J}, {"g", vertcat(g)}};

    // Set options
    Dict opts;

    opts["ipopt.tol"] = 1e-4;
    opts["ipopt.max_iter"] = 20;
//    opts["ipopt.hessian_approximation"] = "limited-memory";
//    opts["ipopt.hessian_constant"] = "yes";
    opts["ipopt.linear_solver"] = "ma27";
    opts["ipopt.print_level"] = 0;
    opts["print_time"] = true;
    opts["expand"] = true;  // Removes overhead

    // Create an NLP solver and buffers
    Function solver = nlpsol("nlpsol", "ipopt", nlp, opts);

    std::map<std::string, DM> arg, res;

    // Bounds and initial guess
    arg["lbx"] = v_min;
    arg["ubx"] = v_max;
    arg["lbg"] = 0;
    arg["ubg"] = 0;
    arg["x0"] = v_init;

    // Solve the problem
    res = solver(arg);

    // Optimal solution of the NLP
    vector<double> V_opt(res.at("x"));
    vector<double> J_opt(res.at("f"));

    // Get the optimal control
    Eigen::VectorXd q_dot = Eigen::VectorXd::Zero(state_dim_);
//    Eigen::VectorXd x_new = Eigen::VectorXd::Zero(state_dim_);
    vector<double> x_new;
    SX sx_x_new;
    u_open_loop_.clear();
    x_open_loop_.clear();
    x_new.clear();

    for(int i=0; i<1; ++i)  // Copy only the first optimal control sequence
    {
        for(int j=0; j<control_dim_; ++j)
        {
            q_dot[j] = V_opt.at(state_dim_ + j);
            x_new.push_back(V_opt.at(j));
        }
    }
    sx_x_new = SX::vertcat({x_new});

    vector<double> tmp;

//    // Prepare state guess for next loop
    for(int k=1; k < num_shooting_nodes_; k++)
    {
        tmp.clear();

        for(int i=0; i < state_dim_; i++)
        {
            tmp.push_back((double)V_opt.at((k) * (state_dim_+control_dim_) + i));
        }
        x_open_loop_.push_back(tmp);

        if(k == num_shooting_nodes_-1)
        {
            x_open_loop_.push_back(tmp);
        }
    }

    // Prepare control input guess for next loop
    for(int k=1; k < num_shooting_nodes_; k++)
    {
        tmp.clear();

        for(int i=0; i < control_dim_; i++)
        {
            tmp.push_back((double)V_opt.at((k+1) * (state_dim_+control_dim_) + (i-control_dim_)));
        }
        u_open_loop_.push_back(tmp);

        if(k == num_shooting_nodes_-1)
        {
            u_open_loop_.push_back(tmp);
        }
    }


    // Plot bounding volumes
    geometry_msgs::Point point;
    point.x = 0;
    point.y = 0;
    point.z = 0;

    SX result;
    for( it_scm = self_collision_map_.begin(); it_scm != self_collision_map_.end(); it_scm++)
    {
        vector<string> tmp = it_scm->second;

        for(int i=0; i<tmp.size();i++)
        {
            vector<vector<SX>> SX_vec = bvh_matrix[tmp.at(i)];
            for(int k=0; k<SX_vec.size(); k++)
            {
                SX test = SX::horzcat({SX_vec.at(k).at(0)});
                Function tesssst = Function("test", {x_}, {test});
                result = tesssst(sx_x_new).at(0);
                point.x = (double)result(0);
                point.y = (double)result(1);
                point.z = (double)result(2);

                bv_radius = 0.1;
                visualizeBVH(point, bv_radius, i+i*tmp.size());
            }
        }


        vector<vector<SX>> SX_vec = bvh_matrix[it_scm->first];
        for(int k=0; k<SX_vec.size(); k++)
        {
            SX test = SX::horzcat({SX_vec.at(k).at(0)});
            Function tesssst = Function("test", {x_}, {test});
            result = tesssst(sx_x_new).at(0);
            point.x = (double)result(0);
            point.y = (double)result(1);
            point.z = (double)result(2);

            if(it_scm->first == "body")
            {
                bv_radius = bvb_radius_.at(k);
            }
            else
            {
                bv_radius = 0.1;
            }
            visualizeBVH(point, bv_radius, k+tmp.size()+SX_vec.size());
        }
    }
    return q_dot;
}

Function CobNonlinearMPC::create_integrator(const unsigned int state_dim, const unsigned int control_dim, const double T,
                                            const unsigned int N, SX ode, SX x, SX u, SX L)
{
    // Euler discretize
    double dt = T/((double)N);

    Function f = Function("f", {x, u}, {ode, L});
//
//    f.generate("f");
//
//    // Compile the C-code to a shared library
//    string compile_command = "gcc -fPIC -shared -O3 f.c -o f.so";
//    int flag = system(compile_command.c_str());
//    casadi_assert_message(flag==0, "Compilation failed");
//
//    Function f_ext = external("f");

    MX X0 = MX::sym("X0", state_dim);
    MX U_ = MX::sym("U",control_dim);
    MX X_ = X0;
    MX Q = 0;

    vector<MX> input(2);
    input[0] = X_;
    input[1] = U_;
    MX qdot_new = f(input).at(0);
    MX Q_new = f(input).at(1);

    X_= X_+ dt * qdot_new;
    Q = Q + dt * Q_new;

    Function F = Function("F", {X0, U_}, {X_, Q}, {"x0","p"}, {"xf", "qf"});
    return F.expand("F");   // Remove overhead
}


SX CobNonlinearMPC::dual_quaternion_product(SX q1, SX q2)
{
    SX q1_real = SX::vertcat({q1(0),q1(1),q1(2),q1(3)});
    SX q1_dual = SX::vertcat({q1(4),q1(5),q1(6),q1(7)});
    SX q2_real = SX::vertcat({q2(0),q2(1),q2(2),q2(3)});
    SX q2_dual = SX::vertcat({q2(4),q2(5),q2(6),q2(7)});

    SX q1q2_real = quaternion_product(q1_real,q2_real);
    SX q1_real_q2_dual = quaternion_product(q1_real,q2_dual);
    SX q1_dual_q2_real = quaternion_product(q1_dual,q2_real);

    SX q_prod = SX::vertcat({
        q1q2_real,
        q1_real_q2_dual + q1_dual_q2_real
    });

    return q_prod;
}

SX CobNonlinearMPC::quaternion_product(SX q1, SX q2)
{
    SX q1_v = SX::vertcat({q1(1),q1(2),q1(3)});
    SX q2_v = SX::vertcat({q2(1),q2(2),q2(3)});

    SX c = SX::cross(q1_v,q2_v);

    SX q_new = SX::vertcat({
        q1(0) * q2(0) - dot(q1_v,q2_v),
        q1(0) * q2(1) + q2(0) * q1(1) + c(0),
        q1(0) * q2(2) + q2(0) * q1(2) + c(1),
        q1(0) * q2(3) + q2(0) * q1(3) + c(2)
    });

    return q_new;
}


void CobNonlinearMPC::visualizeBVH(const geometry_msgs::Point point, double radius, int id)
{
    visualization_msgs::Marker marker;
    marker.type = visualization_msgs::Marker::SPHERE;
    marker.lifetime = ros::Duration();
    marker.action = visualization_msgs::Marker::ADD;
    marker.ns = "preview";
    marker.header.frame_id = "odom_combined";


    marker.scale.x = 2*radius;
    marker.scale.y = 2*radius;
    marker.scale.z = 2*radius;

    marker.color.r = 1.0;
    marker.color.g = 0.0;
    marker.color.b = 0.0;
    marker.color.a = 0.1;

    marker_array_.markers.clear();

    marker.id = id;
    marker.pose.position.x = point.x;
    marker.pose.position.y = point.y;
    marker.pose.position.z = point.z;
    marker_array_.markers.push_back(marker);

    marker_pub_.publish(marker_array_);
}
