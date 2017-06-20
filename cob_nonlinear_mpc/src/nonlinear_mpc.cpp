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

#include <cob_nonlinear_mpc/nonlinear_mpc.h>

void MPC::init()
{
    // Total number of NLP variables
    NV = state_dim_*(num_shooting_nodes_+1) +control_dim_*num_shooting_nodes_;

    V = MX::sym("V",NV);
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

    for(int i = 0; i < control_dim_; i++)
    {
        u_init_.push_back(0);
    }

    x_ = _fk_.getX();
    u_ = _fk_.getU();

    // Get end-effector fk
    fk_ = _fk_.getFkVector().at(_fk_.getFkVector().size()-1);

    ROS_WARN_STREAM("MPC initialized");
}

int MPC::get_num_shooting_nodes(){
    return this->num_shooting_nodes_;
}

double MPC::get_time_horizon(){
    return this->time_horizon_;
}

int MPC::get_state_dim(){
    return this->state_dim_;
}

int MPC::get_control_dim(){
    return this->control_dim_;
}

void MPC::set_path_constraints(vector<double> state_path_constraints_min,vector<double> state_path_constraints_max){
    this->state_path_constraints_min_=state_path_constraints_min;
    this->state_path_constraints_max_=state_path_constraints_max;
}

void MPC::set_state_constraints(vector<double> state_terminal_constraints_min,vector<double> state_terminal_constraints_max){
    this->state_terminal_constraints_min_=state_terminal_constraints_min;
    this->state_terminal_constraints_max_=state_terminal_constraints_max;
}

void MPC::set_input_constraints(vector<double> input_constraints_min,vector<double> input_constraints_max){
    this->input_constraints_min_=input_constraints_min;
    this->input_constraints_max_=input_constraints_max;


    // ToDo
    u_min =  this->input_constraints_min_;
    u_max  = this->input_constraints_max_;

    ROS_INFO("Bounds and initial guess for the state");

    x_min  = this->state_path_constraints_min_;
    x_max  = this->state_path_constraints_max_;
    xf_min = this->state_terminal_constraints_min_;
    xf_max = this->state_terminal_constraints_max_;
}

Eigen::MatrixXd MPC::mpc_step(const geometry_msgs::Pose pose, const KDL::JntArray& state)
{
    // Bounds and initial guess for the control
#ifdef __DEBUG__
    ROS_INFO_STREAM("input_constraints_min_: " <<this->input_constraints_min_.size());
    ROS_INFO_STREAM("input_constraints_max_: " <<this->input_constraints_max_.size());
#endif

    //ROS_INFO("ODE right hand side and quadrature");
    SX qdot = SX::vertcat({u_});

    //ROS_INFO("Current Quaternion and Position Vector.");

    double kappa = 0.001; // Small regulation term for numerical stability for the NLP

     SX q_c = SX::vertcat({
         0.5 * sqrt(fk_(0,0) + fk_(1,1) + fk_(2,2) + 1.0 + kappa),
         0.5 * (sign((fk_(2,1) - fk_(1,2)))) * sqrt(fk_(0,0) - fk_(1,1) - fk_(2,2) + 1.0 + kappa),
         0.5 * (sign((fk_(0,2) - fk_(2,0)))) * sqrt(fk_(1,1) - fk_(2,2) - fk_(0,0) + 1.0 + kappa),
         0.5 * (sign((fk_(1,0) - fk_(0,1)))) * sqrt(fk_(2,2) - fk_(0,0) - fk_(1,1) + 1.0 + kappa)
 });

    pos_c = SX::vertcat({fk_(0,3), fk_(1,3), fk_(2,3)}); //current state

    //ROS_INFO("Desired Goal-pose");
    pos_target = SX::vertcat({pose.position.x, pose.position.y, pose.position.z});
    q_target = SX::vertcat({pose.orientation.w, pose.orientation.x, pose.orientation.y, pose.orientation.z});

    // Prevent collision with Base_link
    SX barrier;
    SX dist;

    // Get orientation error
    SX q_c_inverse = SX::vertcat({q_c(0), -q_c(1), -q_c(2), -q_c(3)});
    SX e_quat= quaternion_product(q_c_inverse,q_target);
    SX error_attitute = SX::vertcat({ e_quat(1), e_quat(2), e_quat(3)});

    // L2 norm of the control signal
//    SX R = 1*SX::vertcat({100, 100, 100, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1});
//    SX R = 1*SX::vertcat({0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1});
    SX R = 0.1;
    SX energy = dot(sqrt(R)*u_,sqrt(R)*u_);

    // L2 norm of the states
//    std::vector<int> state_convariance(state_dim_,1);
//    SX S = 0.01*SX::scalar_matrix(state_dim_,1,1);
//    SX motion = dot(sqrt(S)*x_,sqrt(S)*x_);

    //ROS_INFO("Objective");
    SX error=pos_c-pos_target;

    barrier = bv_.getOutputConstraints();

    SX L = 10*dot(pos_c-pos_target,pos_c-pos_target) + energy + 10 * dot(error_attitute,error_attitute) + barrier;
    SX phi = 10*dot(pos_c-pos_target,pos_c-pos_target) + 10 * dot(error_attitute,error_attitute) + barrier;

    //ROS_INFO("Create Euler integrator function");

    Function F = create_integrator(state_dim_, control_dim_, time_horizon_, num_shooting_nodes_, qdot, x_, u_, L);

    Function F_terminal = create_integrator(state_dim_, control_dim_, time_horizon_, num_shooting_nodes_, qdot, x_, u_, phi);


    // Offset in V
    int offset=this->init_shooting_node();

    //ROS_INFO("(Make sure that the size of the variable vector is consistent with the number of variables that we have referenced");
    casadi_assert(offset==NV);

    // Objective function
    MX J = 0;

    //Constraint function and bounds
    vector<MX> g;

    //ROS_INFO("Loop over shooting nodes");
    for(unsigned int k=0; k<num_shooting_nodes_; ++k)
    {
        // Create an evaluation node
        MXDict I_out = F( MXDict{ {"x0", X[k]}, {"p", U[k]} });

        // Save continuity constraints
        g.push_back( I_out.at("xf") - X[k+1] );

        // Add objective function contribution
        J += I_out.at("qf");
    }
    // Terminal cost
    // Create an evaluation node
    MXDict I_term = F_terminal( MXDict{ {"x0", X[num_shooting_nodes_-1]}, {"p", U[num_shooting_nodes_-1]} });
    J += I_term.at("qf");

    MXDict nlp = {{"x", V}, {"f", J}, {"g", vertcat(g)}};

    // Set options
    Dict opts;

    opts["ipopt.tol"] = 1e-4;
    opts["ipopt.max_iter"] = 20;
//    opts["ipopt.hessian_approximation"] = "limited-memory";
//    opts["ipopt.hessian_constant"] = "yes";
    opts["ipopt.linear_solver"] = "ma27";
    opts["ipopt.print_level"] = 0;
    opts["print_time"] = false;
    opts["expand"] = false;  // Removes overhead

    ros::Time time = ros::Time::now();

    Function solver = nlpsol("nlpsol", "ipopt", nlp, opts);

    std::map<std::string, DM> arg, res;

    arg["lbx"] = min_state;
    arg["ubx"] = max_state;
    arg["lbg"] = 0;
    arg["ubg"] = 0;
    arg["x0"] = init_state;


    res = solver(arg);

    ros::Time time_new = ros::Time::now();

    ROS_INFO_STREAM("NLP time: " << (time_new - time).toSec());
    // Optimal solution of the NLP
    vector<double> V_opt(res.at("x"));
    vector<double> J_opt(res.at("f"));

    Eigen::VectorXd q_dot = Eigen::VectorXd::Zero(state_dim_);

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

    // Safe optimal control sequence at time t_k and take it as inital guess at t_k+1

    bv_.plotBoundingVolumes(sx_x_new);

    KDL::Frame ef_pos = forward_kinematics(state);

    return q_dot;
}
int MPC::init_shooting_node()
{
    // Offset in V
    int offset=0;
    X.clear();
    U.clear();
    min_state.clear();
    max_state.clear();
    init_state.clear();

    //THIS CAN BE SERIOULSY OPTIMISED
    for(unsigned int k=0; k<num_shooting_nodes_; ++k)
    {
        //ROS_INFO("Local state");
        X.push_back( V.nz(Slice(offset,offset+state_dim_)));

        if(k==0)
        {
            min_state.insert(min_state.end(), x0_min.begin(), x0_min.end());
            max_state.insert(max_state.end(), x0_max.begin(), x0_max.end());
        }
        else
        {
            min_state.insert(min_state.end(), x_min.begin(), x_min.end());
            max_state.insert(max_state.end(), x_max.begin(), x_max.end());
        }
        init_state.insert(init_state.end(), x_init.begin(), x_init.end());
        offset += state_dim_;

        //ROS_INFO("Local control via shift initialization");
        U.push_back( V.nz(Slice(offset,offset+control_dim_)));
        min_state.insert(min_state.end(), u_min.begin(), u_min.end());
        max_state.insert(max_state.end(), u_max.begin(), u_max.end());

        init_state.insert(init_state.end(), u_init_.begin(), u_init_.end());
        offset += control_dim_;
    }

    //ROS_INFO("State at end");
    X.push_back(V.nz(Slice(offset,offset+state_dim_)));
    min_state.insert(min_state.end(), xf_min.begin(), xf_min.end());
    max_state.insert(max_state.end(), xf_max.begin(), xf_max.end());
    init_state.insert(init_state.end(), u_init_.begin(), u_init_.end());

    x0_min.clear();
    x0_max.clear();
    x_init.clear();
    return offset += state_dim_;
}

KDL::Frame MPC::forward_kinematics(const KDL::JntArray& state){

    KDL::Frame ef_pos; //POsition of the end effector

    //SX p_base = SX::vertcat({fk_base_(0,3), fk_base_(1,3), fk_base_(2,3)});

    Function fk = Function("fk_", {x_}, {pos_c});
    //Function fk_base = Function("fk_base_", {x_}, {p_base});
    vector<double> x;
    for(unsigned int i=0; i < state.rows();i++)
    {
        x.push_back(state(i));
    }
    SX test_v = fk(SX::vertcat({x})).at(0);
    //SX test_base = fk_base(SX::vertcat({x_new})).at(0);

    ef_pos.p.x((double)test_v(0));
    ef_pos.p.y((double)test_v(1));
    ef_pos.p.z((double)test_v(2));
    ROS_WARN_STREAM("Current Position: \n" << ef_pos.p.x() << " "<< ef_pos.p.y() << " "<< ef_pos.p.z() << " ");
    return ef_pos;
}

Function MPC::create_integrator(const unsigned int state_dim, const unsigned int control_dim, const double T,
                                            const unsigned int N, SX ode, SX x, SX u, SX L)
{
    // Euler discretize
    double dt = T/((double)N);

    Function f = Function("f", {x, u}, {ode, L});

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

SX MPC::quaternion_product(SX q1, SX q2)
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

void MPC::setBoundingVolumes(BoundingVolume &bv)
{
    bv_ = bv;
}

void MPC::setForwardKinematics(ForwardKinematics &fk)
{
    _fk_ = fk;
}
