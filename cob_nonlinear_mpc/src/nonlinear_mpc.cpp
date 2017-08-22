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

    for(int i=0; i<control_dim_; ++i)
    {
        u_apply_.push_back(0);
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

    previous_command.push_back(Eigen::VectorXd::Zero(state_dim_));
    previous_command.push_back(Eigen::VectorXd::Zero(state_dim_));

    x_ = _fk_.getX();
    u_ = _fk_.getU();

    // Get end-effector fk
    fk_ = _fk_.getFkVector().at(_fk_.getFkVector().size()-1);


    // Degree of interpolating polynomial
        int d = 2;
        p_order_ = d;

        // Size of the finite elements
        h_ = time_horizon_/(double)num_shooting_nodes_;

        // Choose collocation points
        vector<double> tau_root = collocation_points(d, "radau");
        tau_root.insert(tau_root.begin(), 0);

        // All Time points
        vector<vector<double> > T_mat((double)num_shooting_nodes_,vector<double>(d+1));

        for(int k=0; k<num_shooting_nodes_; k++)
        {
            for(int j=0; j<d+1; j++)
            {
                T_mat[k][j] = h_*(k + tau_root[j]);
            }
        }

        // Coefficients of the quadrature function
        vector<double> B(d+1);

        // Coefficients of the collocation equation
        vector<vector<double> > C(d+1,vector<double>(d+1));
        vector<vector<double> > C2(d+1,vector<double>(d+1));

        // Coefficients of the continuity equation
        vector<double> D(d+1);

        // For all collocation points
        for(int j=0; j<d+1; ++j)
        {
            // Construct Lagrange polynomials to get the polynomial basis at the collocation point
            Polynomial p = 1;
            for(int r=0; r<d+1; ++r)
            {
                if(r!=j)
                {
                    p *= Polynomial(-tau_root[r],1)/(tau_root[j]-tau_root[r]);
                }
            }

            // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
            D[j] = p(1.0);

            // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
            Polynomial dp = p.derivative();
            Polynomial dpp = dp.derivative();
            for(int r=0; r<d+1; ++r)
            {
                C[j][r] = dp(tau_root[r]);
                C2[j][r] = dpp(tau_root[r]);
            }
            Polynomial p_int = p.anti_derivative();
            B[j] = p_int(1.0);
        }

        C_ = C;
        C2_ = C2;
        B_ = B;
        D_ = D;
        T_mat_= T_mat;

        // Total number of NLP variables
        int nxd = num_shooting_nodes_ * (p_order_+1)*state_dim_;
        int nu = num_shooting_nodes_ * control_dim_;
        int nxf = state_dim_;
//        NV = control_dim_ + nxd + nu + nxf;
        NV = nxd + nu + nxf;
        V = MX::sym("V",NV);

    ROS_WARN_STREAM("MPC initialized");
}

void MPC::set_coordination_weights(vector<double> masses){
    ROS_INFO("Setting weights");
    ROS_INFO_STREAM("Mass size"<<masses.size()<< " control size " <<u_.size());
    R.resize(masses.size(),1);

    for(int i=0;i<masses.size();i++)
    {
        R(i)=u_(i)*masses.at(i);
    }
    ROS_INFO("Done Setting weights");
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
    SX t_sym = SX::sym("t", 1);
    SX u1_sym = SX::sym("U1",control_dim_);
    SX u2_sym = SX::sym("U2",control_dim_);

    SX a1_sym = SX::sym("A1",control_dim_);

    // Approximate acceleration
    SX accel = (u2_sym-u1_sym)/h_;
    Function accl = Function("accl", {u1_sym, u2_sym}, {accel});

    SX L_acc = 0.02*dot(accel, accel);

    SX jerk = (a1_sym-accel)/h_;
    SX L_jerk = 0.0003*dot(jerk, jerk);


    // Bounds and initial guess for the control
#ifdef __DEBUG__
    ROS_INFO_STREAM("input_constraints_min_: " <<this->input_constraints_min_.size());
    ROS_INFO_STREAM("input_constraints_max_: " <<this->input_constraints_max_.size());
#endif

    //ROS_INFO("ODE right hand side and quadrature");
    SX qdot = SX::vertcat({u_});

    //ROS_INFO("Current Quaternion and Position Vector.");

    q_c = SX::vertcat({ //current quaternion
        0.5 * sqrt(fk_(0,0) + fk_(1,1) + fk_(2,2) + 1.0),
        0.5 * ((fk_(2,1) - fk_(1,2))) / sqrt(fk_(0,0) + fk_(1,1) + fk_(2,2) + 1.0),
        0.5 * ((fk_(0,2) - fk_(2,0))) / sqrt(fk_(1,1) + fk_(2,2) + fk_(0,0) + 1.0),
        0.5 * ((fk_(1,0) - fk_(0,1))) / sqrt(fk_(2,2) + fk_(0,0) + fk_(1,1) + 1.0)
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
    SX q = SX::vertcat({e_quat(0), e_quat(1), e_quat(2),e_quat(3)});
    Function qk = Function("e_quat", {x_}, {q});
    vector<double> x;
    for(unsigned int i=0; i < state.rows();i++)
    {
        x.push_back(state(i));
    }
    SX test_v = qk(SX::vertcat({x})).at(0);

    //ROS_INFO("L2 norm of the control signal");
    this->acceleration_coordination(state);

    SX energy = SX::dot(R,u_);
    energy.print(std::cout);

    SX R2 = 0.005*SX::vertcat({100, 100, 100, 1,1,1,1,1,1,1});
    SX energy2 = dot(sqrt(R2)*u_,sqrt(R2)*u_);
//    SX energy2 = dot(sqrt(R2)*u_,sqrt(R2)*u_);

    SX error=pos_c-pos_target;

    barrier = bv_.getOutputConstraints();
//    ROS_INFO("Objective");
    SX L = energy2 + 10*dot(pos_c-pos_target,pos_c-pos_target) + 10 * dot(error_attitute,error_attitute) + barrier + L_jerk;

    SX phi = 0*dot(pos_c-pos_target,pos_c-pos_target) +   0*dot(error_attitute,error_attitute)+ 0*barrier;
    // Offset in V
    int offset=this->init_shooting_node();

    ROS_WARN("TEST");

//    ROS_INFO("(Make sure that the size of the variable vector is consistent with the number of variables that we have referenced");
    casadi_assert(offset==NV);

    Function F = Function("F", {t_sym, x_, u_, u1_sym ,u2_sym, a1_sym}, {qdot,L, accel}, {"t","x0","u", "u1", "u2", "a1"}, {"qdot","cost", "accel"});

    Function Ff = Function("Ff", {x_}, {phi});

    // Objective function
    MX J = 0;

    //Constraint function and bounds
    vector<MX> g;

    MX X_f;

    MX U0 = MX::vertcat({u_apply_});
    MX A0 = MX::vertcat({a0_});

    for(unsigned int k=0; k<num_shooting_nodes_; ++k)
    {
        for(unsigned int j=1; j<p_order_+1; j++)
        {
            // Expression for the state derivative at the collocation point
            MX xp_jk = 0;
            MX xpp_jk = 0;

            for(int r=0; r<p_order_+1; ++r)
            {
                xp_jk += C_[r][j]*X_[k][r];
                xpp_jk += C2_[r][j]*X_[k][r];
            }

            // Create an evaluation node
            MXDict I_out = F( MXDict{ {"t", T_mat_[k][j]}, {"x0", X_[k][j]}, {"u", U_[k]} , {"u1", U0}, {"u2", U_[k]} , {"a1", A0} });

            A0 = I_out.at("accel");

            // Save continuity constraints
            g.push_back( h_*I_out.at("qdot") - xp_jk );
            // Add objective function contribution
            J += B_[j] * I_out.at("cost")*h_;

        }

        if(k==num_shooting_nodes_-1)
        {
            X_f = X_[k][p_order_];
        }

        MX xf_k = 0;
        for (unsigned int r=0; r<p_order_+1; r++)
        {
            xf_k += D_[r]*X_[k][r];
        }
        g.push_back( X_[k+1][0] - xf_k );

        U0 = U_[k];
    }
//    MX Terminal_cost = Ff( X_f).at(0);
//    J += Terminal_cost;

//    ROS_INFO("NLP");
    MXDict nlp = {{"x", V}, {"f", J}, {"g", vertcat(g)}};

    // Set options
    Dict opts;

    opts["ipopt.tol"] = 1e-2;
    opts["ipopt.max_iter"] = 10;
//    opts["ipopt.hessian_approximation"] = "limited-memory";
//    opts["ipopt.hessian_constant"] = "yes";
    opts["ipopt.linear_solver"] = "ma27";
    opts["ipopt.print_level"] = 1;
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

//    ROS_INFO_STREAM("V_opt: " << V_opt);

    Eigen::VectorXd q_dot = Eigen::VectorXd::Zero(state_dim_);
    Eigen::VectorXd q_dot2 = Eigen::VectorXd::Zero(state_dim_);
    SX sx_x_new;
    u_open_loop_.clear();
    x_open_loop_.clear();
    x_new.clear();
    u_apply_.clear();

    for(int i=0; i<1; ++i)  // Copy only the first optimal control sequence
    {
        for(int j=0; j<control_dim_; ++j)
        {
            q_dot[j] = V_opt.at((p_order_+1)*state_dim_ + j);

            x_new.push_back(V_opt.at(j));
            u_apply_.push_back(q_dot[j]);
        }
    }



    a0_.clear();
    vector<double> temp;

    for(int i=0; i<control_dim_; ++i)
    {
        q_dot2[i] = V_opt.at((p_order_+1)*state_dim_*2 + i);
    }

    Eigen::VectorXd acc = Eigen::VectorXd::Zero(state_dim_);

    acc = (q_dot2 - q_dot)/h_;

    for(int i=0; i<control_dim_; ++i)
    {
        a0_.push_back(acc[i]);
    }
    sx_x_new = SX::vertcat({x_new});

    // Safe optimal control sequence at time t_k and take it as inital guess at t_k+1

    bv_.plotBoundingVolumes(sx_x_new);

    //KDL::Frame ef_pos = forward_kinematics(state);

    previous_command.push_back(q_dot);
    return q_dot;
}
int MPC::init_shooting_node()
{
//    ROS_INFO_STREAM("Init Nodes");
    // Offset in V
    int offset=0;
    X_.clear();
    U_.clear();
    min_state.clear();
    max_state.clear();
    init_state.clear();

    vector<vector<MX> > X(num_shooting_nodes_+1,vector<MX>(p_order_+1));
    vector<MX> U(p_order_+1);


    // U0 - for acceleration
//    U.push_back( V.nz(Slice(offset,offset+control_dim_)));
//    min_state.insert(min_state.end(), u_apply_.begin(), u_apply_.end());
//    max_state.insert(max_state.end(), u_apply_.begin(), u_apply_.end());
//    init_state.insert(init_state.end(), u_apply_.begin(), u_apply_.end());
//    offset += control_dim_;
    for(unsigned int k=0; k<num_shooting_nodes_; ++k)
    {
        for(unsigned int j=0; j<p_order_+1; j++)
        {
            // Local state
            X[k][j] =  V.nz(Slice(offset,offset+state_dim_));

            init_state.insert(init_state.end(), x_init.begin(), x_init.end());

            // Bounds
            if(k==0 && j==0)
            {
                min_state.insert(min_state.end(), x_init.begin(), x_init.end());
                max_state.insert(max_state.end(), x_init.begin(), x_init.end());
            }
            else    // Path constraints from now on
            {
                min_state.insert(min_state.end(), x_min.begin(), x_min.end());
                max_state.insert(max_state.end(), x_max.begin(), x_max.end());
            }

            offset += state_dim_;
        }

        // Local control via shift initialization
        U.push_back( V.nz(Slice(offset,offset+control_dim_)));
        min_state.insert(min_state.end(), u_min.begin(), u_min.end());
        max_state.insert(max_state.end(), u_max.begin(), u_max.end());
        init_state.insert(init_state.end(), u_init_.begin(), u_init_.end());

        offset += control_dim_;
    }

    ROS_INFO("State at end");

    // State at end
    X[num_shooting_nodes_][0] = V.nz(Slice(offset,offset+state_dim_));
    min_state.insert(min_state.end(), xf_min.begin(), xf_min.end());
    max_state.insert(max_state.end(), xf_max.begin(), xf_max.end());
    init_state.insert(init_state.end(), x_init.begin(), x_init.end());

    x0_min.clear();
    x0_max.clear();
    x_init.clear();

    U_ = U;
    X_ = X;

    return offset += state_dim_;
}

void MPC::acceleration_coordination(const KDL::JntArray& state){
    Eigen::VectorXd acc = previous_command.at(1)-previous_command.at(0);
    for(int i=0;i<control_dim_;i++){
        //R(i)=R(i)*acc[i]/time_horizon_;
    }
    ROS_INFO_STREAM("Acceleration: "<< acc);
    previous_command.pop_back();
}

KDL::Frame MPC::forward_kinematics(const KDL::JntArray& state){

    KDL::Frame ef_pos; //POsition of the end effector

    Function fk = Function("fk_", {x_}, {pos_c});

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
    ROS_INFO_STREAM("Joint values:" << x);
    ROS_WARN_STREAM("Current Position: \n" << ef_pos.p.x() << " "<< ef_pos.p.y() << " "<< ef_pos.p.z() << " ");
    //ROS_WARN_STREAM("Base Position: \n" << (double)test_base(0) << " "<< (double)test_base(1) << " "<< (double)test_base(2) << " ");
    ROS_WARN_STREAM("Target Position: \n" << pos_target);
    return ef_pos;
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
