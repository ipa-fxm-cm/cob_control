#include <vector>
#include <ros/ros.h>
#include <std_msgs/Float64.h>
#include <kdl/jntarray.hpp>
#include <cob_twist_controller/utils/simpson_integrator.h>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include <termios.h>
#include <signal.h>

class SimpsonIntegratorTester
{
public:
    SimpsonIntegratorTester()
    {
        dof_ = 1;
        q_.resize(dof_);
        KDL::SetToZero(q_);
        q_dot_.resize(dof_);
        KDL::SetToZero(q_dot_);
        simpson_integrated_q_dot_.resize(dof_);
        KDL::SetToZero(simpson_integrated_q_dot_);
        euler_integrated_q_dot_.resize(dof_);
        KDL::SetToZero(euler_integrated_q_dot_);
        simpson_derived_q_dot_.resize(dof_);
        KDL::SetToZero(simpson_derived_q_dot_);



        integrator_.reset(new SimpsonIntegrator(dof_));
        output_sin_pub_ = nh_.advertise<std_msgs::Float64>("sin", 1);
        output_derived_simpson_pub_ = nh_.advertise<std_msgs::Float64>("derived_simpson", 1);
        output_simpson_sin_pub_ = nh_.advertise<std_msgs::Float64>("simpson", 1);
        output_cos_pub_ = nh_.advertise<std_msgs::Float64>("cos", 1);
        output_euler_sin_pub_ = nh_.advertise<std_msgs::Float64>("euler", 1);
        stop_=false;
        time_ = 0;
    }


    ~SimpsonIntegratorTester()
    {}

    void run()
    {
        ros::Rate r(100.0);

//        boost::thread start_thread;
//        start_thread = boost::thread(boost::bind(&SimpsonIntegratorTester::stopIntegration, this));
//        ros::AsyncSpinner spinner(0);
//        spinner.start();
        ROS_INFO("Start integration \n Enter any key to stop it.");

        ros::Time time = ros::Time::now();
        ros::Time last_update_time = time;
        ros::Duration period = time - last_update_time;
        double old_pos = -99;
        euler_integrated_q_dot_(0) = -cos(time_);
        while(ros::ok())
        {
            time = ros::Time::now();
            period = time - last_update_time;

            std::vector<double> next_q;
            std::vector<double> next_q_dot;

            q_dot_(0) = sin(time_);
            q_(0)     = -cos(time_);

            if (integrator_->updateIntegration(q_dot_, q_, next_q, next_q_dot))
            {
                simpson_integrated_q_dot_(0) = next_q[0];
                q_dot_(0) = next_q_dot[0];
            }

            euler_integrated_q_dot_(0) += q_dot_(0) * period.toSec();

            if(old_pos != -99)
            {
                simpson_derived_q_dot_(0) = (euler_integrated_q_dot_(0) - old_pos)/period.toSec();
            }


            std_msgs::Float64 simpson_q_;
            simpson_q_.data = simpson_integrated_q_dot_(0);
            std_msgs::Float64 euler_q_;
            euler_q_.data = euler_integrated_q_dot_(0);
            std_msgs::Float64 real_q;
            real_q.data = q_(0);
            std_msgs::Float64 derived_q_;
            derived_q_.data = simpson_derived_q_dot_(0);
            std_msgs::Float64 real_sin;
            real_sin.data = q_dot_(0);


            output_simpson_sin_pub_.publish(simpson_q_);
            output_cos_pub_.publish(real_q);
            output_euler_sin_pub_.publish(euler_q_);
            output_euler_sin_pub_.publish(euler_q_);
            output_derived_simpson_pub_.publish(derived_q_);
            output_sin_pub_.publish(real_sin);

            time_+=period.toSec();
            last_update_time = time;
            old_pos = euler_integrated_q_dot_(0);
            ros::spinOnce();
            r.sleep();
        }
    }

//        tcsetattr(kfd, TCSANOW, &cooked);

    void stopIntegration()
    {
        c = 0x0;
        // get the console in raw mode
        tcgetattr(kfd, &cooked);
        memcpy(&raw, &cooked, sizeof(struct termios));
        raw.c_lflag &=~ (ICANON | ECHO);
        // Setting a new line, then end of file
        raw.c_cc[VEOL] = 1;
        raw.c_cc[VEOF] = 2;
        tcsetattr(kfd, TCSANOW, &raw);

        while(ros::ok())
        {
            if(read(kfd, &c, 1) < 0)
            {
                perror("read():");
                exit(-1);
            }
            if(c == 0x61)
            {
                stop_ = true;
                break;
            }
        }
    }


    ros::NodeHandle nh_;
    ros::Publisher output_sin_pub_;
    ros::Publisher output_derived_simpson_pub_;

    ros::Publisher output_simpson_sin_pub_;
    ros::Publisher output_cos_pub_;
    ros::Publisher output_euler_sin_pub_;

    KDL::JntArray q_;
    KDL::JntArray q_dot_;
    KDL::JntArray simpson_integrated_q_dot_;
    KDL::JntArray euler_integrated_q_dot_;
    KDL::JntArray simpson_derived_q_dot_;


    boost::shared_ptr<SimpsonIntegrator> integrator_;
    bool stop_;
    double time_;
    unsigned int dof_;

    /// For Keyboard commands
    char c;
    int kfd;
    struct termios cooked, raw;

};



int main(int argc, char **argv)
{
    ros::init(argc, argv, "test_simpson_integrator_node");

    SimpsonIntegratorTester sit;
    sit.run();
    ros::spin();
    return 0;
}
