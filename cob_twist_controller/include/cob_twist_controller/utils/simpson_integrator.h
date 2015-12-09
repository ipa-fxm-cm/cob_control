/*!
 *****************************************************************
 * \file
 *
 * \note
 *   Copyright (c) 2015 \n
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
 *   Author: Felix Messmer, email: felix.messmer@ipa.fraunhofer.de
 *
 * \date Date of creation: November, 2015
 *
 * \brief
 *   This Class provides a helper for integrating velocities using Simpson method
 *
 ****************************************************************/

#ifndef COB_TWIST_CONTROLLER_UTILS_SIMPSON_INTEGRATOR_H
#define COB_TWIST_CONTROLLER_UTILS_SIMPSON_INTEGRATOR_H

#include <vector>

#include <ros/ros.h>
#include <kdl/jntarray.hpp>

#include "cob_twist_controller/utils/moving_average.h"

class SimpsonIntegrator
{
    public:
        explicit SimpsonIntegrator(const uint8_t dof)
        {
            dof_ = dof;
            for (uint8_t i = 0; i < dof_; i++)
            {
                // ma_.push_back(new MovingAvgSimple_double_t(3));
                ma_measured_vel_.push_back(new MovingAvgExponential_double_t(0.4));
                ma_measured_pos_.push_back(new MovingAvgExponential_double_t(0.4));
            }

            last_update_time_ = ros::Time(0.0);
            last_period_ = ros::Duration(0.0);
        }
        ~SimpsonIntegrator()
        {}


        void resetIntegration()
        {
            // resetting outdated values
            vel_last_.clear();
            vel_before_last_.clear();

            // resetting moving average
            for (unsigned int i = 0; i < dof_; ++i)
            {
                ma_measured_vel_[i]->reset();
                ma_measured_pos_[i]->reset();
            }

        }

        bool updateIntegration(const KDL::JntArray& q_dot_ik,
                               const KDL::JntArray& current_q,
                               std::vector<double>& pos,
                               std::vector<double>& vel)
        {
            ros::Time now = ros::Time::now();
            ros::Duration period = now - last_update_time_;

            std::vector<double> pos_filtered, vel_filtered;


            bool value_valid = false;
            pos.clear();
            vel.clear();
            pos_filtered.clear();
            vel_filtered.clear();

            // ToDo: Test these conditions and find good thresholds
            // if (period.toSec() > 2*last_period_.toSec())  // missed about a cycle
            if (period.toSec() > ros::Duration(0.5).toSec())  // missed about 'max_command_silence'
            {
                ROS_WARN_STREAM("reset Integration: " << period.toSec());
                resetIntegration();
            }


            for (unsigned int i = 0; i < dof_; ++i)
            {
                double filtered_q_dot_ik = 0.0;
                double filtered_q_ik = 0.0;

                ma_measured_vel_[i]->addElement(q_dot_ik(i));
                ma_measured_pos_[i]->addElement(current_q(i));

                if (ma_measured_vel_[i]->calcMovingAverage(filtered_q_dot_ik))
                {
                    vel_filtered.push_back(filtered_q_dot_ik);
                }

                if (ma_measured_pos_[i]->calcMovingAverage(filtered_q_ik))
                {
                    pos_filtered.push_back(filtered_q_ik);
                }
            }

            if (!vel_before_last_.empty())
            {
                for (unsigned int i = 0; i < dof_; ++i)
                {
                    // Simpson
                    double integration_value = static_cast<double>(period.toSec() / 6.0 * (vel_before_last_[i] + 4.0 * (vel_before_last_[i] + vel_last_[i]) + vel_before_last_[i] + vel_last_[i] + vel_filtered[i]) + pos_filtered[i]);

                    pos.push_back(integration_value);
                    vel.push_back(vel_filtered[i]);

                }
                value_valid = true;
            }

            // Continuously shift the vectors for simpson integration
            vel_before_last_.clear();
            for (unsigned int i=0; i < vel_last_.size(); ++i)
            {
                vel_before_last_.push_back(vel_last_[i]);
            }

            vel_last_.clear();
            for (unsigned int i=0; i < vel_filtered.size(); ++i)
            {
                vel_last_.push_back(vel_filtered[i]);
            }

            last_update_time_ = now;
            last_period_ = period;

            return value_valid;
        }

    private:
        std::vector<MovingAvgBase_double_t*> ma_measured_vel_;
        std::vector<MovingAvgBase_double_t*> ma_measured_pos_;
        uint8_t dof_;
        std::vector<double> vel_last_, vel_before_last_;
        ros::Time last_update_time_;
        ros::Duration last_period_;
};

#endif  // COB_TWIST_CONTROLLER_UTILS_SIMPSON_INTEGRATOR_H
