#include <cob_nonlinear_mpc/bounding_volumes.h>

void BoundingVolume::addCollisionBall(const geometry_msgs::Point point, double radius, int id)
{
    visualization_msgs::Marker marker;
    marker.type = visualization_msgs::Marker::SPHERE;
    marker.lifetime = ros::Duration();
    marker.action = visualization_msgs::Marker::ADD;
    marker.ns = "preview";
    marker.header.frame_id = robot_.root_frame;

    marker.scale.x = 2*radius;
    marker.scale.y = 2*radius;
    marker.scale.z = 2*radius;

    marker.color.r = 1.0;
    marker.color.g = 0.0;
    marker.color.b = 0.0;
    marker.color.a = 0.1;

    marker.id = id;
    marker.pose.position.x = point.x;
    marker.pose.position.y = point.y;
    marker.pose.position.z = point.z;
    marker_array_.markers.push_back(marker);

}

void BoundingVolume::plotBoundingVolumes(SX x_current)
{

    SX x = getForwardKinematic().getX();
    SX result;
    std::unordered_map<std::string, std::vector<std::string> >::iterator it_scm;
    geometry_msgs::Point point;
    double bv_radius = 0.1;
    double id = 0;
    marker_array_.markers.clear();

    for( it_scm = robot_.self_collision_map_.begin(); it_scm != robot_.self_collision_map_.end(); it_scm++)
    {
        vector<string> tmp = it_scm->second;

        for(int i=0; i<tmp.size();i++)
        {
            vector<vector<SX>> SX_vec = bv_mat[tmp.at(i)];
            for(int k=0; k<SX_vec.size(); k++)
            {
                SX sym_value = SX::horzcat({SX_vec.at(k).at(0)});
                Function sym_function = Function("sym_value", {x}, {sym_value});
                result = sym_function(x_current).at(0);
                point.x = (double)result(0);
                point.y = (double)result(1);
                point.z = (double)result(2);

                bv_radius = 0.12;
                addCollisionBall(point, bv_radius, id);
                id++;
            }
        }


        vector<vector<SX>> SX_vec = bv_mat[it_scm->first];
        for(int k=0; k<SX_vec.size(); k++)
        {
            SX sym_value = SX::horzcat({SX_vec.at(k).at(0)});
            Function sym_function = Function("sym_value", {x}, {sym_value});
            result = sym_function(x_current).at(0);
            point.x = (double)result(0);
            point.y = (double)result(1);
            point.z = (double)result(2);

            if(it_scm->first == "body")
            {
                bv_radius = bvb_radius_.at(k);
            }
            else
            {
                bv_radius = 0.12;
            }
            addCollisionBall(point, bv_radius, id);
            id++;
        }
    }

    marker_pub_.publish(marker_array_);
}


void BoundingVolume::generate_bounding_volumes()
{
    std::vector<SX> fk_vector = fk_.getFkVector();
    std::vector<T_BVH> transform_vector = fk_.getTransformVector();
    SX fk_base = fk_.getFkBase();

    // Get bounding volume forward kinematics
    for(int i=0; i<transform_vector.size(); i++)
    {
        T_BVH bvh = transform_vector.at(i);
        std::vector<SX> bvh_arm;
        if(i-1<0)
        {
            SX transform = mtimes(fk_vector.at(i),bvh.T);
            SX tmp = SX::vertcat({transform(0,3), transform(1,3), transform(2,3)});
            bvh_arm.push_back(tmp);
            bv_mat[bvh.link].push_back(bvh_arm);

            if(bvh.constraint)
            {
                bvh_arm.clear();
                transform = mtimes(fk_vector.at(i),bvh.BVH_p);
                tmp = SX::vertcat({transform(0,3), transform(1,3), transform(2,3)});
                bvh_arm.push_back(tmp);
                bv_mat[bvh.link].push_back(bvh_arm);
            }
        }
        else
        {
            bvh_arm.clear();
            SX transform = mtimes(fk_vector.at(i-1),bvh.T);
            SX tmp = SX::vertcat({transform(0,3), transform(1,3), transform(2,3)});
            bvh_arm.push_back(tmp);
            bv_mat[bvh.link].push_back(bvh_arm);
            bvh_arm.clear();

            if(bvh.constraint)
            {
                tmp.clear();
                transform = mtimes(fk_vector.at(i-1),bvh.BVH_p);
                tmp = SX::vertcat({transform(0,3), transform(1,3), transform(2,3)});
                bvh_arm.push_back(tmp);
                bv_mat[bvh.link].push_back(bvh_arm);
            }
        }
    }
    if(robot_.base_active_)
    {
        for(int i=0; i<bvb_positions_.size(); i++)
        {
            std::vector<SX> base_bvh;
            SX tmp = SX::vertcat({fk_base(0,3), fk_base(1,3), fk_base(2,3)+bvb_positions_.at(i)});
            base_bvh.push_back(tmp);
            bv_mat["body"].push_back(base_bvh);
        }
    }
}

SX BoundingVolume::getOutputConstraints()
{
    std::unordered_map<std::string, std::vector<std::string> >::iterator it_scm;
    SX dist;
    SX barrier;
    int counter = 0;
    double bv_radius;


    for( it_scm = robot_.self_collision_map_.begin(); it_scm != robot_.self_collision_map_.end(); it_scm++)
    {
        std::vector<string> scm_collision_links = it_scm->second;
        for(int i=0; i<scm_collision_links.size(); i++)
        {

            vector<vector<SX>> p1_mat = bv_mat[it_scm->first];
            vector<vector<SX>> p2_mat = bv_mat[scm_collision_links.at(i)];


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

    return barrier;
}


bool BoundingVolume::setBVBpositions(vector<double> bvb_pos)
{
    bvb_positions_ = bvb_pos;
    return true;
}

bool BoundingVolume::setBVBradius(vector<double> bvb_rad)
{
    bvb_radius_ = bvb_rad;
    return true;
}

void BoundingVolume::setRobot(Robot robot)
{
    robot_ = robot;
}

Robot BoundingVolume::getRobot()
{
    return robot_;
}

void BoundingVolume::setForwardKinematic(ForwardKinematics fk)
{
    fk_ = fk;
}

ForwardKinematics BoundingVolume::getForwardKinematic()
{
    return fk_;
}
