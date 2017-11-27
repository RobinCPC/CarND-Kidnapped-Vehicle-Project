/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    // set number of particles
    this->num_particles = 25;

    // use random std generator
    default_random_engine gen;

    //create normal Gaussian distribution for each particle
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    // initialize particles and their weights
    for (int i=0; i < num_particles; ++i){
        Particle p;
        p.id = i, p.x=dist_x(gen), p.y=dist_y(gen), p.theta=dist_theta(gen);
        p.weight = 1.;
        this->particles.push_back(p);
        weights.push_back(1.0);
    }

    this->is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    // use random std generator
    default_random_engine gen;

    for(int i=0; i < this->num_particles; ++i){
        double x = particles[i].x;
        double y = particles[i].y;
        double theta = particles[i].theta;
        if(fabs(yaw_rate) > 0){
            x = x + (velocity/yaw_rate)*( sin(theta + yaw_rate*delta_t) - sin(theta));
            y = y + (velocity/yaw_rate)*(-cos(theta + yaw_rate*delta_t) + cos(theta));
            theta += yaw_rate * delta_t;
        }else{
            x = x + velocity * delta_t * cos(theta);
            y = y + velocity * delta_t * sin(theta);
        }
            
        //create normal Gaussian distribution for each particle
        normal_distribution<double> noise_x(0, std_pos[0]);
        normal_distribution<double> noise_y(0, std_pos[1]);
        normal_distribution<double> noise_theta(0, std_pos[2]);

        // add measurement back with noise TODO: add noise
        particles[i].x = x + noise_x(gen);
        particles[i].y = y + noise_y(gen);
        particles[i].theta = theta + noise_theta(gen);
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    for(size_t i=0; i < this->particles.size(); ++i){
        vector<int> assoc;
        vector<double> s_x;
        vector<double> s_y;
        // position particle w.r.t  map
        double x_p = particles[i].x;
        double y_p = particles[i].y;
        double theta_p = particles[i].theta;
        double weight = 1.0;
        double gauss_norm = (1. / (2 * M_PI * std_landmark[0] * std_landmark[1]));

        // Data Transformation and Associations
        for(size_t j=0; j < observations.size(); ++j){
            // position of observation w.r.t. vehicle 
            double x_c = observations[j].x;
            double y_c = observations[j].y;
            // position of observation w.r.t. map
            double x_m = x_p + (cos(theta_p)*x_c) - (sin(theta_p) * y_c);
            double y_m = y_p + (sin(theta_p)*x_c) + (cos(theta_p) * y_c);

            // TODO: rewrite with more simplier way (should get landmark index not id_i)
            vector< pair<double, int>> dist_mark;
            // user nearest neighbor to find associations landmark in the map
            for(auto mark : map_landmarks.landmark_list){
                double distance = dist(x_m, y_m, mark.x_f, mark.y_f);
                if (distance < sensor_range)
                    dist_mark.push_back(make_pair(distance, mark.id_i));
            }
            vector<pair<double, int>>::iterator res_it = std::min_element(dist_mark.begin(), dist_mark.end());
            int cls_mark_id = dist_mark[ std::distance(dist_mark.begin(), res_it) ].second;
            
            // assign mark_id and sens_xy
            assoc.push_back(cls_mark_id);
            s_x.push_back(x_m);
            s_y.push_back(y_m);
            
            // calculate weight of each obs (Multivariable-Gaussian Probobility)
            double mu_x = map_landmarks.landmark_list[cls_mark_id].x_f;
            double mu_y = map_landmarks.landmark_list[cls_mark_id].y_f;
            double exponent = (pow(x_m - mu_x,2)) / (2 * pow(std_landmark[0],2)) + 
                              (pow(y_m - mu_y,2)) / (2 * pow(std_landmark[1],2));
            weight *= gauss_norm * exp(-exponent);
            
        }
        // update particle element (weight)
        particles[i].weight = weight;
        particles[i].associations = assoc;
        particles[i].sense_x = s_x;
        particles[i].sense_y = s_y;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    double total_sum = 0.;
    // clear weights
    this->weights.empty();
    for (auto p : particles)
        total_sum += p.weight;
    std::vector<int> weights_list;
    for(auto p : particles){
        this->weights.push_back(p.weight);
        weights_list.push_back( (int)(100*(p.weight)/total_sum) );
    }

    // use random std generator
    default_random_engine gen;
    //std::random_device rd;
    //std::mt19937 gen(rd());

    std::discrete_distribution<int> dis_d( weights.begin(), weights.end());

    vector<Particle> neo_ps;
    for(int i = 0; i < this->num_particles; ++i){
        neo_ps.push_back( particles[dis_d(gen)] );
    }
    this->particles = neo_ps;


}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
