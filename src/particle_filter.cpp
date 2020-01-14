/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

using namespace std;

static default_random_engine generator;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // number of particles
    num_particles = 100;
    
    // initialization of normal distributions
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
      
    // initialization of particles
    for (int i = 0; i < num_particles; i++) {
        Particle particle;
        particle.id = i;
        particle.x = dist_x(generator);
        particle.y = dist_y(generator);
        particle.theta = dist_theta(generator);
        particle.weight = 1.0;

        particles.push_back(particle);
        weights.push_back(1.0);
    }
    
    // set particle filtes initialized
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // initialization variables
    double p_x, p_y, theta, std_x, std_y, std_theta, new_x, new_y, new_theta;
    std_x = std_pos[0];
    std_y = std_pos[1];
    std_theta = std_pos[2];
    
    // Predict each particle
    for (int i = 0; i < num_particles; i++) {
        p_x = particles[i].x;
        p_y = particles[i].y;
        theta = particles[i].theta;
        
        // if yaw rate is almost 0, ignore it
        if (fabs(yaw_rate) < 0.00001) {
            new_x = p_x + (velocity * delta_t * cos(theta));
            new_y = p_y + (velocity * delta_t * sin(theta));
        } else {
            new_x = p_x + (velocity/yaw_rate) * (sin(theta + (yaw_rate*delta_t)) - sin(theta));
            new_y = p_y + (velocity/yaw_rate) * (cos(theta) - cos(theta + (yaw_rate*delta_t)));
            new_theta = theta + yaw_rate*delta_t;
        }
        
        // calculate predicted location
        normal_distribution<double> dist_x(new_x, std_x);
        normal_distribution<double> dist_y(new_y, std_y);
        normal_distribution<double> dist_theta(new_theta, std_theta);
        
        // add noise
        particles[i].x = dist_x(generator);
        particles[i].y = dist_y(generator);
        particles[i].theta = dist_theta(generator);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
    // for each observation find predicted landmark
    for (int i = 0; i < observations.size(); i++) {
        // init some variables
        LandmarkObs obs = observations[i];
        int map_id = -1;
        double min_distance = 1000000000000000000;
        
        // find nearest landmark (nearest neighbor)
        for (int j = 0; j < predicted.size(); j++) {
            LandmarkObs pred = predicted[j];
            double distance = dist(obs.x,obs.y,pred.x,pred.y);
            
            if (min_distance > distance) {
                min_distance = distance;
                map_id = pred.id;
            }
        }
        observations[i].id = map_id;
    }

}

double calculateMultivariateGaussian(double x, double y, double mu_x, double mu_y, double sigma_x, double sigma_y) {
    // Calculate multivariate 2D gaussian
    return (1/(2*M_PI*sigma_x*sigma_y)) * exp((-1.0/2) * (pow(x-mu_x,2) / (2*pow(sigma_x, 2)) + (pow(y-mu_y,2) / (2*pow(sigma_y, 2)))));
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], const vector<LandmarkObs> &observations, const Map &map_landmarks) {
    // update each particles weight
    for (int i = 0; i < num_particles; i++) {
        Particle p = particles[i];
        vector<LandmarkObs> predicted;
        
        // select landmarks within sensor range
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            double landmark_x = map_landmarks.landmark_list[j].x_f;
            double landmark_y = map_landmarks.landmark_list[j].y_f;
            int landmark_id = map_landmarks.landmark_list[j].id_i;
            
            double distance_landmark_particle = dist(p.x, p.y, landmark_x, landmark_y);
            
            if (distance_landmark_particle <= sensor_range) {
                predicted.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
            }
        }
        
        // transform observations to world coordinates
        vector<LandmarkObs> transformed_obs;
        for (int j = 0; j < observations.size(); j++) {
            // calculate x,y
            double trans_x = observations[j].x * cos(p.theta) - observations[j].y * sin(p.theta) + p.x;
            double trans_y = observations[j].x * sin(p.theta) + observations[j].y * cos(p.theta) + p.y;
            transformed_obs.push_back(LandmarkObs{observations[j].id, trans_x, trans_y});
        }
        
        // associate closest predicted landmark to observation
        dataAssociation(predicted, transformed_obs);
        
        //reset weight
        particles[i].weight = 1.0;
        
        // calculate weight of particle
        for (int j = 0; j < transformed_obs.size(); j++) {
            double obs_x, obs_y, pred_x, pred_y;
            obs_x = transformed_obs[j].x;
            obs_y = transformed_obs[j].y;
            
            // find associated landmark x and y
            int associated_pred_id = transformed_obs[j].id;
            for (int k = 0; k < predicted.size(); k++) {
                if (predicted[k].id == associated_pred_id) {
                    pred_x = predicted[k].x;
                    pred_y = predicted[k].y;
                    break;
                }
            }
            particles[i].weight *= calculateMultivariateGaussian(pred_x, pred_y, obs_x, obs_y, std_landmark[0], std_landmark[1]);
        }
        // update weight in weights array
        weights[i] = particles[i].weight;
    }
}

void ParticleFilter::resample() {
    // get maximal weight
    double max_weight = *max_element(weights.begin(), weights.end());
    
    // init uniform distributions
    uniform_int_distribution<> uniform_int(0,num_particles-1);
    uniform_real_distribution<double> uniform(0.0, 2*max_weight);
    
    // init variables
    vector<Particle> resambled_particles;
    double beta = 0;
    
    // select index by random
    int index = uniform_int(generator);
    
    // for each particle...
    for (int i = 0; i < num_particles; i++) {
        // choose particles randomly but weighted from particle array (particle with higher weight is more likely than a particle with lower weight)
        beta += uniform(generator);
        while (weights[index] < beta) {
            beta -= weights[index];
            index = (index+1) % num_particles;
        }
        resambled_particles.push_back(particles[index]);
    }
    // update particles array
    particles = resambled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, const vector<int>& associations, const vector<double>& sense_x, const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
