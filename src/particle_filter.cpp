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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    /*
     TODO: Set the number of particles. Initialize all particles to first position (based on estimates of x, y, theta and their uncertainties from GPS) and all weights to 1.
     TODO: Add random Gaussian noise to each particle.
     NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    */
    num_particles = 1000;  // TODO: Set the number of particles
    
    // initialization variables
    std::default_random_engine gen;
    double std_x = std[0];
    double std_y = std[1];
    double std_theta = std[2];
    
    // initialization of normal distributions
    std::normal_distribution<double> dist_x(x, std_x);
    std::normal_distribution<double> dist_y(y, std_y);
    std::normal_distribution<double> dist_theta(theta, std_theta);
      
    // initialization of particles
    for (int i = 0; i < num_particles; i++) {
        Particle particle = Particle();
        particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1;

        particles.push_back(particle);
    }
    
    // set particle filtes initialized
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    
    // initialization variables
    std::default_random_engine gen;
    double std_x = std_pos[0];
    double std_y = std_pos[1];
    double std_theta = std_pos[2];
    
    // Predict each particle
    for (int i = 0; i < num_particles; i++) {
        // calculate predicted location
        double theta = particles[i].theta;
        double x_predicted = particles[i].x + velocity/yaw_rate * (sin(theta + yaw_rate*delta_t) - sin(theta));
        double y_predicted = particles[i].y + velocity/yaw_rate * (cos(theta) - cos(theta + yaw_rate*delta_t));
        double theta_predicted = theta + yaw_rate*delta_t;
        
        // add noise
        std::normal_distribution<double> dist_x(x_predicted, std_x);
        std::normal_distribution<double> dist_y(y_predicted, std_y);
        std::normal_distribution<double> dist_theta(theta_predicted, std_theta);
        
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    for (int i = 0; i < observations.size(); i++) {
        LandmarkObs obs = observation[i];
        int map_id = 0;
        double min_distance = -1;
        
        for (int j = 0; j < predicted.size(); j++) {
            LandmarkObs pred = predicted[i];
            double distance = dist(obs.x,obs.y,pred.x,pred.y);
            
            if (min_distance < 0 or min_distance > distance) {
                min_distance = distance;
                map_id = pred.id;
            }
        }
        observations[i].id = map_id;
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    
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
        
        for (int j = 0; j < transformed_obs.size(); j++) {
            double obs_x, obs_y, pred_x, pred_y;
            obs_x = transformed_obs[j].x;
            obs_y = transformed_obs[j].y;
            
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
        weights[i] = particles[i].weight;
    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
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
