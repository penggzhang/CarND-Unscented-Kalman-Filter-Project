#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 6;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  // initially set to false
  is_initialized_ = false;
  
  // state dimension
  n_x_ = 5;
  
  // augmented state dimension
  n_aug_ = 7;

  // sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
    
  // weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
    
  // Set weights
  double w = 0.5 / (lambda_ + n_aug_);
  weights_.fill(w);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
    
  // Laser noise covariance matrix
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;
    
  // Radar noise covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
    
  if (!is_initialized_) {
    cout << "EKF: " << endl;
    
    // Set px, py, v, yaw and yawd as zeros
    x_ << 0, 0, 0, 0, 0;
      
    // Initialize state variables px and py with first measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates
      float ro = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      x_(0) = ro * cos(phi);
      x_(1) = ro * sin(phi);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
      
    // Initialize covariance
    P_ << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 25, 0, 0,
      0, 0, 0, 10, 0,
      0, 0, 0, 0, 1;

    // Set flag
    is_initialized_ = true;
      
    // Set timestamp
    time_us_ = meas_package.timestamp_;
      
    // Done initializing, no need to predict or update
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
    
  // Compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
    
  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
    
  // Use the sensor type to perform the update step.
  // Update the state and covariance matrices.
  if (use_radar_ && (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
    // Radar update
    UpdateRadar(meas_package);
  } else if (use_laser_ && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
    // Laser update
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    
  /*******************************************
    Augment state and covariance for noise
  ********************************************/

  // Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
    
  // Create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_aug_ - 2) = 0;
  x_aug(n_aug_ - 1) = 0;
    
  // Create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_aug_-2, n_aug_-2) = std_a_ * std_a_;
  P_aug(n_aug_-1, n_aug_-1) = std_yawdd_ * std_yawdd_;

  /*******************************************
    Generate augmented sigma points
  ********************************************/

  // Create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
  // Calculate square root matrix
  MatrixXd A = P_aug.llt().matrixL();
    
  // Generate augmented sigma points
  Xsig_aug.col(0) = x_aug;
  double spreading_factor = sqrt(lambda_ + n_aug_);
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + spreading_factor * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - spreading_factor * A.col(i);
  }

  /*******************************************
    Predict sigma points
  ********************************************/
  
  // Create vector for predicted sigma point
  VectorXd pred_pt = VectorXd(n_x_);
    
  // Iterate over all sigma points
  // and make prediction for each point
  for (int i = 0; i < 2*n_aug_+1; i++) {
        
    // Extract values of sigma point for better readability
    const double px = Xsig_aug(0, i);
    const double py = Xsig_aug(1, i);
    const double v = Xsig_aug(2, i);
    const double yaw = Xsig_aug(3, i);
    const double yawd = Xsig_aug(4, i);
    const double nu_a = Xsig_aug(5, i);
    const double nu_yawdd = Xsig_aug(6, i);
        
    // Pre-compute to avoid repeated calculations
    const double delta_t_2 = delta_t * delta_t;
    const double sin_yaw = sin(yaw);
    const double cos_yaw = cos(yaw);
    const double noise_px = 0.5 * nu_a * cos_yaw * delta_t_2;
    const double noise_py = 0.5 * nu_a * sin_yaw * delta_t_2;
        
    // Predict state vector and add noise
    // Avoid division by zero
    if (fabs(yawd) > 0.001) {
        pred_pt(0) = px + v / yawd * (sin(yaw + yawd * delta_t) - sin_yaw) + noise_px;
        pred_pt(1) = py + v / yawd * (-cos(yaw + yawd * delta_t) + cos_yaw) + noise_py;
    } else {
        pred_pt(0) = px + v * cos_yaw * delta_t + noise_px;
        pred_pt(1) = py + v * sin_yaw * delta_t + noise_py;
    }
        pred_pt(2) = v + nu_a * delta_t;
        pred_pt(3) = yaw + yawd * delta_t + 0.5 * nu_yawdd * delta_t_2;
        pred_pt(4) = yawd + nu_yawdd * delta_t;
        
    // Write predicted sigma point into right column
    Xsig_pred_.col(i) = pred_pt;
  }

  /*******************************************
    Calculate posterior state and covariance
  ********************************************/
    
  // Calculate state mean
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
    
  // Calculate covariance
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //iterate over sigma points
    // Residual
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
      
    // Angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    P_ += weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    
  /***************************************
    Predict measurement
  ****************************************/
    
  // Measurement dimension
  int n_z = 2;
    
  // Create vector for predicted measurement
  VectorXd z_pred = VectorXd(n_z);
    
  // Create measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
    
  // Create transforming matrix
  MatrixXd H = MatrixXd(n_z, n_x_);
  H << 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0;
    
  // Map posterior prediction into measurement space
  z_pred = H * x_;
    
  // Calculate measurement covariance
  S = H * P_ * H.transpose();
    
  // Add measurement noise
  S += R_laser_;
    
  /************************************
    Update
  *************************************/
    
  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
    
  // Calculate cross correlation matrix
  Tc = P_ * H.transpose();

  // Calculate Kalman gain;
  MatrixXd K = Tc * S.inverse();
    
  // Residual between actual and predicted measurement
  VectorXd z_error = meas_package.raw_measurements_ - z_pred;
    
  // Angle normalization
  while (z_error(1)> M_PI) z_error(1)-=2.*M_PI;
  while (z_error(1)<-M_PI) z_error(1)+=2.*M_PI;
    
  // Update state mean and covariance matrix
  x_ += K * z_error;
  P_ -= K * S * K.transpose();
    
  /************************************
     Calculate the radar NIS
  *************************************/
  
  NIS_laser_ = z_error.transpose() * S.inverse() * z_error;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /***************************************
    Predict measurement
  ****************************************/

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(3, 2*n_aug_+1);
    
  // Create vector for predicted measurement
  VectorXd z_pred = VectorXd(3);
    
  // Create measurement covariance matrix S
  MatrixXd S = MatrixXd(3, 3);
    
  // Define a small value to avoid division by zero
  const double eps = 0.0001;
    
  // Transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    // Create variables for better readability
    const double px = Xsig_pred_(0, i);
    const double py = Xsig_pred_(1, i);
    const double v = Xsig_pred_(2, i);
    const double yaw = Xsig_pred_(3, i);
        
    // Predict measurement
    const double radr = sqrt(px * px + py * py);
    double radphi;
    if (fabs(px)>=eps && fabs(py)>=eps) {
      radphi = atan2(py, px);
    } else {
      radphi = 0;
    }
    const double radrd = (px * v * cos(yaw) + py * v * sin(yaw)) / (std::max(eps, radr));
        
    // Set corresponding values of Zsig
    Zsig(0, i) = radr;
    Zsig(1, i) = radphi;
    Zsig(2, i) = radrd;
  }
    
  // Calculate predicted measurement mean
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //iterate over sigma points
    z_pred += weights_(i) * Zsig.col(i);
  }
    
  // Calculate measurement covariance
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //iterate over sigma points
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
      
    // Angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    S += weights_(i) * z_diff * z_diff.transpose();
  }
    
  // Add measurement noise
  S += R_radar_;
    
  /************************************
    Update
  *************************************/
    
  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 3);
    
  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) {
    // Residuals
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;
        
    // Angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Calculate Kalman gain;
  MatrixXd K = Tc * S.inverse();
    
  // Residual between actual and predicted measurement
  VectorXd z_error = meas_package.raw_measurements_ - z_pred;
    
  // Angle normalization
  while (z_error(1)> M_PI) z_error(1)-=2.*M_PI;
  while (z_error(1)<-M_PI) z_error(1)+=2.*M_PI;
    
  // Update state mean and covariance matrix
  x_ += K * z_error;
  P_ -= K * S * K.transpose();
    
  /************************************
     Calculate the radar NIS
  *************************************/
  
  NIS_radar_ = z_error.transpose() * S.inverse() * z_error;
}