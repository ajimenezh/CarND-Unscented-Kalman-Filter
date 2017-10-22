#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <mutex>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// Uncomment to generate debug info
#define _DEBUG 1
#define PI acos(-1.0)

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

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(5, 2 * 7 + 1);

  //predicted state covariance matrix
  P_cov_ = MatrixXd(5, 5);

  //predicted state mean
  x_mean_ = VectorXd(5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

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

  previous_timestamp_ = 0;

  is_initialized_ = false;

  isInCalc = false;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

std::mutex g_mutex;

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    if (!use_laser_ &&
            meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
        return;
    }

    if (!use_radar_ &&
            meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        return;
    }

    /*****************************************************************************
    *  Initialization
    ****************************************************************************/

    if (!is_initialized_) {

        // first measurement
        x_ = VectorXd(5);
        x_ << 1, 1, 0, 0, 0;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            float ro = meas_package.raw_measurements_(0);
            float phi = meas_package.raw_measurements_(1);
            float ro_dot = meas_package.raw_measurements_(2);
            x_(0) = ro*cos(phi);
            x_(1) = ro*sin(phi);
            x_(2) = ro_dot;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);
        }

        P_.fill(0.0);
        for (int i=0; i<5; i++) P_(i, i) = 1.0;

        previous_timestamp_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    g_mutex.lock();
    if (isInCalc) return;
    isInCalc = true;
    g_mutex.unlock();

    delta_t_ = (double)(meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = meas_package.timestamp_;

    if (_DEBUG)
    {
        std::cout<<"Measurement"<<std::endl;
        std::cout<<meas_package.sensor_type_<<std::endl;
        std::cout<<meas_package.timestamp_<<std::endl;
        std::cout<<delta_t_<<std::endl;
    }

    Prediction(delta_t_);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        UpdateRadar(meas_package);
    }
    else {
        UpdateLidar(meas_package);
    }

    g_mutex.lock();
    isInCalc = false;
    g_mutex.unlock();

    if (_DEBUG)
    {
        cout << "x_ = " << x_ << endl;
        cout << "P_ = " << P_ << endl;
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

    //set state dimension
    int n_x = 5;

    //set augmented dimension
    int n_aug = 7;

    //define spreading parameter
    double lambda = 3 - n_aug;

    SigmaPointPrediction(&Xsig_pred_);

    VectorXd weights = VectorXd(2*n_aug+1);

    weights(0) = lambda/(lambda+n_aug);
    for (int i=1; i<2*n_aug+1; i++) {
      weights(i) = 0.5/(n_aug+lambda);
    }

    //predicted state mean
    x_mean_.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
      x_mean_ = x_mean_ + weights(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P_cov_.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_mean_;

      //angle normalization
      while (x_diff(3) > PI) x_diff(3) -= 2*PI;
      while (x_diff(3) < -PI) x_diff(3) += 2*PI;

      P_cov_ = P_cov_ + weights(i) * x_diff * x_diff.transpose() ;
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    UpdateStateLidar(meas_package.raw_measurements_, &x_, &P_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    UpdateStateRadar(meas_package.raw_measurements_, &x_, &P_);
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //create augmented mean vector
  static VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  static MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  static MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(n_x) = x_;
  x_aug[n_x] = 0;
  x_aug[n_x+1] = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x, n_x) = P_;
  P_aug(n_x, n_x) = std_a_*std_a_;
  P_aug(n_x+1, n_x+1) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();

  //create augmented sigma points
    //set first column of sigma point matrix
  Xsig_aug.fill(0.0);
  Xsig_aug.col(0)  = x_aug;

  //set remaining sigma points
  for (int i = 0; i < n_aug; i++)
  {
    Xsig_aug.col(i+1)     = x_aug + sqrt(lambda+n_aug) * A_aug.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * A_aug.col(i);
  }

  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //create example sigma point matrix
  static MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  AugmentedSigmaPoints(&Xsig_aug);

  //create matrix with predicted sigma points as columns
  static MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  double delta_t = delta_t_; //time diff in sec

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  for (int i=0; i<2*n_aug+1; i++) {
      if (fabs(Xsig_aug(4, i)) < 1.0e-6) {
          Xsig_pred(0, i) = Xsig_aug(0, i) + Xsig_aug(2, i)*cos(Xsig_aug(3, i))*delta_t;
          Xsig_pred(1, i) = Xsig_aug(1, i) + Xsig_aug(2, i)*sin(Xsig_aug(3, i))*delta_t;
          Xsig_pred(2, i) = Xsig_aug(2, i) + 0;
          Xsig_pred(3, i) = Xsig_aug(3, i) + 0;
          Xsig_pred(4, i) = Xsig_aug(4, i) + 0;

          Xsig_pred(0, i) += 0.5*delta_t*delta_t*cos(Xsig_aug(3, i))*Xsig_aug(5, i);
          Xsig_pred(1, i) += 0.5*delta_t*delta_t*sin(Xsig_aug(3, i))*Xsig_aug(5, i);
          Xsig_pred(2, i) += delta_t*Xsig_aug(5, i);
          Xsig_pred(3, i) += 0.5*delta_t*delta_t*Xsig_aug(6, i);
          Xsig_pred(4, i) += delta_t*Xsig_aug(6, i);
      }
      else {
          Xsig_pred(0, i) = Xsig_aug(0, i) + Xsig_aug(2, i)/Xsig_aug(4, i)*
          (sin(Xsig_aug(3, i) + Xsig_aug(4, i)*delta_t) - sin(Xsig_aug(3, i)));
          Xsig_pred(1, i) = Xsig_aug(1, i) + Xsig_aug(2, i)/Xsig_aug(4, i)*
          (-cos(Xsig_aug(3, i) + Xsig_aug(4, i)*delta_t) + cos(Xsig_aug(3, i)));
          Xsig_pred(2, i) = Xsig_aug(2, i) + 0;
          Xsig_pred(3, i) = Xsig_aug(3, i) + Xsig_aug(4, i)*delta_t;
          Xsig_pred(4, i) = Xsig_aug(4, i) + 0;

          Xsig_pred(0, i) += 0.5*delta_t*delta_t*cos(Xsig_aug(3, i))*Xsig_aug(5, i);
          Xsig_pred(1, i) += 0.5*delta_t*delta_t*sin(Xsig_aug(3, i))*Xsig_aug(5, i);
          Xsig_pred(2, i) += delta_t*Xsig_aug(5, i);
          Xsig_pred(3, i) += 0.5*delta_t*delta_t*Xsig_aug(6, i);
          Xsig_pred(4, i) += delta_t*Xsig_aug(6, i);
      }

  }

  //write result
  *Xsig_out = Xsig_pred;

}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Z_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  weights(0) = lambda/(lambda+n_aug);
  for (int i=1; i<2*n_aug+1; i++) {
    weights(i) = 0.5/(n_aug+lambda);
  }

  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

   SigmaPointPrediction(&Xsig_pred);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Z_out = Zsig;
}

void UKF::UpdateStateRadar(VectorXd z, VectorXd* x_out, MatrixXd* P_out) {

  std::cout<<"Radar"<<std::endl;

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  weights(0) = lambda/(lambda+n_aug);
  for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
    weights(i) = 0.5/(n_aug+lambda);
  }

  //create example matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //create example vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //create example matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);

  PredictRadarMeasurement(&z_pred, &S, &Zsig);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  // initial state vector
  VectorXd x = x_mean_;

  // initial covariance matrix
  MatrixXd P = P_cov_;

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();

  while (x(3)> M_PI) x(3)-=2.*PI;
  while (x(3)<-M_PI) x(3)+=2.*PI;

  //write result
  *x_out = x;
  *P_out = P;

}

void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Z_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimensiont
  int n_z = 2;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  weights(0) = lambda/(lambda+n_aug);
  for (int i=1; i<2*n_aug+1; i++) {
    weights(i) = 0.5/(n_aug+lambda);
  }

  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
   MatrixXd P = MatrixXd(n_x, n_x);

   SigmaPointPrediction(&Xsig_pred);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Z_out = Zsig;
}

void UKF::UpdateStateLidar(VectorXd z, VectorXd* x_out, MatrixXd* P_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  weights(0) = lambda/(lambda+n_aug);
  for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
    weights(i) = 0.5/(n_aug+lambda);
  }

  //create example matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //create example vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //create example matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);

  PredictLidarMeasurement(&z_pred, &S, &Zsig);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  // initial state vector
  VectorXd x = x_mean_;

  // initial covariance matrix
  MatrixXd P = P_cov_;

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points


    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();

  while (x(3)> M_PI) x(3)-=2.*PI;
  while (x(3)<-M_PI) x(3)+=2.*PI;

  //write result
  *x_out = x;
  *P_out = P;

}
