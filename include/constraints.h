#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum ConstraintType{DISTANCE, COLLISION} ConstraintType;   //You can expand it for more constraints
typedef enum ConstraintEqualityType{EQUALITY, INEQUALITY} ConstraintEqualityType;

//there is such constraints per two variables that are equal. That is, for every attached vertex there are three such constraints for (x,y,z);
class Constraint{
public:
  
  int m1, m2;                     //Two participating meshes (can be the same)  - auxiliary data for users (constraint class shouldn't use that)
  int v1, v2;                     //Two vertices from the respective meshes - auxiliary data for users (constraint class shouldn't use that)
  double invMass1, invMass2;       //inverse masses of two bodies
  double refValue;                //Reference values to use in the constraint, when needed (like distance)
  bool isUpper;                   //in case this is an inequality constraints, whether it's an upper or a lower bound
  RowVector3d refVector;             //Reference vector when needed (like vector)
  double CRCoeff;                 //velocity bias
  ConstraintType constraintType;  //The type of the constraint, and will affect the value and the gradient. This SHOULD NOT change after initialization!
  ConstraintEqualityType constraintEqualityType;  //whether the constraint is an equality or an inequality
  bool isResolved; // Marks whether the constraint has been resolved in the current iteration
  bool isNotResolved;   // Marks whether the constraint is active (needs to be resolved)



  Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const bool _isUpper, const int& _m1, const int& _v1, const int& _m2, const int& _v2, const double& _invMass1, const double& _invMass2, const RowVector3d& _refVector, const double& _refValue, const double& _CRCoeff):constraintType(_constraintType), constraintEqualityType(_constraintEqualityType), isUpper(_isUpper), m1(_m1), v1(_v1), m2(_m2), v2(_v2), invMass1(_invMass1), invMass2(_invMass2),  refValue(_refValue), CRCoeff(_CRCoeff){
    refVector=_refVector;
  }
  
  ~Constraint(){}
  
  
  
  //computes the impulse needed for all particles to resolve the velocity constraint, and corrects the velocities accordingly.
  //The velocities are a vector (vCOM1, w1, vCOM2, w2) in both input and output.
  //returns true if constraint was already valid with "currVelocities", and false otherwise (false means there was a correction done)
  bool resolve_velocity_constraint(const MatrixXd& currCOMPositions, const MatrixXd& currVertexPositions, const MatrixXd& currCOMVelocities, const MatrixXd& currAngVelocities, const Matrix3d& invInertiaTensor1, const Matrix3d& invInertiaTensor2, MatrixXd& correctedCOMVelocities, MatrixXd& correctedAngVelocities, double tolerance){

    //get positions
    RowVector3d COM1 = currCOMPositions.row(0);
    RowVector3d COM2 = currCOMPositions.row(1);


    RowVector3d p1 = currVertexPositions.row(0);
    RowVector3d p2 = currVertexPositions.row(1);

    //position vectors relative to COM
    RowVector3d r1 = p1 - COM1;
    RowVector3d r2 = p2 - COM2;

    RowVector3d vCOM1 = currCOMVelocities.row(0);
    RowVector3d vCOM2 = currCOMVelocities.row(1);
    RowVector3d w1 = currAngVelocities.row(0);
    RowVector3d w2 = currAngVelocities.row(1);

    //get constraint normal 
    RowVector3d n = (p1 - p2).normalized();

    //get the angular components
    RowVector3d r1CrossN = r1.cross(n);
    RowVector3d r2CrossN = r2.cross(n);

    //construct the Jacobian matrix()
    VectorXd J(12);
    J << n.transpose(), r1CrossN.transpose(), -n.transpose(), -r2CrossN.transpose();
    //cout << "J matrix: \n" << J << endl;

    //construct Velocity vector
    VectorXd velocitiesVector(12);
    velocitiesVector << vCOM1.transpose(), w1.transpose(), vCOM2.transpose(), w2.transpose();

    //cout << "J size: " << J.rows() << " x " << J.cols() << endl;
   //cout << "inverseMass size: " << inverseMass.rows() << " x " << inverseMass.cols() << endl;

    //construct inverse mass matrix
    //check if both meshes are fixed (have zero inverse mass)
    if (invMass1 == 0.0 && invMass2 == 0.0) {
  
        return true;  //constraint is already valid
    }

    MatrixXd inverseMass(12, 12);
    inverseMass.setZero();

    //cout << "inverseMass matrix:\n" << inverseMass << endl;
    Matrix3d invMassMatrix1 = invMass1 * Matrix3d::Identity();
    Matrix3d invMassMatrix2 = invMass2 * Matrix3d::Identity();

    inverseMass.block<3, 3>(0, 0) = invMassMatrix1;  
    inverseMass.block<3, 3>(3, 3) = invInertiaTensor1;  
    inverseMass.block<3, 3>(6, 6) = invMassMatrix2;  
    inverseMass.block<3, 3>(9, 9) = invInertiaTensor2;  

    if (abs((J.transpose() * velocitiesVector)(0, 0)) < tolerance)
    {
      correctedCOMVelocities=currCOMVelocities;
      correctedAngVelocities=currAngVelocities;
      return true;
    }

    //compute the lagrange multiplyer
    MatrixXd temp = J.transpose() * inverseMass;
    MatrixXd denominatorM = temp * J;
    //MatrixXd denominatorM = J  * inverseMass * J.transpose();
    double denominator = denominatorM(0, 0);

    if (abs(denominator) < 1e-8) 
    {
      correctedCOMVelocities = currCOMVelocities;
      correctedAngVelocities = currAngVelocities;
      return true;
    }

    //apply the impulse correction to velocities
    double impulseMagnitude = -((1.0 + CRCoeff) * J.dot(velocitiesVector)) / denominator;
    // cout << " lambda " << impulseMagnitude << endl;
    // cout << " J: " << J << endl;
    // cout << "velocity change "  << inverseMass * J.transpose() * impulseMagnitude << endl;
    
    // velocities
    VectorXd temp2 = J.transpose() * impulseMagnitude;  //1 x 12 row vector multiplied by impulse
    VectorXd velocityCorrection = inverseMass * temp2;  
    //VectorXd velocityCorrection = inverseMass * J.transpose() * impulseMagnitude;
    velocitiesVector += velocityCorrection;

    correctedCOMVelocities = MatrixXd(2, 3);
    correctedAngVelocities = MatrixXd(2, 3);

    // Assign values directly without resizing
    correctedCOMVelocities.row(0) = velocitiesVector.segment(0, 3).transpose();
    correctedCOMVelocities.row(1) = velocitiesVector.segment(6, 3).transpose();
    correctedAngVelocities.row(0) = velocitiesVector.segment(3, 3).transpose();
    correctedAngVelocities.row(1) = velocitiesVector.segment(9, 3).transpose();

    return false;  //   
  }


  //projects the position unto the constraint
  //returns true if constraint was already good
  bool resolve_position_constraint(const MatrixXd& currCOMPositions, const MatrixXd& currConstPositions, MatrixXd& correctedCOMPositions, double tolerance){
    
    //current positions
    RowVector3d p1 = currConstPositions.row(0);
    RowVector3d p2 = currConstPositions.row(1);


    //compute current distance
    double currentD = (p2 - p1).norm();
    double CiP = currentD - refValue;

    if (constraintEqualityType == INEQUALITY) 
    {
      if (isUpper) {
        // distance shouldn;t be greater than refValue
        if (CiP <= tolerance)
        {
          correctedCOMPositions = currCOMPositions;
          return true;
        }
    } else {
        
        if (CiP >= -tolerance)
        {
          correctedCOMPositions = currCOMPositions;
          return true;
        }
      }
    }
    else if (constraintEqualityType == EQUALITY)
    {
      if (abs(CiP) <= tolerance)
      {
        correctedCOMPositions=currCOMPositions;
        return true;
      }
    }
    

    Vector3d n = (p2 - p1).normalized();           //direction vector
    Vector3d correction = (CiP / (invMass1 + invMass2)) * n;


    correctedCOMPositions = currCOMPositions;
    correctedCOMPositions.row(0) += invMass1 * correction;
    correctedCOMPositions.row(1) -= invMass2 * correction;


   // cout << "Correction applied: " << correction.transpose() << endl;

    return false;


   // correctedCOMPositions=currCOMPositions;
 
   // return true;
  
  }



};
#endif /* constraints_h */
