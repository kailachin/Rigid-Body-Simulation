#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <fstream>
#include "ccd.h"
#include "volInt.h"
#include "auxfunctions.h"
#include "readMESH.h"
#include "mesh.h"
#include "constraints.h"
#include "Grid.h"
#include <unordered_map>


using namespace Eigen;
using namespace std;





//This class contains the entire scene operations, and the engine time loop.
class Scene{
public:
  double currTime;
  vector<Mesh> meshes;
  vector<Constraint> constraints;
  Mesh groundMesh;
  
  //Mostly for visualization
  MatrixXi allF, constEdges;
  MatrixXd currV, currConstVertices;

  //for the grid (broad-phase space subdivision)
  Grid spatialGrid;

  
  //BVHNode* rootBVH = nullptr;

  //for constraint scheduling
  unordered_map<int, vector<int>> mConstraints;
  
  //adding an objects. You do not need to update this generally
  void add_mesh(const MatrixXd& V, const MatrixXi& F, const MatrixXi& T, const double density, const bool isFixed, const RowVector3d& COM, const RowVector4d& orientation){
    
    Mesh m(V,F, T, density, isFixed, COM, orientation);
    meshes.push_back(m);
    //cout<<"m.origV.row(0): "<<m.origV.row(0)<<endl;
    //cout<<"m.currV.row(0): "<<m.currV.row(0)<<endl;
    
    MatrixXi newAllF(allF.rows()+F.rows(),3);
    newAllF<<allF, (F.array()+currV.rows()).matrix();
    allF = newAllF;
    MatrixXd newCurrV(currV.rows()+V.rows(),3);
    newCurrV<<currV, m.currV;
    currV = newCurrV;
  }
  
  /*********************************************************************
   This function handles a collision between objects ro1 and ro2 when found, by assigning impulses to both objects.
   Input: RigidObjects m1, m2
   depth: the depth of penetration
   contactNormal: the normal of the conact measured m1->m2
   penPosition: a point on m2 such that if m2 <= m2 + depth*contactNormal, then penPosition+depth*contactNormal is the common contact point
   CRCoeff: the coefficient of restitution
   *********************************************************************/
  void handle_collision(Mesh& m1, Mesh& m2,const double& depth, const RowVector3d& contactNormal,const RowVector3d& penPosition, const double CRCoeff){

    RowVector3d contactPosition;

    
    if (m1.isFixed)
    {
      m2.currV.rowwise() += depth * contactNormal;
      m2.COM += depth * contactNormal;
      contactPosition = penPosition + depth * contactNormal;
    } 
    else if (m2.isFixed)
    {
      m1.currV.rowwise() -= depth * contactNormal;
      m1.COM -= depth * contactNormal;
      contactPosition = penPosition - depth * contactNormal;
    } 
    else 
    { 
      double totalWeight = m1.totalInvMass + m2.totalInvMass;
      double weight2 = m2.totalInvMass / totalWeight;
      double weight1 = m1.totalInvMass / totalWeight;
      m1.currV.rowwise() -= weight1 * depth * contactNormal;
      m1.COM -= weight1 * depth * contactNormal;
      m2.currV.rowwise() += weight2 * depth * contactNormal;
      m2.COM += weight2 * depth * contactNormal;
      
      contactPosition = penPosition + weight2 * depth * contactNormal;
    }
    
    //calculate r1 and r2 (distance from COM of each mesh to p12)
    RowVector3d r1 = contactPosition - m1.COM;
    RowVector3d r2 = contactPosition - m2.COM;
    
    //compute velocities at p12
    RowVector3d v1 = m1.comVelocity + m1.angVelocity.cross(r1);
    RowVector3d v2 = m2.comVelocity + m2.angVelocity.cross(r2);

    RowVector3d v12 = v1 - v2;                  //v1 - v2
    double vn = contactNormal.dot(v12);
    
    //get inverse inertia tensors
    Matrix3d invIT1 = m1.get_curr_inv_IT();
    Matrix3d invIT2 = m2.get_curr_inv_IT();


    Vector3d r1_cross_n = r1.cross(contactNormal);  // cross product, result is a Vector3d
    double r1Inverse = r1_cross_n.dot(invIT1 * r1_cross_n);  
   // Vector3d term1_cross_r1 = term1.cross(r1);  

    Vector3d r2_cross_n = r2.cross(contactNormal);  
    double r2Inverse = r2_cross_n.dot(invIT2 * r2_cross_n);  
    //Vector3d term2_cross_r2 = term2.cross(r2); 

    
    //double denominator = totalInvMass + contactNormal.dot(term1_cross_r1 + term2_cross_r2);
    double denominator = m1.totalInvMass + m2.totalInvMass + r1Inverse + r2Inverse;

    double numerator = -(1.00 + CRCoeff) * vn;
    
    double j = numerator / denominator;

    RowVector3d impulse = j * contactNormal;
    
    
    if (impulse.norm() >10e-6)
    {

    // if (!m1.isFixed) {
    //   //m1.comVelocity += invMass1 * j * contactNormal;
    m1.currImpulses.push_back(Impulse(contactPosition, impulse));
       
      
      // m2.currImpulses.push_baxk(impulse);
    // }
    // if (!m2.isFixed) {
    //   //m2.comVelocity -= invMass2 * j * contactNormal;

    m2.currImpulses.push_back(Impulse(contactPosition, -impulse));

      
     // m2.currImpulses.push_back(impulse);
    m1.update_impulse_velocities();
    m2.update_impulse_velocities();
    }
}


  
  /*********************************************************************
   This function handles a single time step by:
   1. Integrating velocities, positions, and orientations by the timeStep
   2. detecting and handling collisions with the coefficient of restitutation CRCoeff
   3. updating the visual scene in fullV and fullT
   *********************************************************************/
  void update_scene(double timeStep, double CRCoeff, int maxIterations, double tolerance){

    //integrating velocity, position and orientation from forces and previous states
    for (int i=0;i<meshes.size();i++)
      meshes[i].integrate(timeStep);
    
  //  spatialGrid.cellObjects.clear();

  //   for (int i = 0; i < meshes.size(); i++)
  //   {
      
  //     spatialGrid.insert(i, meshes[i].currV, meshes[i].isFixed);     
  //   }

  //   // for (int i = 0; i < meshes.size(); i++) 
  //   // {
  //   //   if (!meshes[i].isFixed) {
  //   //     spatialGrid.insert(i, meshes[i].currV);
  //   // }
  //   // }


  //   double depth;
  //   RowVector3d contactNormal, penPosition;

  // unordered_map<int, vector<int>>::iterator i;
  //   for (i = spatialGrid.cellObjects.begin(); i != spatialGrid.cellObjects.end(); i++)
  //   {
  //     const vector<int>& objectsInCell = i->second;

  //     for (int i = 0; i < objectsInCell.size(); ++i) {
  //       for (int j = i + 1; j < objectsInCell.size(); ++j) {
  //           int object1 = objectsInCell[i];
  //           int object2 = objectsInCell[j];
        
  //           if (meshes[object1].is_collide(meshes[object2], depth, contactNormal, penPosition)) {
  //               if (meshes[object1].isFixed || meshes[object2].isFixed) {
  //                   //if there is a collision with fixed
  //                   if (meshes[object1].isFixed) {
                      
  //                       handle_collision(meshes[object2], meshes[object1], depth, contactNormal, penPosition, CRCoeff);
  //                   } else if (meshes[object2].isFixed) {
                        
  //                       handle_collision(meshes[object1], meshes[object2], depth, contactNormal, penPosition, CRCoeff);
  //                   }
  //               } else {
  //                   //not fixed
  //                   handle_collision(meshes[object1], meshes[object2], depth, contactNormal, penPosition, CRCoeff);
  //               }
  //           }
  //           // // Check if the two objects are colliding
  //           // if (meshes[object1].is_collide(meshes[object2], depth, contactNormal, penPosition)) {
  //           //     handle_collision(meshes[object1], meshes[object2], depth, contactNormal, penPosition, CRCoeff);
  //           // }
  //       }
  //     }
  //   }


 

    
    //detecting and handling collisions when found
    //This is done exhaustively: checking every two objects in the scene.
    double depth;
    RowVector3d contactNormal, penPosition;
    for (int i=0;i<meshes.size();i++)
      for (int j=i+1;j<meshes.size();j++)
        if (meshes[i].is_collide(meshes[j],depth, contactNormal, penPosition))
          handle_collision(meshes[i], meshes[j],depth, contactNormal, penPosition,CRCoeff);    



    //colliding with the pseudo-mesh of the ground
    for (int i=0;i<meshes.size();i++){
      int minyIndex;
      double minY = meshes[i].currV.col(1).minCoeff(&minyIndex);
      //linear resolution
      if (minY<=0.0)
        handle_collision(meshes[i], groundMesh, minY, {0.0,1.0,0.0},meshes[i].currV.row(minyIndex),CRCoeff);
    }
    
    //Resolving constraints
    //unordered_set<int> affectedMeshes;

    int currIteration=0;
    int zeroStreak=0;  //how many consecutive constraints are already below tolerance without any change; the algorithm stops if all are.
    int currConstIndex=0;

    //unordered_set<int> flaggedConstraints;

    while ((zeroStreak<constraints.size())&&(currIteration*constraints.size()<maxIterations)){
      

      Constraint currConstraint=constraints[currConstIndex];
      
      RowVector3d origConstPos1=meshes[currConstraint.m1].origV.row(currConstraint.v1);
      RowVector3d origConstPos2=meshes[currConstraint.m2].origV.row(currConstraint.v2);
      
      RowVector3d currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation)+meshes[currConstraint.m1].COM;
      RowVector3d currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation)+meshes[currConstraint.m2].COM;
      
      MatrixXd currCOMPositions(2,3); currCOMPositions<<meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
      MatrixXd currConstPositions(2,3); currConstPositions<<currConstPos1, currConstPos2;
      
      MatrixXd correctedCOMPositions;
      
      bool positionWasValid=currConstraint.resolve_position_constraint(currCOMPositions, currConstPositions,correctedCOMPositions, tolerance);
      
      if (positionWasValid){
        zeroStreak++;
      }else{
        //only update the COM and angular velocity, don't both updating all currV because it might change again during this loop!
        zeroStreak=0;
        
        meshes[currConstraint.m1].COM = correctedCOMPositions.row(0);
        meshes[currConstraint.m2].COM = correctedCOMPositions.row(1);
        
        //resolving velocity
        currConstPos1 = QRot(origConstPos1, meshes[currConstraint.m1].orientation)+meshes[currConstraint.m1].COM;
        currConstPos2 = QRot(origConstPos2, meshes[currConstraint.m2].orientation)+meshes[currConstraint.m2].COM;
        //cout<<"(currConstPos1-currConstPos2).norm(): "<<(currConstPos1-currConstPos2).norm()<<endl;
        //cout<<"(meshes[currConstraint.m1].currV.row(currConstraint.v1)-meshes[currConstraint.m2].currV.row(currConstraint.v2)).norm(): "<<(meshes[currConstraint.m1].currV.row(currConstraint.v1)-meshes[currConstraint.m2].currV.row(currConstraint.v2)).norm()<<endl;
        currCOMPositions<<meshes[currConstraint.m1].COM, meshes[currConstraint.m2].COM;
        currConstPositions<<currConstPos1, currConstPos2;
        MatrixXd currCOMVelocities(2,3); currCOMVelocities<<meshes[currConstraint.m1].comVelocity, meshes[currConstraint.m2].comVelocity;
        MatrixXd currAngVelocities(2,3); currAngVelocities<<meshes[currConstraint.m1].angVelocity, meshes[currConstraint.m2].angVelocity;
        
        Matrix3d invInertiaTensor1=meshes[currConstraint.m1].get_curr_inv_IT();
        Matrix3d invInertiaTensor2=meshes[currConstraint.m2].get_curr_inv_IT();
        MatrixXd correctedCOMVelocities, correctedAngVelocities, correctedCOMPositions;
        
        bool velocityWasValid=currConstraint.resolve_velocity_constraint(currCOMPositions, currConstPositions, currCOMVelocities, currAngVelocities, invInertiaTensor1, invInertiaTensor2, correctedCOMVelocities,correctedAngVelocities, tolerance);
        
        if (!velocityWasValid){
          meshes[currConstraint.m1].comVelocity =correctedCOMVelocities.row(0);
          meshes[currConstraint.m2].comVelocity =correctedCOMVelocities.row(1);
          
          meshes[currConstraint.m1].angVelocity =correctedAngVelocities.row(0);
          meshes[currConstraint.m2].angVelocity =correctedAngVelocities.row(1);
        }
      }
      
      currIteration++;
      currConstIndex=(currConstIndex+1)%(constraints.size());
    }
    
    if (currIteration*constraints.size()>=maxIterations)
      cout<<"Constraint resolution reached maxIterations without resolving!"<<endl;
    
    currTime+=timeStep;
    
    //updating meshes and visualization
    for (int i=0;i<meshes.size();i++)
      for (int j=0;j<meshes[i].currV.rows();j++)
        meshes[i].currV.row(j)<<QRot(meshes[i].origV.row(j), meshes[i].orientation)+meshes[i].COM;
    
    int currVOffset=0;
    for (int i=0;i<meshes.size();i++){
      currV.block(currVOffset, 0, meshes[i].currV.rows(), 3) = meshes[i].currV;
      currVOffset+=meshes[i].currV.rows();
    }
    for (int i=0;i<constraints.size();i+=2){   //jumping bc we have constraint pairs
      currConstVertices.row(i) = meshes[constraints[i].m1].currV.row(constraints[i].v1);
      currConstVertices.row(i+1) = meshes[constraints[i].m2].currV.row(constraints[i].v2);
    }
  }
  
  //loading a scene from the scene .txt files
  //you do not need to update this function
  bool load_scene(const std::string sceneFileName, const std::string constraintFileName){
    
    ifstream sceneFileHandle, constraintFileHandle;
    sceneFileHandle.open(DATA_PATH "/" + sceneFileName);
    if (!sceneFileHandle.is_open())
      return false;
    int numofObjects;
    vector<MatrixXd> objectVertices;
    vector<int> objectIDs;

    currTime=0;
    sceneFileHandle>>numofObjects;
    for (int i=0;i<numofObjects;i++){
      MatrixXi objT, objF;
      MatrixXd objV;
      std::string MESHFileName;
      bool isFixed; 
      double density;
      RowVector3d userCOM;
      RowVector4d userOrientation;
      sceneFileHandle>>MESHFileName>>density>>isFixed>>userCOM(0)>>userCOM(1)>>userCOM(2)>>userOrientation(0)>>userOrientation(1)>>userOrientation(2)>>userOrientation(3);
      userOrientation.normalize();
      readMESH(DATA_PATH "/" + MESHFileName,objV,objF, objT);
      objectVertices.push_back(objV);
      objectIDs.push_back(i);

      //fixing weird orientation problem
      MatrixXi tempF(objF.rows(),3);
      tempF<<objF.col(2), objF.col(1), objF.col(0);
      objF=tempF;
      
      add_mesh(objV,objF, objT,density, isFixed, userCOM, userOrientation);

    }

    //build_bvh();
    spatialGrid.initialize(objectVertices);

    for (int i = 0; i < objectVertices.size(); i++)
    {
      spatialGrid.insert(objectIDs[i], objectVertices[i]);    //populate the grid initially
    }
    
    //adding ground mesh artifically
    groundMesh = Mesh(MatrixXd(0,3), MatrixXi(0,3), MatrixXi(0,4), 0.0, true, RowVector3d::Zero(), RowVector4d::Zero());
    
    //Loading constraints
    int numofConstraints;
    constraintFileHandle.open(DATA_PATH "/" + constraintFileName);
    if (!constraintFileHandle.is_open())
      return false;
    constraintFileHandle>>numofConstraints;
    currConstVertices.resize(numofConstraints*2,3);
    constEdges.resize(numofConstraints,2);

    mConstraints.clear();

    for (int i=0;i<numofConstraints;i++){
      int attachM1, attachM2, attachV1, attachV2;
      double lowerBound, upperBound;
      constraintFileHandle>>attachM1>>attachV1>>attachM2>>attachV2>>lowerBound>>upperBound;
      //cout<<"Constraints: "<<attachM1<<","<<attachV1<<","<<attachM2<<","<<attachV2<<","<<lowerBound<<","<<upperBound<<endl;
      
      double initDist=(meshes[attachM1].currV.row(attachV1)-meshes[attachM2].currV.row(attachV2)).norm();
      //cout<<"initDist: "<<initDist<<endl;
      double invMass1 = (meshes[attachM1].isFixed ? 0.0 : meshes[attachM1].totalInvMass);  //fixed meshes have infinite mass
      double invMass2 = (meshes[attachM2].isFixed ? 0.0 : meshes[attachM2].totalInvMass);
      constraints.push_back(Constraint(DISTANCE, INEQUALITY,false, attachM1, attachV1, attachM2, attachV2, invMass1,invMass2,RowVector3d::Zero(), lowerBound*initDist, 0.0));
      constraints.push_back(Constraint(DISTANCE, INEQUALITY,true, attachM1, attachV1, attachM2, attachV2, invMass1,invMass2,RowVector3d::Zero(), upperBound*initDist, 0.0));
      currConstVertices.row(2*i) = meshes[attachM1].currV.row(attachV1);
      currConstVertices.row(2*i+1) = meshes[attachM2].currV.row(attachV2);
      constEdges.row(i)<<2*i, 2*i+1;
    }
    
    for (int i = 0; i < constraints.size(); i++)
    {
      const Constraint& constraint = constraints[i];
      mConstraints[constraint.m1].push_back(i);
      mConstraints[constraint.m2].push_back(i);
    }

    return true;
  }
  
  
  Scene(){allF.resize(0,3); currV.resize(0,3);}
  ~Scene(){}






};



#endif
