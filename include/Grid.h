#ifndef GRID_HEADER_FILE
#define GRID_HEADER_FILE

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <unordered_set>
#include <unordered_map>

using namespace std;
using namespace Eigen;


class Grid {
public:

    double gridSize;
    //min and max boundig boxes         
    Vector3d minBounds;      
    Vector3d maxBounds;      
    int cellX, cellY, cellZ;        // grid resolution

    unordered_map<int, vector<int>> cellObjects;        //hashmap cell index→list of objects
    unordered_set<int> fixedObjects;

    Grid(){}
    ~Grid(){}


    // //constructor
    // Grid(const vector<MatrixXd>& objectVertices, double scaleFactor = 2.0) {
    //     getBounds(objectVertices);
    //     setGridSize(scaleFactor);
    // }

    void initialize(const vector<MatrixXd>& objectVertices, double scaleFactor = 1.0){

        getBounds(objectVertices);
        setGridSize(scaleFactor);
    }

    // // Get the total number of cells in the grid
    //int getNumberOfCells() const {
    //     return numCellsX * numCellsY * numCellsZ;
    // }

    //compute min/max boundsfrom teh objects
    void getBounds(const vector<MatrixXd>& objectVertices) {


        minBounds = Vector3d::Constant(numeric_limits<double>::max());
        maxBounds = Vector3d::Constant(numeric_limits<double>::lowest());

        for (int i = 0; i < objectVertices.size(); ++i) 
        {
            const MatrixXd& vertices = objectVertices[i];

            Vector3d min = vertices.colwise().minCoeff();
            Vector3d max = vertices.colwise().maxCoeff();

            minBounds = minBounds.cwiseMin(min);
            maxBounds = maxBounds.cwiseMax(max);
        }
        
    }

    //set the grid resolution based on the bounding box and scale factor
    void setGridSize(double scaleFactor) {

        Vector3d boxSize = maxBounds - minBounds;


        gridSize = boxSize.maxCoeff() / scaleFactor;         

        cellX = ceil(boxSize.x() / gridSize);
        cellY = ceil(boxSize.y() / gridSize);
        cellZ = ceil(boxSize.z() / gridSize);

    }


    // changing general position to grid pos
    int getCellIndex(const Vector3d& position) const {

        // Calculate grid coordinates
        int x = static_cast<int>((position.x() - minBounds.x()) / gridSize);
        int y = static_cast<int>((position.y() - minBounds.y()) / gridSize);
        int z = static_cast<int>((position.z() - minBounds.z()) / gridSize);

        
        int cellIndex = x + cellX * (y + cellY * z);

        return cellIndex;

    }


    //insert object into grid
    void insert(int object, const MatrixXd& vertices, bool isFixed = false) {

        unordered_set<int> visitedCells; 

        if (isFixed)
        {
            fixedObjects.insert(object);
            
        }

        for (int i = 0; i < vertices.rows(); i++) 
        {
            int cellIndex = getCellIndex(vertices.row(i));

            if (visitedCells.find(cellIndex) == visitedCells.end()) 
            {   
                cellObjects[cellIndex].push_back(object);
                visitedCells.insert(cellIndex);
                //cout << "Object " << objectID << " inserted into cell " << cellIndex << endl;
            }


        }

    }

    

        // void insertObject(int objectID, const MatrixXd& vertices) {
        //     for (int i = 0; i < vertices.rows(); i++) {
        //         int cellIndex = getCellIndex(vertices.row(i));
        //         cellObjects[cellIndex].push_back(objectID);
        //     }
        // }



    };



#endif


