#ifndef CLASS_HEADER_FILE
#define CLASS_HEADER_FILE

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


class Grid {
public:

    double gridSize;         // grid cell size
    Vector3d minBounds;      // minimum bounding box
    Vector3d maxBounds;      // maximum bounding box
    int numCellsX, numCellsY, numcellsZ;        // grid resolution

    unordered_map<int, vector<int>> cellObjects;        // Hashmap: cell index → list of objects

    // Constructor: Computes grid from objects
    Grid(const vector<MatrixXd>& objectVertices, double scaleFactor = 2.0) {
        computeBounds(objectVertices);
        setGridSize(scaleFactor);
    }

    // Compute global min/max bounds from all objects
    void computeBounds(const vector<MatrixXd>& objectVertices) {
        minBounds = Vector3d::Constant(numeric_limits<double>::max());
        maxBounds = Vector3d::Constant(numeric_limits<double>::lowest());

        for (int i = 0; i < objectVertices.size(); ++i) {
            const MatrixXd& vertices = objectVertices[i];
            Vector3d objMin = vertices.colwise().minCoeff();
            Vector3d objMax = vertices.colwise().maxCoeff();

            minBounds = minBounds.cwiseMin(objMin);
            maxBounds = maxBounds.cwiseMax(objMax);
        }
    }

    // Set grid resolution based on bounding box
    void setGridSize(double scaleFactor) {
        Vector3d boxSize = maxBounds - minBounds;
        h = boxSize.maxCoeff() / scaleFactor; // Set cell size
        numCellsX = ceil(boxSize.x() / h);
        numCellsY = ceil(boxSize.y() / h);
        numCellsZ = ceil(boxSize.z() / h);
    }

    // Convert world position to grid index
    int getCellIndex(const Vector3d& position) const {
        int x = floor((position.x() - minBounds.x()) / h);
        int y = floor((position.y() - minBounds.y()) / h);
        int z = floor((position.z() - minBounds.z()) / h);
        return x + numCellsX * (y + numCellsY * z);
    }

    // Insert an object into the grid
    void insertObject(int objectID, const MatrixXd& vertices) {
        for (int i = 0; i < vertices.rows(); i++) {
            int cellIndex = getCellIndex(vertices.row(i));
            cellObjects[cellIndex].push_back(objectID);
        }
    }

};



#endif


