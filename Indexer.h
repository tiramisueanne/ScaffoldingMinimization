#pragma once
#include <map>
#include <vector>
using namespace std;
// A class to get the index into a matrix of only internal nodes
class Indexer {
    int numVert;
    vector<int> indexes;
    map<pair<int, int>, int> edgeIndex;

   public:
    Indexer() : numVert(-1) {};
    Indexer(int numVertices, set<int> internalNodes,
            set<pair<int, int>> allEdges);

    inline int indexVert(int vertex) { return indexes[vertex]; }

    // If the edge is v1 - v2
    inline int indexEdge(int v1, int v2) {
        return edgeIndex[pair<int, int>(v1, v2)];
    }

    inline const map<pair<int, int>, int> &edgeMap() { return edgeIndex; }
};

inline Indexer::Indexer(int numVertices, set<int> internalNodes,
                 set<pair<int, int>> allEdges)
    : numVert(numVertices), indexes(numVert) {
    int matrixIndex = 0;
    for (int internal : internalNodes) {
        indexes[internal] = matrixIndex;
        matrixIndex++;
    }
    int matIndexEdge = 0;
    for (const auto edge : allEdges) {
        if (edgeIndex.find(edge) == edgeIndex.end()) {
            edgeIndex[edge] = matIndexEdge;
            edgeIndex[pair<int, int>(edge.second, edge.first)] = matIndexEdge;
            matIndexEdge++;
        }
    }
}
