#ifndef CREATE_ANALYTIC_H
#define CREATE_ANALYTIC_H

#include <PCU.h>
#include <lionPrint.h>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include <gmi.h>
#include <gmi_analytic.h>
#include <gmi_mesh.h>
#include <gmi_null.h>

#include <ma.h>

#include <cassert>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <math.h>
#include <pcu_util.h>


// model inspection tools
char const* const typeName[4] =
{"vertex",
 "  edge",
 "  face",
 "region"};

void printInfo(gmi_model* model, int dim);
void visualizeFace(gmi_model* model, gmi_ent* entity, int n, int m, const char* fileName);

agm_bdry add_bdry(gmi_model* m, gmi_ent* e);
agm_use add_adj(gmi_model* m, agm_bdry b, int tag);
agm_use add_adj(gmi_model* m, agm_bdry b, int dim, int tag);

extern double boxLength;
extern double boxWidth;
extern double boxHeight;
extern double sphereRadius;
extern double xyz_offset[3];
extern apf::Mesh2* m;

struct Enclosure{
    
    std::vector<int> vertexMap;
    double vertRan[1][2]={{0.0,0.0}};
    int vertPer=0;

    std::vector<int> edgeMap;
    int edgePer = 0;
    double edgeRan[1][2] = {{0.0,1.0}};

    std::vector<int> faceMap;
    int faPer[2] = {0, 0};
    double faRan[2][2] = {{0,1},{0,1}};

    int regionID = 92; // fixed ID
    void makeBox2D(gmi_model* model); 
    void makeBox3D(gmi_model* model); 
};

class Sphere{
    public:
    int faceID; //values are dictated by spatial tools
    double radius;
    double offset[3];
    int dim; 
    int faPer[2] = {1, 0};
    double faRan[2][2] = {{0,6.28318530718},{0.0,apf::pi}};

    Sphere(int x){
        dim = x;
        if(dim==2){
            faRan[1][1] = 0.0; 
            faceID=6;
        }
        if(dim==3)
            faceID=9;
    }
    Sphere(){}
    void makeSphere(gmi_model* model);
};

class PiercingCylinder{
    public:
    int faceID; //values are dictated by spatial tools
    double radius;
    double offset[3];
    int dim; 
    int faPer[2] = {1, 0};
    double faRan[2][2] = {{0,6.28318530718},{0.0,apf::pi}};

    PiercingCylinder(){
        dim = 3;
        faceID=11;
    }
    void makePiercingCylinder(gmi_model* model);
};


namespace Reparam{
    void reparameterizeEntities(gmi_model*model,apf::Mesh2*m,Enclosure box, Sphere sphere);
    void reparameterizeEntities(gmi_model*model,apf::Mesh2*m,Enclosure box, PiercingCylinder cylinder, Sphere circle1, Sphere circle2);
    void reparameterize2D(gmi_model*model,apf::Mesh2*m,Enclosure box, Sphere sphere);
    void reparameterize3D(gmi_model*model,apf::Mesh2*m,Enclosure box, Sphere sphere);
}
void freeField(apf::Field*& );
#endif
