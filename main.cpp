#include "main.h"
#include <ma.h>
#include <lionPrint.h>
#include <PCU.h>
#include <gmi.h>
#include <gmi_mesh.h>
#include <gmi_null.h>

agm_bdry add_bdry(gmi_model* m, gmi_ent* e)
{
  return agm_add_bdry(gmi_analytic_topo(m), agm_from_gmi(e));
}

agm_use add_adj(gmi_model* m, agm_bdry b, int tag)
{
  agm* topo = gmi_analytic_topo(m);
  int dim = agm_dim_from_type(agm_bounds(topo, b).type);
  gmi_ent* de = gmi_find(m, dim - 1, tag);
  return agm_add_use(topo, b, agm_from_gmi(de));
}

//anonymous namespace to make all functions local to this file scope

double xyz_offset[3];
double sphereRadius;
double boxLength;
double boxWidth;
double boxHeight;
int geomDim;
Enclosure modelBox;
Sphere modelSphere;
Sphere modelCircle1;
Sphere modelCircle2;
PiercingCylinder modelPiercingCylinder;
apf::Mesh2* m;

namespace 
{
    void vert0(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = 0.0;

    }
    void vert1(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = boxLength;
      x[1] = 0.0;
      x[2] = 0.0;
    }

    void vert2(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = boxLength;
      x[1] = boxWidth;
      x[2] = 0.0;
    }

    void vert3(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = 0.0;
      x[1] = boxWidth;
      x[2] = 0.0;
    }

    void vert4(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = boxHeight;

    }
    void vert5(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = boxLength;
      x[1] = 0.0;
      x[2] = boxHeight;
    }

    void vert6(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = boxLength;
      x[1] = boxWidth;
      x[2] = boxHeight;
    }

    void vert7(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = 0.0;
      x[1] = boxWidth;
      x[2] = boxHeight;
    }

    void edge0(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = 0.0;
      x[2] = 0.0;
    }

    void edge1(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = boxWidth*p[0];
      x[2] = 0.0;
    }

    void edge2(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth;
      x[2] = 0.0;
    }

    void edge3(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = boxWidth*p[0];
      x[2] = 0.0;
    }

    void edge4(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = 0.0;
      x[2] = boxHeight;
    }

    void edge5(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = boxWidth*p[0];
      x[2] = boxHeight;
    }

    void edge6(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth;
      x[2] = boxHeight;
    }

    void edge7(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = boxWidth*p[0];
      x[2] = boxHeight;
    }

    void edge8(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = boxHeight*p[0];
    }

    void edge9(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = 0.0;
      x[2] = boxHeight*p[0];
    }

    void edge10(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = boxWidth;
      x[2] = boxHeight*p[0];
    }

    void edge11(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = boxWidth;
      x[2] = boxHeight*p[0];
    }

    void face0(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = 0.0;
      x[2] = boxHeight*p[1];
    }

    void face1(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = boxWidth*p[0];
      x[2] = boxHeight*p[1];
    }

    void face2(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth;
      x[2] = boxHeight*p[1];
    }

    void face3(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = boxWidth*p[0];
      x[2] = boxHeight*p[1];
    }
    void face4(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth*p[1];
      x[2] = 0.0;
    }
    void face5(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth*p[1];
      x[2] = boxHeight;
    }

    void reparamVert_zero(double const from[2], double to[2], void*)
    {
      (void)from;
      to[0] = 0;
      to[1] = 0;
    }
    void reparamVert_one(double const from[2], double to[2], void*)
    {
      (void)from;
      to[0] = 1;
      to[1] = 0;
    }

    void reparamREdge_0(double const from[2], double to[2], void*)
    {
      to[0] = from[0];
      to[1] = 0.0;
    }

    void reparamREdge_1(double const from[2], double to[2], void*)
    {
      to[0] = 0.0;
      to[1] = from[0];
    }
    void reparamREdge_2(double const from[2], double to[2], void*)
    {
      to[0] = from[0];
      to[1] = 1.0;
    }

    void reparamREdge_3(double const from[2], double to[2], void*)
    {
      to[0] = 1.0; 
      to[1] = from[0];
    }

    void regionFunction(double const p[2], double x[3], void*)
    {
      (void)p;
      (void)x;
    }

    //from circle to box
/*
    void reparam_Circle(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        double x = sphereRadius*cos(from[0])+xyz_offset[0];
        double y = sphereRadius*sin(from[0])+xyz_offset[1];
        to[0] = x/boxLength;
        to[1] = y/boxWidth;
    }
*/
    void reparam_Circle0(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        double x = sphereRadius*cos(from[0])+xyz_offset[0];
        double y = sphereRadius*sin(from[0])+xyz_offset[1];
        to[0] = x/boxLength;
        to[1] = y/boxWidth;
        to[2] = 0.0;
    }
    void reparam_Circle1(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        double x = sphereRadius*cos(from[0])+xyz_offset[0];
        double y = sphereRadius*sin(from[0])+xyz_offset[1];
        to[0] = x/boxLength;
        to[1] = y/boxWidth;
        to[2] = boxHeight;
    }

    void reparam_CircleCylinder0(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        //double x = sphereRadius*cos(from[0])+xyz_offset[0];
        //double y = sphereRadius*sin(from[0])+xyz_offset[1];
        //to[0] = x/boxLength;
        //to[1] = y/boxWidth;
        to[0] = from[0];
        to[1] = 0.0;//from[1];
        //to[2] = 0;
    }
    void reparam_CircleCylinder1(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        //double x = sphereRadius*cos(from[0])+xyz_offset[0];
        //double y = sphereRadius*sin(from[0])+xyz_offset[1];
        //to[0] = x/boxLength;
        //to[1] = y/boxWidth;
        to[0] = from[0];
        to[1] = 1.0;//from[1];
        //to[2] = 1.0;//boxHeight;
    }


//need to set these functions separately from a member because the function signatures are otherwise modified.
//Likewise the internals need to use static variables to still connect them with other objects

    void sphereFace(double const p[2], double x[3], void*)
    {
      if(geomDim == 2){
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]);
        x[2] = 0.0;
      }
      else if(geomDim == 3){
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]) * sin(p[1]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]) * sin(p[1]);
        x[2] = xyz_offset[2]+sphereRadius*cos(p[1]);
      }
    }
    void circleFace0(double const p[2], double x[3], void*)
    {
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]);
        x[2] = 0.0;
    }

    void circleFace1(double const p[2], double x[3], void*)
    {
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]);
        x[2] = boxHeight;
    }

    void cylinderFace(double const p[2], double x[3], void*)
    {
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]);
        //x[1] = xyz_offset[1]+sphereRadius*sin(p[0]) * sin(p[1]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]);
        x[2] = boxHeight*p[1]; //need to double check this
    }


}

typedef void (*EntityMapArray) (double const p[2], double x[3], void*);
typedef void (*ParametricFunctionArray) (double const from[2], double to[2], void*);

void Enclosure::makeBox3D(gmi_model* model)
{
  //making a box

  vertexMap = {58,56,54,60,5,10,15,2};
  gmi_ent* g_vert[vertexMap.size()];
  EntityMapArray vertexPoints[] = {
        vert0,
        vert1,
        vert2,
        vert3,
        vert4,
        vert5,
        vert6,
        vert7
    };

  for(auto i=0; i<vertexMap.size();i++)
    g_vert[i] = gmi_add_analytic(model, 0, vertexMap[i], vertexPoints[i], &vertPer, vertRan, 0); 

  edgeMap = {50,48,46,52,11,16,20,6,73,72,71,74};
  gmi_ent* g_edge[edgeMap.size()];
  EntityMapArray edgeEntities[] = {
        edge0,
        edge1,
        edge2,
        edge3,
        edge4,
        edge5,
        edge6,
        edge7,
        edge8,
        edge9,
        edge10,
        edge11
  };
  for(int i=0;i<edgeMap.size();i++)
    g_edge[i] = gmi_add_analytic(model, 1, edgeMap[i], edgeEntities[i], &edgePer, edgeRan, 0);

  //reparameterize vertices on edges
  agm_bdry b;
  std::vector<std::pair<int,int>> indexPairs = {{0,1}, {1,2}, {2,3}, {3,0},
            {4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}  };


  //std::vector<std::pair<int,int>> indexPairs = {{0,1}, {1,2}, {3,2}, {0,3},
  //          {4,5},{5,6},{7,6},{4,7},{0,4},{1,5},{2,6},{3,7}  };

/*
  for(int i=0;i<edgeMap.size();i++){

    b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[i]));
    agm_use edgeUse0 = add_adj(model, b, vertexMap[indexPairs[i].first]);
    agm_use edgeUse0_1 = add_adj(model,b,vertexMap[indexPairs[i].second]);
    gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
    gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);
  }
*/

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[0]));
  agm_use edgeUse0 = add_adj(model, b, vertexMap[0]);
  agm_use edgeUse0_1 = add_adj(model,b,vertexMap[1]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[1]));
  edgeUse0 = add_adj(model, b, vertexMap[1]);
  edgeUse0_1 = add_adj(model,b,vertexMap[2]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[2]));
  edgeUse0 = add_adj(model, b, vertexMap[2]);
  edgeUse0_1 = add_adj(model,b,vertexMap[3]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_one, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_zero, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[3]));
  edgeUse0 = add_adj(model, b, vertexMap[3]);
  edgeUse0_1 = add_adj(model,b,vertexMap[0]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_one, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_zero, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[4]));
  edgeUse0 = add_adj(model, b, vertexMap[4]);
  edgeUse0_1 = add_adj(model,b,vertexMap[5]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[5]));
  edgeUse0 = add_adj(model, b, vertexMap[5]);
  edgeUse0_1 = add_adj(model,b,vertexMap[6]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[6]));
  edgeUse0 = add_adj(model, b, vertexMap[6]);
  edgeUse0_1 = add_adj(model,b,vertexMap[7]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_one, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_zero, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[7]));
  edgeUse0 = add_adj(model, b, vertexMap[7]);
  edgeUse0_1 = add_adj(model,b,vertexMap[4]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_one, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_zero, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[8]));
  edgeUse0 = add_adj(model, b, vertexMap[0]);
  edgeUse0_1 = add_adj(model,b,vertexMap[4]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[9]));
  edgeUse0 = add_adj(model, b, vertexMap[1]);
  edgeUse0_1 = add_adj(model,b,vertexMap[5]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[10]));
  edgeUse0 = add_adj(model, b, vertexMap[2]);
  edgeUse0_1 = add_adj(model,b,vertexMap[6]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[11]));
  edgeUse0 = add_adj(model, b, vertexMap[3]);
  edgeUse0_1 = add_adj(model,b,vertexMap[7]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);


  //reparameterize vertices on edges
  //make faces

  EntityMapArray faceEntities[] = {
        face0,
        face1,
        face2,
        face3,
        face4,
        face5,
  };

  //faceMap = {80,78,76,82,42,24};
  faceMap = {5,4,3,2,1,6};
  gmi_ent* g_face[faceMap.size()];
  for(int i=0;i<edgeMap.size();i++){
    g_face[i] = gmi_add_analytic(model, 2, faceMap[i], faceEntities[i], faPer, faRan, 0);
  }

  //reparam edges onto face
  int numFaces = faceMap.size();
  int numEdgesFace = 4;

  int edgeLoop[numFaces][numEdgesFace] = {{50,72,11,73},{72,48,71,16},{46,74,20,71},{52,73,6,74},{50,52,46,48},{11,16,20,6}}; 
  int edgeReparamLoop[numFaces][numEdgesFace] = {{0,3,2,1},{1,0,3,2},{0,1,2,3},{0,1,2,3},{0,1,2,3},{0,3,2,1}}; 
  
  //typedef void (*ParametricFunctionArray) (double const from[2], double to[2], void*);
  ParametricFunctionArray edgeFaceFunction[] = 
    {
      reparamREdge_0,
      reparamREdge_1,
      reparamREdge_2,
      reparamREdge_3,
    };

  for(int i=0; i<numFaces;i++)
  {
    b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_face[i]));
    for(int j=0; j<numEdgesFace;j++)
    {
      agm_use faceUse = add_adj(model, b, edgeLoop[i][j]);
      gmi_add_analytic_reparam(model, faceUse, edgeFaceFunction[edgeReparamLoop[i][j]], 0);
    }
  }


  gmi_add_analytic_cell(model,3,regionID);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,3,regionID)));
  for(int i=0; i<numFaces;i++)
  {
    agm_use regionUse = add_adj(model, b, faceMap[i]);
    gmi_add_analytic_reparam(model, regionUse, regionFunction, 0);
  }

  return;
}

void PiercingCylinder::makePiercingCylinder(gmi_model* model)
{ 

  sphereRadius = radius;
  gmi_add_analytic(model, dim-1, faceID, cylinderFace, faPer, faRan, 0);
}


void Sphere::makeSphere(gmi_model* model)
{ 

  sphereRadius = radius;
  gmi_add_analytic(model, dim-1, faceID, sphereFace, faPer, faRan, 0);
}


void setParameterizationCylinder(gmi_model* model,apf::Mesh2* m, Enclosure box, PiercingCylinder cylinder,Sphere circle1, Sphere circle2 )
{
  //Get the classification of each entity in the SimMesh
  apf::MeshEntity* ent;
  for(int i =0;i<4;i++)
  {
    apf::MeshIterator* it = m->begin(i);
    while( (ent = m->iterate(it)))
    {
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      if(modelTag > 139 && modelTag< 148)
      {
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,cylinder.faceID));
      }
      else if(modelTag >= 200){
        //std::cout<<circle1.faceID<<" "<<circle2.faceID<<" "<<modelTag<<std::endl;
        if(modelTag > 200 && modelTag <= 204)
            m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,circle1.faceID));
        else if(modelTag >= 209 && modelTag <= 212)
            m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,circle2.faceID));
      }
      else{
        //std::cout<<"belongs nowhere? "<<modelType<<" "<<modelTag<<std::endl;
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,modelTag));
      }

    }
    m->end(it);
  }
  m->setModel(model);
  m->acceptChanges();

  Reparam::reparameterizeEntities(model,m,box,cylinder,circle1,circle2);
 
}

void Reparam::reparameterizeEntities(gmi_model*model,apf::Mesh2*m,Enclosure box, PiercingCylinder cylinder, Sphere circle1, Sphere circle2){
    //reparameterizeCylinder(model,m,box, cylinder,circle1,circle2);

    //Need to set the parametric coordinates of each of the boundary vertices
    std::map<int,int> edgeParam;
    int edgeScales[12] = {0,1,0,1,0,1,0,1,2,2,2,2};
    double edgeLengths[3] = {boxLength,boxWidth,boxHeight};
    for(int i=0;i<12;i++)
    {
        edgeParam[box.edgeMap[i]] = edgeScales[i];
    }

    std::map<int,int(*)[2]> faceParam;
    int faceScales[6][2] = {{0,2},{1,2},{0,2},{1,2},{0,1},{0,1}};
    for(int i = 0; i<6;i++)
    {
        faceParam[box.faceMap[i]] = &(faceScales[i]); 
    }
  
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* ent;
    while( (ent = m->iterate(it)) )
    {
        apf::ModelEntity* g_ent = m->toModel(ent);
    
        apf::MeshEntity* ev[2];
        m->getDownward(ent,0,ev);
        int modelTag = m->getModelTag(g_ent);
        int modelType = m->getModelType(g_ent);
        if(modelType<3 && modelType!=0)
        {
            apf::Vector3 pt;
            apf::Vector3 newParam;
            m->getPoint(ent,0,pt);
            if(modelType==1)
            {
                if(modelTag == circle1.faceID){
                    double argy = (pt[1]-xyz_offset[1]);
                    double argx = (pt[0]-xyz_offset[0]);
                    if(argx == 0 && argy ==0)
                        newParam[0] = 0.0; // not sure if this will happen or if this is right
                    else 
                        newParam[0] = atan2(argy,argx);
                    m->setParam(ent,newParam);
                }
                else if(modelTag == circle2.faceID){
                    double argy = (pt[1]-xyz_offset[1]);
                    double argx = (pt[0]-xyz_offset[0]);
                    if(argx == 0 && argy ==0)
                        newParam[0] = 0.0; // not sure if this will happen or if this is right
                    else 
                        newParam[0] = atan2(argy,argx);
                    m->setParam(ent,newParam);
                }
                else{
                    int relevantIndex = edgeParam[modelTag];
                    newParam[0]=pt[relevantIndex]/edgeLengths[relevantIndex];
                    m->setParam(ent,newParam);
                }
            }
            else if (modelType==2 && modelTag!=cylinder.faceID)
            { 
                int* relevantIndex = faceParam[modelTag][0]; //size is 2
                newParam[0] = pt[relevantIndex[0]]/edgeLengths[relevantIndex[0]];
                newParam[1] = pt[relevantIndex[1]]/edgeLengths[relevantIndex[1]];
                m->setParam(ent,newParam);
            }
            else if (modelType==2 && modelTag == cylinder.faceID)
            {
                double argy = (pt[1]-xyz_offset[1]);
                double argx = (pt[0]-xyz_offset[0]);
                if(argx == 0 && argy ==0)
                    newParam[0] = 0.0; // not sure if this will happen or if this is right
                else 
                    newParam[0] = atan2(argy,argx);
                newParam[1] = pt[2]/boxHeight;

                m->setParam(ent,newParam);
            }

        } //end if
  } //end while
  m->end(it);
  m->acceptChanges();

}

void initialAdapt_analytic(){
  //at this point, hmin and hmax haven't been set yet
  //apf::Field* size_initial = apf::createLagrangeField(m,"size_initial",apf::SCALAR,1);
  lion_set_verbosity(1);
  apf::Field* size_init = apf::createLagrangeField(m,"proteus_init",apf::SCALAR,1);
  //sizeFieldList.push(size_init);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* ent;
  //double hmin = 0.0125;
  //double hmax = 0.2;
  while( (ent = m->iterate(it)) )
  {
    apf::Vector3 pt;
    m->getPoint(ent,0,pt);
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      if(modelTag == modelPiercingCylinder.faceID && modelType==2 || modelType==1 && (modelTag == modelCircle1.faceID || modelTag == modelCircle2.faceID))
          apf::setScalar(size_init,ent,0,0.4);
      else
          apf::setScalar(size_init,ent,0,0.4);
  }
  m->end(it);

  //gradeMesh(1.5);
  //isotropicIntersect();
  
  apf::writeVtkFiles("initialProteus",m);
  ma::Input* in = ma::configure(m,size_init);
  in->maximumIterations = 10;
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  in->shouldFixShape = true;
  in->debugFolder="./debug_fine";
  ma::adaptVerbose(in,true);
  m->verify();
  apf::writeVtkFiles("middleProteus",m);
  
  freeField(size_init);

  size_init = apf::createLagrangeField(m,"proteus_initial",apf::SCALAR,1);
  //sizeFieldList.push(size_init);
  it = m->begin(0);
  while( (ent = m->iterate(it)) )
  {
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      if(modelTag == modelPiercingCylinder.faceID && modelType==2 || modelType==1 && (modelTag == modelCircle1.faceID || modelTag == modelCircle2.faceID))
          apf::setScalar(size_init,ent,0,0.1025);
      else
          apf::setScalar(size_init,ent,0,0.1025);
  }
  m->end(it);

  //gradeMesh(1.5);
   
  in = ma::configure(m,size_init);
  in->maximumIterations = 10;
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  in->shouldFixShape = true;
  in->debugFolder="./debug_fine2";
  ma::adaptVerbose(in,true);
  m->verify();
  
  apf::writeVtkFiles("finalProteus",m);
  freeField(size_init);
}


void updateSphereCoordinates(double*sphereCenter)
{
  xyz_offset[0] = sphereCenter[0];
  xyz_offset[1] = sphereCenter[1];
  xyz_offset[2] = sphereCenter[2];

  char buffer[100];
  sprintf(buffer,"Checking coordinates at update %f %f %f",xyz_offset[0],xyz_offset[1],xyz_offset[2]);

}

int main(int argc, char** argv)
{
  //int comm_size = PCU_Comm_Peers();
  //int comm_rank = PCU_Comm_Self();
  //m = apf::loadMdsMesh(".null", meshFile);
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  gmi_register_null();

  m = apf::loadMdsFromGmsh(gmi_load(".null"), "CylinderBox.msh");
  m->verify();

  boxLength = 2.5;
  boxWidth = 0.41;
  boxHeight = 0.41;
  int dim = 3;
  double radius = 0.05;
  double  coordinates[] = {0.5,0.2,0.0};
  updateSphereCoordinates(&(coordinates[0]));

  //create analytic model
  gmi_model* model = gmi_make_analytic();
  
  //add sphere

  PiercingCylinder cylinder = PiercingCylinder();
  cylinder.radius = radius;

  cylinder.makePiercingCylinder(model);
  //int edgeID[2] = {200,201};
  int edgeID[2] = {9,10};
  Sphere circle1 = Sphere(2);
  circle1.faceID = edgeID[0];
  circle1.radius = radius;
  //can't use makeSphere because of paramterization function is not general
  gmi_add_analytic(model, 1, circle1.faceID, circleFace0, circle1.faPer, circle1.faRan, 0);

  Sphere circle2 = Sphere(2);
  circle2.faceID = edgeID[1];
  circle2.radius = radius;
  //circle2.makeSphere(model);
  gmi_add_analytic(model, 1, circle2.faceID, circleFace1, circle2.faPer, circle2.faRan, 0);

  //add the box
  geomDim = dim;
  Enclosure box;
  
  //add special flag to makeBox to include circles? 
  //need to classify vertices on these circles independently from the cylindrical face
  
  box.makeBox3D(model);
  //create 2 circles, parameterize it as part of the two faces 1 and 6, faces 4+5

  //face ids come from faceMap
  agm_bdry b;
  agm_use regionUse;

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,2,1)));
  regionUse = add_adj(model, b, circle1.faceID);
  gmi_add_analytic_reparam(model, regionUse, reparam_Circle0, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,2,6)));
  regionUse = add_adj(model, b, circle2.faceID);
  gmi_add_analytic_reparam(model, regionUse, reparam_Circle1, 0);

  //parameterize the cylinder to the region
  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,dim,box.regionID)));
  regionUse = add_adj(model, b, cylinder.faceID);
  gmi_add_analytic_reparam(model, regionUse, regionFunction, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,2,cylinder.faceID)));
  regionUse = add_adj(model, b, circle1.faceID);
  gmi_add_analytic_reparam(model, regionUse, reparam_CircleCylinder0, 0);
  regionUse = add_adj(model, b, circle2.faceID);
  gmi_add_analytic_reparam(model, regionUse, reparam_CircleCylinder1, 0);

  setParameterizationCylinder(model,m,box,cylinder,circle1,circle2);

  modelBox = box;
  modelCircle1 = circle1;
  modelCircle2 = circle2;
  modelPiercingCylinder = cylinder;
  m->verify();

  initialAdapt_analytic();
  return 0;

}

void freeField(apf::Field*& f)
/**
 * @brief Destroy field to avoid memory leak
 * 
 * Every field needs to eventually be destroyed. Returns nothing.
 */
{
  if (f) {
    apf::destroyField(f);
    f = 0;
  }
}



