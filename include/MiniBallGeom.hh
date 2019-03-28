/************************************************************************
 * \file MiniBallGeom.hh
 *
 * \brief simple geometrical objects
 *
 * \author hans.boie@mpi-hd.mpg.de
 *
 ************************************************************************/

#ifndef MiniBallGeom_H
#define MiniBallGeom_H 1

#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#define TRAP_THRESH 1e-10
#define MINIBALL_GEOM_ZERO 1e-10
#define GET_TRAP_DL 0

/************************************************************************
 * \class MiniBallLine
 *
 * \brief mathematical straight line
 *
 ************************************************************************/
class MiniBallLine
{

public:
  MiniBallLine();
  MiniBallLine(G4ThreeVector position, G4ThreeVector direction);
  ~MiniBallLine();

public:
  void SetPosition(G4ThreeVector p_position) { position = p_position; };
  void SetDirection(G4ThreeVector p_direction) { direction = p_direction; };

  G4ThreeVector GetPosition() const { return(position); };
  G4ThreeVector GetDirection() const { return(direction); };

  void Normalize();

  friend std::ostream & operator << (std::ostream & os,
				     const MiniBallLine & line);

private:
  G4ThreeVector position;
  G4ThreeVector direction;

};

/************************************************************************
 * \class MiniBallLine
 *
 * \brief mathematical plane
 *
 ************************************************************************/
class MiniBallPlane 
{
  
public:
  MiniBallPlane();
  MiniBallPlane(G4ThreeVector position, G4ThreeVector direction1,
		G4ThreeVector direction2);
  MiniBallPlane(G4ThreeVector position, G4ThreeVector normal);

  ~MiniBallPlane();

public:
  void SetPosition(G4ThreeVector p_position) { position = p_position; };
  void SetDirection1(G4ThreeVector p_direction1) 
                              { direction1 = p_direction1; };
  void SetDirection2(G4ThreeVector p_direction2) 
                              { direction2 = p_direction2; };

  G4ThreeVector GetPosition() const { return(position); };
  G4ThreeVector GetDirection1() const { return(direction1); };
  G4ThreeVector GetDirection2() const { return(direction2); };

  MiniBallPlane operator * (const G4RotationMatrix &);
  MiniBallPlane& operator *= (const G4RotationMatrix &);
  MiniBallPlane& transform(const G4RotationMatrix &);

  G4ThreeVector GetNormal();
  G4double GetDistance(G4ThreeVector point);
  
  MiniBallLine Intersect(MiniBallPlane second_plane);
  G4ThreeVector Intersect(MiniBallPlane second_plane,
			  MiniBallPlane thrid_plane);
  G4ThreeVector Intersect(MiniBallLine line);

  friend std::ostream & operator << (std::ostream & os,
				     const MiniBallPlane & plane);

private:
  G4ThreeVector position;
  G4ThreeVector direction1;
  G4ThreeVector direction2;

};

/************************************************************************
 * \class MiniBallGeom
 *
 * \brief construction of a G4Trap solid from intersecting planes
 *
 ************************************************************************/
class MiniBallGeom
{

public:
  MiniBallGeom();
  ~MiniBallGeom();

public:
  static G4Transform3D GetTrap(G4ThreeVector x[8], G4double* pDz,
			       G4double* pTheta, G4double* pPhi,
			       G4double* pDy1, G4double* pDx1, G4double* pDx2,
			       G4double* pAlp1,
			       G4double* pDy2, G4double* pDx3, G4double* pDx4,
			       G4double* pAlp2);

  static G4Transform3D GetTrap(MiniBallPlane minus_x, MiniBallPlane plus_x,
			       MiniBallPlane minus_y, MiniBallPlane plus_y,
			       MiniBallPlane minus_z, MiniBallPlane plus_z,
			       G4double* pDz,
			       G4double* pTheta, G4double* pPhi,
			       G4double* pDy1, G4double* pDx1, G4double* pDx2,
			       G4double* pAlp1,
			       G4double* pDy2, G4double* pDx3, G4double* pDx4,
			       G4double* pAlp2);

  static void CleanVector(G4ThreeVector* vec);
};

#endif
