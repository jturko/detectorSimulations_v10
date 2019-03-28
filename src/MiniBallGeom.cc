/////////////////////////////////////////////////////////////////////////////
// simple geometrical objects
////////////////////////////////////////////////////////////////////////////

#include "MiniBallGeom.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

#include <iostream>

MiniBallGeom::MiniBallGeom()
{
}

MiniBallGeom::~MiniBallGeom()
{
}

G4Transform3D MiniBallGeom::GetTrap(G4ThreeVector* x, G4double* pDz,
				    G4double* pTheta, G4double* pPhi,
				    G4double* pDy1, G4double* pDx1,
				    G4double* pDx2,
				    G4double* pAlp1,
				    G4double* pDy2, G4double* pDx3,
				    G4double* pDx4,
				    G4double* pAlp2)
{
  G4ThreeVector tr[8];
  
  G4ThreeVector midpoint
    = (x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7])/8;

  // move to origin
  for(int q=0; q<8; q++)
    tr[q] = x[q]-midpoint;

  if(GET_TRAP_DL > 0)
    {
      G4cout << "moved to origin" << G4endl;
      G4cout << tr[0] << ", " << tr[1] << ", "
	     << tr[2] << ", " << tr[3] << G4endl;
      G4cout << tr[4] << ", " << tr[5] << ", "
	     << tr[6] << ", " << tr[7] << G4endl;
    }

  MiniBallPlane lower_plane(x[0], x[1]-x[0], x[2]-x[0]);
  G4ThreeVector normal_vec = lower_plane.GetNormal();

  // move upper and lower plane parallel to z-plane
  for(int q=0; q<8; q++)
    {
    tr[q].rotateZ(-normal_vec.phi());
    tr[q].rotateY(-normal_vec.theta());
    }

  if(GET_TRAP_DL > 0)
    {
      G4cout << "parallel to z-plane" << G4endl;
      G4cout << tr[0] << ", " << tr[1] << ", "
	     << tr[2] << ", " << tr[3] << G4endl;
      G4cout << tr[4] << ", " << tr[5] << ", "
	     << tr[6] << ", " << tr[7] << G4endl;
    }

  G4ThreeVector y_direction = tr[1] - tr[0];

  // rotate y-edges to -y-direction
  for(int q=0; q<8; q++)
    tr[q].rotateZ(-y_direction.phi());

  if(GET_TRAP_DL > 0)
    {
      G4cout << "rotate to -y-direction" << G4endl;
      G4cout << tr[0] << ", " << tr[1] << ", "
	     << tr[2] << ", " << tr[3] << G4endl;
      G4cout << tr[4] << ", " << tr[5] << ", "
	     << tr[6] << ", " << tr[7] << G4endl;
    }

  G4bool spin_180 = false;
  if(tr[0].z() > 0)
    {
      spin_180 = true;
      for(int q=0; q<8; q++)
	tr[q].rotateX(M_PI);
    }

  if(GET_TRAP_DL > 0)
    {
      G4cout << "spin 180" << G4endl;
      G4cout << tr[0] << ", " << tr[1] << ", "
	     << tr[2] << ", " << tr[3] << G4endl;
      G4cout << tr[4] << ", " << tr[5] << ", "
	     << tr[6] << ", " << tr[7] << G4endl;
    }

  G4ThreeVector y_displacement(0,
		     (tr[0].y()+tr[2].y()+tr[4].y()+tr[6].y())/4, 0);

  // move to origin
  for(int q=0; q<8; q++)
    tr[q] -= y_displacement;

  G4ThreeVector x_displacement((tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
				+tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x())/8,
			       0, 0);

  // move to origin
  for(int q=0; q<8; q++)
    tr[q] -= x_displacement;


  if(GET_TRAP_DL > 0)
    {
      G4cout << "prefinal trap" << G4endl;
      G4cout << tr[0] << ", " << tr[1] << ", "
	     << tr[2] << ", " << tr[3] << G4endl;
      G4cout << tr[4] << ", " << tr[5] << ", "
	     << tr[6] << ", " << tr[7] << G4endl;
    }

  if( fabs(tr[0].z()-tr[1].z()) < TRAP_THRESH )
    tr[1].setZ(tr[0].z());
  if( fabs(tr[0].z()-tr[2].z()) < TRAP_THRESH )
    tr[2].setZ(tr[0].z());
  if( fabs(tr[0].z()-tr[3].z()) < TRAP_THRESH )
    tr[3].setZ(tr[0].z());
  if( fabs(tr[0].z()+tr[4].z()) < TRAP_THRESH )
    tr[4].setZ(-tr[0].z());
  if( fabs(tr[4].z()-tr[5].z()) < TRAP_THRESH )
    tr[5].setZ(tr[4].z());
  if( fabs(tr[4].z()-tr[6].z()) < TRAP_THRESH )
    tr[6].setZ(tr[4].z());
  if( fabs(tr[4].z()-tr[7].z()) < TRAP_THRESH )
    tr[7].setZ(tr[4].z());  
  if( fabs(tr[0].y()+tr[2].y()+tr[4].y()+tr[6].y()) < TRAP_THRESH )
    {
      tr[0].setY(-tr[2].y()-tr[4].y()-tr[6].y());
      if( fabs(tr[0].y()+tr[2].y()+tr[4].y()+tr[6].y()) > 0 )
	tr[2].setY(-tr[0].y()-tr[4].y()-tr[6].y());
      if( fabs(tr[0].y()+tr[2].y()+tr[4].y()+tr[6].y()) > 0 )
	tr[4].setY(-tr[0].y()-tr[2].y()-tr[6].y());
      if( fabs(tr[0].y()+tr[2].y()+tr[4].y()+tr[6].y()) > 0 )
	tr[6].setY(-tr[0].y()-tr[2].y()-tr[4].y());
    }
  if( fabs(tr[0].y()-tr[1].y()) < TRAP_THRESH )
    tr[1].setY(tr[0].y());
  if( fabs(tr[2].y()-tr[3].y()) < TRAP_THRESH )
    tr[3].setY(tr[2].y());
  if( fabs(tr[4].y()-tr[5].y()) < TRAP_THRESH )
    tr[5].setY(tr[4].y());
  if( fabs(tr[6].y()-tr[7].y()) < TRAP_THRESH )
    tr[7].setY(tr[6].y());
  if( fabs(tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
	   +tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x()) < TRAP_THRESH )
    {
      tr[0].setX(-tr[1].x()-tr[2].x()-tr[3].x()
		 -tr[4].x()-tr[5].x()-tr[6].x()-tr[7].x());
      if( fabs(tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
	       +tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x()) > 0 )
	tr[1].setX(-tr[0].x()-tr[2].x()-tr[3].x()
		   -tr[4].x()-tr[5].x()-tr[6].x()-tr[7].x());
      if( fabs(tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
	       +tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x()) > 0 )
	tr[2].setX(-tr[0].x()-tr[1].x()-tr[3].x()
		   -tr[4].x()-tr[5].x()-tr[6].x()-tr[7].x());
      if( fabs(tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
	       +tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x()) > 0 )
	tr[3].setX(-tr[0].x()-tr[1].x()-tr[2].x()
		   -tr[4].x()-tr[5].x()-tr[6].x()-tr[7].x());
      if( fabs(tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
	       +tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x()) > 0 )
	tr[4].setX(-tr[0].x()-tr[1].x()-tr[2].x()
		   -tr[3].x()-tr[5].x()-tr[6].x()-tr[7].x());
      if( fabs(tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
	       +tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x()) > 0 )
	tr[5].setX(-tr[0].x()-tr[1].x()-tr[2].x()
		   -tr[3].x()-tr[4].x()-tr[6].x()-tr[7].x());
      if( fabs(tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
	       +tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x()) > 0 )
	tr[6].setX(-tr[0].x()-tr[1].x()-tr[2].x()
		   -tr[3].x()-tr[4].x()-tr[5].x()-tr[7].x());
      if( fabs(tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
	       +tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x()) > 0 )
	tr[7].setX(-tr[0].x()-tr[1].x()-tr[2].x()
		   -tr[3].x()-tr[4].x()-tr[5].x()-tr[6].x());
    }

  if(GET_TRAP_DL > 0)
    {
      G4cout << "final trap" << G4endl;
      G4cout << tr[0] << ", " << tr[1] << ", "
	     << tr[2] << ", " << tr[3] << G4endl;
      G4cout << tr[4] << ", " << tr[5] << ", "
	     << tr[6] << ", " << tr[7] << G4endl;
      G4cout << "1: " << tr[0].z()-tr[1].z() << G4endl;
      G4cout << "2: " << tr[0].z()-tr[2].z() << G4endl;
      G4cout << "3: " << tr[0].z()-tr[3].z() << G4endl;
      G4cout << "4: " << tr[4].z()-tr[5].z() << G4endl;
      G4cout << "5: " << tr[4].z()-tr[6].z() << G4endl;
      G4cout << "6: " << tr[4].z()-tr[7].z() << G4endl;
      G4cout << "7: " << tr[0].z()+tr[4].z() << G4endl;
      G4cout << "8: " << tr[0].y()-tr[1].y() << G4endl;
      G4cout << "9: " << tr[2].y()-tr[3].y() << G4endl;
      G4cout << "10: " << tr[4].y()-tr[5].y() << G4endl;
      G4cout << "11: " << tr[6].y()-tr[7].y() << G4endl;
      G4cout << "12: " << tr[0].y()+tr[2].y()+tr[4].y()+tr[6].y() << G4endl;
      G4cout << "13: " << tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x()
	     +tr[4].x()+tr[5].x()+tr[6].x()+tr[7].x()<< G4endl;
    }

  G4double Talpha1, Talpha2, TthetaCphi, TthetaSphi;

  *pDz = (tr[7]).z() ;
      
  *pDy1 = ((tr[2]).y()-(tr[1]).y())*0.5;
  *pDx1 = ((tr[1]).x()-(tr[0]).x())*0.5;
  *pDx2 = ((tr[3]).x()-(tr[2]).x())*0.5;
  Talpha1 = ((tr[2]).x()+(tr[3]).x()-(tr[1]).x()-(tr[0]).x())*0.25/(*pDy1);
  *pAlp1 = atan(Talpha1);

  *pDy2 = ((tr[6]).y()-(tr[5]).y())*0.5;
  *pDx3 = ((tr[5]).x()-(tr[4]).x())*0.5;
  *pDx4 = ((tr[7]).x()-(tr[6]).x())*0.5;
  Talpha2 = ((tr[6]).x()+(tr[7]).x()-(tr[5]).x()-(tr[4]).x())*0.25/(*pDy2);
  *pAlp2 = atan(Talpha2);
  
  TthetaCphi = ((tr[4]).x()+(*pDy2)*Talpha2+(*pDx3))/(*pDz);
  TthetaSphi = ((tr[4]).y()+(*pDy2))/(*pDz);

  if(TthetaCphi!=0 && TthetaSphi!=0)
    {
      *pPhi = atan(TthetaSphi/TthetaCphi);
      *pTheta = atan(TthetaSphi/sin(*pPhi));
    }
  else if(TthetaCphi==0 && TthetaSphi==0)
    {
      *pPhi = 0;
      *pTheta = 0;
    }
  else if(TthetaCphi==0)
    {
      *pPhi = M_PI/2;
      *pTheta = atan(TthetaSphi);
    }
  else
    {
      *pPhi = 0;
      *pTheta = atan(TthetaCphi);
    }

  if(GET_TRAP_DL > 0)
    {
      G4cout << "Dz: " << *pDz << G4endl;
      G4cout << "Theta: " << *pTheta << G4endl;
      G4cout << "Phi: " << *pPhi << G4endl;
      G4cout << "Dy1: " << *pDy1 << G4endl;
      G4cout << "Dx1: " << *pDx1 << G4endl;
      G4cout << "Dx2: " << *pDx2 << G4endl;
      G4cout << "Alp1: " << *pAlp1 << G4endl;
      G4cout << "Dy2: " << *pDy2 << G4endl;
      G4cout << "Dx3: " << *pDx3 << G4endl;
      G4cout << "Dx4: " << *pDx4 << G4endl;
      G4cout << "Alp2: " << *pAlp2 << G4endl;
    }

  G4RotationMatrix rot1;
  G4ThreeVector transl1 = y_displacement+x_displacement;
  G4Transform3D transf1(rot1, transl1);

  G4RotationMatrix rot1a;
  if(spin_180)
    rot1a.rotateX(-M_PI);
  G4ThreeVector transl1a(0, 0, 0);
  G4Transform3D transf1a(rot1a, transl1a);

  G4RotationMatrix rot2;
  rot2.rotateZ(y_direction.phi());
  rot2.rotateY(normal_vec.theta());
  rot2.rotateZ(normal_vec.phi());
  G4ThreeVector transl2(0, 0, 0);
  G4Transform3D transf2(rot2, transl2);

  G4RotationMatrix rot3;
  G4ThreeVector transl3 = midpoint;
  G4Transform3D transf3(rot3, transl3);

  G4Transform3D ret = transf3*transf2*transf1a*transf1;

  return(ret);
}

G4Transform3D MiniBallGeom::GetTrap(MiniBallPlane minus_x,
				    MiniBallPlane plus_x, 
				    MiniBallPlane minus_y,
				    MiniBallPlane plus_y,
				    MiniBallPlane minus_z,
				    MiniBallPlane plus_z,
				    G4double* pDz,
				    G4double* pTheta, G4double* pPhi,
				    G4double* pDy1, G4double* pDx1,
				    G4double* pDx2,
				    G4double* pAlp1,
				    G4double* pDy2, G4double* pDx3,
				    G4double* pDx4,
				    G4double* pAlp2)
{
  G4ThreeVector x[8];

  x[0]=minus_z.Intersect(minus_y, minus_x);
  x[1]=minus_z.Intersect(minus_y, plus_x);
  x[2]=minus_z.Intersect(plus_y, minus_x);
  x[3]=minus_z.Intersect(plus_y, plus_x);
  x[4]=plus_z.Intersect(minus_y, minus_x);
  x[5]=plus_z.Intersect(minus_y, plus_x);
  x[6]=plus_z.Intersect(plus_y, minus_x);
  x[7]=plus_z.Intersect(plus_y, plus_x);

  return(GetTrap(x, pDz, pTheta, pPhi,
		 pDy1, pDx1, pDx2, pAlp1,
		 pDy2, pDx3, pDx4, pAlp2));
}

void MiniBallGeom::CleanVector(G4ThreeVector* vec)
{
  if(fabs(vec->x()) < MINIBALL_GEOM_ZERO)
    vec->setX(0.);
  if(fabs(vec->y()) < MINIBALL_GEOM_ZERO)
    vec->setY(0.);
  if(fabs(vec->z()) < MINIBALL_GEOM_ZERO)
    vec->setZ(0.);
}

MiniBallLine::MiniBallLine(G4ThreeVector p_position, G4ThreeVector p_direction)
  : position(p_position), direction(p_direction)
{
}

MiniBallLine::MiniBallLine()
{
  G4ThreeVector p_position(0, 0, 0);
  position = p_position;
  G4ThreeVector p_direction(1, 1, 1);
  direction = p_direction;
}

MiniBallLine::~MiniBallLine()
{
}

void MiniBallLine::Normalize()
{
  G4double sum_sqr = direction.x() * direction.x()
                         + direction.y() * direction.y()
                         + direction.z() * direction.z();

  G4double scalar = direction.x()*position.x()
                         + direction.y()*position.y()
                         + direction.z()*position.z();

  position = -scalar/sum_sqr*direction + position;
}

std::ostream & operator<< (std::ostream & os, const MiniBallLine & line) 
{
  G4ThreeVector pos = line.GetPosition();
  G4ThreeVector dir = line.GetDirection();
  return os << "s*" << dir << " + " << pos; 
}

MiniBallPlane::MiniBallPlane(G4ThreeVector p_position, 
    G4ThreeVector p_direction1, G4ThreeVector p_direction2)
  : position(p_position), direction1(p_direction1), direction2(p_direction2)
{
}

MiniBallPlane::MiniBallPlane(G4ThreeVector p_position,
			     G4ThreeVector p_normal)
  : position(p_position)
{
  G4RotationMatrix transf;
  transf.rotateY(p_normal.theta());
  transf.rotateZ(p_normal.phi());

  direction1=transf*G4ThreeVector(1, 0, 0);
  direction2=transf*G4ThreeVector(0, 1, 0);

  MiniBallGeom::CleanVector(&direction1);
  MiniBallGeom::CleanVector(&direction2);
}

MiniBallPlane::MiniBallPlane()
{
  G4ThreeVector p_position(0, 0, 0);
  position = p_position;
  G4ThreeVector p_direction1(1, 0, 0);
  direction1 = p_direction1;
  G4ThreeVector p_direction2(0, 1, 0);
  direction2 = p_direction2;
}

MiniBallPlane::~MiniBallPlane()
{
}

MiniBallPlane MiniBallPlane::operator * (const G4RotationMatrix& rot)
{
  G4ThreeVector new_position=rot*position;
  MiniBallGeom::CleanVector(&new_position);
  G4ThreeVector new_direction1=rot*direction1;
  MiniBallGeom::CleanVector(&new_direction1);
  G4ThreeVector new_direction2=rot*direction2;
  MiniBallGeom::CleanVector(&new_direction2);

  return(MiniBallPlane(new_position, new_direction1, new_direction2));
}

MiniBallPlane & MiniBallPlane::operator *= (const G4RotationMatrix& rot)
{
  position=rot*position;
  direction1=rot*direction1;
  direction2=rot*direction2;
  return *this;
}

MiniBallPlane & MiniBallPlane::transform(const G4RotationMatrix & rot)
{
  position=rot*position;
  direction1=rot*direction1;
  direction2=rot*direction2;
  return *this;
}

MiniBallLine MiniBallPlane::Intersect(MiniBallPlane second_plane)
{
  G4double a[5][3];
  MiniBallLine ret;
  
  a[0][0] = direction1.x();
  a[0][1] = direction1.y();
  a[0][2] = direction1.z();
  a[1][0] = direction2.x();
  a[1][1] = direction2.y();
  a[1][2] = direction2.z();

  G4ThreeVector scnd_dir1 = second_plane.GetDirection1();
  a[2][0] = -scnd_dir1.x();
  a[2][1] = -scnd_dir1.y();
  a[2][2] = -scnd_dir1.z();

  G4ThreeVector scnd_dir2 = second_plane.GetDirection2();
  a[3][0] = -scnd_dir2.x();
  a[3][1] = -scnd_dir2.y();
  a[3][2] = -scnd_dir2.z();

  G4ThreeVector scnd_pos = second_plane.GetPosition();
  a[4][0] = position.x() - scnd_pos.x();
  a[4][1] = position.y() - scnd_pos.y();
  a[4][2] = position.z() - scnd_pos.z();

  if(a[0][0] == 0)
    {
      G4int swap_row;
      if(a[0][1] != 0)
	swap_row = 1;
      else if(a[0][2] != 0)
	swap_row = 2;
      else
	{
	  G4cerr << "MiniBallPlane::Intersect: Illegal plane structure"
		 << G4endl;
	  return(ret);
	}
      
      G4double tmp;
      for(G4int q=0; q<5; q++)
	{
	  tmp = a[q][0];
	  a[q][0] = a[q][swap_row];
	  a[q][swap_row] = tmp;
	}
    }

  G4double factor = a[0][0];
  for(G4int q=0; q<5; q++)
    a[q][0] /= factor;

  if(a[0][1] != 0)
    {
      factor = a[0][1];
      for(G4int q=0; q<5; q++)
	a[q][1] = a[q][0] - a[q][1]/factor;
    }

  if(a[0][2] != 0)
    {
      factor = a[0][2];
      for(G4int q=0; q<5; q++)
	a[q][2] = a[q][0] - a[q][2]/factor;
    }

  if(a[1][1] == 0)
    {
      if(a[1][2] == 0)
	{
	  G4cerr << "MiniBallPlane::Intersect: calculation error" << G4endl;
	  return(ret);
	}
      G4double tmp;
      for(G4int q=0; q<5; q++)
	{
	  tmp = a[q][1];
	  a[q][1] = a[q][2];
	  a[q][2] = tmp;
	}
    }

  factor = a[1][1];
  for(G4int q=1; q<5; q++)
    a[q][1] /= factor;

  if(a[1][2] != 0)
    {
      factor = a[1][2];
      for(G4int q=1; q<5; q++)
	a[q][2] = a[q][1] - a[q][2]/factor;
    }

  if(a[2][2] == 0 && a[3][2] == 0)
    {
      G4cerr << "MiniBallPlane::Intersect: calculation error 2" << G4endl;
      return(ret);
    }

  if(fabs(a[2][2]) > fabs(a[3][2]) && a[2][2] != 0)
    {
      factor = a[2][2];
      for(G4int q=2; q<5; q++)
	a[q][2] /= factor;

      ret.SetPosition(scnd_pos-a[4][2]*scnd_dir1);
      ret.SetDirection(scnd_dir2-a[3][2]*scnd_dir1);
    }
  else
    {
      factor = a[3][2];
      for(G4int q=2; q<5; q++)
	a[q][2] /= factor;

      ret.SetPosition(scnd_pos-a[4][2]*scnd_dir2);
      ret.SetDirection(scnd_dir1-a[2][2]*scnd_dir2);
    }
    
  return(ret);
}

G4ThreeVector MiniBallPlane::Intersect(MiniBallPlane plane1,
				       MiniBallPlane plane2)
{
  MiniBallLine intersection_line = this->Intersect(plane1);
  return(plane2.Intersect(intersection_line));
}

G4ThreeVector MiniBallPlane::Intersect(MiniBallLine line)
{
  G4double a[4][3];
  G4ThreeVector ret(0, 0, 0);

  a[0][0] = direction1.x();
  a[0][1] = direction1.y();
  a[0][2] = direction1.z();
  a[1][0] = direction2.x();
  a[1][1] = direction2.y();
  a[1][2] = direction2.z();

  G4ThreeVector line_dir = line.GetDirection();
  a[2][0] = -line_dir.x();
  a[2][1] = -line_dir.y();
  a[2][2] = -line_dir.z();

  G4ThreeVector line_pos = line.GetPosition();
  a[3][0] = position.x() - line_pos.x();
  a[3][1] = position.y() - line_pos.y();
  a[3][2] = position.z() - line_pos.z();

  if(a[0][0] == 0)
    {
      G4int swap_row;
      if(a[0][1] != 0)
	swap_row = 1;
      else if(a[0][2] != 0)
	swap_row = 2;
      else
	{
	  G4cerr << "MiniBallPlane::Intersect: Illegal plane structure"
		 << G4endl;
	  return(ret);
	}
      
      G4double tmp;
      for(G4int q=0; q<4; q++)
	{
	  tmp = a[q][0];
	  a[q][0] = a[q][swap_row];
	  a[q][swap_row] = tmp;
	}
    }

  G4double factor = a[0][0];
  for(G4int q=0; q<4; q++)
    a[q][0] /= factor;

  if(a[0][1] != 0)
    {
      factor = a[0][1];
      for(G4int q=0; q<4; q++)
	a[q][1] = a[q][0] - a[q][1]/factor;
    }

  if(a[0][2] != 0)
    {
      factor = a[0][2];
      for(G4int q=0; q<4; q++)
	a[q][2] = a[q][0] - a[q][2]/factor;
    }

  if(a[1][1] == 0)
    {
      if(a[1][2] == 0)
	{
	  G4cerr << "MiniBallPlane::Intersect: calculation error" << G4endl;
	  return(ret);
	}
      G4double tmp;
      for(G4int q=0; q<4; q++)
	{
	  tmp = a[q][1];
	  a[q][1] = a[q][2];
	  a[q][2] = tmp;
	}
    }

  factor = a[1][1];
  for(G4int q=1; q<4; q++)
    a[q][1] /= factor;

  if(a[1][2] != 0)
    {
      factor = a[1][2];
      for(G4int q=1; q<4; q++)
	a[q][2] = a[q][1] - a[q][2]/factor;
    }

  if(a[2][2] == 0)
    {
      G4cerr << "MiniBallPlane::Intersect: calculation error" << G4endl;
      return(ret);
    }

  ret = line_pos - a[3][2]/a[2][2]*line_dir;
  return(ret);
}

G4ThreeVector MiniBallPlane::GetNormal()
{
  G4double a[3][2];

  G4double c1, c2, c3;
  
  a[0][0] = direction1.x();
  a[1][0] = direction1.y();
  a[2][0] = direction1.z();

  a[0][1] = direction2.x();
  a[1][1] = direction2.y();
  a[2][1] = direction2.z();

  if(a[0][0] == 0 && a[0][1] == 0)
    return(G4ThreeVector(1, 0, 0));

  if(a[0][0] == 0)
    {
      G4double tmp;
      for(int q=0; q<3; q++)
	{
	  tmp = a[q][0];
	  a[q][0] = a[q][1];
	  a[q][1] = tmp;
	}
    }

  G4double factor = a[0][0];
  for(int q=0; q<3; q++)
    a[q][0] /= factor;

  if(a[0][1] != 0)
    {
      factor = a[0][1];
      for(int q=0; q<3; q++)
	a[q][1] = a[q][0] - a[q][1]/factor;
    }

  if(a[1][1] == 0)
    {
      c2 = 1;
      c3 = 0;
    }
  else
    {
      c3 = 1;
      c2 = -a[2][1]/a[1][1];
    }

  c1 = -a[1][0]*c2 - a[2][0]*c3;

  G4double mag = sqrt(c1*c1 + c2*c2 + c3*c3);
  G4ThreeVector ret(c1/mag, c2/mag, c3/mag);

  return(ret);
}

G4double MiniBallPlane::GetDistance(G4ThreeVector point)
{
  MiniBallLine normal_line(point, GetNormal());

  G4ThreeVector lot = Intersect(normal_line);

  G4ThreeVector dist = point - lot;

  return(dist.mag());
}
  
std::ostream & operator<< (std::ostream & os, const MiniBallPlane & plane) 
{
  G4ThreeVector pos = plane.GetPosition();
  G4ThreeVector dir1 = plane.GetDirection1();
  G4ThreeVector dir2 = plane.GetDirection2();
  return os << "s*" << dir1 << " + r*" << dir2 << " + " << pos; 
}

