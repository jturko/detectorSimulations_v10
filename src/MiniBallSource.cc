/************************************************************************
 * \file MiniBallSource.hh
 *
 * \class MiniBallSource
 * \brief mother class for referrence sources
 *
 * \author hans.boie@mpi-hd.mpg.de
 *
 ************************************************************************/

#include "MiniBallSource.hh"

#include "G4Gamma.hh"
#include "Randomize.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


#define THETA_MAX M_PI

MiniBallSource::MiniBallSource(G4double x_pos, G4double y_pos, G4double z_pos)
  : source_x(x_pos), source_y(y_pos), source_z(z_pos),
    cleared_statistics(false), checked_settings(false),
    angular_distribution(false), events(0L), ang_err_msg(false)
{
  G4int n_particle = 1;

  particleGun = new G4ParticleGun(n_particle);

  particleGun->SetParticleDefinition(G4Gamma::GammaDefinition());
}

MiniBallSource::~MiniBallSource()
{}

void MiniBallSource::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->SetParticlePosition(G4ThreeVector(source_x, source_y,
						 source_z));

  GenerateDecay();

  for(G4int gamma_nb=0; gamma_nb<GetNbOfGammas(); gamma_nb++)
    {
      particleGun->SetParticleEnergy(GetGammaEnergy(gamma_nb));
      particleGun->SetParticleMomentumDirection(GetGammaDirection(gamma_nb));
      particleGun->GeneratePrimaryVertex(anEvent);
    }
}

void MiniBallSource::GenerateDecay()
{
  if(!cleared_statistics)
    ClearStatistics();

  if(!checked_settings)
    CheckSettings();

  // clear old values
  gamma_energy.erase(gamma_energy.begin(),gamma_energy.end());
  gamma_direction.erase(gamma_direction.begin(),gamma_direction.end());
  gamma_transform.erase(gamma_transform.begin(),gamma_transform.end());

  decay_level_history.erase(decay_level_history.begin(),
			    decay_level_history.end());
  decay_type_history.erase(decay_type_history.begin(),
			    decay_type_history.end());

  // get populated level
  size_t nb_of_feeded_levels = feeding->decay.size();
  G4double summed_prob=0;

  for(size_t q=0; q<nb_of_feeded_levels; q++)
    summed_prob += feeding->prob[q];

  G4double decay = G4RandFlat::shoot(0., summed_prob);

  summed_prob=feeding->prob[0];
  size_t q=0;
  while(q<(nb_of_feeded_levels-1) && decay>summed_prob)
    {
      q++;
      summed_prob += feeding->prob[q];
    }

  G4int populated_level = feeding->decay[q];

  if(populated_level<0 || populated_level>=(G4int)level_structure.size())
    {
      G4cerr << "trying to populate non existing level " << populated_level
	     << G4endl;
      return;
    }

  level_structure[populated_level]->pop_stat++;
  level_structure[populated_level]->feed_stat++;

  if(MB_SOURCE_DL>0)
    std::cout << "populated level id: " << populated_level
	 << " energy: " << level_structure[populated_level]->energy/keV
	 << "keV" << G4endl;


  while(level_structure[populated_level]->energy > 0.*keV)
    {
      populated_level = DecayLevel(populated_level);
      if(populated_level<0 || populated_level>=(G4int)level_structure.size())
	{
	  G4cerr << "trying to populate non existing level " << populated_level
		 << G4endl;
	  return;
	}
    }

  events++;
}

void MiniBallSource::GenerateDecay(G4int populated_level)
{
  if(!cleared_statistics)
    ClearStatistics();

  if(!checked_settings)
    CheckSettings();

  // clear old values
  gamma_energy.erase(gamma_energy.begin(),gamma_energy.end());
  gamma_direction.erase(gamma_direction.begin(),gamma_direction.end());
  gamma_transform.erase(gamma_transform.begin(),gamma_transform.end());

  decay_level_history.erase(decay_level_history.begin(),
			    decay_level_history.end());
  decay_type_history.erase(decay_type_history.begin(),
			    decay_type_history.end());

  if(populated_level<0 || populated_level>=(G4int)level_structure.size())
    {
      G4cerr << "trying to populate non existing level " << populated_level
	     << G4endl;
      return;
    }

  level_structure[populated_level]->pop_stat++;
  level_structure[populated_level]->feed_stat++;

  if(MB_SOURCE_DL>0)
    std::cout << "populated level id: " << populated_level
	 << " energy: " << level_structure[populated_level]->energy/keV
	 << "keV" << G4endl;


  while(level_structure[populated_level]->energy > 0.*keV)
    {
      populated_level = DecayLevel(populated_level);
      if(populated_level<0 || populated_level>=(G4int)level_structure.size())
	{
	  G4cerr << "trying to populate non existing level " << populated_level
		 << G4endl;
	  return;
	}
    }

  events++;
}

G4int MiniBallSource::DecayLevel(G4int level_id)
{
  if(level_id<0 || level_id>=(G4int)level_structure.size())
    {
      G4cerr << "trying to decay from non existing level " << level_id
	     << G4endl;
      return(0);
    }

  size_t nb_of_decay_levels = level_structure[level_id]->decay.size();
  G4double summed_prob=0;

  for(size_t q=0; q<nb_of_decay_levels; q++)
    summed_prob += level_structure[level_id]->prob[q];

  G4double decay = G4RandFlat::shoot(0., summed_prob);

  if(MB_SOURCE_DL>0)
    std::cout<<"DecayLevel("<<level_id<<"): nb_of_decay_levels = "<<nb_of_decay_levels<<", decay = "<<decay<<" (0-"<<summed_prob<<")"<<std::endl;

  summed_prob=level_structure[level_id]->prob[0];
  size_t q=0;
  while(q<(nb_of_decay_levels-1) && decay>summed_prob)
    {
      if(MB_SOURCE_DL>0)
	std::cout<<"q = "<<q<<", summed_prob = "<<summed_prob<<std::endl;
      q++;
      summed_prob += level_structure[level_id]->prob[q];
    }

  G4int new_level_id = level_structure[level_id]->decay[q];

  if(MB_SOURCE_DL>0)
    std::cout<<"chose decay["<<q<<"] for level_structure["<<level_id<<"] (= "<<level_structure[level_id]->decay[q]<<")"<<std::endl;

  if(new_level_id<0 || new_level_id>=(G4int)level_structure.size())
    {
      G4cerr << "trying to decay to non existing level " << new_level_id
	     << G4endl;
      return(0);
    }
  
  level_structure[new_level_id]->feed_stat++;

  G4double decay_energy
    = level_structure[level_id]->energy-level_structure[new_level_id]->energy;

  if(decay_energy<0)
    {
      G4cerr << "negative decay energy from level " << level_id
	     << " to level " << new_level_id << G4endl;
      return(0);
    }

  G4double internal_conversion=0;
  if(level_structure[level_id]->decay_type[q]==E1)
    internal_conversion = GetInternalConversionCoeffE1(decay_energy);
  else if(level_structure[level_id]->decay_type[q]==E2)
    internal_conversion = GetInternalConversionCoeffE2(decay_energy);

  G4bool generate_gamma = true;
  if(internal_conversion>0)
    {
      G4double conv_prob = G4RandFlat::shoot(0., 1+internal_conversion);
      if(conv_prob>1)
	generate_gamma = false;
    }

  G4int decay_type;
  if(generate_gamma)
    {
      if(MB_SOURCE_DL>0)
	std::cout << "decayed to level id: " << new_level_id
	     << " energy: " << level_structure[new_level_id]->energy/keV
	     << "keV decay energy: " << decay_energy/keV << "keV" << G4endl;
      
      decay_type = GetDecayType(level_structure[level_id]->decay_type[q],
				level_structure[level_id]->delta[q]);

      G4double phi, theta;

      // get random phi  
      phi = G4RandFlat::shoot(-M_PI, M_PI);

      if(!angular_distribution || gamma_transform.size()==0)
	{
	  //generate isotropic direction
	  G4double cos_theta = G4RandFlat::shoot(-1., 1.);
	  theta = acos(cos_theta);
	}
      else
	theta
	  = GetAngularDistribution(level_id, new_level_id,
				   decay_type);

      G4RotationMatrix transform;
      transform.rotateY(theta);
      transform.rotateZ(phi);
  
      G4ThreeVector direction(0., 0., 1);
      direction.transform(transform);

      for(G4int trans_id=(G4int)gamma_transform.size()-1;
	  trans_id>=0;trans_id--)
	direction.transform(gamma_transform[trans_id]);

      gamma_energy.push_back(decay_energy);
      gamma_direction.push_back(direction);
      gamma_transform.push_back(transform);

      level_structure[level_id]->decay_stat[q]++;      
    }
  else
    {
      level_structure[level_id]->conv_stat[q]++;      

      if(MB_SOURCE_DL>0)
	std::cout << "decayed to level id: " << new_level_id
	     << " energy: " << level_structure[new_level_id]->energy/keV
	     << "keV internal conversion" << G4endl;

      decay_type = INT_CONV;
    }

  //fill history
  decay_level_history.push_back(level_id);    
  decay_type_history.push_back(decay_type);

  return(new_level_id);
}

G4double MiniBallSource::GetAngularDistribution(G4int level_id,
						G4int new_level_id,
						G4int decay_type)
{
  size_t history_last = decay_level_history.size()-1;

  if(level_structure[new_level_id]->energy > 0.*keV)
    {
      if(!ang_err_msg)
	{
	  G4cout << "gamma-gamma angular distribution for more than two "
		 << "gammas not implemented, using isotropic distribution"
		 << G4endl;
	  ang_err_msg=true;
	}
      G4double cos_theta = G4RandFlat::shoot(-1., 1.);
      G4double theta = acos(cos_theta);
      return(theta);
    }

  if(decay_type == UNKNOWN || decay_type_history[history_last] == UNKNOWN
     || decay_type_history[history_last] == INT_CONV)
    {
      G4double cos_theta = G4RandFlat::shoot(-1., 1.);
      G4double theta = acos(cos_theta);
      return(theta);
    }      

  G4double I_i = level_structure[decay_level_history[history_last]]->l;
  G4double I = level_structure[level_id]->l;
  G4double I_f = level_structure[new_level_id]->l;

  G4double P_i = level_structure[decay_level_history[history_last]]->p;
  G4double P = level_structure[level_id]->p;
  G4double P_f = level_structure[new_level_id]->p;

  G4double L1, L2;
  G4bool P1, P2;
  GetDecayValues(decay_type_history[history_last], &L1, &P1);
  GetDecayValues(decay_type, &L2, &P2);

  if(I_i > I + L1 || I_i < fabs(I - L1) )
    {
      G4cerr << "inconsistent decay values (I_i=" << I_i
	     << ", I=" << I << ", L1=" << L1 << ")" << G4endl;
      return(0.);
    }
  if(I > I_f + L2 || I < fabs(I_f - L2) )
    {
      G4cerr << "inconsistent decay values (I=" << I
	     << ", I_f=" << I_f << ", L2=" << L2 << ")" << G4endl;
      return(0.);
    }
  if(P1 && P_i==P)
    {
      G4cerr << "inconsistent decay values (P_i=" << P_i
	     << ", P=" << P << ", P1=true)" << G4endl;
      return(0.);
    }
  if(!P1 && P_i!=P)
    {
      G4cerr << "inconsistent decay values (P_i=" << P_i
	     << ", P=" << P << ", P1=false)" << G4endl;
      return(0.);
    }
  if(P2 && P==P_f)
    {
      G4cerr << "inconsistent decay values (P=" << P
	     << ", P_f=" << P_f << ", P2=true)" << G4endl;
      return(0.);
    }
  if(!P2 && P!=P_f)
    {
      G4cerr << "inconsistent decay values (P=" << P
	     << ", P_f=" << P_f << ", P2=false)" << G4endl;
      return(0.);
    }

  std::vector<MiniBallSourceDecayBranch> branches;

  for(G4double m_i=-I_i; m_i<=I_i; m_i+=1.)
    for(G4double m1=-1; m1<=1; m1+=2)
      for(G4double m_f=-I_f; m_f<=I_f; m_f+=1.)
	  {
	    G4double m = m_i-m1;
	    G4double m2 = m-m_f;

	    if(m<-I || m>I)
	      continue;
	    if(m2<-L2 || m2>L2)
	      continue;

	    G4double CG1 = GetClebschGordan(I, m, L1, m1, I_i, m_i);
	    G4double CG2 = GetClebschGordan(I_f, m_f, L2, m2, I, m);
	    G4double F_L1 = GetF((G4int)L1, (G4int)m1, 0.);
	    G4double prob = CG1*CG1*F_L1*CG2*CG2;

	    branches.push_back(MiniBallSourceDecayBranch(m_i, m_f, m,
							 m1, m2, prob));

	    if(MB_SOURCE_DL>1)
	      std::cout << "I_i=" << I_i << ", m_i=" << m_i
		   << ", I=" << I << ", m=" << m
		   << ", I_f=" << I_f << ", m_f=" << m_f
		   << ", L1=" << L1 << ", m1=" << m1
		   << ", L2=" << L2 << ", m2=" << m2
		   << ", CG1=" << CG1 << ", CG2=" << CG2
		   << ", F_L1=" << F_L1 << ", value=" << prob << G4endl;
	  }

  G4double summed_prob=0.;
  for(size_t branch_nb=0; branch_nb<branches.size(); branch_nb++)
    summed_prob += branches[branch_nb].prob;

  G4double rand_prob = G4RandFlat::shoot(0., summed_prob);

  summed_prob=branches[0].prob;
  size_t q=0;
  while(q<(branches.size()-1) && rand_prob>summed_prob)
    {
      q++;
      summed_prob += branches[q].prob;
    }
  
  G4double theta = RandTheta((G4int)L2, (G4int)branches[q].m2);

  if(MB_SOURCE_DL>1)
    std::cout << "decay: I_i=" << I_i << ", m_i=" << branches[q].m_i
	 << ", I=" << I << ", m=" << branches[q].m
	 << ", I_f=" << I_f << ", m_f=" << branches[q].m_f
	 << ", L1=" << L1 << ", m1=" << branches[q].m1
	 << ", L2=" << L2 << ", m2=" << branches[q].m2
	 << ", theta=" << theta << G4endl;
  
  return(theta);
}

G4int MiniBallSource::GetDecayType(G4int decay_type, G4double delta)
{
  if(decay_type == M1_E2)
    {
      G4double prob = G4RandFlat::shoot(0., 1.);
      if( prob < delta*delta/(delta*delta+1) )
	return(E2);
      else
	return(M1);
    }
  else
    return(decay_type);
}

void MiniBallSource::GetDecayValues(G4int decay_type, 
				    G4double* Ldiff, G4bool* Pdiff)
{
  if(decay_type == E1)
    {
      *Ldiff = 1.;
      *Pdiff = true;
    }
  else if(decay_type == E2)
    {
      *Ldiff = 2.;
      *Pdiff = false;
    }
  else if(decay_type == M1)
    {
      *Ldiff = 1.;
      *Pdiff = false;
    }
  else
    {
      G4cout << "unimplemented decay_type: " << decay_type << G4endl;
      return;
    }
}

G4int MiniBallSource::GetNbOfGammas()
{
  return((G4int)gamma_energy.size());
}

G4ThreeVector MiniBallSource::GetGammaDirection(G4int gamma_nb)
{
  if(gamma_nb<0 || gamma_nb>=(G4int)gamma_direction.size())
    return(G4ThreeVector(0., 0., 1.));
  else
    return(gamma_direction[gamma_nb]);
}

G4double MiniBallSource::GetGammaEnergy(G4int gamma_nb)
{
  if(gamma_nb<0 || gamma_nb>=(G4int)gamma_energy.size())
    return(0.*keV);
  else
    return(gamma_energy[gamma_nb]);
}

G4double MiniBallSource::RandTheta(G4int l, G4int m)
{
  G4double theta;
  G4double func_theta;
  G4double func_rand;
  G4double func_max = GetFmax(l, m);

  do
    {
      theta = G4RandFlat::shoot(0.,THETA_MAX);
      func_theta = GetF(l, m, theta)*sin(theta);
      func_rand = G4RandFlat::shoot(0., func_max);
    }
  while(func_rand > func_theta);

  return(theta);
}

G4double MiniBallSource::GetF(G4int l, G4int m, G4double theta)
{
  G4double sin_theta = sin(theta);
  G4double sin_theta2 = sin_theta*sin_theta;
  G4double cos_theta = cos(theta);
  G4double cos_theta2 = cos_theta*cos_theta;
  G4double cos_theta4 = cos_theta2*cos_theta2;

  if(l == 1)
    {
      //dipole radiation
      if(m==0)
	return(3*sin_theta2);
      else if(m==-1 || m==1)
	return(3./2.*(1+cos_theta2));
      else
	{
	  G4cerr << "MiniBallSource::RandTheta: undefined m-value"
		 << " (l=" << l << ", m=" << m << G4endl;
	  return(0.);
	}
    }
  else if(l == 2)
    {
      //quadropole radiation
      if(m==0)
	return(15.*sin_theta2*cos_theta2);
      else if(m==-1 || m==1)
	return(5./2.*(1-3*cos_theta2+4*cos_theta4));
      else if(m==-2 || m==2)
	return(5./2.*(1-cos_theta4));
      else
	{
	  G4cerr << "MiniBallSource::RandTheta: undefined m-value"
		 << " (l=" << l << ", m=" << m << G4endl;
	  return(0.);
	}
    }
  else
    {
      G4cerr << "MiniBallSource::RandTheta: l-value: " << l
	     << " not implemented" << G4endl;
      return(0.);
    }
}

G4double MiniBallSource::GetFmax(G4int l, G4int m)
{
  if(l == 1)
    {
      //dipole radiation
      if(m==0)
	return(3.);
      else if(m==-1 || m==1)
	return(2.*sin(acos(sqrt(1/3.))));
      else
	{
	  G4cerr << "MiniBallSource::RandTheta: undefined m-value"
		 << " (l=" << l << ", m=" << m << G4endl;
	  return(0.);
	}
    }
  else if(l == 2)
    {
      //quadropole radiation
      if(m==0)
	{
	  G4double sin_acos_s_2_5 = sin(acos(sqrt(2./5.)));
	  G4double sin_acos_s_2_5_3
	    = sin_acos_s_2_5*sin_acos_s_2_5*sin_acos_s_2_5;
	  return(6.*sin_acos_s_2_5_3);
	}
      else if(m==-1 || m==1)
	return(5.);
      else if(m==-2 || m==2)
	return(5./2.);
      else
	{
	  G4cerr << "MiniBallSource::RandTheta: undefined m-value"
		 << " (l=" << l << ", m=" << m << G4endl;
	  return(0.);
	}
    }
  else
    {
      G4cerr << "MiniBallSource::RandTheta: l-value: " << l
	     << " not implemented" << G4endl;
      return(0.);
    }
}

G4double MiniBallSource::GetClebschGordan(G4double j1, G4double m1,
					  G4double j2, G4double m2,
					  G4double J, G4double M)
{
  if(j1 < 0)
    {
      G4cerr << "MiniBallSource::GetClebschGordon: j1 < 0 ( j1="
	     << j1 << " )" << G4endl;
      return(0.);
    }

  if(fabs(m1) > j1)
    {
      G4cerr << "MiniBallSource::GetClebschGordon: |m1| > j1 ( m1="
	     << m1 << ", j1=" << j1 << " )" << G4endl;
      return(0.);
    }

  if(j2 < 0)
    {
      G4cerr << "MiniBallSource::GetClebschGordon: j2 < 0 ( j1="
	     << j2 << " )" << G4endl;
      return(0.);
    }

  if(fabs(m2) > j2)
    {
      G4cerr << "MiniBallSource::GetClebschGordon: |m2| > j2 ( m2="
	     << m2 << ", j2=" << j2 << " )" << G4endl;
      return(0.);
    }

  if(J < 0)
    {
      G4cerr << "MiniBallSource::GetClebschGordon: J < 0 ( J="
	     << J << " )" << G4endl;
      return(0.);
    }

  if(fabs(M) > J)
    {
      G4cerr << "MiniBallSource::GetClebschGordon: |M| > J ( M="
	     << M << ", J=" << J << " )" << G4endl;
      return(0.);
    }

  if(M != m1+m2)
    {
      G4cerr << "MiniBallSource::GetClebschGordon: M != m1+m2 ( M="
	     << M << ", m1=" << m1 << ", m2=" << m2 << " )" << G4endl;
      return(0.);
    }

  if(J < fabs(j1-j2)) 
    {
      G4cerr << "MiniBallSource::GetClebschGordon: J < |j1-j2| ( J="
	     << J << ", j1=" << j1 << ", j2=" << j2 << " )" << G4endl;
      return(0.);
    }

  if(J > j1+j2)
    {
      G4cerr << "MiniBallSource::GetClebschGordon: J > j1+j2 ( J="
	     << J << ", j1=" << j1 << ", j2=" << j2 << " )" << G4endl;
      return(0.);
    }

  G4double ret_value=1.;

  if(M<0)
    {
      G4int exp = (G4int)(j1+j2-J);
      if(exp%2 == 1)
	ret_value=-1;

      M=-M;
      m1=-m1;
      m2=-m2;
    }

  if(j2 == 1)
    {
      if(j1 == 0 && J == 1)
	return(1.);
      else if(j1 == 1 && J == 0 && M == 0)
	{
	  if(m2 == -1)
	    ret_value *= sqrt(1./3.);
	  else if(m2 == 0)
	    ret_value *= -sqrt(1./3.);
	  else if(m2 == 1)
	    ret_value *= sqrt(1./3.);
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 3./2. && J == 1./2. && M == 1./2.)
	{
	  if(m2 == -1)
	    ret_value *= sqrt(1./2.);
	  else if(m2 == 0)
	    ret_value *= -sqrt(1./3.);
	  else if(m2 == 1)
	    ret_value *= sqrt(1./6.);
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 1 && J == 1)
	{
	  if(M == 0)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./2.);
	      else if(m2 == 0)
		ret_value *= 0;
	      else if(m2 == 1)
		ret_value *= -sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(1./2.);
	      else if(m2 == 1)
		ret_value *= -sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 2 && J == 1)
	{
	  if(M == 0)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(3./10.);
	      else if(m2 == 0)
		ret_value *= -sqrt(2./5.);
	      else if(m2 == 1)
		ret_value *= sqrt(3./10.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(3./5.);
	      else if(m2 == 0)
		ret_value *= -sqrt(3./10.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./10.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 3./2. && J == 3./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(2./5.);
	      else if(m2 == 0)
		ret_value *= sqrt(1./15.);
	      else if(m2 == 1)
		ret_value *= -2*sqrt(2./15.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(3./5.);
	      else if(m2 == 1)
		ret_value *= -sqrt(2./5.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 5./2. && J == 3./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(2./5.);
	      else if(m2 == 0)
		ret_value *= -sqrt(2./5.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./5.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(2./3.);
	      else if(m2 == 0)
		ret_value *= -2*sqrt(1./15.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./15.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 1. && J == 2.)
	{
	  if(M == 0)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./6.);
	      else if(m2 == 0)
		ret_value *= sqrt(2./3.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./6.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(1./2.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == 1)
		ret_value *= 1.;
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 2. && J == 2.)
	{
	  if(M == 0)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./2.);
	      else if(m2 == 0)
		ret_value *= 0.;
	      else if(m2 == 1)
		ret_value *= -sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./3.);
	      else if(m2 == 0)
		ret_value *= sqrt(1./6.);
	      else if(m2 == 1)
		ret_value *= -sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(2./3.);
	      else if(m2 == 1)
		ret_value *= -sqrt(1./3.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 3. && J == 2.)
	{
	  if(M == 0)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(2./7.);
	      else if(m2 == 0)
		ret_value *= -sqrt(3./7.);
	      else if(m2 == 1)
		ret_value *= -sqrt(2./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(10./21.);
	      else if(m2 == 0)
		ret_value *= -2*sqrt(2./21.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(5./7.);
	      else if(m2 == 0)
		ret_value *= -sqrt(5./21.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./21.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 3./2. && J == 5./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./10.);
	      else if(m2 == 0)
		ret_value *= sqrt(3./5.);
	      else if(m2 == 1)
		ret_value *= sqrt(3./10.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(2./5.);
	      else if(m2 == 1)
		ret_value *= sqrt(3./5.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 5./2.)
	    {
	      if(m2 == 1)
		ret_value *= 1.;
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 5./2. && J == 5./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -1)
		ret_value *= 4*sqrt(1./35.);
	      else if(m2 == 0)
		ret_value *= sqrt(1./35.);
	      else if(m2 == 1)
		ret_value *= -3*sqrt(2./35.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(2./7.);
	      else if(m2 == 0)
		ret_value *= 3*sqrt(1./35.);
	      else if(m2 == 1)
		ret_value *= -4*sqrt(1./35.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 5./2.)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(5./7.);
	      else if(m2 == 1)
		ret_value *= -sqrt(2./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 7./2. && J == 5./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(5./14.);
	      else if(m2 == 0)
		ret_value *= -sqrt(3./7.);
	      else if(m2 == 1)
		ret_value *= sqrt(3./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == -1)
		ret_value *= 1./2.*sqrt(15./7.);
	      else if(m2 == 0)
		ret_value *= -sqrt(5./14.);
	      else if(m2 == 1)
		ret_value *= 1./2.*sqrt(3./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 5./2.)
	    {
	      if(m2 == -1)
		ret_value *= 1./2.*sqrt(3.);
	      else if(m2 == 0)
		ret_value *= -sqrt(3./14.);
	      else if(m2 == 1)
		ret_value *= 1./2.*sqrt(1./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 2. && J == 3.)
	{
	  if(M == 0)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./5.);
	      else if(m2 == 0)
		ret_value *= sqrt(3./5.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./5.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./15.);
	      else if(m2 == 0)
		ret_value *= 2*sqrt(2./15.);
	      else if(m2 == 1)
		ret_value *= sqrt(2./5.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(1./3.);
	      else if(m2 == 1)
		ret_value *= sqrt(2./3.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3)
	    {
	      if(m2 == 1)
		ret_value *= 1.;
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 3. && J == 3.)
	{
	  if(M == 0)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./2.);
	      else if(m2 == 0)
		ret_value *= 0.;
	      else if(m2 == 1)
		ret_value *= -sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -1)
		ret_value *= 1./2.*sqrt(5./3.);
	      else if(m2 == 0)
		ret_value *= 1./2.*sqrt(1./3.);
	      else if(m2 == 1)
		ret_value *= -sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == -1)
		ret_value *= 1./2.;
	      else if(m2 == 0)
		ret_value *= sqrt(1./3.);
	      else if(m2 == 1)
		ret_value *= -1./2.*sqrt(5./3.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3)
	    {
	      if(m2 == 0)
		ret_value *= 1./2.*sqrt(3.);
	      else if(m2 == 1)
		ret_value *= -1./2.;
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else
	{
	  G4cerr << "MiniBallSource::GetClebschGordon: unimplemented case"
		 << "(j1=" << j1 << ", m1=" << m1
		 << ", j2=" << j2 << ", m2=" << m2
		 << ", J=" << J << ", M=" << M << ")" << G4endl;
	  return(0.);
	}
    }
  else if(j2 == 2)
    {
      if(j1 == 0 && J == 2)
	ret_value *= 1;
      else if(j1 == 1 && J == 2)
	{
	  if(M == 0)
	    {
	      if(m2 == -1)
		ret_value *= 1./2.*sqrt(2.);
	      else if(m2 == 0)
		ret_value *= 0.;
	      else if(m2 == 1)
		ret_value *= -1./2.*sqrt(2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == 0)
		ret_value *= 1./2.*sqrt(2.);
	      else if(m2 == 1)
		ret_value *= -1./3.*sqrt(3./2.);
	      else if(m2 == 2)
		ret_value *= -1./3.*sqrt(3.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == 1)
		ret_value *= 1./3.*sqrt(3.);
	      else if(m2 == 2)
		ret_value *= -1./3.*sqrt(6.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 2 && J == 0 && M == 0)
	{
	  if(m2 == -2)
	    ret_value *= sqrt(1./5.);
	  else if(m2 == -1)
	    ret_value *= -sqrt(1./5.);
	  else if(m2 == 0)
	    ret_value *= sqrt(1./5.);
	  else if(m2 == 1)
	    ret_value *= -sqrt(1./5.);
	  else if(m2 == 2)
	    ret_value *= sqrt(1./5.);
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 5./2. && J == 1./2. && M == 1./2.)
	{
	  if(m2 == -2)
	    ret_value *= sqrt(1./3.);
	  else if(m2 == -1)
	    ret_value *= -2*sqrt(1./15.);
	  else if(m2 == 0)
	    ret_value *= sqrt(1./5.);
	  else if(m2 == 1)
	    ret_value *= -sqrt(2./15.);
	  else if(m2 == 2)
	    ret_value *= sqrt(1./15.);
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 2 && J == 1)
	{
	  if(M == 0)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(2./5.);
	      else if(m2 == -1)
		ret_value *= -sqrt(1./10.);
	      else if(m2 == 0)
		ret_value *= 0;
	      else if(m2 == 1)
		ret_value *= sqrt(1./10.);
	      else if(m2 == 2)
		ret_value *= -sqrt(2./5.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./5.);
	      else if(m2 == 0)
		ret_value *= -sqrt(3./10.);
	      else if(m2 == 1)
		ret_value *= sqrt(3./10.);
	      else if(m2 == 2)
		ret_value *= -sqrt(1./5.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 3 && J == 1)
	{
	  if(M == 0)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./7.);
	      else if(m2 == -1)
		ret_value *= -2*sqrt(2./35.);
	      else if(m2 == 0)
		ret_value *= 3*sqrt(1./35.);
	      else if(m2 == 1)
		ret_value *= -2*sqrt(2./35.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(3./7.);
	      else if(m2 == -1)
		ret_value *= -sqrt(2./7.);
	      else if(m2 == 0)
		ret_value *= sqrt(6./35.);
	      else if(m2 == 1)
		ret_value *= -sqrt(3./35.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./35.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 5./2. && J == 3./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -2)
		ret_value *= 2*sqrt(2./21.);
	      else if(m2 == -1)
		ret_value *= -sqrt(2./105.);
	      else if(m2 == 0)
		ret_value *= -sqrt(2./35.);
	      else if(m2 == 1)
		ret_value *= sqrt(5./21.);
	      else if(m2 == 2)
		ret_value *= -4*sqrt(2./105.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(2./7.);
	      else if(m2 == 0)
		ret_value *= -2*sqrt(3./35.);
	      else if(m2 == 1)
		ret_value *= 3*sqrt(1./35.);
	      else if(m2 == 2)
		ret_value *= -2*sqrt(1./35.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 7./2. && J == 3./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(3./14.);
	      else if(m2 == -1)
		ret_value *= -sqrt(2./7.);
	      else if(m2 == 0)
		ret_value *= 3*sqrt(1./35.);
	      else if(m2 == 1)
		ret_value *= -sqrt(6./35.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./2.);
	      else if(m2 == -1)
		ret_value *= -sqrt(2./7.);
	      else if(m2 == 0)
		ret_value *= sqrt(1./7.);
	      else if(m2 == 1)
		ret_value *= -sqrt(2./35.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./70.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 2 && J == 2)
	{
	  if(M == 0)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(2./7.);
	      else if(m2 == -1)
		ret_value *= sqrt(1./14.);
	      else if(m2 == 0)
		ret_value *= -sqrt(2./7.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./14.);
	      else if(m2 == 2)
		ret_value *= sqrt(2./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(3./7.);
	      else if(m2 == 0)
		ret_value *= -sqrt(1./14.);
	      else if(m2 == 1)
		ret_value *= -sqrt(1./14.);
	      else if(m2 == 2)
		ret_value *= sqrt(3./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(2./7.);
	      else if(m2 == 1)
		ret_value *= -sqrt(3./7.);
	      else if(m2 == 2)
		ret_value *= sqrt(2./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 3 && J == 2)
	{
	  if(M == 0)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(5./14.);
	      else if(m2 == -1)
		ret_value *= -sqrt(1./7.);
	      else if(m2 == 0)
		ret_value *= 0.;
	      else if(m2 == 1)
		ret_value *= sqrt(1./7.);
	      else if(m2 == 2)
		ret_value *= -sqrt(5./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(5./14.);
	      else if(m2 == -1)
		ret_value *= 0.;
	      else if(m2 == 0)
		ret_value *= -sqrt(1./7.);
	      else if(m2 == 1)
		ret_value *= sqrt(2./7.);
	      else if(m2 == 2)
		ret_value *= -sqrt(3./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(5./14.);
	      else if(m2 == 0)
		ret_value *= -sqrt(5./14.);
	      else if(m2 == 1)
		ret_value *= sqrt(3./14.);
	      else if(m2 == 2)
		ret_value *= -sqrt(1./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 4 && J == 2)
	{
	  if(M == 0)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(5./42.);
	      else if(m2 == -1)
		ret_value *= -sqrt(5./21.);
	      else if(m2 == 0)
		ret_value *= sqrt(2./7.);
	      else if(m2 == 1)
		ret_value *= -sqrt(5./21.);
	      else if(m2 == 2)
		ret_value *= sqrt(5./42.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -2)
		ret_value *= 1./3.*sqrt(5./2.);
	      else if(m2 == -1)
		ret_value *= -2./3.*sqrt(5./7.);
	      else if(m2 == 0)
		ret_value *= sqrt(5./21.);
	      else if(m2 == 1)
		ret_value *= -2./3.*sqrt(2./7.);
	      else if(m2 == 2)
		ret_value *= 1./3.*sqrt(5./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == -2)
		ret_value *= 1./3.*sqrt(5.);
	      else if(m2 == -1)
		ret_value *= -1./3.*sqrt(5./2.);
	      else if(m2 == 0)
		ret_value *= sqrt(5./42.);
	      else if(m2 == 1)
		ret_value *= -1./3.*sqrt(5./14.);
	      else if(m2 == 2)
		ret_value *= 1./3.*sqrt(1./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 5./2. && J == 5./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(3./14.);
	      else if(m2 == -1)
		ret_value *= sqrt(6./35.);
	      else if(m2 == 0)
		ret_value *= -2*sqrt(2./35.);
	      else if(m2 == 1)
		ret_value *= 0.;
	      else if(m2 == 2)
		ret_value *= 3*sqrt(3./70.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(3./7.);
	      else if(m2 == 0)
		ret_value *= -sqrt(1./70.);
	      else if(m2 == 1)
		ret_value *= -sqrt(6./35.);
	      else if(m2 == 2)
		ret_value *= 3.*sqrt(3./70.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 5./2.)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(5./14.);
	      else if(m2 == 1)
		ret_value *= -sqrt(3./7.);
	      else if(m2 == 2)
		ret_value *= sqrt(3./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 7./2. && J == 5./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -2)
		ret_value *= 2*sqrt(2./21.);
	      else if(m2 == -1)
		ret_value *= -sqrt(1./14.);
	      else if(m2 == 0)
		ret_value *= -sqrt(1./35.);
	      else if(m2 == 1)
		ret_value *= sqrt(7./30.);
	      else if(m2 == 2)
		ret_value *= -sqrt(2./7.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./3.);
	      else if(m2 == -1)
		ret_value *= 1./2.*sqrt(1./21.);
	      else if(m2 == 0)
		ret_value *= -sqrt(3./14.);
	      else if(m2 == 1)
		ret_value *= 11./2.*sqrt(1./105.);
	      else if(m2 == 2)
		ret_value *= -4.*sqrt(1./105.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 5./2.)
	    {
	      if(m2 == -1)
		ret_value *= 1./2.*sqrt(5./3.);
	      else if(m2 == 0)
		ret_value *= -sqrt(5./14.);
	      else if(m2 == 1)
		ret_value *= 1./2.*sqrt(5./7.);
	      else if(m2 == 2)
		ret_value *= -sqrt(1./21.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 9./2. && J == 5./2.)
	{
	  if(M == 1./2.)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./6.);
	      else if(m2 == -1)
		ret_value *= -sqrt(2./7.);
	      else if(m2 == 0)
		ret_value *= sqrt(2./7.);
	      else if(m2 == 1)
		ret_value *= -2*sqrt(1./21.);
	      else if(m2 == 2)
		ret_value *= sqrt(2./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3./2.)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./3.);
	      else if(m2 == -1)
		ret_value *= -sqrt(1./3.);
	      else if(m2 == 0)
		ret_value *= sqrt(3./14.);
	      else if(m2 == 1)
		ret_value *= -sqrt(2./21.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./42.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 5./2.)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(3./5.);
	      else if(m2 == -1)
		ret_value *= -2.*sqrt(1./15.);
	      else if(m2 == 0)
		ret_value *= sqrt(1./10.);
	      else if(m2 == 1)
		ret_value *= -sqrt(1./35.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./210.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 2 && J == 3)
	{
	  if(M == 0)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./10.);
	      else if(m2 == -1)
		ret_value *= sqrt(2./5.);
	      else if(m2 == 0)
		ret_value *= 0.;
	      else if(m2 == 1)
		ret_value *= -sqrt(2./5.);
	      else if(m2 == 2)
		ret_value *= -sqrt(1./10.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(3./10.);
	      else if(m2 == 0)
		ret_value *= sqrt(1./5.);
	      else if(m2 == 1)
		ret_value *= -sqrt(1./5.);
	      else if(m2 == 2)
		ret_value *= -sqrt(3./10.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(1./2.);
	      else if(m2 == 1)
		ret_value *= 0.;
	      else if(m2 == 2)
		ret_value *= -sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3)
	    {
	      if(m2 == 1)
		ret_value *= sqrt(1./2.);
	      else if(m2 == 2)
		ret_value *= -sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 3 && J == 3)
	{
	  if(M == 0)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./3.);
	      else if(m2 == -1)
		ret_value *= sqrt(2./30.);
	      else if(m2 == 0)
		ret_value *= -2*sqrt(1./15.);
	      else if(m2 == 1)
		ret_value *= sqrt(1./30.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./3.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./6.);
	      else if(m2 == -1)
		ret_value *= 1./2.;
	      else if(m2 == 0)
		ret_value *= -1./2.*sqrt(3./5.);
	      else if(m2 == 1)
		ret_value *= -sqrt(1./30.);
	      else if(m2 == 2)
		ret_value *= sqrt(2./5.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == -1)
		ret_value *= 1./2.*sqrt(5./3.);
	      else if(m2 == 0)
		ret_value *= 0.;
	      else if(m2 == 1)
		ret_value *= -1./2.;
	      else if(m2 == 2)
		ret_value *= sqrt(1./3.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3)
	    {
	      if(m2 == 0)
		ret_value *= 1./2.*sqrt(5./3.);
	      else if(m2 == 1)
		ret_value *= -1./2.*sqrt(5./3.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./6.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 4 && J == 3)
	{
	  if(M == 0)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./3.);
	      else if(m2 == -1)
		ret_value *= -sqrt(1./6.);
	      else if(m2 == 0)
		ret_value *= 0;
	      else if(m2 == 1)
		ret_value *= sqrt(1./6.);
	      else if(m2 == 2)
		ret_value *= -sqrt(1./3.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -2)
		ret_value *= 1./3.*sqrt(7./2.);
	      else if(m2 == -1)
		ret_value *= -1./6.;
	      else if(m2 == 0)
		ret_value *= -1./2.*sqrt(1./3.);
	      else if(m2 == 1)
		ret_value *= 1./3.*sqrt(5./2.);
	      else if(m2 == 2)
		ret_value *= -1./3.*sqrt(2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == -2)
		ret_value *= 1./3.*sqrt(14./5.);
	      else if(m2 == -1)
		ret_value *= 1./6.*sqrt(7./5.);
	      else if(m2 == 0)
		ret_value *= -2.*sqrt(1./15.);
	      else if(m2 == 1)
		ret_value *= 7./6.*sqrt(1./5.);
	      else if(m2 == 2)
		ret_value *= -1./3.;
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(7./15.);
	      else if(m2 == 0)
		ret_value *= -1./2.*sqrt(7./5.);
	      else if(m2 == 1)
		ret_value *= 1./2.*sqrt(3./5.);
	      else if(m2 == 2)
		ret_value *= -sqrt(1./30.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else if(j1 == 2 && J == 4)
	{
	  if(M == 0)
	    {
	      if(m2 == -2)
		ret_value *= sqrt(1./70.);
	      else if(m2 == -1)
		ret_value *= 2*sqrt(2./35.);
	      else if(m2 == 0)
		ret_value *= 3*sqrt(2./35.);
	      else if(m2 == 1)
		ret_value *= 2*sqrt(2./35.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./70.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 1)
	    {
	      if(m2 == -1)
		ret_value *= sqrt(1./14.);
	      else if(m2 == 0)
		ret_value *= sqrt(3./7.);
	      else if(m2 == 1)
		ret_value *= sqrt(3./7.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 2)
	    {
	      if(m2 == 0)
		ret_value *= sqrt(3./14.);
	      else if(m2 == 1)
		ret_value *= 2.*sqrt(1./7.);
	      else if(m2 == 2)
		ret_value *= sqrt(3./14.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 3)
	    {
	      if(m2 == 1)
		ret_value *= sqrt(1./2.);
	      else if(m2 == 2)
		ret_value *= sqrt(1./2.);
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else if(M == 4)
	    {
	      if(m2 == 2)
		ret_value *= 1.;
	      else
		{
		  G4cerr << "MiniBallSource::GetClebschGordon: error"
			 << G4endl;
		  return(0.);
		}
	    }
	  else
	    {
	      G4cerr << "MiniBallSource::GetClebschGordon: error" << G4endl;
	      return(0.);
	    }
	}
      else
	{
	  G4cerr << "MiniBallSource::GetClebschGordon: unimplemented case"
		 << "(j1=" << j1 << ", m1=" << m1
		 << ", j2=" << j2 << ", m2=" << m2
		 << ", J=" << J << ", M=" << M << ")" << G4endl;
	  return(0.);
	}
    }
  else
    {
      G4cerr << "MiniBallSource::GetClebschGordon: unimplemented case"
		 << "(j1=" << j1 << ", m1=" << m1
		 << ", j2=" << j2 << ", m2=" << m2
		 << ", J=" << J << ", M=" << M << ")" << G4endl;
      return(0.);
    }

  return(ret_value);
}

void MiniBallSource::SetInternalConversionCoeffE1(G4double c, G4double en)
{
  conv_coeff_E1.push_back(c);
  conv_en_E1.push_back(en/keV);
}

void MiniBallSource::SetInternalConversionCoeffE2(G4double c, G4double en)
{
  conv_coeff_E2.push_back(c);
  conv_en_E2.push_back(en/keV);
}

G4double MiniBallSource::GetInternalConversionCoeffE1(G4double en)
{
  if(conv_coeff_E1.size() == 0)
    return(0.);

  return(GetInternalConversionCoeff(E1, en));
}

G4double MiniBallSource::GetInternalConversionCoeffE2(G4double en)
{
  if(conv_coeff_E2.size() == 0)
    return(0.);

  return(GetInternalConversionCoeff(E2, en));
}

G4double MiniBallSource::GetInternalConversionCoeff(G4int type, G4double en)
{
  G4double decay_en = en/keV;

  std::vector<G4double>* coeff_p;
  std::vector<G4double>* en_p;

  if(type==E1)
    {
      coeff_p = &conv_coeff_E1;
      en_p = &conv_en_E1;
    }
  else if(type==E2)
    {
      coeff_p = &conv_coeff_E2;
      en_p = &conv_en_E2;
    }
  else
    {
      G4cerr << "MiniBallSource::GetInternalConversionCoeff: type not known"
	     << G4endl;
      return(0.);
    }

  G4double prev_en=-1;
  G4double next_en=9e9;
  G4double prev_coeff=0;
  G4double next_coeff=0;

  for(size_t entry=0; entry<coeff_p->size(); entry++)
    {
      if((*en_p)[entry] <= decay_en && (*en_p)[entry]>prev_en)
	{
	  prev_en = (*en_p)[entry];
	  prev_coeff = (*coeff_p)[entry];
	}
      if((*en_p)[entry] > decay_en && (*en_p)[entry]<next_en)
	{
	  next_en = (*en_p)[entry];
	  next_coeff = (*coeff_p)[entry];
	}
    }
    
  G4double a, b;
  if(prev_en != -1 && next_en != 9e9)
    {
      b = log(next_coeff/prev_coeff)/log(next_en/prev_en);
      a = prev_coeff/pow(prev_en, b);
    }
  else if(prev_en != -1)
    {
      b = -2.5;
      a = prev_coeff/pow(prev_en, b);
    }
  else
    {
      b = -2.5;
      a = next_coeff/pow(next_en, b);
    }

  return(a*pow(decay_en, b));
}

void MiniBallSource::InternalConversionCorrectProbabilities()
{
  for(size_t level_nb=0; level_nb<level_structure.size(); level_nb++)
    {
      MiniBallSourceLevel* level = level_structure[level_nb];

      for(size_t decay_nb=0; decay_nb<level->decay.size(); decay_nb++)
	{
	  G4int dest_level = level->decay[decay_nb];
	  G4double decay_energy
	    = level->energy-level_structure[dest_level]->energy;
	  G4double conv_coeff=0;

	  if(level->decay_type[decay_nb] == E1)
	    conv_coeff = GetInternalConversionCoeffE1(decay_energy);
	  else if(level->decay_type[decay_nb] == E2)
	    conv_coeff = GetInternalConversionCoeffE2(decay_energy);

	  level->prob[decay_nb] *= (1+conv_coeff);
	}
    }
}
	  
void MiniBallSource::ClearStatistics()
{
  for(size_t level=0; level<level_structure.size(); level++)
    {
      level_structure[level]->pop_stat=0L;
      level_structure[level]->feed_stat=0L;
      level_structure[level]->decay_stat.erase(
		      level_structure[level]->decay_stat.begin(),
		      level_structure[level]->decay_stat.end());
      level_structure[level]->conv_stat.erase(
		      level_structure[level]->conv_stat.begin(),
		      level_structure[level]->conv_stat.end());
      for(size_t dec_nb=0;dec_nb<level_structure[level]->decay.size();dec_nb++)
	{
	  level_structure[level]->decay_stat.push_back(0L);
	  level_structure[level]->conv_stat.push_back(0L);
	}
    }
  cleared_statistics=true;
}

void MiniBallSource::PrintStatistics()
{
  G4cout << "MiniBallSource::PrintStatistics:" << G4endl << G4endl;
  G4cout << "   generated events: " << events << G4endl;

  if(events==0 || !cleared_statistics)
    return;

  G4cout << G4endl << "   feeding:" << G4endl;

  for(size_t level=0; level<level_structure.size(); level++)
    G4cout << "      level id: " << level
	   << " energy: " << level_structure[level]->energy/keV
	   << "keV - directly feeded: " << level_structure[level]->pop_stat
	   << " (="
	   << (G4double)level_structure[level]->pop_stat/(G4double)events*100
	   << "%) - feeded: " << level_structure[level]->feed_stat
	   << " (="
	   << (G4double)level_structure[level]->feed_stat/(G4double)events*100
	   << "%)" << G4endl;

  G4cout << G4endl << "   decays:" << G4endl;

  for(size_t level=0; level<level_structure.size(); level++)
    for(size_t dec_nb=0;dec_nb<level_structure[level]->decay.size();dec_nb++)
      {
	G4int dec_level = level_structure[level]->decay[dec_nb];
	G4cout << "      from " << level << " to " << dec_level
	       << " energy: " 
	       << (level_structure[level]->energy
		   -level_structure[dec_level]->energy)/keV
	       << "keV - decays: "
	       << level_structure[level]->decay_stat[dec_nb]
	       << " (="
	       << (G4double)level_structure[level]->decay_stat[dec_nb]
	                                            /(G4double)events*100
	       << "%)" << G4endl;
	if(level_structure[level]->conv_stat[dec_nb]>0L)
	  G4cout << "      from " << level << " to " << dec_level
		 << " internal conversion - decays: "
		 << level_structure[level]->conv_stat[dec_nb]
		 << " (="
		 << (G4double)level_structure[level]->conv_stat[dec_nb]
	                                              /(G4double)events*100
		 << "%)" << G4endl;
      }

  G4cout << G4endl << G4endl;
}

void MiniBallSource::CheckSettings()
{
  for(size_t level=0; level<level_structure.size(); level++)
    {
      size_t nb_of_decay = level_structure[level]->decay.size();
      if(level_structure[level]->decay_type.size() != nb_of_decay
	 || level_structure[level]->prob.size() != nb_of_decay)
	{
	  G4cerr << "MiniBallSource: inconsistent setting for level "
		 << "id: " << level << " energy: "
		 << level_structure[level]->energy/keV << "keV" << G4endl;
	  exit(-1);
	}
      if(level_structure[level]->energy != 0. && nb_of_decay == 0)
	{
	  G4cerr << "MiniBallSource: level without decay information "
		 << "id: " << level << " energy: "
		 << level_structure[level]->energy/keV << "keV" << G4endl;
	  exit(-1);
	}
    }
  checked_settings=true;
}

void MiniBallSourceLevel::AddDecay(G4int p_level_id, G4double p_prob,
				   G4int p_decay_type, G4double p_delta)
{
  decay.push_back(p_level_id);
  prob.push_back(p_prob);
  decay_type.push_back(p_decay_type);
  delta.push_back(p_delta);
}

MiniBallSourceDecayBranch::MiniBallSourceDecayBranch(G4double p_m_i,
						     G4double p_m_f,
						     G4double p_m,
						     G4double p_m1,
						     G4double p_m2,
						     G4double p_prob)
{
  m_i = p_m_i;
  m_f = p_m_f;
  m = p_m;
  m1 = p_m1;
  m2 = p_m2;
  prob = p_prob;
}
