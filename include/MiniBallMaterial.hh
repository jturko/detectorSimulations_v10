///////////////////////////////////////////////////////////////////////////////
/** \file 
 *  \brief  Materials used for the miniball detectors
 *  \author Hans Boie */
///////////////////////////////////////////////////////////////////////////////

#ifndef MiniBallMaterial_H
#define MiniBallMaterial_H 1

class G4Material;
class G4Element;

class MiniBallMaterial
{
	public:
		static MiniBallMaterial* Get();
		~MiniBallMaterial();

	public:
		G4Material* GetMaterialGe() { return(Ge); };
		G4Material* GetMaterialAl() { return(Al); };
		G4Material* GetMaterialVacuum() { return(Vacuum); };

	private:
		MiniBallMaterial();
		static MiniBallMaterial* fMiniBallMaterial;

		G4Material* Ge;
		G4Material* Al;
		G4Element* N;
		G4Element* O;
		G4Material* Air;
		G4Material* Vacuum;
};

#endif
