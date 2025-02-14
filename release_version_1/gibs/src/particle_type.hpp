////////////////////////////////////////////////////////////////////////////////
// Copyright � 2017, Battelle Memorial Institute
// All rights reserved.
// 1.	Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any person or 
// 	entity lawfully obtaining a copy of this software and associated documentation files 
// 	(hereinafter �the Software�) to redistribute and use the Software in source and binary 
// 	forms, with or without modification.  Such person or entity may use, copy, modify, 
// 	merge, publish, distribute, sublicense, and/or sell copies of the Software, and may 
// 	permit others to do so, subject to the following conditions:
// �	Redistributions of source code must retain the above copyright notice, this list of 
// 	conditions and the following disclaimers. 
// �	Redistributions in binary form must reproduce the above copyright notice, this list of 
// 	conditions and the following disclaimer in the documentation and/or other materials 
// 	provided with the distribution. 
// �	Other than as used herein, neither the name Battelle Memorial Institute or Battelle may 
// 	be used in any form whatsoever without the express written consent of Battelle.  
// 2.	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
// 	CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
// 	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
// 	MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// 	DISCLAIMED. IN NO EVENT SHALL BATTELLE OR CONTRIBUTORS BE 
// 	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
// 	OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// 	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
// 	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
// 	AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
// 	LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
// 	IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
// 	THE POSSIBILITY OF SUCH DAMAGE.
////////////////////////////////////////////////////////////////////////////////
/*
 * particle_type.hpp
 *
 *  Created on: Apr 25, 2015
 *      Author:  Dennis G. Thomas
 *
 *  @file particle_type.hpp
 *  @author Dennis G. Thomas
 *
 *  @brief Class definitions and function prototypes for particle type class
 */

#ifndef PARTICLE_TYPE_HPP_
#define PARTICLE_TYPE_HPP_

#include "generic.hpp"
#include "ion_type.hpp"
#include "solvent_type.hpp"
#include "constants.hpp"
#include "cavity_grid.hpp"
#include "particle_structs.hpp"
#include "input_parameters.hpp"
#include "particletype_structs.hpp"
#include "convert_functions.hpp"


/**
 * @brief Class for storing and accessing the parameters of each particle type
 *
 * Use defined type CParticleType_t for class declarations
 */

class CParticleType{

	/*Vector of struct_particletype type to store the parameters of the particle types  */
	std::vector<struct_particletype> particle_types;

	/*Integer variable to store the number of particle types */
	int nparticle_types;

public:


	// constructor
	CParticleType(int nparticle_types,const CCavityGrid &cavity_grid);

	// destructor
	~CParticleType(){}

	// read and set methods
	void setParameters(const CIonType_t &ion,const CSolvent_t &solvent,
			int niontype,int nsolvent,double temp,double box_vol);

	void readMolarMassAndSetParameters(const char* fname,int nparticle_types);
	void readAndSetLennardJonesParameters(const char* fname,int nparticle_types);
	void setParticleTypeAvgNum(int i,double num_avg){
		particle_types[i].num_avg = num_avg;
	}


	void readAndSetPIDInitValues(const char* fname);

	void initializeLJCellCavityIndices(int nparticle_types,const CCavityGrid &cavity_grid);
	void addLJCellCavityIndices(int ptcltype_index,int lj_cellindex1D,int cavity_cellindex1D);

	void removeCavityCellIndexFromLJCell(int ptcltype_index,int lj_cellindex1D,int cavity_cellindex1D);

	void setLJSwitch(int ptcltype_index,bool lj){

		this->particle_types[ptcltype_index].lj_switch = lj;
	}

	// get methods

	bool getLJSwitch(int ptcltype_index){
		return this->particle_types[ptcltype_index].lj_switch;
	}

	struct_ljcell_cavityindices getLJCellCavityCellIndices(int ptcltype_index, int ljcellindex_1D) const {
		return this->particle_types[ptcltype_index].ljcell_cavityindices.find(ljcellindex_1D)->second;
	}
	struct_particletype getParticleType (int i) const {return particle_types[i];}
	int getNumParticleTypes()const {return nparticle_types;}
	int getNumCavities(int i)const {return particle_types[i].num_cavities;}
	int getNum(int i) const {return particle_types[i].num;}
	int getNumCavitiesInSegment(int itype, int iseg)const
	{return particle_types[itype].num_cavities_in_segment[iseg];}

	double getMuOld(int i) const {return particle_types[i].mu_old;}
	double getMuEx(int i) const {return particle_types[i].mu_ex;}
	double getBip1(int i)const {return particle_types[i].Bip1;}
	double getEpsi(int i) const {return particle_types[i].epsi;}

	std::string getLabel(int i) const {return particle_types[i].label;}
	double getRadius(int i) const {return particle_types[i].radius;}
	double getConc(int i) const {return particle_types[i].conc;}
	double getCharge(int i) const {return particle_types[i].charge;}

	// update methods
	void updateNumCavities(int i, int num_cavities){particle_types[i].num_cavities = num_cavities;}
	void updateNumTypes(int i, int num){particle_types[i].num = num;}
	void updateNumAvg(int i,double num_avg){particle_types[i].num_avg = num_avg;}
	void updateConcAvg(int i,double conc_avg){particle_types[i].conc_avg = conc_avg;}
	void assignIonTypeParametersToParticle(struct_particle &ptcl, struct_iontype iontype);
	void assignParticleTypeParametersToParticle(struct_particle &ptcl, int i);
	void incrementParticleTypeCountByOne(int itype){ ++particle_types[itype].num;}
	void decrementParticleTypeCountByOne(int itype){ --particle_types[itype].num;}

	void updateChemicalPotential(int iter,const CInputParameters_t &parameters,const CBox_t box);

	void updateChemPotentialPID(int ind1,double vol,double temp);
	void targetChemPotentialPID(int ind1,double vol,double temp);
	void replaceOldChemicalPotentialWithNew();

	void updateChemPotentialPIDwSlothCorr(int ind1,double vol,double temp,double solvdiel,double boxlen);
	void updateChemPotentialPIDtargetNum(int ind1,double vol,double temp);

	void updateNumCavitiesInSegments(int itype,int iseg, int diff){
		particle_types[itype].num_cavities_in_segment[iseg] += diff;
	}
	// write methods
	void writeChemicalPotentialToFile(int iter);
	void writePIDStateParameters(int istep1);

	void resetNumConcAvg(){
		for(int i=0;i<nparticle_types;i++){
			particle_types[i].num_avg = 0;
			particle_types[i].conc_avg = 0;
		}
	}
};

typedef class CParticleType CParticleType_t;


#endif /* PARTICLE_TYPE_HPP_ */
