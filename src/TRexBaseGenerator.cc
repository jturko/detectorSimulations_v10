/*
 * TRexBaseGenerator.cc
 *
 *  Created on: Jun 16, 2014
 *      Author: sklupp
 */

#include "TRexBaseGenerator.hh"

TRexBaseGenerator::TRexBaseGenerator() {
	// TODO Auto-generated constructor stub
}

TRexBaseGenerator::~TRexBaseGenerator() {
	// TODO Auto-generated destructor stub
}

void TRexBaseGenerator::FillTree() {
	fTree->Fill();
}

void TRexBaseGenerator::SetTree(TTree *tree){
	fTree = tree;
}
