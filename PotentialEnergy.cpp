/*
 * PotentialEnergy.cpp
 *
 *  Created on: 12 cze 2018
 *      Author: oramus
 */

#include "PotentialEnergy.h"
#include <mpi.h>

PotentialEnergy::PotentialEnergy() {
    numberOfCallsPE = 0;
	numberOfCallsPESQ = 0;
}

PotentialEnergy::~PotentialEnergy() {
}

long PotentialEnergy::shareCounter( long value ) {
    long result;
    MPI_Reduce( &value, &result, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
    return result;
}

void PotentialEnergy::shareNumberOfCalls() {
    numberOfCallsPE = shareCounter( numberOfCallsPE );
    numberOfCallsPESQ = shareCounter( numberOfCallsPESQ );
}