/*
 * MonteCarlo.cpp
 *
 *  Created on: 13 cze 2018
 *      Author: oramus
 */

#include "MonteCarlo.h"
#include "Consts.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <omp.h>

using namespace std;

MonteCarlo::MonteCarlo() : MAX_RANDOM( 1.0 / ( 1.0 + RAND_MAX )) {
}

MonteCarlo::~MonteCarlo() {
	// TODO Auto-generated destructor stub
}

void MonteCarlo::calcInitialDr() {
	double drMin = calcAvrMinDistance();
	dr = DR_INITIAL_RATIO * drMin;
}

double MonteCarlo::calcAvrMinDistance() {
	double drMinSQ = 100000.0;
	double tmp;
	for ( int i = 0; i < particles->getNumberOfParticles(); i++ ) {
		tmp = particles->getDistanceSQToClosest(i);
		if ( tmp < drMinSQ )
			drMinSQ = tmp;
	}
	return sqrt( drMinSQ );
}

void MonteCarlo::setParticles( Particles *particles ) {
	this->particles = particles;
	calcInitialDr();
}

void MonteCarlo::setPotential( PotentialEnergy *energy ) {
	this->energy = energy;
}

double MonteCarlo::calcContribution( int idx, double xx, double yy ) {
	double sum = 0;
	for ( int i = 0; i < idx; i++ ) {
		sum += energy->getPotentialEnergyDistanceSQ( particles->getDistanceBetweenSQ( i, xx, yy ));
	}
	for ( int i = idx+1; i < particles->getNumberOfParticles(); i++ ) {
		sum += energy->getPotentialEnergyDistanceSQ( particles->getDistanceBetweenSQ( i, xx, yy ));
	}
	
	return sum;
}

double MonteCarlo::calcContribution2( int idx, double xx, double yy ) {
	int rank;
	int global_size;
	myMPI->MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	myMPI->MPI_Comm_size(MPI_COMM_WORLD, &global_size);

	int separator = this->particles->getNumberOfParticles() / global_size;

	double sum = 0;

	if (rank != global_size-1) {
		for (int i =rank * separator; i < (rank+1)*separator; i++ ) {
			if (i != idx)
				sum += energy->getPotentialEnergyDistanceSQ( particles->getDistanceBetweenSQ( i, xx, yy ));
		}
	} else {
		for (int i = rank * separator; i < this->particles->getNumberOfParticles(); i++ ) {
			if (i != idx)
				sum += energy->getPotentialEnergyDistanceSQ( particles->getDistanceBetweenSQ( i, xx, yy ));
		}
 	}
    return sum;
}

double MonteCarlo::calcTotalPotentialEnergy() {
	double tmp = 0;
	for ( int i = 0; i < particles->getNumberOfParticles(); i++ )
		tmp += calcContribution( i, particles->getX( i ), particles->getY( i ) );

	totalEp = tmp * 0.5;

	return totalEp;
}

double MonteCarlo::deltaEp(int idx, double oldX, double oldY, double newX, double newY ) {
    int rank;
    myMPI->MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
	double contribution = (calcContribution2(idx, newX, newY) - calcContribution2(idx, oldX, oldY));
	double full_contribution;
	myMPI->MPI_Reduce(&contribution, &full_contribution, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	return full_contribution;
}

void MonteCarlo::shareParticles() {
	int rank;
	myMPI->MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double x_buff[this->particles->getNumberOfParticles()];
	double y_buff[this->particles->getNumberOfParticles()];

	//copy data to buffers
	if (rank == 0)
	{
		for ( int i = 0; i < particles->getNumberOfParticles(); i++ ) {
			x_buff[i] = particles->getX(i);
			y_buff[i] = particles->getY(i);
		}
	}

	myMPI->MPI_Bcast(&x_buff,this->particles->getNumberOfParticles(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	myMPI->MPI_Bcast(&y_buff,this->particles->getNumberOfParticles(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	if (rank != 0)
	{
		for ( int i = 0; i < particles->getNumberOfParticles(); i++ ) {
			this->particles->setXY(i, x_buff[i], y_buff[i]);
		}
	}
}

// proces o rank=0 po zakoĹczeniu tej metody musi zawieraÄ
// zaktualizowane pozycje czÄstek
void MonteCarlo::gatherParticles() {
}

void MonteCarlo::calcMC( int draws ) {
	int accepted = 0;
	int idx;
	double xnew, ynew, xold, yold, dE, prob;
	int proc;
	myMPI->MPI_Comm_rank(MPI_COMM_WORLD, &proc);

	for ( int i = 0; i < draws; i++ ) {

// którą z cząstek będzemy próbowali przestawić
		if (proc == 0)
		{
			idx = (int)( particles->getNumberOfParticles() * rnd() );
	// stara pozycja dla czastki
			xold = particles->getX( idx );
			yold = particles->getY( idx );
	// nowa pozycja dla czastki
			xnew = xold + dr * ( rnd() - 0.5 );
			ynew = yold + dr * ( rnd() - 0.5 );
		}

        myMPI->MPI_Bcast(&idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
        myMPI->MPI_Bcast(&xold, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        myMPI->MPI_Bcast(&yold, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        myMPI->MPI_Bcast(&xnew, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        myMPI->MPI_Bcast(&ynew, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		dE = deltaEp( idx, xold, yold, xnew, ynew );

		int setNewPosition = 0;
		if (proc == 0)
		{
			prob = exp( dE * kBTinv );
			if ( rnd() < prob ) {
				setNewPosition = 1;
			} else {
				setNewPosition = 0;
			}
		}		
		myMPI->MPI_Bcast(&setNewPosition, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if ( setNewPosition == 1) {
            particles->setXY(idx, xnew, ynew);
			
			totalEp += dE;
			accepted++;
		}
	}

// zmiana dr jeĹli zmian byĹo ponado 50%, to
// dr roĹnie, jeĹli byĹo mniej, to dr maleje.
	if (proc == 0) {
        if (accepted * 2 > draws) {
            dr *= (1.0 + DR_CORRECTION);
        }
        else {
            dr *= (1.0 - DR_CORRECTION);
        }
        if (dr > DR_MAX)
            dr = DR_MAX;

        if (dr < DR_MIN)
            dr = DR_MIN;
    }
}
