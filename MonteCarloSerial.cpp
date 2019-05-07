/*
 * MonteCarloSerial.cpp
 *
 *  Created on: 13 cze 2018
 *      Author: oramus
 */

#include "MonteCarloSerial.h"
#include "Consts.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <omp.h>

using namespace std;

MonteCarloSerial::MonteCarloSerial() : MAX_RANDOM( 1.0 / ( 1.0 + RAND_MAX )) {
}

MonteCarloSerial::~MonteCarloSerial() {
	// TODO Auto-generated destructor stub
}

void MonteCarloSerial::calcInitialDr() {
	double drMin = calcAvrMinDistance();
	dr = DR_INITIAL_RATIO * drMin;
}

double MonteCarloSerial::calcAvrMinDistance() {
	double drMinSQ = 100000.0;
	double tmp;
	for ( int i = 0; i < particles->getNumberOfParticles(); i++ ) {
		tmp = particles->getDistanceSQToClosest(i);
		if ( tmp < drMinSQ )
			drMinSQ = tmp;
	}
	return sqrt( drMinSQ );
}

void MonteCarloSerial::setParticles( Particles *particles ) {
	this->particles = particles;
	calcInitialDr();
}

void MonteCarloSerial::setPotential( PotentialEnergy *energy ) {
	this->energy = energy;
}

double MonteCarloSerial::calcContribution( int idx, double xx, double yy ) {
	double sum = 0;
	for ( int i = 0; i < idx; i++ ) {
		sum += energy->getPotentialEnergyDistanceSQ( particles->getDistanceBetweenSQ( i, xx, yy ));
	}
	for ( int i = idx+1; i < particles->getNumberOfParticles(); i++ ) {
		sum += energy->getPotentialEnergyDistanceSQ( particles->getDistanceBetweenSQ( i, xx, yy ));
	}
	return sum;
}

double MonteCarloSerial::calcTotalPotentialEnergy() {
	double tmp = 0;

	for ( int i = 0; i < particles->getNumberOfParticles(); i++ )
		tmp += calcContribution( i, particles->getX( i ), particles->getY( i ) );

	totalEp = tmp * 0.5;

	return totalEp;
}

double MonteCarloSerial::deltaEp(int idx, double oldX, double oldY, double newX, double newY ) {
	return calcContribution( idx, newX, newY ) - calcContribution( idx, oldX, oldY );
}

// rozesłanie położeń cząstek z procesu o rank=0 do pozostałych
void MonteCarloSerial::shareParticles() {}

// proces o rank=0 po zakończeniu tej metody musi zawierać
// zaktualizowane pozycje cząstek
void MonteCarloSerial::gatherParticles() {}

void MonteCarloSerial::calcMC( int draws ) {
	int accepted = 0;
	int idx;
	double xnew, ynew, xold, yold, dE, prob;

	for ( int i = 0; i < draws; i++ ) {
// którą z cząstek będzemy próbowali przestawić
		idx = (int)( particles->getNumberOfParticles() * rnd() );
// stara pozycja dla czastki
		xold = particles->getX( idx );
		yold = particles->getY( idx );
// nowa pozycja dla czastki
		xnew = xold + dr * ( rnd() - 0.5 );
		ynew = yold + dr * ( rnd() - 0.5 );

// wyliczamy zmianę energii potencjalnej gdy cząstka idx
// przestawiana jest z pozycji old na new
		dE = deltaEp( idx, xold, yold, xnew, ynew );
// pradopodobieństwo zależy od temperatury
		prob = exp( dE * kBTinv );
// czy zaakceptowano zmianę położenia ?
		if ( rnd() < prob ) {
// tak zaakceptowano -> zmiana położenia i energii
			particles->setXY( idx, xnew, ynew );
			totalEp += dE;
			accepted++;
		}
	}

// zmiana dr jeśli zmian było ponado 50%, to
// dr rośnie, jeśli było mniej, to dr maleje.
	if ( accepted * 2 > draws ) {
		dr *= ( 1.0 + DR_CORRECTION );
	} else {
		dr *= ( 1.0 - DR_CORRECTION );
	}
	if ( dr > DR_MAX )
		dr = DR_MAX;

	if ( dr < DR_MIN ) 
		dr = DR_MIN;
}

