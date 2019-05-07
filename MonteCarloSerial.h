/*
 * MonteCarloSerial.h
 *
 *  Created on: 13 cze 2018
 *      Author: oramus
 */

#ifndef MonteCarloSerial_H_
#define MonteCarloSerial_H_

#include"Particles.h"
#include "PotentialEnergy.h"
#include "MyMPI.h"

class MonteCarloSerial {
private:
	double dx, dy, dr;
	double kBTinv;
	double MAX_RANDOM;
	Particles *particles;
	PotentialEnergy *energy;
	MyMPI *myMPI;
	double totalEp;
	void calcInitialDr();
	double calcContribution( int idx, double xx, double yy );
	double deltaEp( int idx, double oldX, double oldY, double newX, double newY );
	double rnd() {
		return random() * MAX_RANDOM;
	}
public:
	MonteCarloSerial();
	virtual ~MonteCarloSerial();
	void setParticles( Particles *particles );
	void setPotential( PotentialEnergy *energy );
	void calcMC( int draws );
	double calcAvrMinDistance();
	double calcTotalPotentialEnergy();

	double getTotalPotentialEnergy() {
		return totalEp;
	}
	void setMyMPI( MyMPI *myMPI ) {
		this->myMPI = myMPI;
	}
	void setKBTinv( double kBTinv ) {
		this->kBTinv = -kBTinv;
	}

	void shareParticles();
	void gatherParticles();
};

#endif /* MonteCarloSerial_H_ */
