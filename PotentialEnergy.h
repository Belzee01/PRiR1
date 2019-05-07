/*
 * PotentialEnergy.h
 *
 *  Created on: 12 cze 2018
 *      Author: oramus
 */

#ifndef POTENTIALENERGY_H_
#define POTENTIALENERGY_H_

class PotentialEnergy {
private:
	long shareCounter( long value );
protected:
	long numberOfCallsPE;
	long numberOfCallsPESQ;
public:
	PotentialEnergy();
	virtual ~PotentialEnergy();

	virtual double getPotentialEnergy( double distance ) = 0;
	virtual double getPotentialEnergyDistanceSQ( double distanceSQ ) = 0;
	virtual double getForce( double distance ) = 0;
	virtual double getForceDistanceSQ( double distanceSQ ) = 0;

	void shareNumberOfCalls();
	long getNumberOfCallsPE() {
		return numberOfCallsPE;
	}
	long getNumberOfCallsPESQ() {
		return numberOfCallsPESQ;
	}
};

#endif /* POTENTIALENERGY_H_ */
