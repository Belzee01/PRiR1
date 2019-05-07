#include<iostream>
#include<time.h>
#include<math.h>
#include"LennardJonesPotential.h"
#include"MonteCarloSerial.h"
#include"MonteCarlo.h"
#include"Particles.h"
#include"MyMPI.h"

using namespace std;

struct Result {
	double calcTime;
	double avrMinDist;
	double EtotInitial;
	double EtotGet;
	double EtotCalc;
	int processes;
	long numberOfCallsPE;
	long numberOfCallsPESQ;
};

void showResult( string txt, Result *result ) {
	cout << txt << "Liczba procesów: ......................... " << result->processes << endl;
    cout << txt << "Etot initial: ............................ " << result->EtotInitial << endl;
	cout << txt << "Etot get: ................................ " << result->EtotGet << endl;
	cout << txt << "Etot calc: ............................... " << result->EtotCalc << endl;
	cout << txt << "Srednia odleglosc do najblizszego sasiada: " << result->avrMinDist << endl;
	cout << txt << "Obliczenia wykonano w czasie: ............ " << result->calcTime << " sekund" << endl;
	cout << txt << "Liczba wywołan getPotentialEnergy: ....... " << result->numberOfCallsPE << endl;
	cout << txt << "Liczba wywołan getPotentialEnergySQ: ..... " << result->numberOfCallsPESQ << endl;
}

Result* calcSerial( Particles *particles, MyMPI *mmpi ) {
	Result *result = new Result();

	LennardJonesPotential *p = new LennardJonesPotential();
	MonteCarloSerial *mc = new MonteCarloSerial();

	mc->setParticles(particles);
	mc->setPotential(p);

    result->EtotInitial = mc->calcTotalPotentialEnergy();

	double kBT = 0.9;

// start pomiaru czasu
    double tStart = mmpi->MPI_Wtime();

	for (int i = 0; i < TEMPERATURES; i++) {
		mc->setKBTinv(kBT); // ustalenie parametrów temperatury
		mc->calcMC( PARTICLES * PARTICLES_PER_TEMP_MULTI );
		kBT += 0.1;
	}

    double tStop = mmpi->MPI_Wtime();
// koniec pomiaru czasu 

    result->EtotGet = mc->getTotalPotentialEnergy();
	result->EtotCalc = mc->calcTotalPotentialEnergy();
	result->avrMinDist = mc->calcAvrMinDistance();
	result->calcTime = tStop - tStart;
	result->processes = 1;
	result->numberOfCallsPE = p->getNumberOfCallsPE();
	result->numberOfCallsPESQ= p->getNumberOfCallsPESQ();

    return result;
}

Result* calc( Particles *particles, int myRank, MyMPI *mmpi ) {
	LennardJonesPotential *p = new LennardJonesPotential();
	MonteCarlo *mc = new MonteCarlo();
	mc->setParticles(particles);
	mc->setPotential(p);
	mc->setMyMPI( mmpi );

	Result *result = new Result();

	if ( myRank == 0 ) {
  	    result->EtotInitial = mc->calcTotalPotentialEnergy();
	}
     	
	double kBT = 0.9;

// start pomiaru czasu
    double tStart = mmpi->MPI_Wtime();

	mc->shareParticles();
	for (int i = 0; i < TEMPERATURES; i++) {
		mc->setKBTinv(kBT); // ustalenie parametrów temperatury
		mc->calcMC( PARTICLES * PARTICLES_PER_TEMP_MULTI );
		kBT += 0.1;
	}
	mc->gatherParticles();

    double tStop = mmpi->MPI_Wtime();
// koniec pomiaru czasu 

	p->shareNumberOfCalls();

	if ( myRank == 0 ) {
 	    mmpi->MPI_Comm_size( MPI_COMM_WORLD, &result->processes );
		result->EtotGet = mc->getTotalPotentialEnergy();
		result->EtotCalc = mc->calcTotalPotentialEnergy();
		result->avrMinDist = mc->calcAvrMinDistance();
		result->calcTime = tStop - tStart;
  		result->numberOfCallsPE = p->getNumberOfCallsPE();
		result->numberOfCallsPESQ= p->getNumberOfCallsPESQ();
	}

	return result;
}

bool equalsAbs( double v1, double v2, double acc ) {
	return fabs(v1-v2) < acc;
}

bool equalsRel( double v1, double v2, double acc ) {

    acc = acc * fabs( v1 + v2 ) / 2.0;

	return fabs(v1-v2) < acc;
}

double efficiency( Result *serial, Result *parallel ) {
	return serial->calcTime / ( parallel->calcTime * parallel->processes );
}

void test( Result *serial, Result *parallel ) {
	// sprawdzamy czy etotcalc zgadza sie z etotget

    double ENERGY_CALC_GET_MAX_RELATIVE_ERROR = 0.00001;
    double ENERGY_CALC_MAX_RELATIVE_ERROR = 0.1;
	double EFFICIENCY = 0.7;

    bool resultOK = true;

	double eff = efficiency( serial, parallel ); 

    cout << "      --------------------------------------------" << endl;
	cout << "      ---- Efektywność obliczeń: " << eff * 100.0 << "%" << endl;
    cout << "      --------------------------------------------" << endl;

	if ( ! equalsRel( serial->EtotCalc, serial->EtotGet, ENERGY_CALC_GET_MAX_RELATIVE_ERROR ) ) {
		cout << "BŁĄD: test należy ponowić z innym seed" << endl;
		return;
	}

	if ( ! equalsRel( parallel->EtotCalc, parallel->EtotGet, ENERGY_CALC_GET_MAX_RELATIVE_ERROR ) ) {
		cout << "BŁĄD: zbyt duża rozbieżność energi get i calc" << endl;
		resultOK = false;
	}

    if ( ! equalsRel(parallel->EtotCalc, serial->EtotCalc, ENERGY_CALC_MAX_RELATIVE_ERROR ) ) {
		cout << "BŁĄD: zbyt duża rozbieżność pomiędzy wynikiem sekwencyjnych i równoległych obliczeń energii" << endl;
		resultOK = false;
	}

	if ( ! equalsAbs( parallel->avrMinDist, serial->avrMinDist, 0.05) ) {
		cout << "BŁĄD: zbyt duża rozbieżność w ułożeniu cząstek" << endl;
		resultOK = false;
	}

	if ( eff < EFFICIENCY ) {
		cout << "BŁĄD: Uzyskana efektywność obliczeń jest niezadawalająca. Jest " << eff <<
		     " limit " << EFFICIENCY << endl;
		resultOK = false;
	}

    double paralleCalcTimeExpected = serial->calcTime / ( EFFICIENCY * parallel->processes );

	if ( parallel->calcTime > paralleCalcTimeExpected ) {
		cout << "BŁĄD: oczekiwano krótszego czasu obliczeń. Jest " << parallel->calcTime 
		<< " limit " << paralleCalcTimeExpected << endl; 
		resultOK = false;
	}

	if ( parallel->numberOfCallsPE != serial->numberOfCallsPE ) {
		cout << "BŁĄD: Liczba wywołan getPotentialEnergy w wersji sekwencyjnej i równoległej jest inna" << endl;
		resultOK = false;
	}

	if ( parallel->numberOfCallsPESQ != serial->numberOfCallsPESQ ) {
		cout << "BŁĄD: Liczba wywołan getPotentialEnergySQ w wersji sekwencyjnej i równoległej jest inna" << endl;
		resultOK = false;
	}

	cout << endl << endl;
	if ( resultOK ) {
		cout << "------------------- OK --- OK --- OK ------------------";
	} else {
		cout << "- BŁĄD - BŁĄD - BŁĄD - BŁĄD - BŁĄD - BŁĄD - BŁĄD - BŁĄD -";
	}
	cout << endl << endl << endl;
}


int main(int argc, char **argv) {
//	srandom( time( NULL ) );
	srandom( SEED );

	MyMPI *mmpi = new MyMPI();
	mmpi->MPI_Init( &argc, &argv );

	int myRank;
	mmpi->MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

	Particles *pa = new Particles(PARTICLES);
	Particles *paBackup;
	pa->setBoxSize( BOX_SIZE );

    Result *serialCalculations;

	if ( myRank == 0 ) {
		// inicjalizacja położeń tylko w jednym z procesów
		pa->initializePositions(BOX_SIZE, DISTANCE);
		paBackup = new Particles( pa ); // tu tworzona jest kopia położeń cząstek
		// przyda się w testach...
		cout << "- s - t - a - r - t -" << endl;
    	serialCalculations = calcSerial( paBackup, mmpi );
	}

	mmpi->MPI_Barrier( MPI_COMM_WORLD );

	cout << "- s - t - a - r - t - - - p - a - r - a - l - l - e - l -" << endl;
	Result *parallelCalculations = calc( pa, myRank, mmpi );

	mmpi->MPI_Barrier( MPI_COMM_WORLD );
	mmpi->MPI_Finalize();

	if ( myRank == 0 ) {
	   cout << endl;
  	   showResult( "Serial>   ", serialCalculations );
	   showResult( "Parallel> ", parallelCalculations );
	   test( serialCalculations, parallelCalculations );
	}
}
