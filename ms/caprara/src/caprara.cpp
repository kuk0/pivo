#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits.h>

#include "caprara.hpp"

#include "inversion_median_alberto.h"

using namespace std;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void caprara(const vector<vector<int> > &genomes,
		int lower_bound, int upper_bound,
		int allmed, int lazy, int &mediancnt,
		vector<vector<int> > &medians){

	int **medians_ptr=NULL,	// medians returned by caprara
		n = 0;				// number of genes
	struct genome_struct *gs[3];
	int *genes;

	if(genomes.size() != 3){
		cout << "caprara: needs exactly 3 genomes!"<<endl;
		exit(1);
	}

	n = genomes[0].size();

	mediancnt = 0;
	medians.clear();

	genes = (int *) malloc(n * sizeof(int));
	for(unsigned i=0; i<3; i++){
		gs[i] = (genome_struct *) malloc(sizeof(genome_struct));
	}

		// fill with data (this is a dirty hack that is not guaranteed to work)
	for(unsigned i=0; i<genomes.size(); i++){
		gs[i]->genes = (int*) &(genomes[i][0]);
	}

	if(upper_bound == INT_MAX)
		upper_bound = 0;
	if(lower_bound == INT_MAX)
		lower_bound = 0;

		// call the median solver
	albert_inversion_median_noncircular_mb ( gs, n, genes, &mediancnt,
		&medians_ptr, lower_bound, upper_bound, 
//dynamic_restrictions,
//		circular, 
allmed, lazy);

	for(unsigned i = 0 ; i < (unsigned)mediancnt ; i++){
		vector<int> temp(n, 0);
		for( unsigned j=0; j<n; j++ ){
			temp[j] = medians_ptr[i][j];
		}
		medians.push_back(temp);
		free(medians_ptr[i]);
	}


	sort(medians.begin(), medians.end());
	unique(medians.begin(), medians.end());
	mediancnt = medians.size();

	for(unsigned i=0; i<3; i++){
		free( gs[i] );
	}
	free(genes);
	free(medians_ptr);
	return;
}

int main(int argc, char *argv[] ){
	vector<vector<int> > genomes;
	vector<vector<int> > medians; 
	int mediancnt;

	srand(time(NULL));

	for( unsigned i=0; i<3; i++ ){
		genomes.push_back(vector<int>(10, 0));
		for( unsigned j=0; j<10; j++ ){
			genomes[i][j]=j+1;
		}
		random_shuffle(genomes[i].begin(), genomes[i].end());
		cout << i<< " : ";
		for( unsigned j=0; j<genomes[i].size(); j++ ){
			cout << genomes[i][j]<<" " ;
		}
		cout << endl;
	}

	caprara(genomes, 0, INT_MAX, ALL_MEDIAN, NOTLAZY, mediancnt,
		medians);


	for( unsigned i=0; i<medians.size(); i++ ){
		cout <<"median "<< i<< " : ";
		for( unsigned j=0; j<medians[i].size(); j++ ){
			cout << medians[i][j]<<" " ;
		}
		cout << endl;
	}

}
