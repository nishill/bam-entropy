// main.cpp

/**************************************************/
// Nicholas Hill								  //
// University of California, Santa Cruz		      //
// 												  //
// This program computes the entropy of a given	  //
// input bam file based on the formula from 	  //
// information theory.							  //
//												  //
// Research done as part of the SSRC/CRSS and 	  //
// The Genome Data Engine project.				  //
//												  //	
/**************************************************/

#include <map> 
#include <cstdio>
#include <cstring>
#include <cctype>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <stdexcept>
#include <cmath>
using namespace std;

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "api/SamSequence.h"
#include "api/SamSequenceDictionary.h"
using namespace BamTools;

#include "BamAlignmentIterator.h"
#include "CompressionIterator.h"

typedef map <string, size_t> G_Map;
typedef map <string, size_t>::const_iterator MapIterator;
typedef pair<char,char> BaseQualPair;

// measure_entropy()
// compute probability then finds 
// the information theory entropy value
double measure_entropy(const G_Map& prob_map, double base_count ) {
	MapIterator iter = prob_map.begin();
	double entropy = 0.0;
	for (; iter != prob_map.end(); iter++) {
		double prob = iter->second/base_count;
		entropy += (prob*log2(prob));   
	}
	return -(entropy);
} 


// slide_window()
// slide across alignments with ConstBamAlignmentIterator
void slide_window(G_Map& lkup, BamAlignmentReader bar, size_t k) { 
	ConstBamAlignmentIterator iter = bar.begin();
	ConstBamAlignmentIterator end  = bar.end();
    while ( iter != end ) {
		string token;
		for (size_t i = 0; i < k; i++ ) {
			ConstBamAlignmentIterator::Value v = *iter;
			token += v.first;
			token += v.second;
			iter++;
			if ( iter == end ) break;	
		}
        G_Map::iterator p_entry = lkup.find(token);
        if ( p_entry == lkup.end() ) {   
        	lkup.insert(pair<string, size_t>(token, 1u));
        	lkup.insert(pair<string, size_t>(token, 2u));
		}   
	 	else {
			p_entry->second++;		
		}
		
	}
	
    //
	// print out map of <base,quality> keys and occurence counts 
	//
    /*for (MapIterator iter = lkup.begin(); iter != lkup.end(); ++iter) {
    	cout << iter->first << " => " << iter->second << endl;
    } */  

}

// sumQualScores()
size_t sumQualScores(const quality_map& qm ) {
    size_t qual_scores = 0;
    quality_map::const_iterator qual_iter = qm.begin();
    while ( qual_iter != qm.end()){ 
        qual_scores += qual_iter->second;
        qual_iter++;
    }
    return qual_scores;
}


// qualToProb()
double qualToProb( char q_ascii ) {
    if ( q_ascii < '!') q_ascii = '!';
    if ( q_ascii > 'K') q_ascii = 'K';
    double q_value = q_ascii - '!';
    double p_value = pow(10.0, -q_value / 10.0);
    return p_value;
}

// probToQual()
char probToQual( double p_value ) {
    assert(p_value >= 0.0 and p_value <= 1.0 ); 
    double q_value = -10.0* log10(p_value);
    // add 0.5 to account for integer rounding
    assert(q_value >= 0.0 and q_value < 41.5 );
    char q_int_value = static_cast<char>(q_value + 0.5);
    char result = '!' + q_int_value;
    return result;
}

// my_func()
BaseQualPair my_func(const base_map& bm) {
    // find base in quality_map most populated entry
    base_map::const_iterator base_iter = bm.begin();
    base_map::const_iterator maxBase_iter = bm.end();
    size_t max_entries = 0;
    while ( base_iter != bm.end()) {
        size_t qual_scores = sumQualScores(base_iter->second);
        if ( qual_scores > max_entries){
            max_entries = qual_scores;
            maxBase_iter = base_iter;
        }
        base_iter++;
    }

    if ( base_iter == maxBase_iter ) return BaseQualPair(NULL, NULL );
   
    size_t total_entries = 0u;
    double sum = 0.0;
    double avg = 0.0;
    quality_map::const_iterator qual_iter = 
        maxBase_iter->second.begin();
    while ( qual_iter != maxBase_iter->second.end() ) {    
        total_entries += qual_iter->second;
        sum += qualToProb(qual_iter->first)*qual_iter->second;
        qual_iter++;
    }
    if ( total_entries > 0 ) {
        avg = sum / total_entries;
    }
    char q_ascii = probToQual(avg);
    return BaseQualPair(maxBase_iter->first, q_ascii);
}

// slide_windowTwo()
// slide across alignments with ConstBamAlignmentIterator
void slide_windowTwo(G_Map& lkup, BA_Reader bar, size_t k) { 
	BA_Reader::ListIterator iter = bar.lbegin(my_func);
	BA_Reader::ListIterator end  = bar.lend();
    while ( iter != end ) {
		string token;
		for (size_t i = 0; i < k; i++ ) {
			BaseQualPair value = *iter;
			token += value.first;
			token += value.second;
			iter++;
			if ( iter == end ) break;	
		}
        G_Map::iterator p_entry = lkup.find(token);
        if ( p_entry == lkup.end() ) {   
        	lkup.insert(pair<string, size_t>(token, 1u));
        	lkup.insert(pair<string, size_t>(token, 2u));
		}   
	 	else {
            p_entry->second++;
        }
    }
    //
	// print out map of <base,quality> keys and occurence counts 
	//
    for (MapIterator iter = lkup.begin(); iter != lkup.end(); ++iter) {
    	cout << iter->first << " => " << iter->second << endl;
    }  
    
}


// main()
int main(int argc, char** argv) {
	string bam_file;
	if ( argc == 2 ) {
		bam_file = argv[1];
	}
	else {
		cout << "Usage: incompatible number of arguments." << endl;
		return EXIT_FAILURE;
	}
   
	BamReader reader;
	reader.Open(bam_file);
	if(!reader.Open(bam_file)){
		cout << "Usage: unable to open BAM file." << endl;	
		return EXIT_FAILURE;
	}
    cout << endl;
	cout << "File: *********" <<bam_file
	<<"*************" <<endl;
	 
    BamAlignmentReader bar ( bam_file );
	G_Map m;
	//slide_window(m, bar, 1);
    ConstBamAlignmentIterator iter = bar.begin ();
	int count = 0;
	while ( iter != bar.end() ) {
    	iter++;		
		count++;
	}	
	
    //
    // calculate and print out entropy value
	//
    /*double chaos = measure_entropy(m, count);
	cout << "The measured entropy value of " << bam_file
		 <<" is " << chaos << "." << endl;
	cout << endl;*/

	//
    // print out alignments w/ qualities
	//
    /*BamAlignment alignment;
	while ( reader.GetNextAlignment(alignment)) {
		cout << alignment.AlignedBases << endl;
		cout << alignment.Qualities << end;
		cout << endl;
	}*/
	
	BA_Reader ba_reader (bam_file);
    //ba_reader.print_tree();


    //ba_reader.summarizeBases(my_func);

    G_Map map;
    slide_windowTwo(map, ba_reader,9);
    
    int countTwo = 0;
    BA_Reader::ListIterator li = ba_reader.lbegin(my_func);
    for (; li != ba_reader.lend(); li++) {
        pair<char, char> iter_pair = *li;       
        //cout << "base= " << iter_pair.first 
        //    << "quality= " << iter_pair.second;
        //li++;
        countTwo++;
    }
    
    double chaosTwo = measure_entropy(map, countTwo);
    cout << chaosTwo << endl;

    reader.Close();
	return EXIT_SUCCESS;
}
