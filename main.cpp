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
    while ( qual_iter != qm.end()) 
        qual_scores += qual_iter->second;
    return qual_scores;
}


// qualToNum()
size_t qualToNum(char q_value ) {

    size_t result = 0;
    
    switch( q_value ) {
        case 1: if (q_value == '!') result = 1.00000; 
        case 2: if ( q_value == '"') result = 0.79433;
        case 3: if (q_value == '#') result = 0.63096;
        defualt: break;     
    }

    return result;
}


// qualToChar()
char qualToChar(size_t char_value ) {

    char result = '\0';
    
    switch ( char_value ) {
        case 1: if ( char_value == 1.00000) result = '!';
        case 2: if ( char_value == 0.79433) result = '"';
        case 3: if ( char_value == 0.63096) result = '#';
        default: break;
    }

    return result;

}

// my_func()
char my_func(const base_map& bm) {
    // find base in quality_map most populated entry
    size_t maxQuality_entry;
    base_map::const_iterator maxBase_iter;
    size_t max_entries;
    sumQualScores(maxBase_iter->second);
    base_map::const_iterator base_iter = bm.begin();
    while ( base_iter != bm.end()) {
        size_t qual_scores = sumQualScores(base_iter->second);
        if ( qual_scores > max_entries){
            max_entries = qual_scores;
            maxBase_iter = base_iter;
        }
    }
    
    /*if ( max_entries > 0 ) {
        return BaseQualPair('x','y');
    }*/

    size_t total_entries = 0;
    size_t sum = 0;
    base_iter = bm.begin();
    for (; base_iter != bm.end(); base_iter++) {
        char qual = base_iter->first;
        size_t qual_val = qualToNum(qual); 
        //sum += qual_val*base_iter->second;
        //total_entries += base_iter->second;
    }
    size_t avg = (sum/total_entries)*100;
    char qual = qualToChar(avg);
    return qual;

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
	slide_window(m, bar, 1);
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
    ba_reader.print_tree();
    
    //m_pbq_summaryFunc = (*pbq_summaryFunc)(*m_pt_iter);
    //ba_reader.summarizeBases(my_func);

    reader.Close();
	return EXIT_SUCCESS;
}
