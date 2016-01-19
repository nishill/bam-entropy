// BamAlignmentIterator.h
// The header interface for the ConstBamAlignmentIterator class
// and the BamAlignmentReader class

#ifndef __BAMALIGNMENTITERATOR_H_
#define __BAMALIGNMENTITERATOR_H_

#include <cstdio>
#include <string>
#include <stdexcept>
#include <utility>
#include <iostream>
using namespace std;

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
using namespace BamTools;

class ConstBamAlignmentIterator  {
	
	public:
		typedef pair<char,char> Value ;
		bool operator == ( const ConstBamAlignmentIterator & ) const;
    	bool operator != ( const ConstBamAlignmentIterator & ) const;
    	ConstBamAlignmentIterator & operator ++ ( int );
		Value operator * () const;
		
	private:
    	BamReader m_bamReader;
    	BamAlignment m_alignment;
    	size_t m_alignmentSeqIdx;
    	size_t m_idx;
    	bool m_eoa;
    	ConstBamAlignmentIterator ( const string & fileName );
    	ConstBamAlignmentIterator (); 
    	friend class BamAlignmentReader;
};

class BamAlignmentReader {
	public:
    	BamAlignmentReader ( const string & fileName );
    	ConstBamAlignmentIterator begin () const;
    	ConstBamAlignmentIterator end () const; 
	private:
    	string m_fileName;
};

#endif 
