// BamAlignmnetIterator.cpp
// A class implementation for the ConstBamAlignmentIterator.h

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
#include "BamAlignmentIterator.h"
using namespace std;

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "api/SamSequence.h"
#include "api/SamSequenceDictionary.h"
using namespace BamTools;


// BamAlignmentReader() constructor
BamAlignmentReader :: BamAlignmentReader ( const string & fileName ) :
    m_fileName ( fileName ){}

// BamAlignmentReader::begin()
ConstBamAlignmentIterator BamAlignmentReader::begin() const{
    return ConstBamAlignmentIterator ( m_fileName );
}

// BamAlignmentReader::end()
ConstBamAlignmentIterator BamAlignmentReader::end() const{
    return ConstBamAlignmentIterator ();
}

// ConstBamAlignmentIterator() default constructor
ConstBamAlignmentIterator::ConstBamAlignmentIterator() :
    m_alignmentSeqIdx ( 0 ),
    m_idx ( 0 ),
    m_eoa ( true ){}

// ConstBamAlignmentIterator() constructor w/ fileName
ConstBamAlignmentIterator::ConstBamAlignmentIterator ( const string & fileName ) :
    m_alignmentSeqIdx ( 0 ),
    m_idx ( 0 ),
    m_eoa ( false )
{
    m_bamReader.Open ( fileName );
    bool success = m_bamReader.GetNextAlignment ( m_alignment );
    if ( ! success || m_idx >= m_alignment.AlignedBases.size() ) {
        m_eoa = true;
    }
    m_alignmentSeqIdx++;
}

// ConstBamAlignmentIterator::operator++()
ConstBamAlignmentIterator& ConstBamAlignmentIterator::operator++(int){
    if ( !m_eoa ) {
        if ( m_idx + 1u < m_alignment.AlignedBases.size () ) {
            m_idx++;
        }
        else {
            bool success = m_bamReader.GetNextAlignment ( m_alignment );
            if ( !success || m_idx >= m_alignment.AlignedBases.size () ) {
               m_eoa = true;
            }
            m_alignmentSeqIdx++;
            m_idx = 0;
        }
    }
    return *this;
}

// ConstBamAlignmentIterator::Value::operator*()
ConstBamAlignmentIterator::Value
    ConstBamAlignmentIterator::operator*() const {

    return ConstBamAlignmentIterator::Value (
                m_alignment.AlignedBases[m_idx],
                m_alignment.Qualities[m_idx] );
}

// ConstBamAlignmentIterator::operator==()
bool ConstBamAlignmentIterator::operator==
		( const ConstBamAlignmentIterator & rhs ) const {
    if ( m_eoa != rhs.m_eoa ) {
        return false;
    }
    else if ( !m_eoa) {
        if (m_idx != rhs.m_idx) {
            return false;
        }
        if ( m_alignmentSeqIdx != rhs.m_alignmentSeqIdx ) {
           return false;
        }
    }
    return true;
}

// ConstBamAlignmentIterator::operator!=()
bool ConstBamAlignmentIterator::operator!=
		( const ConstBamAlignmentIterator& rhs) const {
    return ! this->operator == ( rhs);
}
