// CompressionIterator.cpp
// A class implementation of CompressionIterator.h

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
#include "CompressionIterator.h"
#include <cmath>
using namespace std;

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "api/SamSequence.h"
#include "api/SamSequenceDictionary.h"
using namespace BamTools;

// Base character content
typedef list<char> BaseContent;

// det_mag() - determine magnitude of quality score
char BA_Reader::det_mag( char q_ascii ) {

    if ( q_ascii < '!') q_ascii = '!';
    if ( q_ascii > 'K') q_ascii = 'K';
    double q_value = q_ascii - '!';
    double p_value = pow(10.0, -q_value / 10.0);

    char strong = 'S';
    char weak = 'W';

    if ( p_value >= 0.00100 ) return strong;
    else if ( p_value < 0.00100 )return weak;

} 


// print_tree()
void BA_Reader::print_tree() {
    
    int32_t curr_pos = m_first_pos;
    PT::iterator pos_iter = m_pt.begin();
    while ( pos_iter != m_pt.end()) {
        cout << "Genome position " << curr_pos  << endl;
        base_map::iterator base_iter = pos_iter->begin();
        while ( base_iter != pos_iter->end()) {
            char curr_base = base_iter->first;
            cout << "\t Current base is " << curr_base <<endl;
            quality_map::iterator qual_iter = base_iter->second.begin();
            while (qual_iter != base_iter->second.end()) {
                cout << "\t\t Quality scores are " << 
                qual_iter->first<< " => "  << qual_iter->second << endl;
                qual_iter++;
            }
            base_iter++;
        }
        pos_iter++;
        curr_pos++;
    }

} 

// m_insert_base()
void BA_Reader::m_insert_base( base_map & p_map, 
                    const BamAlignment & ba, int32_t base_pos ){ 
    char base_key = ba.AlignedBases[base_pos];
    char quality_key = ba.Qualities[base_pos];
    // indicates inconsistency within the bam file
    if ( !isgraph(base_key) or !isgraph(quality_key)) return;
    
    typedef pair < base_map :: iterator, bool > BaseResult;
    BaseResult base_result =
            p_map.emplace(base_key, quality_map());
    typedef pair < quality_map :: iterator, bool > QualityResult;
    QualityResult qr = 
            base_result.first->second.emplace(quality_key, 0u);
    qr.first->second++;
}

// BamAlignmentReader() constructor
BA_Reader::BA_Reader ( const string & fileName ) : 
    m_fileName ( fileName ){

    BamReader bar;
    bar.Open( fileName );
    BamAlignment alignment;
    
    PT::iterator pnext_map = m_pt.end();
    PT::iterator pcurr_map = m_pt.end();
    bool success = bar.GetNextAlignment(alignment);
    if (success ) m_first_pos = alignment.Position;
    while ( success ) {
        int32_t curr_pos = alignment.Position;
        BamAlignment nextAlignment;
        success = bar.GetNextAlignment ( nextAlignment );
        int32_t next_pos;
        if ( success ) { 
            next_pos = nextAlignment.Position;
        }
        else {
            next_pos = -1;
        }
        for ( int32_t i = 0; i < alignment.Length; i++ ) {
            if ( pcurr_map == m_pt.end() ) {
                pcurr_map = m_pt.emplace(m_pt.end());
            } 
            m_insert_base( *pcurr_map, alignment, i );
            if ( i + curr_pos == next_pos ) {
                pnext_map = pcurr_map;   
            }
            pcurr_map++;
        }
        pcurr_map = pnext_map;
        pnext_map = m_pt.end();
        alignment = nextAlignment;
    }

} 

// BamAlignmentReader::begin()
BA_Reader::CompressionIterator BA_Reader::begin() {
    return CompressionIterator (*this);
}

// BamAlignmentReader::end()
BA_Reader::CompressionIterator BA_Reader::end(){
    return CompressionIterator(); 
}

// CompressionIterator::operator==()
bool BA_Reader::CompressionIterator::operator==
        ( const CompressionIterator& rhs) const { 
    if ( m_pba_reader != rhs.m_pba_reader ) {
        return false;
    }
    else if (m_pba_reader == NULL ) {
        return true;
    }
    else {
        return m_pt_iter == rhs.m_pt_iter and
               m_pba_reader == rhs.m_pba_reader and
               m_bm_iter == rhs.m_bm_iter;   

    }
}

// CompressionIterator::operator!=()
bool BA_Reader::CompressionIterator::operator!=
        ( const CompressionIterator& rhs) const {
    return !this->operator== (rhs);
}

// BA_Reader::CompressionIterator() constructor
BA_Reader::CompressionIterator::CompressionIterator(BA_Reader& reader):
    m_pt_iter(reader.m_pt.begin()), m_pba_reader(& reader) {
    while ( m_pt_iter != m_pba_reader->m_pt.end() ) {
        m_bm_iter = m_pt_iter->begin ();
        while ( m_bm_iter != m_pt_iter->end() ) {
            m_qm_iter = m_bm_iter->second.begin ();
            if ( m_qm_iter != m_bm_iter->second.end() ) {
                return;
            }
            m_bm_iter++;
        }
        m_pt_iter++;
    }
    m_pba_reader = NULL;
}
// constructor without arguments 
BA_Reader::CompressionIterator::CompressionIterator():
    m_pba_reader(NULL) {}


// BA_Reader::ComressionIterator::m_next()
void BA_Reader::CompressionIterator::m_next(){
    while( m_pt_iter != m_pba_reader->m_pt.end()) {
        while ( m_bm_iter != m_pt_iter->end()) {
            if ( m_qm_iter != m_bm_iter->second.end()) {
                return;
            }
            m_bm_iter++;
            if (m_bm_iter != m_pt_iter->end()) {
                m_qm_iter = m_bm_iter->second.begin();
            }
        }
        m_pt_iter++;
        if (m_pt_iter != m_pba_reader->m_pt.end()) {
            m_bm_iter = m_pt_iter->begin();
        }
    }
    m_pba_reader = NULL;
}

// operator++(int)
BA_Reader::CompressionIterator & BA_Reader::
    CompressionIterator::operator++( int ){
    m_qm_iter++;
    m_next ();
    return *this;
}
 
// operator*()
BA_Reader::BaseQualInfo BA_Reader::CompressionIterator::operator * (){
    return BaseQualInfo ( m_bm_iter->first, m_qm_iter->first,
                            m_qm_iter->second );
}


// summarizeBases()
/*void BA_Reader::summarizeBases(bqpf_summaryFunc pFunc ) {
    PT::iterator iter = m_pt.begin();
    while ( iter != m_pt.end() ) {
        BaseQualPairFunc bqpf = (*pFunc)(*iter);
        char base = bqpf.first;
        char qual = bqpf.second;
        iter++;
    }
}*/



//////////////////////////////////////
////    ListIterator Class        ////
//////////////////////////////////////

// ListIterator constructor (for begin)
BA_Reader::ListIterator::ListIterator(BA_Reader& read, 
                                        bqpf_summaryFunc sumFunc): 
    m_pt_iter(read.m_pt.begin()), m_pt_iter_end(read.m_pt.end()),
    m_pSummaryFunc(sumFunc), m_baseLookaheadBuf{0,0}, 
    m_qualLookaheadBuf(0), m_inExonRegion(false){}

// ListIterator constructor (for end)
BA_Reader::ListIterator::ListIterator(BA_Reader& read ):
    m_pt_iter(read.m_pt.end()), m_pt_iter_end(read.m_pt.end()),
    m_pSummaryFunc(NULL), m_baseLookaheadBuf{0,0}, 
    m_qualLookaheadBuf(0), m_inExonRegion(false){}  

// ListIterator lbegin()
BA_Reader::ListIterator BA_Reader::lbegin(bqpf_summaryFunc p_func) {
    return ListIterator(*this, p_func);
}

// ListIterator lend()
BA_Reader::ListIterator BA_Reader::lend(){
    return ListIterator(*this);
}

// ListIterator operator++
BA_Reader::ListIterator BA_Reader::ListIterator::operator++ ( int ) {
    BA_Reader::ListIterator temp = *this;
    //m_next();
    m_pt_iter++;
    return temp; 
}

// ListIterator operator++
BA_Reader::ListIterator& BA_Reader::ListIterator::operator++() {
    //m_next();
    m_pt_iter++;
    return *this;
}

// ListIterator operator*
BA_Reader::BaseQualCode BA_Reader::ListIterator::operator*() {
    //BaseQualCode summary = m_current();
    BaseQualCode summary = (*m_pSummaryFunc)(*m_pt_iter);
    //pair < char, char > summary = (*m_pSummaryFunc)(*m_pt_iter);
    assert ( std :: isgraph ( summary.m_base ) );
    assert ( std :: isgraph ( summary.m_quality ) );
    
    return summary;
}

// ListIterator operator==
bool BA_Reader::ListIterator::operator==(const ListIterator& rhs) const{
    return m_pt_iter == rhs.m_pt_iter;
}

// ListIterator operator!=
bool BA_Reader::ListIterator::operator!=(const ListIterator& rhs ) const {
    return !this->operator==(rhs);
}

// m_advance()
void BA_Reader::ListIterator::m_advance() {
    BaseQualCode cur = (*m_pSummaryFunc)(*m_pt_iter);
    m_pt_iter++;
    m_baseLookaheadBuf[0] = m_baseLookaheadBuf[0] << 8u;
    m_qualLookaheadBuf = m_qualLookaheadBuf << 8u;
    m_baseLookaheadBuf[1] = m_baseLookaheadBuf[1] << 8u;
    m_baseLookaheadBuf[0] = m_baseLookaheadBuf[0] | 
                            static_cast<unsigned>(cur.m_base);
    
    m_qualLookaheadBuf = m_qualLookaheadBuf | 
                            static_cast<unsigned>(cur.m_quality);
    const unsigned compare0 = m_baseLookaheadBuf[0] & 0xffffff;
    m_baseLookaheadBuf[1] = m_baseLookaheadBuf[1] | (m_baseLookaheadBuf[0]>>24);
    const unsigned compare1 = m_baseLookaheadBuf[1] & 0xffffff;
    // start codon
    if ( compare0 == m_Methionine ) m_inExonRegion = true;
    // stop codon 1
    if ( compare1 == m_Ochre ) m_inExonRegion = false;
    // stop codon 2
    if ( compare1 == m_Amber ) m_inExonRegion = false;
    // stop codon 3
    if ( compare1 == m_Opal ) m_inExonRegion = false;
    
}

// m_next()
void BA_Reader::ListIterator::m_next() {
    
    while ( m_pt_iter != this->m_pt_iter_end ) {
        m_advance();
    }

}

// m_current()
BA_Reader::BaseQualCode BA_Reader::ListIterator::m_current() {
    //return pair<char,char>(
        //static_cast<char>((m_baseLookaheadBuf[0]>>16) & 0xff),
        //static_cast<char>((m_qualLookaheadBuf>>16) & 0xff));
    return BA_Reader::BaseQualCode(static_cast<char>((m_baseLookaheadBuf[0]>>16) & 0xff),
                        static_cast<char>((m_qualLookaheadBuf>>16) & 0xff));
                        //m_inExonRegion); 


}
